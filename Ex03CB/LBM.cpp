#include "LBM.h"
#include <fstream>
#include <sstream>
#include <iostream>
#include <stdlib.h>
#include <iomanip>

using namespace std;

#define FLUID 1
#define OBSTACLE 0

enum NAVIGATION{C,E,NE,N,NW,W,SW,S,SE};
realtype directionVec9[9][2] = {{0,0},{1,0},{1,1},{0,1},{-1,1},{-1,0},{-1,-1},{0,-1},{1,-1}};
realtype directionOpp9[9][2] = {{0,0},{-1,0},{-1,-1},{0,-1},{1,-1},{1,0},{1,1},{0,1},{-1,1}};
realtype talpha9[9] = {4.0/9, 1.0/9, 1.0/36, 1.0/9, 1.0/36, 1.0/9, 1.0/36, 1.0/9, 1.0/36};


int LBM::GetNumFluidCells(){
    int numFluidCells = 0;
    //Initialize all f to equilibrium state
    for (int y=0; y<sizey; ++y){
        for (int x=0; x<sizex; ++x){
            if (flags(x,y) == FLUID){
                ++numFluidCells;
            }
        }
    }
    return numFluidCells;
}
void LBM::Solve(){
    if (sizex*sizey == 0){
        cerr<<"The grid size is 0. Cannot proceed"<<endl;
        exit(-1);
    }

    if (flags.GetSize()==0){ //if the geometry is not defined
        sizex += 2;
        sizey += 2;
        flags.SetParams(sizex,sizey,1);
        flags.SetValue(FLUID);
    }
    for (int x=0;x<sizex;++x){  //top and bottom boundary
        flags(x,0) = OBSTACLE;
        flags(x,sizey-1) = OBSTACLE;
    }
    for (int y=1;y<=sizey-2;++y){  //left and right boundary
        flags(0,y) = OBSTACLE;
        flags(sizex-1,y) = OBSTACLE;
    }

    f.SetParams(sizex,sizey,9,STREAM_OPTIMAL);
    ftemp.SetParams(sizex,sizey,9,STREAM_OPTIMAL);
    u.SetParams(sizex,sizey,2);
    density.SetParams(sizex,sizey,1);
    //Initialization
    u.SetValue(0.0);
    density.SetValue(1.0);


    //Initialize all f to equilibrium state
    for (int dir = 0; dir < numDirection; ++dir){
        for (int y=0; y<sizey; ++y){
            for (int x=0; x<sizex; ++x){
                if (flags(x,y) == FLUID){
                    f(x,y,dir) = talpha9[dir]*density(x,y);
                }
            }
        }
    }

    stringstream outfilename;
    for (int t=0;t<=timesteps;++t){
        HandleBoundary();
        Stream();
        UpdateDensity(density);
        UpdateVelocity(u);
        Collide();
        if (t % vtk_step == 0 && t>0){
            outfilename.str("");
            outfilename<<vtk_file.substr(0,vtk_file.size()-4)<<t<<".vtk";
            WriteVTK(outfilename.str());
        }
    }

}
void LBM::Stream(){
    ftemp.SetValue(0);
    for (int dir = 0; dir<numDirection; ++dir){
        for (int y = 0; y < sizey; ++y){
            for (int x = 0; x < sizex; ++x){
                if (flags(x,y) == FLUID){
                    ftemp(x,y,dir) = f(x-directionVec9[dir][0],y-directionVec9[dir][1],dir);
                }
            }
        }
    }
    swap(f,ftemp);
}

void LBM::UpdateDensity(LBMGrid<realtype> &density){
    //break the loops for the sake of performance. The code length is longer but faster
    for (int y = 0; y < sizey; ++y){
        for (int x = 0; x < sizex; ++x){
            if (flags(x,y) == FLUID){
                density(x,y) = 0;
            }
        }
    }

    for (int dir = 0; dir < numDirection; ++dir){
        for (int y = 0; y < sizey; ++y){
            for (int x = 0; x < sizex; ++x){
                if (flags(x,y) == FLUID){
                    density(x,y) += f(x,y,dir);
                }
            }
        }
    }
}
void LBM::UpdateVelocity(LBMGrid<realtype> &u){
    //break the loops for the sake of performance. The code length is longer but faster
    for (int dir = 0; dir < 2; ++dir){
        for (int y = 0; y < sizey; ++y){
            for (int x = 0; x < sizex; ++x){
                if (flags(x,y) == FLUID){
                    u(x,y,dir) = 0.0;
                }
            }
        }
    }
    for (int dir = 0; dir < numDirection; ++dir){
        for (int y = 0; y < sizey; ++y){
            for (int x = 0; x < sizex; ++x){
                if (flags(x,y) == FLUID){
                    u(x,y,0) += f(x,y,dir)*directionVec9[dir][0];
                }
            }
        }
    }
    for (int dir = 0; dir < numDirection; ++dir){
        for (int y = 0; y < sizey; ++y){
            for (int x = 0; x < sizex; ++x){
                if (flags(x,y) == FLUID){
                    u(x,y,1) += f(x,y,dir)*directionVec9[dir][1];
                }
            }
        }
    }
}

void LBM::Collide(){
    for (int dir = 0; dir < numDirection; ++dir){
        for (int y = 0; y < sizey; ++y){
            for (int x = 0; x < sizex; ++x){
                if (flags(x,y) == FLUID){
                    realtype cAlphaU = directionVec9[dir][0]*u(x,y,0)+directionVec9[dir][1]*u(x,y,1);
                    realtype u2 = u(x,y,0)*u(x,y,0)+u(x,y,1)*u(x,y,1);
                    realtype fequilibrium = talpha9[dir]*(density(x,y) + 3*cAlphaU+9.0*cAlphaU*cAlphaU/2-3.0*u2/2);
                    //realtype fequilibrium = talpha9[dir]*density(x,y)*(1 + 3*cAlphaU+9.0*cAlphaU*cAlphaU/2-3.0*u2/2);
                    f(x,y,dir) = (1-omega)*f(x,y,dir)+omega*fequilibrium;
                }
            }
        }
    }
}

void LBM::HandleBoundary(){
    for (int dir = 0; dir<numDirection; ++dir){
        for (int y = 0; y < sizey; ++y){
            for (int x = 0; x < sizex; ++x){
                if (flags(x,y) == OBSTACLE){
                    int neighborX = x + directionVec9[dir][0];
                    int neighborY = y + directionVec9[dir][1];
                    if (neighborX >= 0 && neighborX < sizex
                        && neighborY >=0 && neighborY < sizey
                        && flags(neighborX,neighborY) == FLUID)
                    f(x,y,dir) = f(neighborX,neighborY, (dir+3) % 8 + 1);
                }
            }
        }
    }

    //top boundary (moving lid)
    for (int dir = 0; dir<numDirection; ++dir){
        for (int x = 0; x < sizex; ++x){
            int neighborX = x + directionVec9[dir][0];
            int neighborY = sizey-1 + directionVec9[dir][1];
            if (neighborX >= 0 && neighborX < sizex
                && neighborY >=0 && neighborY < sizey
                && flags(neighborX,neighborY) == FLUID){
                //int oppDir = (dir+3) % 8 + 1;
                //f(x,sizey-1,dir) += - 2*talpha9[oppDir]*3*(directionVec9[oppDir][0] * uwx + directionVec9[oppDir][1] * uwy);
                f(x,sizey-1,dir) += - 2*talpha9[dir]*3*(directionOpp9[dir][0] * uwx + directionOpp9[dir][1] * uwy);
            }
        }
    }

}

template<typename TYPE>
void LBM::WriteVTKdata(ostream& outfile, string type, string paramName, string dataType, TYPE &data){
    outfile<<endl;
    if (type == "SCALARS"){
        outfile<<type<<" "<<paramName<<" "<<dataType<<" 1"<<endl;
        outfile<<"LOOKUP_TABLE default"<<endl;
        for (int y=1; y< data.GetSizeY()-1; ++y){
            for (int x=1; x< data.GetSizeX()-1; ++x){
                outfile<<data(x,y)<<endl;
            }
        }
    }
    else if (type == "VECTORS"){
        outfile<<type<<" "<<paramName<<" "<<dataType<<endl;
        for (int y=1; y< data.GetSizeY()-1; ++y){
            for (int x=1; x< data.GetSizeX()-1; ++x){
                outfile<<data(x,y,0)<<" "<<data(x,y,1)<<" 0"<<endl;
            }
        }
    }
}
void LBM::WriteVTK(string filename){
    ofstream outfile(filename.c_str());
    if (outfile.good()){
        outfile<<"# vtk DataFile Version 4.0"<<endl;
        outfile<<"SiwiRVisFile"<<endl;
        outfile<<"ASCII"<<endl;
        outfile<<"DATASET STRUCTURED_POINTS"<<endl;
        outfile<<"DIMENSIONS "<<sizex-2<<" "<<sizey-2<<" 1"<<endl;
        outfile<<"ORIGIN 0 0 0"<<endl;
        outfile<<"SPACING 1 1 1"<<endl;
        outfile<<"POINT_DATA "<<(sizex-2)*(sizey-2)<<endl;
        WriteVTKdata< LBMGrid<realtype> >(outfile,"SCALARS","flags","double", flags);
        WriteVTKdata< LBMGrid<realtype> >(outfile,"SCALARS","density","double",density);
        WriteVTKdata< LBMGrid<realtype> >(outfile,"VECTORS","velocity","double",u);
        outfile.close();
    }
}
void LBM::ReadFile(string filename){
    ifstream infile(filename.c_str());
    if (infile.good()){
        while (!infile.eof()){
            string line, param;
            stringstream ss;
            getline(infile,line);
            ss<<line;
            ss>>param;
            if (param == "sizex"){
                ss>>sizex;
                if (ss.fail()){
                    cerr<<"Error reading 'sizex' in '"<<filename<<"'"<<endl;
                    exit(-1);
                }
            }
            else if (param == "sizey"){
                ss>>sizey;
                if (ss.fail()){
                    cerr<<"Error reading 'sizey' in '"<<filename<<"'"<<endl;
                    return exit(-1);
                }
            }
            else if (param == "timesteps"){
                ss>>timesteps;
                if (ss.fail()){
                    cerr<<"Error reading 'timesteps' in '"<<filename<<"'"<<endl;
                    return exit(-1);
                }
            }
            else if (param == "omega"){
                ss>>omega;
                if (ss.fail()){
                    cerr<<"Error reading 'omega' in '"<<filename<<"'"<<endl;
                    return exit(-1);
                }
            }
            else if (param == "vtk_file"){
                ss>>vtk_file;
                if (ss.fail()){
                    cerr<<"Error reading 'vtk_file' in '"<<filename<<"'"<<endl;
                    return exit(-1);
                }
            }
            else if (param == "vtk_step"){
                ss>>vtk_step;
                if (ss.fail()){
                    cerr<<"Error reading 'vtk_step' in '"<<filename<<"'"<<endl;
                    return exit(-1);
                }
            }
            else if (param == "geometry"){
                ss>>geometry;
                if (ss.fail()){
                    cerr<<"Error reading 'geometry' in '"<<filename<<"'"<<endl;
                    return exit(-1);
                }
                else{
                    ReadPgm(geometry);
                }
            }
        }
        infile.close();
    }
    else{
        cerr<<"Error reading '"<<filename<<"'"<<endl;
        return exit(-1);
    }
}

//Read and write geometry information to grayscale bitmap format
void LBM::ReadPgm(const std::string filename){
    ifstream infile(filename.c_str());
    string line;
    int rows = 0;
    int cols = 0;
    int index = 0;
    if (infile){
        //read 4 lines of header
        getline(infile,line); //get the first line
        getline(infile,line); //get the second line
        infile>>sizey>>sizex;
        infile>>max_gray_level;
        sizex += 2;
        sizey += 2;
        flags.SetParams(sizex,sizey,1);

        for (int y = 1; y <= sizey-2; ++y){
            for (int x = 1; x <= sizex-2; ++x){
                double value = 0;
                infile>>value;
                if (value == 255)
                    flags(x,y) = FLUID;
                else
                    flags(x,y) = OBSTACLE;
            }
        }

        infile.close();
    }
    else{
        cerr<<"error reading .pgm file"<<endl;
        exit(-1);
    }
}
