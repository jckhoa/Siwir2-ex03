#include "LBM.h"
#include <fstream>
#include <sstream>
#include <iostream>
#include <stdlib.h>
#include <iomanip>

enum NAVIGATION{C,E,NE,N,NW,W,SW,S,SE};
Coordinates<inttype> directionVec9[9] ={Coordinates<inttype>(0,0),Coordinates<inttype>(1,0),Coordinates<inttype>(1,1),
        Coordinates<inttype>(0,1),Coordinates<inttype>(-1,1),Coordinates<inttype>(-1,0),
        Coordinates<inttype>(-1,-1),Coordinates<inttype>(0,-1),Coordinates<inttype>(1,-1)};

realtype talpha9[9] = {4.0/9, 1.0/9, 1.0/36,
                        1.0/9, 1.0/36, 1.0/9,
                        1.0/36, 1.0/9, 1.0/36};
using namespace std;

void LBM::Solve(){
    stringstream outfilename;
    f.SetParams(sizex+2,sizey+2,9);
    ftemp.SetParams(sizex+2,sizey+2,9);
    u.SetParams(sizex+2,sizey+2,2);
    flags.SetParams(sizex+2,sizey+2,1);
    density.SetParams(sizex+2,sizey+2,1);

    //Initialization
    u.SetValue(0.0);
    flags.SetValue(1.0);
    density.SetValue(1.0);
    //Initialize all f to equilibrium state
    for (int y=1; y<=sizey; ++y){
        for (int x=1; x<=sizex; ++x){
            f(x,y,C) = 4.0/9*density(x,y);
            f(x,y,E) = f(x,y,N) = f(x,y,W) = f(x,y,S) = 1.0/9*density(x,y);
            f(x,y,NE) = f(x,y,NW) = f(x,y,SW) = f(x,y,SE) = 1.0/36*density(x,y);
        }
    }

    int vtk_count = 0;
    for (int t=0;t<timesteps;++t){
        HandleBoundary();
        Stream();
        UpdateDensity();
        UpdateVelocity();
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
    for (int dir=0; dir<numDirection; ++dir){
        for (int y=1; y<=sizey; ++y){
            for (int x=1; x<=sizex; ++x){
                ftemp(x,y,dir) = f(x-directionVec9[dir].x,y-directionVec9[dir].y,dir);
            }
        }
    }
    swap(f,ftemp);
}

void LBM::UpdateDensity(){
    //break the loops for the sake of performance. The code length is longer but faster
    for (int y=1; y<=sizey; ++y){
        for (int x=1; x<=sizex; ++x){
            density(x,y) = 0;
        }
    }
    for (int dir=0; dir<numDirection; ++dir){
        for (int y=1; y<=sizey; ++y){
            for (int x=1; x<=sizex; ++x){
                density(x,y) += f(x,y,dir);
            }
        }
    }
}
void LBM::UpdateVelocity(){
    //break the loops for the sake of performance. The code length is longer but faster
    for (int dir = 0; dir<2; ++dir){
        for (int y=1; y<=sizey; ++y){
            for (int x=1; x<=sizex; ++x){
                u(x,y,dir) = 0.0;
            }
        }
    }
    for (int dir=0; dir<numDirection; ++dir){
        for (int y=1; y<=sizey; ++y){
            for (int x=1; x<=sizex; ++x){
                u(x,y,0) += f(x,y,dir)*directionVec9[dir].x;
            }
        }
    }
    for (int dir=0; dir<numDirection; ++dir){
        for (int y=1; y<=sizey; ++y){
            for (int x=1; x<=sizex; ++x){
                u(x,y,1) += f(x,y,dir)*directionVec9[dir].y;
            }
        }
    }
    for (int dir = 0; dir<2; ++dir){
        for (int y=1; y<=sizey; ++y){
            for (int x=1; x<=sizex; ++x){
                u(x,y,dir) /= density(x,y);
            }
        }
    }
}

void LBM::Collide(){
    for (int dir=0; dir<numDirection; ++dir){
        for (int y=1; y<=sizey; ++y){
            for (int x=1; x<=sizex; ++x){
                realtype cAlphaU = directionVec9[dir].x*u(x,y,0)+directionVec9[dir].y*u(x,y,1);
                realtype u2 = u(x,y,0)*u(x,y,0)+u(x,y,1)*u(x,y,1);
                realtype fequilibrium = talpha9[dir]*(density(x,y) + 3*cAlphaU
                                            +9.0/2*cAlphaU*cAlphaU-3.0/2*u2);
                f(x,y,dir) = (1-omega)*f(x,y,dir)+omega*fequilibrium;
            }
        }
    }
}
void LBM::HandleBoundary(){

    //left boundary
    for (int y=1; y <= sizey; ++y){
        f(0,y,NE) = f(1,y+1,SW);
        f(0,y,E) = f(1,y,W);
        f(0,y,SE) = f(1,y-1,NW);
    }
    f(0,0,NE) = f(1,1,SW); //left bottom corner

    //right boundary
    for (int y=1; y <= sizey; ++y){
        f(sizex+1,y,NW) = f(sizex,y+1,SE);
        f(sizex+1,y,W) = f(sizex,y,E);
        f(sizex+1,y,SW) = f(sizex,y-1,NE);
    }

    f(sizex+1,0,NW) = f(sizex,1,SE); //right bottom corner

    //bottom boundary
    for (int x=1; x <= sizex; ++x){
        f(x,0,NW) = f(x-1,1,SE);
        f(x,0,N) = f(x,1,S);
        f(x,0,NE) = f(x+1,1,SW);
    }

    //top boundary (moving lid)
    for (int x=1; x<=sizex; ++x){
        f(x,sizey+1,SW) = f(x-1,sizey,NE) - 2*talpha9[NE]*3/c/c*(directionVec9[NE].x * uwx + directionVec9[NE].y * uwy);
        f(x,sizey+1,S) = f(x,sizey,N) - 2*talpha9[N]*3/c/c*(directionVec9[N].x * uwx + directionVec9[N].y * uwy);
        f(x,sizey+1,SE) = f(x+1,sizey,NW) - 2*talpha9[NW]*3/c/c*(directionVec9[NW].x *uwx + directionVec9[NW].y * uwy);

    }

    f(0,sizey+1,SE) = f(1,sizey,NW) - 2*talpha9[NW]*3/c/c*(directionVec9[NW].x*uwx+directionVec9[NW].y*uwy); //left top corner
    f(sizex+1,sizey+1,6) = f(sizex,sizey,2) - 2*talpha9[2]*3/c/c*(directionVec9[2].x*uwx+directionVec9[2].y*uwy); //right top corner

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
                //cout<<"("<<x<<" ,"<<y<<")"<<endl;
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
    cout<<filename<<endl;
    ofstream outfile(filename.c_str());
    if (outfile.good()){
        outfile<<"# vtk DataFile Version 4.0"<<endl;
        outfile<<"SiwiRVisFile"<<endl;
        outfile<<"ASCII"<<endl;
        outfile<<"DATASET STRUCTURED_POINTS"<<endl;
        outfile<<"DIMENSIONS "<<sizex<<" "<<sizey<<" 1"<<endl;
        outfile<<"ORIGIN 0 0 0"<<endl;
        outfile<<"SPACING 1 1 1"<<endl;
        outfile<<"POINT_DATA "<<sizex*sizey<<endl;
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
            }
        }
        infile.close();
    }
    else{
        cerr<<"Error reading '"<<filename<<"'"<<endl;
        return exit(-1);
    }

    cout<<"sizex = "<<sizex<<",sizey = "<<sizey<<endl;
    cout<<"timesteps = "<<timesteps<<endl;
    cout<<"omega = "<<omega<<endl;
    cout<<"vtk_file = "<<vtk_file<<endl;
    cout<<"vtk_step = "<<vtk_step<<endl;
    cout<<"geometry = "<<geometry<<endl;
}
