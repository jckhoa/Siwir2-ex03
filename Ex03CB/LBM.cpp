#include "LBM.h"
#include <fstream>
#include <sstream>
#include <iostream>
#include <stdlib.h>



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
    u.SetParams(sizex+2,sizey+2,2);
    flags.SetParams(sizex+2,sizey+2,1);
    density.SetParams(sizex+2,sizey+2,1);

    //Initialization
    u.SetValue(0.0);
    flags.SetValue(1.0);
    density.SetValue(1.0);
    for (int x=1; x<=sizex; ++x){
        for (int y=1; y<=sizey; ++y){
            f(x,y,0) = 4.0/9*density(x,y);
            f(x,y,1) = f(x,y,3) = f(x,y,5) = f(x,y,7) = 1.0/9*density(x,y);
            f(x,y,2) = f(x,y,4) = f(x,y,6) = f(x,y,8) = 1.0/36*density(x,y);
        }
    }

    int vtk_count = 0;
    for (int t=0;t<timesteps;++t){
        realtype sumf = 0;
        for (int x=1; x<=sizex; ++x){
            for (int y=1; y<=sizey; ++y){
                for (int dir=0; dir<9; ++dir)
                    sumf += f(x,y,dir);
            }
        }
        cout<<"sumf = "<<sumf<<endl;
        //cout<<"Before boundary handling"<<endl;
       // f.Print();
        HandleBoundary();
        sumf = 0;
        for (int x=1; x<=sizex; ++x){
            for (int y=1; y<=sizey; ++y){
                for (int dir=0; dir<9; ++dir)
                    sumf += f(x,y,dir);
            }
        }
        cout<<"sumf = "<<sumf<<endl;
        //cout<<"Before streaming"<<endl;
       // f.Print();
        Stream();
        sumf = 0;
        for (int x=1; x<=sizex; ++x){
            for (int y=1; y<=sizey; ++y){
                for (int dir=0; dir<9; ++dir)
                    sumf += f(x,y,dir);
            }
        }
        cout<<"sumf = "<<sumf<<endl;
        //cout<<"After streaming"<<endl;
        //f.Print();
        UpdateDensity();
        UpdateVelocity();
        Collide();

        if (t % vtk_step == 0 && t>0){
            outfilename.str("");
            outfilename<<vtk_file.substr(0,vtk_file.size()-4)<<t<<".vtk";
            WriteVTK(outfilename.str());
        }
      //  f.Print();
    }

}
void LBM::Stream(){
    LBMGrid<realtype> ftemp(sizex+2,sizey+2,numDirection);
    ftemp.SetValue(0);
    for (int x=1; x<=sizex; ++x){
        for (int y=1; y<=sizey; ++y){
            for (int dir=0; dir<numDirection; ++dir){
                ftemp(x,y,dir) = f(x-directionVec9[dir].x,y-directionVec9[dir].y,dir);
            }
        }
    }
    swap(f,ftemp);
}

void LBM::UpdateDensity(){
    for (int x=1; x<=sizex; ++x){
        for (int y=1; y<=sizey; ++y){
            density(x,y) = 0;
            for (int dir=0; dir<numDirection; ++dir){
                density(x,y) += f(x,y,dir);
            }
        }
    }
}
void LBM::UpdateVelocity(){
    for (int x=1; x<=sizex; ++x){
        for (int y=1; y<=sizey; ++y){
            u(x,y,0) = 0.0;
            u(x,y,1) = 0.0;
            for (int dir=0; dir<numDirection; ++dir){
                u(x,y,0) += f(x,y,dir)*directionVec9[dir].x;
                u(x,y,1) += f(x,y,dir)*directionVec9[dir].y;
            }
            u(x,y,0) = u(x,y,0)/density(x,y);
            u(x,y,1) = u(x,y,1)/density(x,y);
        }
    }
}

void LBM::Collide(){
    LBMGrid<realtype> ftemp(sizex+2,sizey+2,numDirection);
    ftemp.SetValue(0);
    for (int x=1; x<=sizex; ++x){
        for (int y=1; y<=sizey; ++y){
            for (int dir=0; dir<numDirection; ++dir){
                realtype cAlphaU = directionVec9[dir].x*u(x,y,0)+directionVec9[dir].y*u(x,y,1);
                realtype u2 = u(x,y,0)*u(x,y,0)+u(x,y,1)*u(x,y,1);
                realtype fequilibrium = talpha9[dir]*(density(x,y) + 3/c/c*cAlphaU
                                            +9.0/2/c/c/c/c*cAlphaU*cAlphaU-3.0/2*u2/2/c/c);
                ftemp(x,y,dir) = (1-omega)*f(x,y,dir)+omega*fequilibrium;
            }
        }
    }
    swap(f,ftemp);
}
void LBM::HandleBoundary(){
    LBMGrid<realtype> ftemp(sizex+2,sizey+2,numDirection);
    for (int x=1; x<=sizex; ++x){
        for (int y=1; y<=sizey; ++y){
            for (int dir=0; dir<9; ++dir){
                ftemp(x,y,dir) = f(x,y,dir);
            }
        }

    }
    //left boundary
    for (int y=1; y <= sizey; ++y){
        ftemp(0,y,2) = f(1,y-1,4);
        ftemp(0,y,1) = f(1,y,5);
        ftemp(0,y,8) = f(1,y+1,6);
    }
    ftemp(0,0,2) = f(1,1,6); //left bottom corner

    //right boundary
    for (int y=1; y <= sizey; ++y){
        ftemp(sizex+1,y,4) = f(sizex,y-1,2);
        ftemp(sizex+1,y,5) = f(sizex,y,1);
        ftemp(sizex+1,y,6) = f(sizex,y+1,8);
    }

    ftemp(sizex+1,0,4) = f(sizex,1,8); //right bottom corner

    //bottom boundary
    for (int x=1; x <= sizex; ++x){
        ftemp(x,0,4) = f(x+1,1,6);
        ftemp(x,0,3) = f(x,1,7);
        ftemp(x,0,2) = f(x-1,1,8);
    }

    //top boundary (moving lid)
    for (int x=1; x<=sizex; ++x){
        //f(x,sizey+1,6) = f(x+1,sizey,4) + 2*talpha9[6]*3/c/c*(directionVec9[6].x*uwx+directionVec9[6].y*uwy);
        //f(x,sizey+1,7) = f(x,sizey,3) + 2*talpha9[7]*3/c/c*(directionVec9[7].x*uwx+directionVec9[7].y*uwy);
        //f(x,sizey+1,8) = f(x-1,sizey,2) + 2*talpha9[8]*3/c/c*(directionVec9[8].x*uwx+directionVec9[8].y*uwy);
        ftemp(x,sizey+1,6) = f(x+1,sizey,4);
        ftemp(x,sizey+1,7) = f(x,sizey,3);
        ftemp(x,sizey+1,8) = f(x-1,sizey,2);

    }

    //f(0,sizey+1,8) = f(1,sizey,4) + 2*talpha9[8]*3/c/c*(directionVec9[8].x*uwx+directionVec9[8].y*uwy); //left top corner
    //f(sizex+1,sizey+1,6) = f(sizex,sizey,2)+ 2*talpha9[6]*3/c/c*(directionVec9[6].x*uwx+directionVec9[6].y*uwy); //right top corner
    ftemp(0,sizey+1,8) = f(1,sizey,4);
    ftemp(sizex+1,sizey+1,6) = f(sizex,sizey,2);
    swap(f,ftemp);
}

template<typename TYPE>
void LBM::WriteVTKdata(ostream& outfile, string type, string paramName, string dataType, TYPE &data){
    outfile<<endl;
    if (type == "SCALARS"){
        outfile<<type<<" "<<paramName<<" "<<dataType<<" 1"<<endl;
        outfile<<"LOOKUP_TABLE default"<<endl;
        for (int x=1; x< data.GetSizeX()-1; ++x){
            for (int y=1; y< data.GetSizeY()-1; ++y){
                outfile<<data(x,y)<<endl;
                //cout<<"("<<x<<" ,"<<y<<")"<<endl;
            }
        }
    }
    else if (type == "VECTORS"){
        outfile<<type<<" "<<paramName<<" "<<dataType<<endl;
        for (int x=1; x< data.GetSizeX()-1; ++x){
            for (int y=1; y< data.GetSizeY()-1; ++y){
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
