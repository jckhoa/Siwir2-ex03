#include "LBM.h"
#include <fstream>
#include <sstream>
#include <iostream>
#include <stdlib.h>

using namespace std;

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
