#ifndef LBM_H_INCLUDED
#define LBM_H_INCLUDED

#include <iostream>
#include "LBMGrid.h"

using namespace std;

template<typename TYPE>
struct Coordinates{
    Coordinates():x(0),y(0){};
    Coordinates(TYPE _x, TYPE _y):x(_x),y(_y){};
    ~Coordinates(){};
    TYPE x,y;
};

class LBM{
public:
    LBM():sizex(0),sizey(0),numDirection(9),omega(0.0),c(1.0),timesteps(0),uwx(0.0),uwy(0.0)
                    ,vtk_file(""),vtk_step(0),geometry(""){};
    ~LBM(){};
    void ReadFile(string filename);
    void Solve();

private:
    void Stream();
    void Collide();
    void HandleBoundary();

    template<typename TYPE>
    void WriteVTKdata(ostream& outstream, string type, string paramName, string dataType, TYPE &data);
    void WriteVTK(string filename);
    void UpdateDensity();   //update density
    void UpdateVelocity();   //update u

    inttype sizex, sizey, numDirection;
    inttype timesteps, vtk_step;
    realtype omega, c;
    realtype uwx, uwy;
    string vtk_file;
    string geometry;

    LBMGrid<realtype> f;
    LBMGrid<realtype> u;
    LBMGrid<realtype> flags;
    LBMGrid<realtype> density;

};

#endif // LBM_H_INCLUDED
