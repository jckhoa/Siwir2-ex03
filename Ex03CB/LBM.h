#ifndef LBM_H_INCLUDED
#define LBM_H_INCLUDED

#include <iostream>
#include "LBMGrid.h"

using namespace std;

class LBM{
public:
    LBM():sizex(0),sizey(0),numDirection(0),omega(0.0),timesteps(0),vtk_file(""),vtk_step(0),geometry(""){};
    ~LBM(){};
    void ReadFile(string filename);
    void Solve();

private:
    void Stream();
    void Collide();
    void HandleBoundary();

    inttype sizex, sizey, numDirection;
    inttype timesteps, vtk_step;
    realtype omega;
    string vtk_file;
    string geometry;

    LBMGrid<realtype> f, ftilt;
    LBMGrid<realtype> c;
    LBMGrid<inttype> flags;
    LBMGrid<realtype> density;

};

#endif // LBM_H_INCLUDED
