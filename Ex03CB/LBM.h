#ifndef LBM_H_INCLUDED
#define LBM_H_INCLUDED

#include <iostream>
#include "LBMGrid.h"

using namespace std;
/*
template<typename TYPE>
struct Coordinates{
    Coordinates():x(0),y(0){};
    Coordinates(TYPE _x, TYPE _y):x(_x),y(_y){};
    ~Coordinates(){};
    TYPE x,y;
};
*/


class LBM{
public:
    LBM():sizex(0),sizey(0),max_gray_level(0),numDirection(9),omega(0.0),timesteps(0),uwx(0.08),uwy(0.0)
                    ,vtk_file(""),vtk_step(0),geometry(""){};
    ~LBM(){};
    void ReadFile(string filename);
    void Solve();
    int GetNumFluidCells();
private:
    void Stream();
    void Collide();
    void HandleBoundary();

    template<typename TYPE>
    void WriteVTKdata(ostream& outstream, string type, string paramName, string dataType, TYPE &data);
    void WriteVTK(string filename);
    void UpdateDensity(LBMGrid<realtype> &density);   //update density
    void UpdateVelocity(LBMGrid<realtype> &u);   //update u
    realtype GetSum1D(LBMGrid<realtype> &density){
        realtype sum = 0;
        for (int x=1; x<density.GetSizeX()-1; ++x)
            for (int y=1; y<density.GetSizeY()-1; ++y)
                sum += density(x,y);
        return sum;
    }
    realtype GetDiff1D(LBMGrid<realtype> &density, LBMGrid<realtype> &densitytemp){
        realtype sum = 0;
        for (int x=1; x<density.GetSizeX()-1; ++x)
            for (int y=1; y<density.GetSizeY()-1; ++y)
                sum += (density(x,y)-densitytemp(x,y))*(density(x,y)-densitytemp(x,y));
        return sum;
    }
    realtype GetDiff2D(LBMGrid<realtype> &u, LBMGrid<realtype> &utemp){
        realtype sum = 0;
        for (int x=1; x<u.GetSizeX()-1; ++x)
            for (int y=1; y<u.GetSizeY()-1; ++y)
                sum += (u(x,y,0)-utemp(x,y,0))*(u(x,y,1)-utemp(x,y,1));
        return sum;
    }
    realtype GetSumF(){
        realtype sum = 0;
        for (int dir=0; dir< numDirection; ++dir){
            for (int x=1; x<=sizex; ++x)
                for (int y=1; y<=sizey; ++y)
                    sum += f(x,y,dir);
        }
        return sum;
    }
    realtype GetSumVelocity(LBMGrid<realtype> &u){
        realtype sum = 0;
        for (int dir=0; dir< 2; ++dir){
            for (int x=1; x<=sizex; ++x)
                for (int y=1; y<=sizey; ++y)
                    sum += u(x,y,dir)*u(x,y,dir);
        }
        return sum;
    }

    //Read and write geometry information to grayscale bitmap format
    void ReadPgm(const std::string filename);

    inttype sizex, sizey, numDirection;
    inttype timesteps, vtk_step;
    inttype max_gray_level;
    realtype omega;
    realtype uwx, uwy;
    string vtk_file;
    string geometry;

    LBMGrid<realtype> f, ftemp;
    LBMGrid<realtype> u;
    LBMGrid<realtype> flags;
    LBMGrid<realtype> density;

};

#endif // LBM_H_INCLUDED
