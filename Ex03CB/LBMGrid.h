#ifndef LBMGRID_H_INCLUDED
#define LBMGRID_H_INCLUDED
#include <iostream>
#include <iomanip>
#include <stdlib.h>

typedef double realtype;
typedef int inttype;

using namespace std;

enum LbmGridMode{COLLISION_OPTIMAL, STREAM_OPTIMAL};

template<typename TYPE>
class LBMGrid{
public:
    LBMGrid(inttype sizex, inttype sizey, inttype numDirections = 1, LbmGridMode mode = COLLISION_OPTIMAL): sizex(sizex), sizey(sizey), numDirections(numDirections), mode(mode),
                                    dataSize(sizex*sizey*numDirections),data(new TYPE[sizex*sizey*numDirections]),
                                    dummy(0){};
    LBMGrid():sizex(0),sizey(0),numDirections(0),data(0),dataSize(0){};
    void SetParams(inttype _sizex, inttype _sizey, inttype _numDirections = 1, LbmGridMode _mode = COLLISION_OPTIMAL){
        if (_sizex*_sizey*_numDirections >0)
            dataSize = _sizex*_sizey*_numDirections;
        else{
            cerr<<"LBMGrid.SetParams(...) failed the datasize is not a positive number"<<endl;
            exit(-1);
        }
        sizex = _sizex;
        sizey = _sizey;
        mode = _mode;
        numDirections = _numDirections;
        if (data != 0) delete[] data;
        data = new TYPE[dataSize];
    }

    ~LBMGrid(){
        if (data != 0){
            delete[] data;
        }
    };

    TYPE& operator()(inttype x, inttype y, inttype direction = 0){
        if (mode == COLLISION_OPTIMAL){
            int index = y*sizex*numDirections + x*numDirections + direction;
            if (index>=dataSize || index<0){
                cerr<<"The index is out of range "<<endl;
                return dummy;
            }
            else{
                return data[y*sizex*numDirections + x*numDirections + direction];
            }
        }
        else{
            //stream optimal
            int index = direction*sizex*sizey + y*sizex + x;
            if (index>=dataSize || index<0){
                cout<<"data size = "<<dataSize<<endl;
                cerr<<"The index is out of range "<<endl;
                return dummy;
            }
            else{
                return data[direction*sizex*sizey + y*sizex + x]; //
            }
        }
    }

    void SetValue(TYPE value){
        for (int i=0; i<dataSize; ++i)
            data[i] = value;
    }
    inttype GetSize(){return dataSize;}
    inttype GetSizeX(){return sizex;}
    inttype GetSizeY(){return sizey;}

    void swap( LBMGrid<TYPE> & grid )
    {
      std::swap( sizex,grid.sizex );
      std::swap( sizey, grid.sizey );
      std::swap( numDirections, grid.numDirections );
      std::swap( dataSize, grid.dataSize );
      std::swap( data, grid.data );
      std::swap( dummy, grid.dummy );
    }



private:
    inttype sizex,sizey,numDirections;
    inttype dataSize;
    LbmGridMode mode;
    TYPE *data;
    TYPE dummy;
};

template<typename TYPE>
inline void swap( LBMGrid<TYPE>& a, LBMGrid<TYPE>& b )
{
  a.swap(b);
}

#endif // LBMGRID_H_INCLUDED
