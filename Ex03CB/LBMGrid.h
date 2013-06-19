#ifndef LBMGRID_H_INCLUDED
#define LBMGRID_H_INCLUDED
#include <iostream>
#include <iomanip>

typedef double realtype;
typedef int inttype;

using namespace std;

template<typename TYPE>
class LBMGrid{
public:
    LBMGrid(inttype sizex, inttype sizey, inttype numDirections = 1): sizex(sizex), sizey(sizey), numDirections(numDirections),
                                    dataSize(sizex*sizey*numDirections),data(new TYPE[sizex*sizey*numDirections]),
                                    dummy(0){};
    LBMGrid():sizex(0),sizey(0),numDirections(0),data(0){};
    void SetParams(inttype _sizex, inttype _sizey, inttype _numDirections = 1){
        if (_sizex*_sizey*_numDirections >0)
            dataSize = _sizex*_sizey*_numDirections;
        sizex = _sizex;
        sizey = _sizey;
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
        /*
         //collision optimal
        int index = y*sizex*numDirections + x*numDirections + direction;
        if (index>=dataSize || index<0){
            cerr<<"The index is out of range "<<endl;
            return dummy;
        }
        else
            return data[y*sizex*numDirections + x*numDirections + direction];
        */
        //stream optimal
        int index = direction*sizex*sizey + y*sizex + x;
        if (index>=dataSize || index<0){
            cout<<"data size = "<<dataSize<<endl;
            cerr<<"The index is out of range "<<endl;
            return dummy;
        }
        else
            return data[direction*sizex*sizey + y*sizex + x]; //
    }

    void SetValue(TYPE value){
        for (int i=0; i<dataSize; ++i)
            data[i] = value;
    }
    inttype GetSize(){return dataSize;}
    inttype GetSizeX(){return sizex;}
    inttype GetSizeY(){return sizey;}
    void Print(){
        for (int x=0; x<sizex; ++x){
            for (int y=0; y<sizey; ++y){
                if (numDirections>1){
                    cout<<" (";
                    for (int dir = 0; dir<numDirections-1; ++dir)
                        if ((*this)(x,y,dir)<1.0e-10)
                            cout<<"0.00, ";
                        else
                            cout<<setprecision(3)<< (*this)(x,y,dir)<<",";
                    if ((*this)(x,y,numDirections-1)<1.0e-10)
                        cout<<"0.00, ";
                    else
                    cout<<(*this)(x,y,numDirections-1)<<") ";
                }
                else
                    cout<<" "<<(*this)(x,y,0)<<" ";
                cout<<endl;
            }
        }
    }
private:
    inttype sizex,sizey,numDirections;
    inttype dataSize;
    TYPE *data;
    TYPE dummy;
};

#endif // LBMGRID_H_INCLUDED
