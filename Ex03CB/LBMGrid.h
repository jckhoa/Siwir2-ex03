#ifndef LBMGRID_H_INCLUDED
#define LBMGRID_H_INCLUDED

typedef double realtype;
typedef int inttype;

template<typename TYPE>
class LBMGrid{
public:
    LBMGrid(inttype sizex, inttype sizey, inttype numDirections = 1): sizex(sizex), sizey(sizey), numDirections(numDirections),
                                    data(new TYPE[sizex*sizey*numDirections]){};
    LBMGrid():sizex(0),sizey(0),numDirections(0),data(0){};
    void SetParams(inttype _sizex, inttype _sizey, inttype _numDirections = 1){
        sizex = _sizex;
        sizey = _sizey;
        numDirections = _numDirections;
        if (data != 0) delete[] data;
        data = new TYPE[sizex*sizey*numDirections];
    }

    ~LBMGrid(){
        if (data != 0){
            delete[] data;
        }
    };

    TYPE& operator()(inttype x, inttype y, inttype direction = 0){
        return data[y*sizex*numDirections + x*numDirections + direction]; //collision optimal
        //return data[direction*sizex*sizey + y*sizex + x]; //properation optimal
    }

private:
    inttype sizex,sizey,numDirections;
    TYPE *data;
};

#endif // LBMGRID_H_INCLUDED
