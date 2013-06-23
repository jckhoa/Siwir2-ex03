#include <iostream>
#include "LBM.h"
#include <sys/time.h>
using namespace std;
int main(int argc, char** argv)
{
    if (argc!=2){
        cerr<<"There should be only one argument"<<endl;
        return -1;
    }
    string filename = argv[1];

    LBM lbm;
    lbm.ReadFile(filename);
    //initialise timing variables
    timeval start;
    timeval end;
    gettimeofday(&start, NULL);

    lbm.Solve();

    //calculate operation time and print it
    gettimeofday(&end, NULL);
    double operation_time = (end.tv_sec - start.tv_sec)*1000000 + end.tv_usec - start.tv_usec;
    operation_time /= 1000000;
    cout << "Operating time (s): " << operation_time << endl;
    cout<<"Number of fluid cells = "<<lbm.GetNumFluidCells()<<endl;
    cout<<"Mega lattice site updates per second = "<<lbm.GetNumFluidCells()/operation_time/1000<<endl;
    return 0;
}
