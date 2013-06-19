#include <iostream>
#include "LBM.h"

using namespace std;

int main()
{
    LBM lbm;
    lbm.ReadFile("params.dat");
    lbm.Solve();
    return 0;
}
