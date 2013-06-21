#include <iostream>
#include "LBM.h"

using namespace std;

int main()
{
    LBM lbm;
    lbm.ReadFile("params_5x5.dat");
    lbm.Solve();
    return 0;
}
