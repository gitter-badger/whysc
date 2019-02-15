#include <iostream>

#include "AlgebraKernel.h"


typedef WHYSC::AlgebraKernel<double, int> AK;
typedef AK::Array2d Array2d;
int main(int argc, char **argv)
{
    Array2d a{3, 2};
    a.fill();
    std::cout << a;
    return 0;
}
