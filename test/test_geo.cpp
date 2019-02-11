#include <iostream> 

#include "Geometry_kernel.h"

typedef WHYSC::Geometry_kernel<double> GK;
typedef GK::Point_2 Point_2;
typedef GK::Point_3 Point_3;
typedef GK::Vector_2 Vector_2;
typedef GK::Vector_3 Vector_3;

int main(int argc, char **argv)
{
    Point_2 p0 = {0.0, 1.0};
    Point_3 p1 = {0.0, 1.0, 2.0};
    std::cout << p0 << std::endl;
    std::cout << p1 << std::endl;

    Vector_2 v0 = {1.0, 1.0};
    Vector_3 v1 = {3.0, 4.0, 5.0};
    std::cout << v0 << std::endl;
    std::cout << v1 << std::endl;

    p0 += v0;
    std::cout << p0 << std::endl;
    p1 += v1;
    std::cout << p1 << std::endl;

    std::cout << GK::two_thirds() << std::endl;
    std::cout << GK::pi() << std::endl;

    return 0;
}
