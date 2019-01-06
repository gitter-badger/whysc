#include <iostream> 

#include "Point.h"


typedef WHYSC::GeometryObject::Point<double, 2> Point_2;
typedef WHYSC::GeometryObject::Point<double, 3> Point_3;

int main(int argc, char **argv)
{
    Point_2 p0 = {0.0, 1.0};
    Point_3 p1 = {0.0, 1.0, 2.0};
    std::cout << p0 << std::endl;
    std::cout << p1 << std::endl;
    return 0;
}
