#ifndef Geometry_kernel_h
#define Geometry_kernel_h

#include "Point.h"
#include "Vector.h"

namespace WHYSC {

template<typename F=double>
class Geometry_kernel
{
public:
    typedef typename GeometryObject::Point<F, 2> Point_2;
    typedef typename GeometryObject::Point<F, 3> Point_3;
    typedef typename GeometryObject::Vector<F, 2> Vector_2;
    typedef typename GeometryObject::Vector<F, 3> Vector_3;
public:
    static Point_2 point_2(const F * p) { return Point_2{p[0], p[1]};}
    static Point_3 point_3(const F * p) { return Point_3{p[0], p[1], p[2]};}
    static Vector_2 vector_2(const F *v) { return Vector_2{v[0], v[1]};}
    static Vector_3 vector_3(const F *v) { return Vector_3{v[0], v[1], v[2]};}
    static const F two_thirds() { return F(2.0)/F(3.0);}
    static const F one_thirds() { return F(1.0)/F(3.0);}
    static const F pi()  {return F(3.1415926535897931e+0);}
    static const F eps() {return F(1e-12);} 
};

} // end of namespace WHYSC
#endif
