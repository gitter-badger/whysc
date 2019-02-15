#ifndef GeometryKernel_h
#define GeometryKernel_h

#include "Point.h"
#include "Vector.h"
#include "Level_set_function.h"

namespace WHYSC {

template<typename F=double, typename I=int>
class Geometry_kernel_base
{
public:
    typedef F Float;
    typedef I Int;
    typedef typename GeometryObject::Point<F, 2> Point_2;
    typedef typename GeometryObject::Point<F, 3> Point_3;
    typedef typename GeometryObject::Vector<F, 2> Vector_2;
    typedef typename GeometryObject::Vector<F, 3> Vector_3;

    typedef Float (*Curve_2)(const Point_2 & p);
    typedef Float (*Curve_3)(const Point_3 & p);
    typedef Float (*Surface)(const Point_3 & p);

    auto& Level_set_function = LevelSetFunction::Level_set_function<GK> Level_set_function;
    typedef LevelSetFunction::Circle<GK> Level_set_circle;
    typedef LevelSetFunction::Sphere<GK> Level_set_sphere;
    typedef LevelSetFunction::Signed_distance_circle<GK> Signed_distance_circle;
    typedef LevelSetFunction::Signed_distance_sphere<GK> Signed_distance_sphere;

    typedef LevelSetFunction::Union<GK, typename GK::Point_2> Union_2;
    typedef LevelSetFunction::Union<GK, typename GK::Point_3> Union_3;
public:
    static Point_2 point_2(const F * p) { return Point_2{p[0], p[1]};}
    static Point_3 point_3(const F * p) { return Point_3{p[0], p[1], p[2]};}
    static Vector_2 vector_2(const F *v) { return Vector_2{v[0], v[1]};}
    static Vector_3 vector_3(const F *v) { return Vector_3{v[0], v[1], v[2]};}
    static const F two_thirds() { return F(2.0)/F(3.0);}
    static const F one_thirds() { return F(1.0)/F(3.0);}
    static const F pi()  {return F(3.1415926535897931e+0);}
    static const F eps() {return F(1e-12);} 

    static void midpoint_3(const Float * p1, const Float * p2, Float * p)
    {
        p[0] = (p1[0] + p2[0])/2.0;
        p[1] = (p1[1] + p2[1])/2.0;
        p[2] = (p1[2] + p2[2])/2.0;
    }

    static void midpoint_2(const double * p1, const double * p2, double * p)
    {
        p[0] = (p1[0] + p2[0])/2.0;
        p[1] = (p1[1] + p2[1])/2.0;
    }

    static Point_2 midpoint(const Point_2 & p1, const Point_2 & p2)
    {
        Point_2 p{0.0, 0.0};
        p[0] = (p1[0] + p2[0])/2.0;
        p[1] = (p1[1] + p2[1])/2.0;
        return p;
    }

    static Point_3 midpoint(const Point_3 & p1, const Point_3 & p2)
    {
        Point_3 p{0.0, 0.0};
        p[0] = (p1[0] + p2[0])/2.0;
        p[1] = (p1[1] + p2[1])/2.0;
        p[2] = (p1[2] + p2[2])/2.0;
        return p;
    }

    static Float sphere(const Point_3 & p)
    {
        return 
    }

};

} // end of namespace WHYSC
#endif
