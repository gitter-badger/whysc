#ifndef Geometry_kernel_h
#define Geometry_kernel_h

#include "Point_2.h"
#include "Vector_2.h"
#include "Point_3.h"
#include "Vector_3.h"
#include "Level_set_function.h"
#include "Bisection_alg.h"

namespace WHYSC {

template<typename F=double, typename I=int>
class Geometry_kernel_base
{
public:
    typedef F Float;
    typedef I Int;
    typedef typename GeometryObject::Point_2<F> Point_2;
    typedef typename GeometryObject::Point_3<F> Point_3;
    typedef typename GeometryObject::Vector_2<F> Vector_2;
    typedef typename GeometryObject::Vector_3<F> Vector_3;

public:
    static Point_2 point_2(const F * p) { return Point_2(p[0], p[1]);}
    static Point_3 point_3(const F * p) { return Point_3(p[0], p[1], p[2]);}
    static Vector_2 vector_2(const F *v) { return Vector_2(v[0], v[1]);}
    static Vector_3 vector_3(const F *v) { return Vector_3(v[0], v[1], v[2]);}
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

    static void midpoint_2(const Float * p1, const Float * p2, Float * p)
    {
        p[0] = (p1[0] + p2[0])/2.0;
        p[1] = (p1[1] + p2[1])/2.0;
    }

    template<typename Point> 
    static Point midpoint(const Point & p1, const Point & p2)
    {
        Point p;
        Int dim = p.dimension();
        for(auto d = 0; d < dim; d++)
            p[d] = (p1[d] + p2[d])/2.0;
        return p;
    }

    static Int sign(Float val)
    {
        return (Float(0) < val) - (val < Float(0));
    }

};

template<typename F=double, typename I=int>
class Geometry_kernel:public Geometry_kernel_base<F, I>
{
public:
    typedef Geometry_kernel_base<F, I>   GK;
    typedef typename GK::Point_2 Point_2;
    typedef typename GK::Point_3 Point_3;
    typedef typename GK::Float Float;
    typedef typename GK::Int Int;

    typedef Float (*Curve_2)(const Point_2 & p);
    typedef Float (*Curve_3)(const Point_3 & p);
    typedef Float (*Surface)(const Point_3 & p);

    typedef GeoAlg::Bisection_alg<GK>  Bisection_algorithm;

    typedef LevelSetFunction::Circle_2<GK> Circle_2;
    typedef LevelSetFunction::Sphere_3<GK> Sphere_3;
    typedef LevelSetFunction::Signed_distance_circle_2<GK> Signed_distance_circle_2;
    typedef LevelSetFunction::Signed_distance_sphere_3<GK> Signed_distance_sphere_2;
};

} // end of namespace WHYSC
#endif
