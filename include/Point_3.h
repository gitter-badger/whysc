#ifndef Point_3_h
#define Point_3_h

#include <array>
#include <algorithm>
#include <initializer_list>
#include <assert.h>
#include "Vector_3.h"

namespace WHYSC {
namespace GeometryObject {

template<typename F>
class Point_3 : public std::array<F, 3>
{
public:
    typedef typename std::array<F, 3> Base;
public:

    Point_3()
    {
        std::fill_n(this->data(), 3, 0.0);
    }

    Point_3(const std::initializer_list<F> &l)
    { 
        std::copy_n(l.begin(), 3, this->data());
    }

    Point_3(F x, F y, F z)
    {
        this->data()[0] = x;
        this->data()[1] = y;
        this->data()[2] = z;
    }

    Point_3(const Point_3 & p)
    {
        std::copy_n(p.begin(), 3, this->data());
    }

    static int dimension() {return 3;}

    template<class Vector_3>
    Point_3<F> & operator -= (const Vector_3 & rhs)
    {
        for(int d = 0; d < 3; d++)
            this->data()[d] -= rhs[d];
        return *this;
    }

    template<class Vector>
    Point_3<F> & operator += (const Vector & rhs)
    {
        for(int d = 0; d < 3; d++)
            this->data()[d] += rhs[d];
        return *this;
    }
};

template<typename F>
inline Vector_3<F> operator - (const Point_3<F> & p,
                         const Point_3<F> & q)
{
    Vector_3<F> v;
    for(int d = 0; d < 3; d++)
        v[d] = p[d] - q[d];
    return v;
}

template<typename F>
inline Point_3<F> operator + (const Point_3<F> & p,
                        const Vector_3<F> & v)
{
    Point_3<F> q;
    for(int d = 0; d < 3; d++)
        q[d] = p[d] + v[d]; 
    return q;
}

template<typename F>
inline Point_3<F> operator - (const Point_3<F> & p,
                        const Vector_3<F> & v)
{
    Point_3<F> q;
    for(int d = 0; d < 3; d++)
        q[d] = p[d] - v[d]; 
    return q;
}

template<typename F>
std::ostream& operator << (std::ostream & os, const Point_3<F> & p)
{
    return os << "Point_3(" << p[0] << ", " << p[1] << ", " << p[2] << ')';
}

} // end of namespace GeometryObject
} // end of namespace WHYSC
#endif // end of Point_3_h
