#ifndef Point_2_h
#define Point_2_h

#include <array>
#include <algorithm>
#include <initializer_list>
#include <assert.h>
#include "Vector_2.h"

namespace WHYSC {
namespace GeometryObject {

template<typename F>
class Point_2 : public std::array<F, 2>
{
public:
    typedef typename std::array<F, 2> Base;
public:

    Point_2()
    {
        std::fill_n(this->data(), 2, 0.0);
    }

    Point_2(const std::initializer_list<F> &l)
    { 
        std::copy_n(l.begin(), 2, this->data());
    }

    Point_2(F x, F y)
    {
        this->data()[0] = x;
        this->data()[1] = y;
    }

    Point_2(const Point_2 & p)
    {
        std::copy_n(p.begin(), 2, this->data());
    }

    static int dimension() {return 2;}

    template<class Vector_2>
    Point_2<F> & operator -= (const Vector_2 & rhs)
    {
        for(int d = 0; d < 2; d++)
            this->data()[d] -= rhs[d];
        return *this;
    }

    template<class Vector>
    Point_2<F> & operator += (const Vector & rhs)
    {
        for(int d = 0; d < 2; d++)
            this->data()[d] += rhs[d];
        return *this;
    }
};

template<typename F>
inline Vector_2<F> operator - (const Point_2<F> & p,
                         const Point_2<F> & q)
{
    Vector_2<F> v;
    for(int d = 0; d < 2; d++)
        v[d] = p[d] - q[d];
    return v;
}

template<typename F>
inline Point_2<F> operator + (const Point_2<F> & p,
                        const Vector_2<F> & v)
{
    Point_2<F> q;
    for(int d = 0; d < 2; d++)
        q[d] = p[d] + v[d]; 
    return q;
}

template<typename F>
inline Point_2<F> operator - (const Point_2<F> & p,
                        const Vector_2<F> & v)
{
    Point_2<F> q;
    for(int d = 0; d < 2; d++)
        q[d] = p[d] - v[d]; 
    return q;
}

template<typename F>
std::ostream& operator << (std::ostream & os, const Point_2<F> & p)
{
        return os << "Point_2(" << p[0] << ", " <<p[1] <<')';
}

} // end of namespace GeometryObject
} // end of namespace WHYSC
#endif // end of Point_2_h
