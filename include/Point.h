#ifndef Point_h
#define Point_h

#include <array>
#include <algorithm>
#include <initializer_list>
#include <assert.h>
#include "Vector.h"

namespace WHYSC {
namespace GeometryObject {

template<typename F, int DIM>
class Point : public std::array<F, DIM>
{
public:
    typedef typename std::array<F, DIM> Base;
public:

    Point()
    {
        std::fill_n(this->data(), DIM, 0.0);
    }

    Point(const std::initializer_list<F> &l)
    { 
        std::copy_n(l.begin(), DIM, this->data());
    }

    static int dimension() {return DIM;}

    template<class Vector>
    Point<F, DIM> & operator -= (const Vector & rhs)
    {
        for(int d = 0; d < DIM; d++)
            this->data()[d] -= rhs[d];
        return *this;
    }

    template<class Vector>
    Point<F, DIM> & operator += (const Vector & rhs)
    {
        for(int d = 0; d < DIM; d++)
            this->data()[d] += rhs[d];
        return *this;
    }
};

template<typename F, int DIM >
inline Vector<F, DIM> operator - (const Point<F, DIM> & p,
                         const Point<F, DIM> & q)
{
    Vector<F, DIM> v;
    for(int d = 0; d < DIM; d++)
        v[d] = p[d] - q[d];
    return v;
}

template<typename F, int DIM>
inline Point<F, DIM> operator + (const Point<F, DIM> & p,
                        const Vector<F, DIM> & v)
{
    Point<F, DIM> q;
    for(int d = 0; d < DIM; d++)
        q[d] = p[d] + v[d]; 
    return q;
}

template<typename F, int DIM>
inline Point<F, DIM> operator - (const Point<F, DIM> & p,
                        const Vector<F, DIM> & v)
{
    Point<F, DIM> q;
    for(int d = 0; d < DIM; d++)
        q[d] = p[d] - v[d]; 
    return q;
}

template<typename F, int DIM>
std::ostream& operator << (std::ostream & os, const Point<F, DIM> & p)
{
    if( DIM == 2)
        return os << "Point_2(" << p[0] << ", " <<p[1] <<')';
    else if( DIM == 3)
        return os << "Point_3(" << p[0] << ", " << p[1] << ", " << p[2] << ')';
    else
        return os;
}

} // end of namespace GeometryObject
} // end of namespace WHYSC
#endif // end of Point_h
