#ifndef Point_h
#define Point_h

#include <array>
#include <algorithm>
#include <initializer_list>
#include <assert.h>

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

};

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
