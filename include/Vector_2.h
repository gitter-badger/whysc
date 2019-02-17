#ifndef Vector_2_h
#define Vector_2_h

#include <array>
#include <algorithm>
#include <initializer_list>
#include <assert.h>

namespace WHYSC {
namespace GeometryObject {

template<typename F>
class Vector_2 : public std::array<F, 2>
{
public:
    typedef typename std::array<F, 2> Base;

public:

    Vector_2()
    {
        std::fill_n(this->data(), 2, 0.0);
    }

    Vector_2(const std::initializer_list<F> &l)
    { 
        std::copy_n(l.begin(), 2, this->data());
    }

    Vector_2(const Vector_2 &v)
    { 
        std::copy_n(v.begin(), 2, this->data());
    }

    Vector_2(F vx, F vy)
    {
        this->data()[0] = vx;
        this->data()[1] = vy;
    }

    static int dimension() {return 2;}

    double squared_length()
    {
        double sum = 0.0;
        for(int d = 0; d < 2; d++)
            sum += this->data()[d]*this->data()[d];
        return sum;
    }

    template<class RVector_2>
    double operator * (const RVector_2 & w)
    {
        // dot product of vectors
        double sum = 0.0;
        for(int d = 0; d < 2; d++)
            sum += this->data()[d]*w[d];
        return sum;
    }

    Vector_2<F> & operator *= (const F & s)
    {
        for(int d = 0; d < 2; d++)
            this->data()[d] *= s;
        return * this;
    }


    Vector_2<F> & operator /= (const F & s)
    {
        for(int d = 0; d < 2; d++)
            this->data()[d] /= s;
        return * this;
    }

    template<class RVector_2>
    Vector_2<F> & operator += (const RVector_2 & w)
    {
        for(int d = 0; d < 2; d++)
            this->data()[d] += w[d];
        return * this;
    }


    template<class RVector_2>
    Vector_2<F> & operator -= (const RVector_2 & w)
    {
        for(int d = 0; d < 2; d++)
            this->data()[d] -= w[d];
        return * this;
    }

};

template<typename F>
std::ostream& operator << (std::ostream & os, const Vector_2<F> & v)
{
        return os << "Vector_2(" << v[0] << ", " << v[1] <<')';
}

} // end of namespace GeometryObject
} // end of namespace WHYSC
#endif // end of Vector_2_h
