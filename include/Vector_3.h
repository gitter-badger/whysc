#ifndef Vector_3_h
#define Vector_3_h

#include <array>
#include <algorithm>
#include <initializer_list>
#include <assert.h>

namespace WHYSC {
namespace GeometryObject {

template<typename F>
class Vector_3 : public std::array<F, 3>
{
public:
    typedef typename std::array<F, 3> Base;

public:

    Vector_3()
    {
        std::fill_n(this->data(), 3, 0.0);
    }

    Vector_3(const std::initializer_list<F> &l)
    { 
        std::copy_n(l.begin(), 3, this->data());
    }

    Vector_3(const Vector_3 &v)
    { 
        std::copy_n(v.begin(), 3, this->data());
    }

    Vector_3(F vx, F vy, F vz)
    {
        this->data()[0] = vx;
        this->data()[1] = vy;
        this->data()[2] = vz;
    }

    static int dimension() {return 3;}

    double squared_length()
    {
        double sum = 0.0;
        for(int d = 0; d < 3; d++)
            sum += this->data()[d]*this->data()[d];
        return sum;
    }

    template<class RVector_3>
    double operator * (const RVector_3 & w)
    {
        // dot product of vectors
        double sum = 0.0;
        for(int d = 0; d < 3; d++)
            sum += this->data()[d]*w[d];
        return sum;
    }

    Vector_3<F> & operator *= (const F & s)
    {
        for(int d = 0; d < 3; d++)
            this->data()[d] *= s;
        return * this;
    }


    Vector_3<F> & operator /= (const F & s)
    {
        for(int d = 0; d < 3; d++)
            this->data()[d] /= s;
        return * this;
    }

    template<class RVector_3>
    Vector_3<F> & operator += (const RVector_3 & w)
    {
        for(int d = 0; d < 3; d++)
            this->data()[d] += w[d];
        return * this;
    }


    template<class RVector_3>
    Vector_3<F> & operator -= (const RVector_3 & w)
    {
        for(int d = 0; d < 3; d++)
            this->data()[d] -= w[d];
        return * this;
    }

};

template<typename F>
std::ostream& operator << (std::ostream & os, const Vector_3<F> & v)
{
        return os << "Vector_3(" << v[0] << ", " << v[1] <<')';
}

} // end of namespace GeometryObject
} // end of namespace WHYSC
#endif // end of Vector_3_h
