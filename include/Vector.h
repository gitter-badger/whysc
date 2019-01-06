#ifndef Vector_h
#define Vector_h

namespace WHYSC {
namespace GeometryObject {

template<typename F, int DIM>
class Vector : public std::array<F, DIM>
{
public:
    typedef typename std::array<F, DIM> Base;
public:

    Vector()
    {
        std::fill_n(this->data(), DIM, 0.0);
    }

    Vector(const std::initializer_list<F> &l)
    { 
        std::copy_n(l.begin(), DIM, this->data());
    }

    static int dimension() {return DIM;}

};

template<typename F, int DIM>
std::ostream& operator << (std::ostream & os, const Vector<F, DIM> & p)
{
    if( DIM == 2)
        return os << "Vector_2(" << p[0] << ", " <<p[1] <<')';
    else if( DIM == 3)
        return os << "Vector_3(" << p[0] << ", " << p[1] << ", " << p[2] << ')';
    else
        return os;
}

} // end of namespace GeometryObject
} // end of namespace WHYSC
#endif // end of Vector_h
