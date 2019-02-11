#ifndef Matrix_h
#define Matrix_h

#include <array>
#include <algorithm>
#include <initializer_list>
#include <assert.h>

namespace WHYSC {
namespace AlgebraObject {

template<typename F, int ROW, int COL>
class Matrix: public std::array<F, ROW*COL>
{
public:
    Matrix(const std::initializer_list<F> &l)
    { 
        std::copy_n(l.begin(), this->size(), this->data());
    }

};

} // end of namespace AlgebraObject

} // end of namespace WHYSC
#endif // end of Matrix_h
