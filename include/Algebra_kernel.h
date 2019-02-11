#ifndef Algebra_kernel_h
#define Algebra_kernel_h

#include "Matrix.h"
namespace WHYSC {

template<typename F=double>
class Algebra_kernel
{
public:
    typedef typename AlgebraObject::Matrix<F, 2, 2> Matrix22;
    typedef typename AlgebraObject::Matrix<F, 3, 3> Matrix33;
    typedef typename AlgebraObject::Matrix<F, 4, 4> Marrix44;

};

} // end of namespace WHYSC
#endif // end of Algebra_kernel_h
