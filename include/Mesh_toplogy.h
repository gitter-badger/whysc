#ifndef Mesh_toplogy_h
#define Mesh_toplogy_h


#include <vector>
namespace WHYSC {
namespace Mesh {


template<typename I, int DIM>
struct Mesh_toplogy
{
    struct Toplogy 
    {
        std::vector<I> neighbor;
        std::vecotr<I> neighborLocation; // size is N0 + 1
                                         // neighborLocation[-1] == neighbor.size()

        I N0; // the number of class 0 entity
        I N1; // the number of class 1 entity


        /* get the number of neighbors of i-th class 0 entity*/
        I number_of_neighbors(const I i)
        {
            return neighborLocation[i+1] - neighborLocation[i];
        }

        Int & operator()(const Int i, const Int j) 
        {
            return neighbor[neighborLocation[i]+j];
        }

        const Int & operator()(const Int i, const Int j) const
        {
            return neighbor[neighborLocation[i]+j];
        }

        Toplogy()
        {
            N0 = 0;
            N1 = 0;
        }

        Toplogy(I n0, I n1)
        {
            N0 = n0;
            N1 = n1;
            neighborLocation.resize(N0+1);
        }
    };

    Toplogy  toplogy[DiM][DIM]; 
};

} // end of namespace Mesh

} // end of namespace WHYSC


#endif // end of Mesh_toplogy_h
