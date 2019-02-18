#ifndef Quad_h
#define Quad_h

#include <array>
#include <utility>
#include "CellType.h"

namespace WHYSC {
namespace Mesh {

struct Quad
{
    typedef std::array<int, 2> Edge;
    typedef std::pair<int, Edge> IEdge;

    static int dim; // the dimension of cell
    static int NV[3];
    static int ND[3];
    static CellType type;

    static int edge[4][2];
    static int node2edge[4][2];

    int        cell[4];

    int & operator[](const int i) 
    {
        return cell[i];
    }

    const int & operator[](const int i) const
    {
        return cell[i];
    }
};

int Quad::dim = 2;
int Quad::NV[3] = {1, 2, 4};
int Quad::ND[3] = {4, 4, 1};
CellType Quad::type = QUAD;
int Quad::edge[4][2] = {
    {0, 1}, {1, 2}, {2, 3}, {3, 0}
};

int Quad::node2edge[4][2] = {
    {3, 0}, {0, 1}, {1, 2}, {2, 3}
};

std::ostream& operator << (std::ostream & os, const Quad & cell)
{

    auto dim = Quad::dim;
    for(auto i = 0; i < Quad::NV[dim]; i++)
    {
        os << cell[i] << " ";
    }
    return os;
}

} // end of namespace Mesh

} // end of namespace WHYSC

#endif // end of Quad_h
