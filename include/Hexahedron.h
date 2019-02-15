#ifndef cells_h
#define cells_h

#include <array>
#include <utility>
#include "CellType.h"

namespace WHYSC {
namespace Mesh {


struct Hexahedron
{
    typedef std::array<int, 2> Edge;
    typedef std::array<int, 4> Face;
    typedef std::pair<int, Face> IFace;
    typedef std::array<int, 6> Cell2cell;
    typedef std::array<int, 6> Cell2face;
    typedef std::array<int, 12> Cell2edge;
    typedef std::array<int, 4> Face2cell;

    static int dim; // the dimension of cell
    static int NV[4];
    static int ND[4];
    static CellType type;

    static int edge[12][2];
    static int face[6][4];
    int        cell[8];


    static int face2edge[6][4];
    static int edge2face[12][2];
    static int node2edge[8][3];
    static int node2face[8][3];

    int & operator[](const int i) 
    {
        return cell[i];
    }

    const int & operator[](const int i) const
    {
        return cell[i];
    }
};

int Hexahedron::dim = 3;
int Hexahedron::NV[4] = {1, 2, 4, 8}; // 0d entity with 1 vertex; 1d entity with 2 vertex; 2d entity with 4 vertices; 3d entity with 8 vertices
int Hexahedron::ND[4] = {8, 12, 6, 1};// 8 vertices, 12 edges, 6 faces, 1 cell
CellType Hexahedron::type = HEXAHEDRON;

int Hexahedron::edge[12][2] = {
    {0, 1}, {2, 3}, {0, 2}, {1, 3},
    {0, 4}, {1, 5}, {2, 6}, {3, 7},
    {4, 5}, {6, 7}, {4, 6}, {5, 7}
};

int Hexahedron::face[6][4] = {
    {4, 6, 0, 2}, {1, 3, 5, 7},  // left and right faces
    {0, 1, 4, 5}, {6, 7, 2, 3},  // front and back faces
    {4, 5, 6, 7}, {2, 3, 0, 1},  // bottom and top faces
};

int Hexahedron::face2edge[6][4] = {
    {10, 2, 4,  6}, {3, 11, 5, 7},
    { 0, 8, 4,  5}, {9,  1, 6, 7},
    { 8, 9,10, 11}, {1,  0, 2, 3}, 
};

int Hexahedron::edge2face[12][2] = {
    {2, 5}, {5, 3}, {0, 5}, {1, 5},
    {0, 2}, {2, 1}, {3, 0}, {1, 3},
    {4, 2}, {3, 4}, {0, 4}, {4, 1}
};

std::ostream& operator << (std::ostream & os, const Hexahedron & cell)
{

    auto dim = Hexahedron::dim;
    for(auto i = 0; i < Hexahedron::NV[dim]; i++)
    {
        os << cell[i] << " ";
    }
    return os;
}


} // end of namespace Mesh

} // end of namespace WHYSC
#endif // end of cells_h