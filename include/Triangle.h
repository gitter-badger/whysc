#ifndef Triangle_h
#define Triangle_h

namespace WHYSC {

namespace Mesh {

struct Triangle
{
    typedef std::array<int, 2> Edge;
    typedef std::pair<int, Edge> IEdge;

    static int dim; // the dimension of cell
    static int NV[3];
    static int ND[3];
    static CellType type;

    static int edge[3][2];
    static int node2edge[3][2];

    int        cell[3];

    int & operator[](const int i) 
    {
        return cell[i];
    }

    const int & operator[](const int i) const
    {
        return cell[i];
    }
};

int Triangle::dim = 2;
int Triangle::NV[3] = {1, 2, 3};
int Triangle::ND[3] = {3, 2, 1};
CellType Triangle::type = Triangle;
int Triangle::edge[3][2] = {
    {1, 2}, {2, 0}, {1, 0}
};

int Triangle::node2edge[3][2] = {
    {1, 2}, {2, 0}, {0, 1}
};

std::ostream& operator << (std::ostream & os, const Triangle & cell)
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
#endif // end of Triangle_h
