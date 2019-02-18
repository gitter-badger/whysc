#ifndef Polygon_h
#define Polygon_h

namespace WHYSC {

namespace Mesh {

struct Polygon 
{
    typedef std::array<int, 2> Edge;
    typedef std::pair<int, Edge> IEdge;

    static int dim; // the dimension of cell
    static CellType type;

    int NV[3];
    int ND[3];

    int ** edge;
    int ** node2edge;
    int *  cell;

    Polygon(int nv)
    {
        NV[0] = 1;
        NV[1] = 2;
        NV[2] = nv;

        ND[0] = nv;
        ND[1] = nv;
        ND[2] = 1;

        cell = new int[nv];
        edge = new int*[nv];
        node2edge = new int*[nv]
        int i = 0;
        for(auto i = 0; i < nv; i++)
        {
            edge[i] = new int[2];
            node2edge[i] = new int[2];
            edge[i][0] = i;
            edge[i][1] = (i+1)%nv;
            node2edge[i][0] = i-1;
            node2edge[i][1] = i;
        }
        node2edge[0][0] = nv-1;
    }

    int & operator[](const int i) 
    {
        return cell[i];
    }

    const int & operator[](const int i) const
    {
        return cell[i];
    }

    ~Polygon()
    {
        delete[] cell;
        for(auto i = 0; i < ND[2]; i++)
            delete[] edge[i];
        delete[] edge;

        for(auto i = 0; i < ND[1]; i++)
            delete[] node2edge[i];
        delete[] node2edge;
    }

private:
    Polygon() {};
};

int Polygon::dim = 2;
int Polygon::type = POLYGON;

} // end of namespace Mesh

} // end of namespace WHYSC
#endif // end of Polygon_h
