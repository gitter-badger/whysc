#ifndef Polyhedron_h
#define Polyhedron_h

namespace WHYSC {

namespace Mesh {

struct Polygon_face
{
    typedef std::array<int, 2> Edge;
    typedef std::pair<int, Edge> IEdge;

    static int dim; // the dimension of cell
    static CellType type;

    int NV[3];
    int ND[3];

    int ** edge;
    int ** node2edge;
    int *  face;
    int face2cell[2]; // face2cell[0]: the left cell idx of current face 
                      // face2cell[1]: the right cell idx of current face

    Polygon_face(int nv, int lf, int rf)
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
        face2cell[0] = lf;
        face2cell[1] = rf;
    }

    int & operator[](const int i) 
    {
        return face[i];
    }

    const int & operator[](const int i) const
    {
        return face[i];
    }

    ~Polygon_face()
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
    Polygon_face() {};
}

int Polygon_face::dim = 2;
int Polygon_face::type = POLYGON_FACE;

struct Polyhedron
{
    typedef std::array<int, 2> Edge;
    typedef Polygon_face Face; 

    static int dim;
    static CellType type;
    int NV[4]; // NV[3]
    int ND[4];
    std::vector<Face *> face;
    std::vector<bool> direction;

    Polyhedron()
    {
    }

    
};

int Polyhedron::dim = 3;
CellType Polyhedron::type = POLYHEDRON;

} // end of namespace Mesh

} // end of namespace WHYSC
#endif // end of Polyhedron_h
