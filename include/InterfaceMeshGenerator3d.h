#ifndef InterfaceMeshGenerator3d_h
#define InterfaceMeshGenerator3d_h

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <array>
#include <map>
#include <cassert>
#include <cmath>
#include <utility>
#include <sstream>
#include <fstream>        


#include "p8est_interface.h"

// VTK header 
#include <vtkUnstructuredGrid.h>
#include <vtkXMLUnstructuredGridWriter.h>
#include <vtkXMLPUnstructuredGridWriter.h>
#include <vtkSmartPointer.h>
#include <vtkPoints.h>
#include <vtkQuad.h>
#include <vtkTriangle.h>
#include <vtkIntArray.h>
#include <vtkCellData.h>


namespace WHYSC {

namespace MeshAlg {

using namespace P8est;

template<typename GK, typename Mesh>
class InterfaceMeshGenerator3d
{
public:
    typedef typename GK::Float   Float;
    typedef typename GK::Int     Int;
    typedef typename GK::Point_3 Point;
    typedef typename GK::Surface Surface;

    struct CellData {
        Int refine; // 

        CellData()
        {
            refine = 0;
        }
    };

private:
    MPI_Comm  comm;
    Connectivity * conn;
    Forest * forest;
    Nodes * nodes;
    Ghost * ghost;
    PMesh * mesh;

    Surface surface;
    Mesh ds;
public:
    InterfaceMeshGenerator3d(MPI_Comm mpi_comm)
    {
        comm = mpi_comm;
        conn = NULL;
        forest = NULL;
        nodes = NULL;
        ghost = NULL;
        mesh = NULL;
    }

    void create_forest_on_cube_domain(
            Float xmin=0, Float xmax=1, 
            Float ymin=0, Float ymax=1,
            Float zmin=0, Float zmax=1
            )
    { 
        const Topidx num_vertices = 8;
        const Topidx num_trees = 1;
        const Topidx num_ett = 0;
        const Topidx num_ctt = 0;
        const Float       vertices[8 * 3] = {
                  xmin, ymin, zmin,
                  xmax, ymin, zmin,
                  xmin, ymax, zmin,
                  xmax, ymax, zmin,
                  xmin, ymin, zmax,
                  xmax, ymin, zmax,
                  xmin, ymax, zmax,
                  xmax, ymax, zmax
        };
        const Topidx tree_to_vertex[1 * 8] = {0, 1, 2, 3, 4, 5, 6, 7};
        const Topidx tree_to_tree[1 * 6] = {0, 0, 0, 0, 0, 0};
        const int8_t tree_to_face[1 * 6] = {0, 1, 2, 3, 4, 5};

        conn = connectivity_new_copy(num_vertices, num_trees, 0, 0,
                vertices, tree_to_vertex,
                tree_to_tree, tree_to_face,
                NULL, &num_ett, NULL, NULL,
                NULL, &num_ctt, NULL, NULL);
        forest = forest_new_ext(comm, /* communicator */
                conn,                 /* cnnectivity */
                0,                    /* minimum cell per MPI process*/
                0,                    /* minimum level of refinement */
                1,                    /* fill uniform */
                sizeof(CellData),     /* data size */
                NULL, 
                NULL);
    }

    void uniform_refine(int n)
    {
        for(int i = 0; i < n; i++)
            refine(forest, 0, uniform_refine_function, NULL);

    }

    static int uniform_refine_function(Forest * forest, 
            Topidx which_tree, 
            Quadrant * quad)
    {

        return 1;
    }

    void adaptive_refine(int n)
    {
        for(auto i = 0; i < n; i++)
        {
            refine(forest, 0, adaptive_refine_function, NULL);
            partition(forest, 0, NULL);
        }
        balance(forest, P4EST_CONNECT_FULL, NULL);
        partition(forest, 1, NULL);
    }

    static Int adaptive_refine_function(Forest * forest,
            Topidx which_tree,
            Quadrant * quad)
    {

        return 1;
    }

    void to_vtk(const char *file)
    {
        vtk_write_file(forest, NULL, file);
    }

    void destroy_forest()
    {
        /* Destroy the p4est and the connectivity structure. */
        if(forest != NULL)
            forest_destroy(forest);
        if(conn != NULL)
            connectivity_destroy(conn);
        if(nodes != NULL)
            nodes_destroy(nodes);
        if(mesh != NULL)
            mesh_destroy(mesh);
        if(ghost != NULL)
            ghost_destroy(ghost);
    }
};

} // end of namespace MeshAlg

} // end of namespace WHYSC
#endif // end of InterfaceMeshGenerator3d_h
