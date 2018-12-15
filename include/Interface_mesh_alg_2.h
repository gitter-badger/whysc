#ifndef Interface_mesh_alg_2_h
#define Interface_mesh_alg_2_h

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <map>
#include <cmath>

#include "Geometry/Geometry_kernel.h"
#include "p4est_interface.h"

namespace WHYSC {

namespace MeshAlg {

using namespace P4est;

template<class GK = Geometry_kernel<> >
class Interface_mesh_alg_2
{
public:
    typedef typename GK::Point_2 Point_2;
    typedef typename GK::Point_3 Point_3;


public:
    /** Constructor
     *
     * Input:
     *
     *  
     *
     */
    Interface_mesh_alg_2(MPI_Comm mpi_comm) 
    {
        comm = mpi_comm;
        conn = NULL;
        forest = NULL;
    }

    void create_forest_on_rec_domain(double xmin, double xmax, double ymin, double ymax)
    {
        const Topidx num_vertices = 4;
        const Topidx num_trees = 1;
        const Topidx num_ctt = 0;
        const double vertices[4 * 3] = {
            xmin, ymin, 0,
            xmax, ymin, 0,
            xmin, ymax, 0,
            xmax, ymax, 0};
        const Topidx tree_to_vertex[1 * 4] = { 0, 1, 2, 3};
        const Topidx tree_to_tree[1 * 4] = { 0, 0, 0, 0};
        const int8_t tree_to_face[1 * 4] = { 0, 1, 2, 3};
        conn = connectivity_new_copy(
                num_vertices, num_trees, 0,
                vertices, tree_to_vertex,
                tree_to_tree, tree_to_face,
                NULL, &num_ctt, NULL, NULL);
        forest = forest_new(comm, conn, 0, NULL, NULL);
    }

    void uniform_refine(int n)
    {
        for(int i = 0; i < n; i++)
            refine(forest, 0, uniform_refine_function, NULL);

        int partforcoarsen = 0;
        partition(forest, partforcoarsen, NULL);
    }

    static int uniform_refine_function(Forest * forest, 
            Topidx which_tree, 
            Quadrant * quad)
    {
        return 1;
    }

    static int adaptive_refine_function(Forest * forest,
            Topidx which_tree,
            Quadrant * quad)
    {
        return 1;
    }

    void destroy_forest()
    {
        /* Destroy the p4est and the connectivity structure. */
        if(forest != NULL)
            forest_destroy(forest);
        if(conn != NULL)
            connectivity_destroy(conn);
    }

    void get_corners_of_quadrant(Quadrant quad, std::vector<Point_3> & vec)
    {
    }

    void to_vtk(const char *file)
    {
        vtk_write_file(forest, NULL, file);
    }

private:
    MPI_Comm  comm;
    Connectivity * conn;
    Forest * forest;
};

} // end of namespace MeshAlg

} // end of namespace WHYSC
#endif // end of Interface_mesh_alg_2_h
