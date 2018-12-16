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


#define DUMP(a) \
    do { std::cout << #a " is value " << (a) << std::endl; } while(false)

namespace WHYSC {

namespace MeshAlg {

using namespace P4est;

template <typename T> int sign(T val) 
{
    return (T(0) < val) - (val < T(0));
}

template<class GK = Geometry_kernel<> >
class Interface_mesh_alg_2
{
public:
    typedef typename GK::Point_2 Point_2;
    typedef typename GK::Point_3 Point_3;
    typedef typename GK::Delaunay_algorithm_2 Delaunay_algorithm_2;
    typedef typename GK::Bisection_algorithm Bisection_algorithm;
    typedef typename GK::Level_set_function Level_set_function;

public:
    /** Constructor
     *
     * Input:
     *
     *  
     *
     */
    Interface_mesh_alg_2(MPI_Comm mpi_comm, Level_set_function & lfun) 
    {
        comm = mpi_comm;
        conn = NULL;
        forest = NULL;
        nodes = NULL;
        ghost = NULL;
        mesh = NULL;
        fun = lfun;
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


    void adaptive_refine(int n)
    {
        for(int i = 0; i < n; i++)
            refine(forest, 0, adaptive_refine_function, NULL);

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
        int sum  = 0;
        double phi = 0.0;

        auto l = quad->level; // 叶子单元的层数
        auto qx = quad->x; // qx, qy 是逻辑空间中的编号
        auto qy = quad->y; 
        double p[3] = {0.0, 0.0, 0.0}; // 网格所在的空间为 [0, 1]^2, 从逻辑空间实际空间坐标的转换
        auto h = P4EST_QUADRANT_LEN(quad->level); 

        qcoord_to_vertex(forest->connectivity, which_tree, qx, qy, p);
        phi = sqrt(p[0]*p[0] + p[1]*p[1]) - 0.8;
        sum += sign(phi);

        qcoord_to_vertex(forest->connectivity, which_tree, qx+h,   qy, p);
        phi = sqrt(p[0]*p[0] + p[1]*p[1]) - 0.8;
        sum += sign(phi);

        qcoord_to_vertex(forest->connectivity, which_tree, qx,   qy+h, p);
        phi = sqrt(p[0]*p[0] + p[1]*p[1]) - 0.8;
        sum += sign(phi);

        qcoord_to_vertex(forest->connectivity, which_tree, qx+h,   qy+h, p);
        phi = sqrt(p[0]*p[0] + p[1]*p[1]) - 0.8;
        sum += sign(phi);

        if(abs(sum) < 4)
            return 1;
        else 
            return 0;
    }


    void create_ghost()
    {
        ghost = ghost_new(forest, P4EST_CONNECT_FACE);
    }


    void create_mesh()
    {
        mesh = mesh_new(forest, ghost, P4EST_CONNECT_FACE);
        auto q2q = mesh->quad_to_quad;
        auto NC = mesh->local_num_quadrants;
        for(auto i=0; i < NC; i++)
        {
            std::cout << i << ": "
                << q2q[4*i + 0] << ", "
                << q2q[4*i + 1] << ", "
                << q2q[4*i + 2] << ", "
                << q2q[4*i + 3] << std::endl;
        }
    }


    void create_nodes()
    {
        nodes = nodes_new(forest, NULL);

        auto NN = nodes->num_owned_indeps; // NN independent nodes 
        auto NC = forest->local_num_quadrants;

        points.resize(3*NN);

        Locidx * n2c = P4EST_ALLOC(Locidx, NN);
        memset(n2c, -1, NN*sizeof(Locidx));

        auto k = 0;
        for(auto i = 0; i < NC; i++)
        {
            std::cout << i << ": ";
            for(auto j = 0; j < 4; k++, j++)
            {
                auto id = nodes->local_nodes[k];
                std::cout << id << " ";
                if(n2c[id] < 0)
                    n2c[id] = k; // the first appear location in local_nodes!
            }
            std::cout << std::endl;
        }

        // traverse all local trees, one process, one tree ?
        auto tree = tree_array_index(forest->trees, forest->first_local_tree);
        auto quads = &tree->quadrants;

        for(auto i = 0; i < NN; i++)
        {
            auto gidx = n2c[i]/4;
            auto lidx = n2c[i] - 4*gidx;
            auto quad = quadrant_array_index(quads, gidx);

            auto qx = quad->x; // qx, qy 是逻辑空间中的编号
            auto qy = quad->y; 
            auto h = P4EST_QUADRANT_LEN(quad->level); 

            auto i0 = lidx/2;
            auto i1 = lidx - i0*2;
            qcoord_to_vertex(forest->connectivity, forest->first_local_tree, qx+i1*h,   qy+i1*h, points.data()+3*i);

            std::cout << i << ": " 
                << points[3*i + 0] << ", " 
                << points[3*i + 1] << ", " 
                << points[3*i + 2] << ", " 
                << std::endl;
        }

        P4EST_FREE(n2c);

        std::cout << "num_owned_indeps: " << nodes->num_owned_indeps << std::endl;
        std::cout << "num_owned_shared: " << nodes->num_owned_shared << std::endl;
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

    void to_vtk(const char *file)
    {
        vtk_write_file(forest, NULL, file);
    }


    void print()
    {
        std::cout<< "number of vertices: " << conn->num_vertices << std::endl;
        int i = 0;
        for(int i = 0; i < conn->num_vertices; i++)
        {
            std::cout << i << ": " 
                << conn->vertices[3*i+0] << ", "
                << conn->vertices[3*i+1] << ", " 
                << conn->vertices[3*i+2] << std::endl;
        }

        std::cout << "number of corners: " << conn->num_corners << std::endl;

        std::cout << "Information of forest:" << std::endl;
        std::cout << forest->local_num_quadrants << std::endl;
    }

private:
    MPI_Comm  comm;
    Connectivity * conn;
    Forest * forest;
    Nodes * nodes;
    Ghost * ghost;
    Mesh * mesh;
    Level_set_function fun;
    std::vector<double> points;
};

} // end of namespace MeshAlg

} // end of namespace WHYSC
#endif // end of Interface_mesh_alg_2_h
