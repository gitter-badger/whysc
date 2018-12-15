#ifndef Interface_mesh_alg_2_h
#define Interface_mesh_alg_2_h

#include "Geometry/Geometry_kernel.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <map>
#include <cmath>

#include <p4est.h>
#include <p4est_vtk.h>
#include <p4est_mesh.h>

namespace WHYSC {
namespace MeshAlg {

static int
refine_fn (p4est_t * p4est, p4est_topidx_t which_tree,
           p4est_quadrant_t * quadrant)
{

    int sum  = 0;
    double phi = 0.0;

    int8_t l = quadrant->level; // 叶子单元的层数

    if( l > 8) // 如果超过 8 层则停止加密
        return 0;

    p4est_qcoord_t qx = quadrant->x; // qx, qy 是逻辑空间中的编号
    p4est_qcoord_t qy = quadrant->y; 
    double vxy[3] = {0.0, 0.0, 0.0}; // 网格所在的空间为 [0, 1]^2, 从逻辑空间实际空间坐标的转换
    p4est_qcoord_to_vertex(p4est->connectivity, which_tree, qx, qy, vxy);

    // 从 [0, 1]^2 空间到[-1, 1] 空间的转换
    double h = 1.0/pow(2, l); // [0, 1]^2 空间中 第 l 层单元的尺寸
    x = xmin + vxy[0]*(xmax - xmin);
    y = ymin + vxy[1]*(ymax - ymin);
    phi = sqrt(x*x + y*y) - 0.8;
    sum += sign(phi);

    x = xmin + (vxy[0]+h)*(xmax - xmin);
    y = ymin + vxy[1]*(ymax - ymin);
    phi = sqrt(x*x + y*y) - 0.8;
    sum += sign(phi);

    x = xmin + vxy[0]*(xmax - xmin);
    y = ymin + (vxy[1]+h)*(ymax - ymin);
    phi = sqrt(x*x + y*y) - 0.8;
    sum += sign(phi);

    x = xmin + (vxy[0]+h)*(xmax - xmin);
    y = ymin + (vxy[1]+h)*(ymax - ymin);
    phi = sqrt(x*x + y*y) - 0.8;
    sum += sign(phi);

    if(abs(sum) < 4)
        return 1;
    else 
        return 0;
}


template<class GK = Geometry_kernel<> >
class Interface_mesh_alg_2
{
public:
    typedef typename GK::Point_2 Point_2;
    typedef typename GK::Point_3 Point_3;
    typedef p4est_quadrant_t*  Quadrant;
public:
    /** Constructor
     *
     * Input:
     *
     *  
     *
     */
    Interface_mesh_alg_2(int argc, char **argv)
    {
        mpiret = sc_MPI_Init(&argc, &argv);
        SC_CHECK_MPI(mpiret);
        mpicomm = sc_MPI_COMM_WORLD;
    }

    void create_p4est_on_rec_domain(double xmin, double xmax, double ymin, double ymax)
    {
        const p4est_topidx_t num_vertices = 4;
        const p4est_topidx_t num_trees = 1;
        const p4est_topidx_t num_ctt = 0;
        const double        vertices[4 * 3] = {
            xmin, ymin, 0,
            xmax, ymin, 0,
            xmin, ymax, 0,
            xmax, ymax, 0};
        const p4est_topidx_t tree_to_vertex[1 * 4] = { 0, 1, 2, 3};
        const p4est_topidx_t tree_to_tree[1 * 4] = { 0, 0, 0, 0};
        const int8_t tree_to_face[1 * 4] = { 0, 1, 2, 3};
        conn = p4est_connectivity_new_copy(
                num_vertices, num_trees, 0,
                vertices, tree_to_vertex,
                tree_to_tree, tree_to_face,
                NULL, &num_ctt, NULL, NULL);
        p4est = p4est_new(mpicomm, conn, 0, NULL, NULL);
    }

    void get_corners_of_quadrant(p4est_quadrant_t * quadrant, std::vector<Point_3> & vec)
    {
    }

private:
    int mpiret;
    sc_MPI_Comm mpicomm;
    p4est_connectivity_t * conn;
    p4est_t* p4est;

};

} // end of namespace MeshAlg

} // end of namespace WHYSC
#endif // end of Interface_mesh_alg_2_h
