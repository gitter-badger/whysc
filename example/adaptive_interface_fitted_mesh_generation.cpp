#include <p4est.h>
#include <p4est_vtk.h>
#include <p4est_mesh.h>

#include <cmath>
#include <iostream>

using namespace std;

template <typename T> int sign(T val) {
    return (T(0) < val) - (val < T(0));
}

p4est_connectivity_t *
p4est_connectivity_new_rec(double xmin, double xmax, double ymin, double ymax)
{
  const p4est_topidx_t num_vertices = 4;
  const p4est_topidx_t num_trees = 1;
  const p4est_topidx_t num_ctt = 0;
  const double        vertices[4 * 3] = {
    xmin, ymin, 0,
    xmax, ymin, 0,
    xmin, ymax, 0,
    xmax, ymax, 0,
  };
  const p4est_topidx_t tree_to_vertex[1 * 4] = {
    0, 1, 2, 3,
  };
  const p4est_topidx_t tree_to_tree[1 * 4] = {
    0, 0, 0, 0,
  };
  const int8_t        tree_to_face[1 * 4] = {
    0, 1, 2, 3,
  };

  return p4est_connectivity_new_copy (num_vertices, num_trees, 0,
                                      vertices, tree_to_vertex,
                                      tree_to_tree, tree_to_face,
                                      NULL, &num_ctt, NULL, NULL);
}

static int
uniform_refine_fn(p4est_t * p4est, p4est_topidx_t which_tree,
           p4est_quadrant_t * quadrant)

{
    return 1;
}


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

int main(int argc, char **argv)
{
    int mpiret = sc_MPI_Init (&argc, &argv);
    SC_CHECK_MPI(mpiret);
    sc_MPI_Comm mpicomm = sc_MPI_COMM_WORLD;
    p4est_connectivity_t * conn = p4est_connectivity_new_unitsquare();

    /* Create a forest that is not refined; it consists of the root octant. */
    p4est_t* p4est = p4est_new(mpicomm, conn, 0, NULL, NULL);

    for(int i = 0; i < 2; i++)
        p4est_refine(p4est, 0, uniform_refine_fn, NULL);

    //for(int i = 0; i < 10; i++)
    //   p4est_refine(p4est, 0, refine_fn, NULL);

    /* Partition: The quadrants are redistributed for equal element count.  The
     * partition can optionally be modified such that a family of octants, which
     * are possibly ready for coarsening, are never split between processors. */
    int partforcoarsen = 0;
    p4est_partition(p4est, partforcoarsen, NULL);

    /* Write the forest to disk for visualization, one file per processor. */
    p4est_vtk_write_file(p4est, NULL, "test");

    /* Destroy the p4est and the connectivity structure. */
    p4est_destroy(p4est);
    p4est_connectivity_destroy(conn);

    /* Verify that allocations internal to p4est and sc do not leak memory.
     * This should be called if sc_init () has been called earlier. */
    sc_finalize();

    /* This is standard MPI programs.  Without --enable-mpi, this is a dummy. */
    mpiret = sc_MPI_Finalize();
    return 0;
}

