#ifndef p4est_interface_h
#define p4est_interface_h

#include <p4est.h>
#include <p4est_balance.h>
#include <p4est_connectivity.h>
#include <p4est_ghost.h>
#include <p4est_vtk.h>
#include <p4est_mesh.h>

namespace P4est {

typedef p4est_connectivity_t Connectivity;
typedef p4est_t              Forest;
typedef p4est_tree_t         Tree;
typedef p4est_quadrant_t     Quadrant;
typedef p4est_topidx_t       Topidx;
typedef p4est_locidx_t       Locidx;
typedef p4est_ghost_t        Ghost;
typedef sc_MPI_Comm          MPI_Comm;

auto& connectivity_new_copy = p4est_connectivity_new_copy;
auto& connectivity_destroy = p4est_connectivity_destroy;
auto& forest_new = p4est_new;
auto& forest_destroy = p4est_destroy;

auto& refine = p4est_refine;
auto& partition = p4est_partition;
auto& vtk_write_file =  p4est_vtk_write_file;

} // end of namespace P4est

#endif // end of p4est_interface_h
