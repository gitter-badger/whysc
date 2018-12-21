#ifndef p4est_interface_h
#define p4est_interface_h

#include <p4est.h>
#include <p4est_balance.h>
#include <p4est_connectivity.h>
#include <p4est_ghost.h>
#include <p4est_nodes.h>
#include <p4est_vtk.h>
#include <p4est_mesh.h>

namespace P4est {

typedef p4est_connectivity_t Connectivity;
typedef p4est_t              Forest;
typedef p4est_tree_t         Tree;
typedef p4est_quadrant_t     Quadrant;
typedef p4est_topidx_t       Topidx;
typedef p4est_locidx_t       Locidx;
typedef p4est_qcoord_t       Qcoord;
typedef p4est_nodes_t        Nodes;
typedef p4est_mesh_t         Mesh;
typedef p4est_ghost_t        Ghost;
typedef p4est_mesh_face_neighbor_t  MFNeighbor;
typedef sc_MPI_Comm          MPI_Comm;

auto& connectivity_new_copy = p4est_connectivity_new_copy;
auto& connectivity_destroy = p4est_connectivity_destroy;

auto& forest_new = p4est_new;
auto& forest_destroy = p4est_destroy;

auto& nodes_new = p4est_nodes_new;
auto& nodes_destroy = p4est_nodes_destroy;

auto& mesh_new = p4est_mesh_new;
auto& mesh_destroy = p4est_mesh_destroy;

auto& ghost_new = p4est_ghost_new;
auto& ghost_destroy = p4est_ghost_destroy;

auto& tree_array_index = p4est_tree_array_index;
auto& quadrant_array_index = p4est_quadrant_array_index;

auto& refine = p4est_refine;
auto& partition = p4est_partition;
auto& vtk_write_file =  p4est_vtk_write_file;
auto& qcoord_to_vertex = p4est_qcoord_to_vertex;
auto& mesh_face_neighbor_init = p4est_mesh_face_neighbor_init;
auto& mesh_face_neighbor_init2 = p4est_mesh_face_neighbor_init2;
} // end of namespace P4est

#endif // end of p4est_interface_h
