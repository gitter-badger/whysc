#ifndef p8est_interface_h
#define p8est_interface_h

#include <p8est.h>
#include <p8est_balance.h>
#include <p8est_connectivity.h>
#include <p8est_ghost.h>
#include <p8est_nodes.h>
#include <p8est_vtk.h>
#include <p8est_mesh.h>

namespace P8est {

typedef p8est_connectivity_t Connectivity;
typedef p8est_t              Forest;
typedef p8est_tree_t         Tree;
typedef p8est_quadrant_t     Quadrant;
typedef p8est_topidx_t       Topidx;
typedef p8est_locidx_t       Locidx;
typedef p8est_qcoord_t       Qcoord;
typedef p8est_nodes_t        Nodes;
typedef p8est_mesh_t         Mesh;
typedef p8est_ghost_t        Ghost;
typedef p8est_mesh_face_neighbor_t  MFNeighbor;
typedef sc_MPI_Comm          MPI_Comm;

auto& connectivity_new_copy = p8est_connectivity_new_copy;
auto& connectivity_destroy = p8est_connectivity_destroy;

auto& forest_new = p8est_new;
auto& forest_destroy = p8est_destroy;

auto& nodes_new = p8est_nodes_new;
auto& nodes_destroy = p8est_nodes_destroy;

auto& mesh_new = p8est_mesh_new;
auto& mesh_destroy = p8est_mesh_destroy;

auto& ghost_new = p8est_ghost_new;
auto& ghost_destroy = p8est_ghost_destroy;

auto& tree_array_index = p8est_tree_array_index;
auto& quadrant_array_index = p8est_quadrant_array_index;

auto& refine = p8est_refine;
auto& partition = p8est_partition;
auto& balance = p8est_balance;
auto& vtk_write_file =  p8est_vtk_write_file;
auto& qcoord_to_vertex = p8est_qcoord_to_vertex;
auto& mesh_face_neighbor_init = p8est_mesh_face_neighbor_init;
auto& mesh_face_neighbor_init2 = p8est_mesh_face_neighbor_init2;
} // end of namespace P4est
#endif // end of p8est_interface_h
