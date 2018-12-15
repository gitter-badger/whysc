
#include "Interface_mesh_alg_2.h"

#include <iostream>

typedef WHYSC::MeshAlg::Interface_mesh_alg_2<> Interface_mesh_alg_2;

int main(int argc, char **argv)
{
    int mpiret = sc_MPI_Init(&argc, &argv);
    SC_CHECK_MPI(mpiret);
    sc_MPI_Comm mpicomm = sc_MPI_COMM_WORLD;


    auto meshalg = Interface_mesh_alg_2(mpicomm);
    meshalg.create_forest_on_rec_domain(-1, 1, -1, 1);
    meshalg.uniform_refine(2);
    meshalg.to_vtk("test");
    meshalg.destroy_forest();

    sc_finalize();
    mpiret = sc_MPI_Finalize();
    return 0;
}
