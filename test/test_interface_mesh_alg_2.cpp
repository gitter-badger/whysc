
#include "Interface_mesh_alg_2.h"

#include <iostream>

typedef WHYSC::Geometry_kernel<> GK;
typedef WHYSC::MeshAlg::Interface_mesh_alg_2 Interface_mesh_alg_2;
typedef GK::Level_set_circle Level_set_circle;

int main(int argc, char **argv)
{
    int mpiret = sc_MPI_Init(&argc, &argv);
    SC_CHECK_MPI(mpiret);
    sc_MPI_Comm mpicomm = sc_MPI_COMM_WORLD;

    auto circle = Level_set_circle(0, 0, 0.8);
    auto meshalg = Interface_mesh_alg_2(mpicomm, circle);
    meshalg.create_forest_on_rec_domain(-1, 1, -1, 1);

    meshalg.uniform_refine(1);
    meshalg.adaptive_refine(2);
    meshalg.create_data_structure();
    meshalg.print();
    meshalg.to_vtk("test");
    meshalg.destroy_forest();

    sc_finalize();
    mpiret = sc_MPI_Finalize();
    return 0;
}
