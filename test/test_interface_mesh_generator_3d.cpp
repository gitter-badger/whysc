#include <iostream>

#include "Hexahedron.h"
#include "Mesh3d.h"
#include "GeometryKernel.h"
#include "InterfaceMeshGenerator3d.h"

typedef WHYSC::GeometryKernel<double, int> GK;
typedef WHYSC::Mesh::Hexahedron Hex;
typedef WHYSC::Mesh::Mesh3d<GK, Hex> HexMesh;
typedef WHYSC::MeshAlg::InterfaceMeshGenerator3d<GK, HexMesh> IMGenerator;

int main(int argc, char ** argv)
{
    int mpiret = sc_MPI_Init(&argc, &argv);
    SC_CHECK_MPI(mpiret);
    sc_MPI_Comm comm = sc_MPI_COMM_WORLD;

    double start = MPI_Wtime();
    auto alg = IMGenerator(comm);
    alg.create_forest_on_cube_domain(-1, 1, -1, 1, -1, 1);
    alg.uniform_refine(1);
    alg.adaptive_refine(1);
    double end = MPI_Wtime();

    std::cout << "The process took " << end - start << " seconds to run." << std::endl;

    alg.to_vtk("test");
    alg.destroy_forest();

    sc_finalize();
    mpiret = sc_MPI_Finalize();
    return 0;
}
