#include <iostream>

#include "Hexahedron.h"
#include "Mesh_3.h"
#include "Geometry_kernel.h"
#include "Interface_mesh_generator_3.h"

typedef WHYSC::Geometry_kernel<double, int> GK;
typedef WHYSC::Mesh::Hexahedron Hex;
typedef WHYSC::Mesh::Mesh_3<GK, Hex> HexMesh;
typedef GK::Sphere_3  Sphere_3;

typedef WHYSC::MeshAlg::Interface_mesh_generator_3<GK, HexMesh, Sphere_3> IMGenerator;

int main(int argc, char ** argv)
{
    std::string arg = argv[1];
    std::size_t pos;
    int un = std::stoi(arg, &pos);
    std::cout << "Uniform refine: " << un << std::endl;

    arg = argv[2];
    int an = std::stoi(arg, &pos);
    std::cout << "Adaptive refine: " << an << std::endl;
    int mpiret = sc_MPI_Init(&argc, &argv);
    SC_CHECK_MPI(mpiret);
    sc_MPI_Comm comm = sc_MPI_COMM_WORLD;

    double start = MPI_Wtime();

    auto alg = IMGenerator(comm);
    auto sphere = Sphere_3(0.0, 0.0, 0.0, 0.8);
    alg.set_interface(sphere);
    alg.create_forest_on_cube_domain(-1, 1, -1, 1, -1, 1);
    alg.uniform_refine(un);
    alg.adaptive_refine(an);

    double end = MPI_Wtime();
    std::cout << "The process took " << end - start << " seconds to run." << std::endl;

    alg.create_data_structure();
    alg.print();
    alg.to_vtk("test");
    alg.destroy_forest();

    sc_finalize();
    mpiret = sc_MPI_Finalize();
    return 0;
}
