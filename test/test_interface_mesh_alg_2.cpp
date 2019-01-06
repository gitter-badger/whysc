
#include "Interface_mesh_alg_2.h"

#include <iostream>
#include <cmath>
#include <string>

typedef WHYSC::MeshAlg::Interface_mesh_alg_2 Interface_mesh_alg_2;
typedef Interface_mesh_alg_2::Point_2 Point_2;
typedef Interface_mesh_alg_2::Curve Curve;

static double circle(const Point_2 & p)
{
    double val = p[0]*p[0] + p[1]*p[1] - 0.8*0.8;
    return val; 
}

static double petal(const Point_2 & p)
{
    int omega = 3;
    double xc = 0;
    double yc = 0;
    double r0 = 0.5;
    double r1 = 0.4;
    double theta = std::atan2(p[1] - yc, p[0] - xc);
    double r = (p[0] - xc)*(p[0] - xc) + (p[1] - yc)*(p[1] - yc);
    double r2 = r0 + r1*std::sin(omega*theta);
    double val = r - r2;
}

Curve Interface_mesh_alg_2::curve=circle;
//Curve Interface_mesh_alg_2::curve=petal;

int main(int argc, char **argv)
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
    sc_MPI_Comm mpicomm = sc_MPI_COMM_WORLD;

    double start = MPI_Wtime();
    auto meshalg = Interface_mesh_alg_2(mpicomm);
    meshalg.create_forest_on_rec_domain(-1, 1, -1, 1);

    meshalg.uniform_refine(un);
    meshalg.adaptive_refine(an);
    meshalg.create_data_structure();
    meshalg.create_interface_mesh();

    double end = MPI_Wtime();
    std::cout << "The process took " << end - start << " seconds to run." << std::endl;

    meshalg.write_to_vtk();
    //meshalg.to_vtk("test");
    meshalg.destroy_forest();

    sc_finalize();
    mpiret = sc_MPI_Finalize();
    return 0;
}
