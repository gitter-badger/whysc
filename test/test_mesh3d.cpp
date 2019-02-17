#include <iostream>
#include <array> 

#include "Hexahedron.h"
#include "Mesh_3.h"
#include "Geometry_kernel.h"

typedef WHYSC::Geometry_kernel<double, int> GK;
typedef WHYSC::Mesh::Hexahedron Hex;
typedef WHYSC::Mesh::Mesh_3<GK, Hex> Mesh;

int main(int argc, char **argv)
{

    double node[] = {
        0.0, 0.0, 0.0, 
        0.0, 0.0, 0.5,
        0.0, 0.0, 1.0,
        0.0, 0.5, 0.0,
        0.0, 0.5, 0.5,
        0.0, 0.5, 1.0,
        0.0, 1.0, 0.0,
        0.0, 1.0, 0.5,
        0.0, 1.0, 1.0,
        0.5, 0.0, 0.0,
        0.5, 0.0, 0.5,
        0.5, 0.0, 1.0,
        0.5, 0.5, 0.0,
        0.5, 0.5, 0.5,
        0.5, 0.5, 1.0,
        0.5, 1.0, 0.0,
        0.5, 1.0, 0.5,
        0.5, 1.0, 1.0,
        1.0, 0.0, 0.0,
        1.0, 0.0, 0.5,
        1.0, 0.0, 1.0,
        1.0, 0.5, 0.0,
        1.0, 0.5, 0.5,
        1.0, 0.5, 1.0,
        1.0, 1.0, 0.0,
        1.0, 1.0, 0.5,
        1.0, 1.0, 1.0
    };

    int cell[] = {
        0,   9,  12,   3,   1,  10,  13,   4,
        1,  10,  13,   4,   2,  11,  14,   5,
        3,  12,  15,   6,   4,  13,  16,   7,
        4,  13,  16,   7,   5,  14,  17,   8,
        9,  18,  21,  12,  10,  19,  22,  13,
       10,  19,  22,  13,  11,  20,  23,  14,
       12,  21,  24,  15,  13,  22,  25,  16,
       13,  22,  25,  16,  14,  23,  26,  17
    };

    int NN = 27;
    int NC = 8;

    auto mesh = Mesh();
    mesh.init(NN, NC, node, cell);
    mesh.print();
    return 0;
}

