#include "cells.h"

// Tetrahedron

// Hexahedron
int WHYSC::Mesh::Hexahedron::localEdge[12][2] = {
    {0, 1}, {1, 2}, {2, 3}, {0, 3},
    {0, 4}, {1, 5}, {2, 6}, {3, 7},
    {4, 5}, {5, 6}, {6, 7}, {4, 7}
};

int WHYSC::Mesh::Hexahedron::localFace[6][4] = {
    {0, 3, 2, 1}, {4, 5, 6, 7}, // bottom and top faces
    {0, 4, 7, 3}, {1, 2, 6, 5}, // left and right faces
    {0, 1, 5, 4}, {2, 3, 7, 6}  // front and back faces
};

int WHYSC::Mesh::Hexahedron::localFace2edge[6][4] = {
        {3,  2, 1, 0}, {8, 9, 10, 11},
        {4, 11, 7, 3}, {1, 6,  9,  5},
        {0,  5, 8, 4}, {2, 7, 10,  6}
};

int WHYSC::Mesh::Hexahedron::NV = 8;
int WHYSC::Mesh::Hexahedron::NE = 12;
int WHYSC::Mesh::Hexahedron::NF = 6;
