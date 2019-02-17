#ifndef Interface_mesh_generator_3_h
#define Interface_mesh_generator_3_h

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <array>
#include <map>
#include <cassert>
#include <cmath>
#include <utility>
#include <sstream>
#include <fstream>        


#include "p8est_interface.h"

// VTK header 
#include <vtkUnstructuredGrid.h>
#include <vtkXMLUnstructuredGridWriter.h>
#include <vtkXMLPUnstructuredGridWriter.h>
#include <vtkSmartPointer.h>
#include <vtkPoints.h>
#include <vtkQuad.h>
#include <vtkHexahedron.h>
#include <vtkTriangle.h>
#include <vtkIntArray.h>
#include <vtkCellData.h>


namespace WHYSC {

namespace MeshAlg {

using namespace P8est;

template<typename GK, typename Mesh, typename Interface>
class Interface_mesh_generator_3
{
public:
    typedef typename GK::Float   Float;
    typedef typename GK::Int     Int;
    typedef typename GK::Point_3 Point;

    struct CellData {
        Int refine; // 
        CellData()
        {
            refine = 0;
        }
    };

private:
    MPI_Comm  comm;
    Connectivity * conn;
    Forest * forest;
    Nodes * nodes;
    Ghost * ghost;
    PMesh * mesh;
    Interface interface;
    Mesh ds;
public:
    Interface_mesh_generator_3(MPI_Comm mpi_comm)
    {
        comm = mpi_comm;
        conn = NULL;
        forest = NULL;
        nodes = NULL;
        ghost = NULL;
        mesh = NULL;
    }

    void set_interface(Interface inte)
    {
        interface = inte;
    }

    void create_forest_on_cube_domain(
            Float xmin=0, Float xmax=1, 
            Float ymin=0, Float ymax=1,
            Float zmin=0, Float zmax=1
            )
    { 
        const Topidx num_vertices = 8;
        const Topidx num_trees = 1;
        const Topidx num_ett = 0;
        const Topidx num_ctt = 0;
        const Float       vertices[8 * 3] = {
                  xmin, ymin, zmin,
                  xmax, ymin, zmin,
                  xmin, ymax, zmin,
                  xmax, ymax, zmin,
                  xmin, ymin, zmax,
                  xmax, ymin, zmax,
                  xmin, ymax, zmax,
                  xmax, ymax, zmax
        };
        const Topidx tree_to_vertex[1 * 8] = {0, 1, 2, 3, 4, 5, 6, 7};
        const Topidx tree_to_tree[1 * 6] = {0, 0, 0, 0, 0, 0};
        const int8_t tree_to_face[1 * 6] = {0, 1, 2, 3, 4, 5};

        conn = connectivity_new_copy(num_vertices, num_trees, 0, 0,
                vertices, tree_to_vertex,
                tree_to_tree, tree_to_face,
                NULL, &num_ett, NULL, NULL,
                NULL, &num_ctt, NULL, NULL);
        forest = forest_new_ext(comm, /* communicator */
                conn,                 /* cnnectivity */
                0,                    /* minimum cell per MPI process*/
                0,                    /* minimum level of refinement */
                1,                    /* fill uniform */
                sizeof(CellData),     /* data size */
                NULL, 
                NULL);
    }

    void uniform_refine(int n)
    {
        for(int i = 0; i < n; i++)
            refine(forest, 0, uniform_refine_function, NULL);

    }

    static int uniform_refine_function(Forest * forest, 
            Topidx which_tree, 
            Quadrant * quad)
    {
        return 1;
    }

    void adaptive_refine(int n)
    {
        for(auto i = 0; i < n; i++)
        {
            mark_cells();
            refine(forest, 0, adaptive_refine_function, NULL);
            partition(forest, 0, NULL);
        }
        balance(forest, P8EST_CONNECT_FULL, NULL);
        partition(forest, 1, NULL);
    }

    static Int adaptive_refine_function(Forest * forest,
            Topidx which_tree,
            Quadrant * quad)
    {

        auto data = (CellData *)quad->p.user_data;
        return data->refine;
    }

    void mark_cells()
    {
        std::cout << "Entering mark cell!......" << std::endl;
        std::cout << "The number of tree is : " << forest->last_local_tree - forest->first_local_tree + 1 << std::endl;

        for(auto i = forest->first_local_tree; i <= forest->last_local_tree; i++)
        {
            auto tree = tree_array_index(forest->trees, i);
            auto cells = &tree->quadrants;
            auto nc = cells->elem_count;
            std::cout << "The number of quad cells is :" << nc << std::endl;
            for(auto current = 0; current < nc; current++)
            {
                auto cell = quadrant_array_index(cells, current);
                mark_cell(forest, i, cell);
            }
        }
        return;
    }

    void mark_cell(Forest * forest, Topidx which_tree, Quadrant * quad)
    {
        Int sum  = 0;
        Float phi = 0.0;

        auto l = quad->level; // 叶子单元的层数
        auto qx = quad->x; // qx, qy, qz 是逻辑空间中的编号
        auto qy = quad->y; 
        auto qz = quad->z;
        Point p; // 从逻辑空间实际空间坐标的转换
        auto h = P4EST_QUADRANT_LEN(quad->level); 
        for(auto i = 0; i < 8; i++)
        {
            Int i0 = i/4;
            Int i1 = (i - i0*4)/2;
            Int i2 = i - i0*4 - i1*2; 
            qcoord_to_vertex(forest->connectivity, which_tree, qx + i2*h, qy + i1*h, qz + i0*h, p.data());
            phi = interface(p);
            sum += GK::sign(phi);
        }

        auto data = (CellData *)quad->p.user_data;
        if(abs(sum) < 8)
        {
            data->refine = 1;
        }
        else
        {
            data->refine = 0;
        }
    }


    void create_data_structure()
    {
        nodes = nodes_new(forest, NULL);
        ghost = ghost_new(forest, P4EST_CONNECT_FACE);
        mesh = mesh_new(forest, ghost, P4EST_CONNECT_FACE);

        auto NN = nodes->num_owned_indeps; // NN independent nodes 
        auto NC = forest->local_num_quadrants;

        Int idx[8] = {0, 1, 3, 2, 4, 5, 7, 6};
        ds.init(NN, NC, nodes->local_nodes, idx);

        Locidx * n2c = P4EST_ALLOC(Locidx, NN);
        memset(n2c, -1, NN*sizeof(Locidx));

        auto k = 0;
        for(auto i = 0; i < NC; i++)
        {
            for(auto j = 0; j < 8; k++, j++)
            {
                auto id = nodes->local_nodes[k];
                if(n2c[id] < 0)
                    n2c[id] = k; // the first appear location in local_nodes!
            }
        }

        // traverse all local trees, now one process, one tree
        auto tree = tree_array_index(forest->trees, forest->first_local_tree);
        auto quads = &tree->quadrants;

        Int i = 0;
        for(auto& p:ds.nodes)
        {
            auto gidx = n2c[i]/8;
            auto lidx = n2c[i] - 8*gidx;
            auto quad = quadrant_array_index(quads, gidx);

            auto qx = quad->x; 
            auto qy = quad->y; 
            auto qz = quad->z; // qx, qy, qz 是逻辑空间中的编号
            auto h = P4EST_QUADRANT_LEN(quad->level); 

            auto i0 = lidx/4;
            auto i1 = (lidx - i0*4)/2;
            auto i2 = lidx - i0*4 - i1*2;
            qcoord_to_vertex(forest->connectivity, forest->first_local_tree, qx+i2*h, qy+i1*h, qz+i0*h, p.data());
            i++;
        }

        P4EST_FREE(n2c);
    }


    void write_to_vtk(int rank, int nproc, std::ostringstream & fileName)
    {
        using vtkPointsP                    = vtkSmartPointer<vtkPoints>;
        using vtkUnstructuredGridP          = vtkSmartPointer<vtkUnstructuredGrid>;
        using vtkHexahedronP                  = vtkSmartPointer<vtkHexahedron>;
        using vtkXMLUnstructuredGridWriterP = vtkSmartPointer<vtkXMLUnstructuredGridWriter>;
        using vtkIntArrayP                  = vtkSmartPointer<vtkIntArray>;

        auto fname = fileName.str();
        fileName << "_";
        fileName.fill('0');
        fileName.width(4);
        fileName <<std::right << rank;
        auto writer = vtkXMLUnstructuredGridWriterP::New();
        fileName << "." << writer->GetDefaultFileExtension();
        writer->SetFileName((fileName.str()).c_str());

        auto pts = vtkPointsP::New();
        
        auto NN = ds.number_of_nodes();
        pts->SetNumberOfPoints(NN);

        for(const auto & node:ds.nodes)
        {
            pts->SetPoint(i, node[0], node[1], node[2]);
        }

        vtkHexahedronP hex = vtkHexahedronP::New();
        for(const auto & cell : ds.cells)
        {
            int num = 0;
            for(const auto & id : cell)
            {
                (hex->GetPointIds())->SetId(num, id);
                ++num;
            }
            dataSet->InsertNextCell(hex->GetCellType(), hex->GetPointIds());
        }

        return;
    }

    void to_vtk(const char *file)
    {
        vtk_write_file(forest, NULL, file);
    }

    void destroy_forest()
    {
        /* Destroy the p4est and the connectivity structure. */
        if(forest != NULL)
            forest_destroy(forest);
        if(conn != NULL)
            connectivity_destroy(conn);
        if(nodes != NULL)
            nodes_destroy(nodes);
        if(mesh != NULL)
            mesh_destroy(mesh);
        if(ghost != NULL)
            ghost_destroy(ghost);
    }

    void print()
    {
        ds.print();
    }
};

} // end of namespace MeshAlg

} // end of namespace WHYSC
#endif // end of Interface_mesh_generator_3_h
