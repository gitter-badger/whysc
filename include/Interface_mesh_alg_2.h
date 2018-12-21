#ifndef Interface_mesh_alg_2_h
#define Interface_mesh_alg_2_h

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <array>
#include <map>
#include <cassert>
#include <cmath>
#include <utility>

#include "Geometry/Geometry_kernel.h"
#include "p4est_interface.h"


#define DUMP(a) \
    do { std::cout << #a " is value " << (a) << std::endl; } while(false)

namespace WHYSC {

namespace MeshAlg {

using namespace P4est;

template <typename T> int sign(T val) 
{
    return (T(0) < val) - (val < T(0));
}


struct QuadMeshDataStructure
{
    typedef Geometry_kernel<>  GK;
    typedef typename GK::Point_2 Point_2;
    typedef typename GK::Point_3 Point_3;
    typedef typename GK::Delaunay_algorithm_2 Delaunay_algorithm_2;
    typedef typename GK::Bisection_algorithm Bisection_algorithm;
    typedef typename GK::Level_set_function Level_set_function;


    typedef std::array<int, 2> Edge; // the origin edge
    typedef std::array<int, 3> IEdge; // the edge with idx 
    typedef std::array<int, 4> DEdge; // the dual edge
    typedef std::array<int, 4> Neighbor;
    static std::array<int, 8> localEdge;

    typedef std::array<int, 4> Cell;
    

    //data member 
    std::vector<double> points;
    std::vector<Edge> edge;
    std::vector<DEdge> edge2cell;
    std::vector<Cell> cell;
    std::vector<int> celllevel;
    std::vector<int> edgelevel;
    std::vector<double> node;
    std::map<int, Point_2> edge2cnode;
    std::vector<int> nodesign;
    std::vector<double> cnode; // the new cutted node
    int NC;
    int NN;
    int NE;

    //function member
    QuadMeshDataStructure(){}

    void construct(int nn, int nc, Locidx * ln)
    {
        NN = nn;
        NC = nc;

        cell.resize(NC);
        for(auto i = 0; i < NC; i++)
        {
            cell[i][0] = ln[4*i + 0];
            cell[i][1] = ln[4*i + 1];
            cell[i][2] = ln[4*i + 3];
            cell[i][3] = ln[4*i + 2];
        }

        
        std::vector<IEdge> totalEdge(4*NC);
        // initialize total edges 
        int k = 0;
        for(auto i = 0; i < NC; i++)
        {
            totalEdge[k][0] = k;
            totalEdge[k][1] = cell[i][0];
            totalEdge[k][2] = cell[i][1];
            ++k;

            totalEdge[k][0] = k;
            totalEdge[k][1] = cell[i][1];
            totalEdge[k][2] = cell[i][2];
            ++k;

            totalEdge[k][0] = k;
            totalEdge[k][1] = cell[i][2];
            totalEdge[k][2] = cell[i][3];
            ++k;

            totalEdge[k][0] = k;
            totalEdge[k][1] = cell[i][3];
            totalEdge[k][2] = cell[i][0];
            ++k;
        }


        std::cout << std::endl;
        auto sort = [](IEdge & e){std::sort(e.begin()+1, e.end());};
        std::for_each(totalEdge.begin(), totalEdge.end(), sort);
        std::sort(totalEdge.begin(), totalEdge.end(), Less);

        std::list<int> i0;
        std::list<int> i1;
        int i = 0;
        for(; i < totalEdge.size() - 1; i++)
        {
            auto & e0 = totalEdge[i];
            i0.push_back(e0[0]);

            auto & e1 = totalEdge[i+1];
            if(Equal(e0, e1))
            {
                i1.push_back(e1[0]);
                i++;
            }
            else
            {
                i1.push_back(e0[0]);
            }
        }

        // process the last edge of total edges
        if(i == totalEdge.size() - 1)
        {
            i0.push_back(totalEdge[i][0]);
            i1.push_back(totalEdge[i][0]);
        }
        
        
        
        // i0 and i1 should have the same size
        assert(i0.size() == i1.size());
        NE = i0.size();
        edge2cell.resize(NE);
        edge.resize(NE);

        auto it0 = i0.begin();
        auto it1 = i1.begin();
        i = 0;
        for(;it0 != i0.end(); it0++, it1++, i++)
        {
            edge2cell[i][0] = (*it0)/4;
            edge2cell[i][2] = (*it0) - edge2cell[i][0]*4;
            edge2cell[i][1] = (*it1)/4;
            edge2cell[i][3] = (*it1) - edge2cell[i][1]*4;
            auto cidx = edge2cell[i][0];
            auto lidx = edge2cell[i][2]; 
            edge[i][0] = cell[cidx][localEdge[lidx*2 + 0]]; 
            edge[i][1] = cell[cidx][localEdge[lidx*2 + 1]];
        }
    }

    static bool Less(const IEdge& e0, const IEdge& e1)
    {
        if(e0[1] < e1[1])
        {
            return true;
        }
        else if(e0[1] == e1[1] )
        {
            if(e0[2] < e1[2])
                return true;
            else
                return false;
        }
        else
        {
            return false;
        }
    }

    static bool Equal(const IEdge & e0, const IEdge & e1)
    {
        return (e0[1] == e1[1]) && (e0[2] == e1[2]);
    }

    void compute_node_sign(Level_set_function & fun)
    {
        nodesign.resize(NN);
        for(int i=0; i < NN; i++)
        {
            Point_2 p(node[3*i + 0], node[3*i + 1]);
            nodesign[i] = fun.sign(p);
        }
    }

    void find_cutted_edge(Level_set_function & fun)
    {
        // we just consider the finnest edge
        edgelevel.resize(NE);
        int maxlevel = 0;
        for(auto i = 0; i < NE; i++)
        {
            edgelevel[i] = celllevel[edge2cell[i][0]];
            if(edgelevel[i] > maxlevel)
                maxlevel = edgelevel[i];
        }

        for(auto i = 0; i < NE; i++)
        {
            if(edgelevel[i] == maxlevel && nodesign[edge[i][0]]*nodesign[edge[i][1]] < 0)
            {
                Point_2 p0(node[3*edge[i][0] + 0], node[3*edge[i][0] + 1]);
                Point_2 p1(node[3*edge[i][1] + 0], node[3*edge[i][1] + 1]);
                Point_2 m = bisection(fun, p0, p1);
                edge2cnode.insert(std::pair<int, Point_2>(i, m));
            }
        }
    }


    void find_cutted_cell()
    {
    }

    bool is_cutted_cell(int i)
    {
        std::array<int, 4> sign = { 
            nodesign[cell[i][0]],
            nodesign[cell[i][1]],
            nodesign[cell[i][2]],
            nodesign[cell[i][3]]};

        int sum = std::abs(sign[0]) + std::abs(sign[1]) + std::abs(sign[2])
            return (sum < 3) || (sign[0]*sign[1] < 0) ||
            (sign[1]*sign[2] < 0) || (sign[2]*sign[3] < 0) ||  
            (sign[3]*sign[0] < 0);
    }

    int is_special_elem(Interface * mb, const EntityHandle eh)
    {
        std::vector<EntityHandle> conn;
        rval = mb->get_connectivity(&eh, 1, conn, true);
        int sign[conn.size()]; // conn.size() == 4
        rval = mb->tag_get_data(vsign, conn.data(), conn.size(), sign);
        bool flag1 = sign[0] == 0 && sign[2] == 0 && sign[1]*sign[3] < 0;
        bool flag2 = sign[1] == 0 && sign[3] == 0 && sign[0]*sign[2] < 0;

        if(flag1)
            return 1;
        else if(flag2)
            return 2;
        else 
            return 0;
    }

    void print()
    {
        int i = 0;
        for(; i < NN; i++)
        {
            std::cout << i << ": " 
                << node[3*i+0] << " "
                << node[3*i+1] << " "
                << node[3*i+2] << std::endl;
        }

        i = 0;
        for(; i < NC; i++)
        {
            std::cout << i << ": "
                << cell[i][0] << " "
                << cell[i][1] << " "
                << cell[i][2] << " "
                << cell[i][3] << " with level " << celllevel[i] << std::endl;
        }

        for(auto & e : edge)
        {
            std::cout << i++ << ": " << e[0] << ", " << e[1] << std::endl;
        }

        i = 0;
        for(auto & e2c : edge2cell)
        {
            std::cout << i++ << ": " 
                << e2c[0] << ", " << e2c[1] << ", "
                << e2c[2] << ", " << e2c[3] << std::endl;

        }
    }
};

std::array<int, 8> QuadMeshDataStructure::localEdge = {0, 1, 1, 2, 2, 3, 3, 0};

class Interface_mesh_alg_2
{
public:
    typedef Geometry_kernel<>  GK;
    typedef typename GK::Level_set_function Level_set_function;
public:
    /** Constructor
     *
     * Input:
     *
     *  
     *
     */
    Interface_mesh_alg_2(MPI_Comm mpi_comm, Level_set_function & lfun) 
    {
        comm = mpi_comm;
        conn = NULL;
        forest = NULL;
        nodes = NULL;
        ghost = NULL;
        mesh = NULL;
        fun = lfun;
    }

    void create_forest_on_rec_domain(double xmin, double xmax, double ymin, double ymax)
    {
        const Topidx num_vertices = 4;
        const Topidx num_trees = 1;
        const Topidx num_ctt = 0;
        const double vertices[4 * 3] = {
            xmin, ymin, 0,
            xmax, ymin, 0,
            xmin, ymax, 0,
            xmax, ymax, 0};
        const Topidx tree_to_vertex[1 * 4] = { 0, 1, 2, 3};
        const Topidx tree_to_tree[1 * 4] = { 0, 0, 0, 0};
        const int8_t tree_to_face[1 * 4] = { 0, 1, 2, 3};

        conn = connectivity_new_copy(
                num_vertices, num_trees, 0,
                vertices, tree_to_vertex,
                tree_to_tree, tree_to_face,
                NULL, &num_ctt, NULL, NULL);
        forest = forest_new(comm, conn, 0, NULL, NULL);
    }

    void uniform_refine(int n)
    {
        for(int i = 0; i < n; i++)
            refine(forest, 0, uniform_refine_function, NULL);
    }


    void adaptive_refine(int n)
    {
        for(int i = 0; i < n; i++)
            refine(forest, 0, adaptive_refine_function, NULL);

        int partforcoarsen = 0;
        partition(forest, partforcoarsen, NULL);
    }

    static int uniform_refine_function(Forest * forest, 
            Topidx which_tree, 
            Quadrant * quad)
    {

        return 1;
    }

    static int adaptive_refine_function(Forest * forest,
            Topidx which_tree,
            Quadrant * quad)
    {
        int sum  = 0;
        double phi = 0.0;

        auto l = quad->level; // 叶子单元的层数
        auto qx = quad->x; // qx, qy 是逻辑空间中的编号
        auto qy = quad->y; 
        double p[3] = {0.0, 0.0, 0.0}; // 网格所在的空间为 [0, 1]^2, 从逻辑空间实际空间坐标的转换
        auto h = P4EST_QUADRANT_LEN(quad->level); 

        qcoord_to_vertex(forest->connectivity, which_tree, qx, qy, p);
        phi = sqrt(p[0]*p[0] + p[1]*p[1]) - 0.8;
        sum += sign(phi);

        qcoord_to_vertex(forest->connectivity, which_tree, qx+h,   qy, p);
        phi = sqrt(p[0]*p[0] + p[1]*p[1]) - 0.8;
        sum += sign(phi);

        qcoord_to_vertex(forest->connectivity, which_tree, qx,   qy+h, p);
        phi = sqrt(p[0]*p[0] + p[1]*p[1]) - 0.8;
        sum += sign(phi);

        qcoord_to_vertex(forest->connectivity, which_tree, qx+h,   qy+h, p);
        phi = sqrt(p[0]*p[0] + p[1]*p[1]) - 0.8;
        sum += sign(phi);

        if(abs(sum) < 4)
            return 1;
        else 
            return 0;
    }



    void create_data_structure()
    {
        nodes = nodes_new(forest, NULL);
        ghost = ghost_new(forest, P4EST_CONNECT_FACE);
        mesh = mesh_new(forest, ghost, P4EST_CONNECT_FACE);

        auto NN = nodes->num_owned_indeps; // NN independent nodes 
        auto NC = forest->local_num_quadrants;

        ds.construct(NN, NC, nodes->local_nodes);

        ds.node.resize(3*NN);

        Locidx * n2c = P4EST_ALLOC(Locidx, NN);
        memset(n2c, -1, NN*sizeof(Locidx));

        auto k = 0;
        for(auto i = 0; i < NC; i++)
        {
            for(auto j = 0; j < 4; k++, j++)
            {
                auto id = nodes->local_nodes[k];
                if(n2c[id] < 0)
                    n2c[id] = k; // the first appear location in local_nodes!
            }
        }

        // traverse all local trees, one process, one tree ?
        auto tree = tree_array_index(forest->trees, forest->first_local_tree);
        auto quads = &tree->quadrants;

        ds.celllevel.resize(NC);
        for(auto i = 0; i < NC; i++)
        {
            auto quad = quadrant_array_index(quads, i);
            ds.celllevel[i] = quad->level;
        }

        for(auto i = 0; i < NN; i++)
        {
            auto gidx = n2c[i]/4;
            auto lidx = n2c[i] - 4*gidx;
            auto quad = quadrant_array_index(quads, gidx);

            auto qx = quad->x; // qx, qy 是逻辑空间中的编号
            auto qy = quad->y; 
            auto h = P4EST_QUADRANT_LEN(quad->level); 

            auto i0 = lidx/2;
            auto i1 = lidx - i0*2;
            qcoord_to_vertex(forest->connectivity, forest->first_local_tree, qx+i1*h, qy+i0*h, ds.node.data()+3*i);
        }

        P4EST_FREE(n2c);
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


    void to_vtk(const char *file)
    {
        vtk_write_file(forest, NULL, file);
    }


    void print()
    {
        if(forest != NULL)
        {
            std::cout << "Information of forest:" << std::endl;
            std::cout << forest->local_num_quadrants << std::endl;
        }

        if(mesh != NULL)
        {
            auto q2q = mesh->quad_to_quad;
            auto NC = mesh->local_num_quadrants;
            std::cout <<" quad to quad: \n";
            for(auto i = 0; i < NC; i++)
            {
                std::cout << i << ": "
                    << q2q[4*i + 0] << ", "
                    << q2q[4*i + 1] << ", "
                    << q2q[4*i + 2] << ", "
                    << q2q[4*i + 3] << std::endl;
            }

            auto k = 0;
            for(auto i = 0; i < NC; i++)
            {
                std::cout<< i << ": " << std::endl;
                for(auto j = 0; j < 4; k++, j++)
                {
                    auto id = nodes->local_nodes[k];
                    std::cout << id << " ";
                }
                std::cout << std::endl;
            }

            auto NGC = mesh->ghost_num_quadrants;
            std::cout<< "Number of ghost on rank: " << forest->mpirank << " " << NGC << std::endl;
        }

        ds.print();
    }

private:
    MPI_Comm  comm;
    Connectivity * conn;
    Forest * forest;
    Nodes * nodes;
    Ghost * ghost;
    Mesh * mesh;
    Level_set_function fun;
    QuadMeshDataStructure ds;
};

} // end of namespace MeshAlg

} // end of namespace WHYSC
#endif // end of Interface_mesh_alg_2_h
