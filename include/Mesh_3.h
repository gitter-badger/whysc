#ifndef Mesh_3_h
#define Mesh_3_h

#include <vector>
#include <string>
#include <array>
#include <list>
#include <map>
#include <algorithm>
#include <cassert>
#include <cmath>
#include <utility>

#include <sstream>
#include <iostream>
#include <fstream>


#include "mesh_utility.h"

namespace WHYSC {

namespace Mesh {

template<typename GK, typename Cell>
struct Mesh_3
{
    /*typename*/
    typedef typename GK::Float Float;
    typedef typename GK::Int  Int;
    typedef typename GK::Point_3 Point;
    typedef typename Cell::Edge Edge;
    typedef typename Cell::Face Face;
    typedef typename Cell::IFace IFace;
    typedef typename Cell::IEdge IEdge;
    typedef typename Cell::Cell2cell Cell2cell;
    typedef typename Cell::Cell2face Cell2face;
    typedef typename Cell::Cell2edge Cell2edge;
    typedef typename Cell::Face2cell Face2cell;

    /*Data member*/
    std::vector<Point> nodes;
    std::vector<Cell> cells;
    std::vector<Edge> edges;
    std::vector<Face> faces;

    std::vector<Cell2cell> cell2cell;
    std::vector<Cell2face> cell2face;
    std::vector<Cell2edge> cell2edge;
    std::vector<Face2cell> face2cell;

    Int dim;


    Mesh_3()
    {
        dim = 3;
    }

    /*function member*/
    template<typename NodeContainer, typename CellContainer>
    void init(Int NN, Int NC, NodeContainer nc, CellContainer cc, Int * idx =NULL)
    {// 初始化

        if(idx == NULL)
        {
            if(Cell::type == HEXAHEDRON)
            {
                Int i0[8] = {0, 1, 2, 3, 4, 5, 6, 7};
                idx = i0;
            }
            else if(Cell::type == TETRA)
            {
                Int i0[4] = {0, 1, 2, 3};
                idx = i0;
            }
        }

        nodes.resize(NN);
        for(auto i = 0; i < NN; i++)
        {
            for(auto j = 0; j < dim; j++)
                nodes[i][j] = nc[dim*i +j];
        }

        cells.resize(NC);
        for(auto i = 0; i < NC; i++)
        {
            auto V = Cell::NV[dim];
            for(auto j = 0; j < V; j++)
                cells[i][j] = cc[V*i + idx[j]];
        }
        construct_faces();
        construct_edges();
        return;
    }

    template<typename NodeContainer>
    void set_nodes(Int NN, NodeContainer nc)
    {
        nodes.resize(NN);
        for(auto i = 0; i < NN; i++)
        {
            for(auto j = 0; j < dim; j++)
                nodes[i][j] = nc[dim*i +j];
        }
        return;
    }

    template<typename CellContainer>
    void init(Int NN, Int NC, CellContainer cc, Int *idx = NULL)
    {// 初始化

        nodes.resize(NN);
        if(idx == NULL)
        {
            if(Cell::type == HEXAHEDRON)
            {
                Int i0[8] = {0, 1, 2, 3, 4, 5, 6, 7};
                idx = i0;
            }
            else if(Cell::type == TETRA)
            {
                Int i0[4] = {0, 1, 2, 3};
                idx = i0;
            }
        }

        cells.resize(NC);
        for(auto i = 0; i < NC; i++)
        {
            auto V = Cell::NV[dim];
            for(auto j = 0; j < V; j++)
                cells[i][j] = cc[V*i + idx[j]];
        }
        construct_faces();
        construct_edges();
    }

    static bool FaceLess(const Face & f0, const Face& f1)
    {
        for(auto i = 0; i < f0.size(); i++)
        {
            if(f0[i] < f1[i])
            {
                return true;
            }
            else if(f0[i] == f1[i])
            {
                continue;
            }
            else
                return false;
        }
    }

    void construct_faces()
    {// 重建拓扑关系数组
        auto NC = cells.size();
        std::vector<IFace> totalFaces(Cell::ND[2]*NC);
        Int k = 0;
        for(auto i = 0; i < NC; i++)
        {
            auto & cell = cells[i];
            for(auto j = 0; j < Cell::ND[2]; j ++)
            {
                totalFaces[k].first = k;
                for(auto m = 0; m < Cell::NV[2]; m++)
                {
                    auto n = Cell::face[j][m];
                    totalFaces[k].second[m] = cell[n];
                }
                ++k;
            }
        }

        // sort the `totalFaces`
        auto sort = [](IFace & f){std::sort(f.second.begin(), f.second.end());};
        std::for_each(totalFaces.begin(), totalFaces.end(), sort);
        std::sort(totalFaces.begin(), totalFaces.end(), FaceLess);

        std::list<Int> i0;
        std::list<Int> i1;
        Int i = 0;
        for(; i < totalFaces.size() - 1; i++)
        {
            i0.push_back(totalFaces[i].first);

            auto & face0 = totalFaces[i].second;
            auto & face1 = totalFaces[i+1].second;
            if(Equal(face0, face1))
            {
                i1.push_back(totalFaces[i+1].first);
                i++;
            }
            else
            {
                i1.push_back(totalFaces[i].first);
            }
        }
        // process the last face in totalFaces 
        if(i == totalFaces.size() - 1)
        {
            auto idx = totalFaces[i].first;
            i0.push_back(idx);
            i1.push_back(idx);
        }
        // i0 and i1 should have the same size
        assert(i0.size() == i1.size());

        // Get the face, toplogy realitonship of face and cell
        auto NF = i0.size();
        face2cell.resize(NF);
        faces.resize(NF);
        auto it0 = i0.begin();
        auto it1 = i1.begin();
        auto F = Cell::ND[dim-1];
        for(auto i = 0;it0 != i0.end(); it0++, it1++, i++)
        {
            face2cell[i][0] = (*it0)/F;
            face2cell[i][2] = (*it0) - face2cell[i][0]*F;
            face2cell[i][1] = (*it1)/F;
            face2cell[i][3] = (*it1) - face2cell[i][1]*F;
            auto cidx = face2cell[i][0];
            auto lidx = face2cell[i][2]; 
            for(auto j = 0; j < Cell::NV[dim-1]; j++)
            {
                auto k = Cell::face[lidx][j];
                faces[i][j] = cells[cidx][k];
            }
        }

        cell2face.resize(cells.size());
        for(auto i = 0; i < NF; i++)
        {
            cell2face[face2cell[i][0]][face2cell[i][2]] = i;
            cell2face[face2cell[i][1]][face2cell[i][3]] = i;
        }
    }

    static bool EdgeLess(const Edge & e0, const Edge& e1)
    {
        for(auto i = 0; i < e0.size(); i++)
        {
            if(e0[i] < e1[i])
            {
                return true;
            }
            else if(e0[i] == e1[i])
            {
                continue;
            }
            else
                return false;
        }
    }

    void construct_edges()
    {/*construct the unique edges*/
        auto NC = cells.size();
        std::vector<IEdge> totalEdges(Cell::ND[1]*NC);
        Int k = 0;
        for(auto i = 0; i < NC; i++)
        {
            auto & cell = cells[i];
            for(auto j = 0; j < Cell::ND[1]; j++)
            {
                totalEdges[k].first = k;
                
                for(auto m = 0; m < Cell::NV[1]; m++)
                {
                    auto n = Cell::edge[j][m];
                    totalEdges[k].second[m] = cell[n];
                }
                ++k;
            }
        }

        // sort the `totalEdges`
        auto sort = [](IEdge & e){std::sort(e.second.begin(), e.second.end());};
        std::for_each(totalEdges.begin(), totalEdges.end(), sort);
        std::sort(totalEdges.begin(), totalEdges.end(), EdgeLess);

        std::list<IEdge> es; 
        es.push_back(totalEdges[0]);

        Int i = 0;
        Int idx = totalEdges[0].first;
        cell2edge[idx/Cell::ND[1]][idx%Cell::ND[1]] = i;
        for(auto j = 1; j < totalEdges.size(); j++)
        {
            auto & e0 = totalEdges[j-1].second;
            auto & e1 = totalEdges[j].second;
            if(!Equal(e0, e1))
            {
                es.push_back(totalEdges[j]);
                i++;
            }
            idx = totalEdges[j].first;
            cell2edge[idx/Cell::ND[1]][idx%Cell::ND[1]] = i;
        }

        edges.resize(es.size());
        for(auto& e:es)
        {
            edges.push_back(e.second);
        }
    }


    Int number_of_nodes() { return nodes.size();}
    Int number_of_edges() { return edges.size();}
    Int number_of_faces() { return faces.size();}
    Int number_of_cells() { return cells.size();}

    void print()
    {
        std::cout << "Begin output mesh3d datastructure information:" << std::endl;
        std::cout << "The nodes: " << std::endl;
        int i = 0;
        for(auto & node:nodes)
        {
            std::cout << i++ << ": " << node << std::endl;
        }

        std::cout << "The cells:" << std::endl;
        i = 0;
        for(auto & cell:cells)
        {
            std::cout << i++ << ": " << cell << std::endl;
        }

        std::cout << "The faces:" << std::endl;
        i = 0;
        for(auto & face:faces)
        { 
            std::cout << i++ << ": ";
            for(auto & idx:face)
                std::cout << idx << ", ";
            std::cout << std::endl;
        }

        std::cout << "The face2cell:" << std::endl;
        i = 0;
        for(auto & f2c : face2cell)
        {
            std::cout << i++ << ": ";
            for(auto & idx:f2c)
                std::cout << idx << ", ";
            std::cout << std::endl;
        }

        std::cout << "The edges:" << std::endl;
        i = 0;
        for(auto & e : edges)
        {
            std::cout << i++ << ": ";
            for(auto & idx : e)
                std::cout << idx << ", ";
            std::cout << std::endl;
        }
    }
};

} // end of namespace Mesh

} // end of namespace WHYSC
#endif // end of Mesh_3_h
