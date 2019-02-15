#ifndef Mesh3d_h
#define Mesh3d_h

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

namespace WHYSC {

namespace Mesh {

template<typename GK, typename Cell>
struct Mesh3d
{
    /*typename*/
    typedef typename GK::Float Float;
    typedef typename GK::Int  Int;
    typedef typename GK::Point_3 Point;
    typedef typename Cell::Edge Edge;
    typedef typename Cell::Face Face;
    typedef typename Cell::IFace IFace;
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

    /*function member*/
    template<typename NodeContainer, typename CellContainer>
    Mesh3d(Int NN, Int NC, NodeContainer nc, CellContainer cc)
    {// 初始化

        dim = 3;
        nodes.resize(NN);

        for(auto i = 0 ; i < NN; i++)
        { 
            for(auto d = 0; d < dim; d++)
                nodes[i][d] = nc[dim*i+d];
        }

        cells.resize(NC);
        for(auto i = 0; i < NC; i++)
        {
            auto V = Cell::NV[dim];
            for(auto j = 0; j < V; j++)
                cells[i][j] = cc[V*i + j];
        }
        construct();

    }

    void construct()
    {// 重建拓扑关系数组
        auto NC = cells.size();
        std::vector<IFace> totalFaces(Cell::ND[dim-1]*NC);
        Int k = 0;
        for(auto i = 0; i < NC; i++)
        {
            auto & cell = cells[i];
            for(auto j = 0; j < Cell::ND[dim-1]; j ++)
            {
                totalFaces[k].first = k;
                
                for(auto m = 0; m < Cell::NV[dim-1]; m++)
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
        std::sort(totalFaces.begin(), totalFaces.end(), Less);

        std::list<Int> i0;
        std::list<Int> i1;
        int i = 0;
        for(; i < totalFaces.size() - 1; i++)
        {
            auto & idx0 = totalFaces[i].first;
            i0.push_back(idx0);

            auto & face0 = totalFaces[i].second;
            auto & face1 = totalFaces[i+1].second;
            if(Equal(face0, face1))
            {
                auto & idx1 = totalFaces[i+1].first;
                i1.push_back(idx1);
                i++;
            }
            else
            {
                i1.push_back(idx0);
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

    static bool Less(const IFace& f0, const IFace& f1)
    {
        auto & face0 = f0.second;
        auto & face1 = f1.second;

        for(auto i = 0; i < face0.size(); i++)
        {
            if(face0[i] < face1[i])
            {
                return true;
            }
            else if(face0[i] == face1[i])
            {
                continue;
            }
            else
                return false;
        }
    }

    static bool Equal(const Face & f0, const Face & f1)
    {
        bool r = true;
        for(auto i=0; i < f0.size(); i++)
        {
            if(f0[i] != f1[i])
                return false;
            else
                continue;
        }
        return r;  
    }

    Int number_of_cells() { return cells.size();}
    Int number_of_nodes() { return nodes.size();}
    Int number_of_edges() { return edges.size();}
    Int number_of_faces() { return faces.size();}

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
    }
};

} // end of namespace Mesh

} // end of namespace WHYSC
#endif // end of Mesh3d_h
