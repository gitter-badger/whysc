#include <iostream>
#include <vector>
#include <tuple>
#include <algorithm>
#include <iterator>
#include <array>



typedef std::tuple<int, int, int> Edge;

bool Less(const Edge& e0, const Edge& e1)
{
    if(std::get<1>(e0) < std::get<1>(e1))
        return true;
    else if(std::get<1>(e0) == std::get<1>(e1))
    {
        if(std::get<2>(e0) < std::get<2>(e1))
            return true;
        else
            return false;
    }
    else
        return false;
}

bool Equal(const Edge & e0, const Edge & e1)
{
    return (std::get<1>(e0) == std::get<1>(e1)) && (std::get<2>(e0) == std::get<2>(e1));
}

int main(int argc, char **argv)
{

    std::array<int, 3> a;
    a = {0, 1, 2};

    std::copy(a.begin(), a.end(), std::ostream_iterator<int>(std::cout, " "));

    int cell[] = {
        0, 1, 3, 4, 
        1, 2, 4, 5,
        3, 4, 6, 7,
        4, 5, 7, 8};

    std::vector<Edge> totalEdge(16);

    int k = 0;
    for(auto i = 0; i < 4; i++)
    {
        totalEdge[k++] = std::make_tuple(k, cell[4*i + 2], cell[4*i + 0]);
        totalEdge[k++] = std::make_tuple(k, cell[4*i + 1], cell[4*i + 3]);
        totalEdge[k++] = std::make_tuple(k, cell[4*i + 0], cell[4*i + 1]);
        totalEdge[k++] = std::make_tuple(k, cell[4*i + 2], cell[4*i + 3]);
    }

    for(auto & edge : totalEdge)
    {
        if(std::get<1>(edge) > std::get<2>(edge))
        {
            auto a = std::get<1>(edge);
            std::get<1>(edge) = std::get<2>(edge);
            std::get<2>(edge) = a;
        }

    }

    std::sort(totalEdge.begin(), totalEdge.end(), Less);

    for(auto edge : totalEdge)
        std::cout<< std::get<0>(edge) << ":" << std::get<1>(edge) << ", " << std::get<2>(edge)  << std::endl;

    std::vector<int> i0;
    std::vector<int> i1;
    i0.reserve(16);
    i1.reserve(16);
    int i = 0;
    for(; i < totalEdge.size() - 1; i++)
    {
        auto & e0 = totalEdge[i];
        auto & e1 = totalEdge[i+1];
        if(Equal(e0, e1))
        {
            i0.push_back(std::get<0>(e0));
            i1.push_back(std::get<0>(e1));
            i++;
        }
        else
        {
            i0.push_back(std::get<0>(e0));
            i1.push_back(std::get<0>(e0));
        }
    }

    std::cout << "final i: " << i << std::endl;

    if(i == totalEdge.size() - 1)
    {
        i0.push_back(std::get<0>(totalEdge[i]));
        i1.push_back(std::get<0>(totalEdge[i]));
    }

    std::cout << "i0:" << std::endl;
    std::copy(i0.begin(), i0.end(), std::ostream_iterator<int>(std::cout, " "));
    std::cout << std::endl;
    std::cout << "i1:" << std::endl;
    std::copy(i1.begin(), i1.end(), std::ostream_iterator<int>(std::cout, " "));
    std::cout<< std::endl;

    for(auto edge : totalEdge)
        std::cout<< std::get<0>(edge) << ":" << std::get<1>(edge) << ", " << std::get<2>(edge)  << std::endl;
    
    
    return 0;
}
