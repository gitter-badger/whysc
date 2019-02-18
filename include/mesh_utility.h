#ifndef mesh_utility_h
#define mesh_utility_h

namespace WHYSC {
namespace Mesh {

template<typename Entity>
static bool Less(const Entity & e0, const Entity& e1)
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

template<typename Entity>
static bool Equal(const Entity & e0, const Entity & e1)
{
    bool flag = true;
    for(auto i=0; i < e0.size(); i++)
    {
        if(e0[i] != e1[i])
            return false;
        else
            continue;
    }
    return flag;  
}

} // end of namespace Mesh
} // end of namespace WHYSC
#endif // end of mesh_utility_h
