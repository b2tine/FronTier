#include "BV.h"

#ifndef BVH_NODE_H
#define BVH_NODE_H


class BVH_Node
{
    private:


    public:

        BVH_Node() = default;
        ~BVH_Node() = default;

        //delete copy and move ops until there
        //is a good reason not to.
        BVH_Node(const BVH_Node&) = delete;
        BVH_Node& operator=(const BVH_Node&) = delete;
        BVH_Node(BVH_Node&&) = delete;
        BVH_Node& operator=(BVH_Node&&) = delete;

        explicit BVH_Node(FT_HSE*);
};














#endif




