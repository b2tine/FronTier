#include "BV.h"

#ifndef BVH_NODE_H
#define BVH_NODE_H


class BVH_Node
{
    private:

        FT_HSE* hse;

        AABB* getAABB();//temporary to put off introducing template

    public:

        BVH_Node() = default;
        ~BVH_Node() = default;

        //delete copy and move ops until there
        //is a good reason not to.
        BVH_Node(const BVH_Node&) = delete;
        BVH_Node& operator=(const BVH_Node&) = delete;
        BVH_Node(BVH_Node&&) = delete;
        BVH_Node& operator=(BVH_Node&&) = delete;

        //NOTE: hse is a shallow copy of h
        explicit BVH_Node(FT_HSE* h)
            : hse{h}
        {}

        //returns a const pointer to a const FT_HSE
        const FT_HSE* const getHse() const
        {
            return hse;
        }

        /*
        //returns a const pointer to a const BV
        const BV* const getBV() const
        {
            return getAABB();
        }*/
};














#endif




