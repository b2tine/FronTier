#include "BV_Node.h"

#ifndef BVH_H
#define BVH_H



class BVH
{
    private:
        
        BV_Node* root{nullptr};

    public:

        BVH() = default;
        ~BVH() = default;

        //delete copy and move ops until there
        //is a good reason not to.
        BVH(const BVH&) = delete;
        BVH& operator=(const BVH&) = delete;
        BVH(BVH&&) = delete;
        BVH& operator=(BVH&&) = delete;

        const BV_Node* const getRoot() const {return root;}

        static BV_Leaf* createLeaf(FT_HSE* h);
        static BV_Node* createNode(BV_Node* lc, BV_Node* rc);

};




#endif
