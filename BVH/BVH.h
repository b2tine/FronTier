#include "BV_Node.h"

#ifndef BVH_H
#define BVH_H



class BVH
{
    private:
        
        BV_iNode* root{nullptr};

    public:

        BVH() = default;
        ~BVH() = default;

        //delete copy and move ops until there
        //is a good reason not to.
        BVH(const BVH&) = delete;
        BVH& operator=(const BVH&) = delete;
        BVH(BVH&&) = delete;
        BVH& operator=(BVH&&) = delete;

        const BV_iNode* const getRoot() const {return root;}

        static BV_Leaf* createLeaf(FT_HSE* h);
        static BV_iNode* createNode(BV_Node* lc, BV_Node* rc);

};




#endif
