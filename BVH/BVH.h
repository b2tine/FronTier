#include "BVH_Node.h"

#ifndef BVH_H
#define BVH_H



class BVH
{
    private:
        
        RootNode* root{nullptr};

    public:

        BVH() = default;
        ~BVH() {delete root;}

        BVH(const BVH&) = delete;
        BVH& operator=(const BVH&) = delete;
        BVH(BVH&&) = delete;
        BVH& operator=(BVH&&) = delete;

        const RootNode* const getRoot() const {return root;}

        //static LeafNode* createLeaf(FT_HSE* h);
        //static InternalNode* createNode(BV_Node* lc, BV_Node* rc);

};




#endif
