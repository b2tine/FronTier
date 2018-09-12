#include "BVH_Node.h"

#ifndef BVH_H
#define BVH_H



class BVH
{
    private:
        
        std::shared_ptr<RootNode> root;

    public:

        //TODO: Will want to enforce invariant
        //      that root is always initialized,
        //      but requires a build routine
        //      since we want a bottom up construction
              
        BVH() = default;
        ~BVH() = default;

        BVH(const BVH&) = delete;
        BVH& operator=(const BVH&) = delete;
        BVH(BVH&&) = delete;
        BVH& operator=(BVH&&) = delete;

        const std::weak_ptr<const RootNode> getRoot() const
        {
            return std::weak_ptr<RootNode>(root);
        }

        //static std::shared_ptr<LeafNode> createLeaf(FT_HSE* h);
        /*static std::shared_ptr<InternalNode>
            createInternalNode(BV_Node* lc, BV_Node* rc);*/

};




#endif
