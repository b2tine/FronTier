#include "BVH_Node.h"

#ifndef BVH_H
#define BVH_H



class BVH
{
    private:
        
        //std::shared_ptr<InternalNode> root;
        //std::vector<std::shared_ptr<InternalNode>> Leaves;

    public:

        //TODO: Will want to enforce the invariant that
        //      the root is always initialized, but this
        //      requires a build routine since we are
        //      performing a bottom up construction.
              
        BVH() = default;
        ~BVH() = default;

        BVH(const BVH&) = delete;
        BVH& operator=(const BVH&) = delete;
        BVH(BVH&&) = delete;
        BVH& operator=(BVH&&) = delete;

        /*
        const std::weak_ptr<const InternalNode> getRoot() const
        {
            return std::weak_ptr<InternalNode>(root);
        }
        */

        static std::shared_ptr<LeafNode> createLeafNode(Hse* h);
        static std::shared_ptr<InternalNode> createInternalNode(
                std::shared_ptr<BVH_Node> lc, std::shared_ptr<BVH_Node> rc);

};




#endif
