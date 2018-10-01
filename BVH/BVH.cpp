#include "BVH.h"

std::shared_ptr<LeafNode>
BVH::createLeafNode(Hse* h)
{
    return std::make_shared<LeafNode>(h);
}

std::shared_ptr<InternalNode>
BVH::createInternalNode(std::shared_ptr<BVH_Node> lc,
        std::shared_ptr<BVH_Node> rc)
{
    auto node = std::make_shared<InternalNode>(lc,rc);
    node->setChildren(lc,rc);
    return std::move(node);
}
