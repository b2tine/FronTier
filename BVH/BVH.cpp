#include "BVH.h"

BV_Leaf* BVH::createLeaf(FT_HSE* h);
{
    return new BV_Leaf(h);
}

BV_iNode* BVH::createNode(BV_Node* lc, BV_Node* rc);
{
    return new BV_iNode(lc,rc);
}

