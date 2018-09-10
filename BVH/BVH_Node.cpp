#include "BVH_Node.h"

/////////////////////////////////////////////////////////
////////           BVH_Node Methods             ////////
///////////////////////////////////////////////////////


BoundingVolume* BVH_Node::createLeafBV(FT_HSE* h, const BV_Type& bvtype)
{
    switch(bvtype)
    {
        case BV_Type::AABB:
            return new AABB(h);
            break;
    }
}

BoundingVolume* BVH_Node::createParentBV(BVH_Node* lc, BVH_Node* rc)
{
    BoundingVolume* lbv = lc->getBV(); 
    BoundingVolume* rbv = rc->getBV(); 
    assert(lbv->getBvType() == rbv->getBvType());

    BV_Type bvtype = lbv->getBvType();
    switch(bvtype)
    {
        case BV_Type::AABB:
            return new AABB();
            break;
    }
}

/////////////////////////////////////////////////////////
////////           InternalNode Methods             ////////
///////////////////////////////////////////////////////

InternalNode::InternalNode(BVH_Node* lc, BVH_Node* rc)
{
    assert(lc != nullptr && rc != nullptr);

    this->setLeft(lc);
    this->setRight(rc);
    rc->setParent(this);
    lc->setParent(this);
    //TODO: need factory function for BoundingVolume class
    bv = BVH_Node::createParentBV(lc,rc);
}

const InternalNode* const InternalNode::getSibling() const
{
    auto p = this->getParent();
    if( p != nullptr )
    {
        return p->getLeft() == this ?
            p->getRight() : p->getLeft();
    }
}


/////////////////////////////////////////////////////////
////////            LeafNode Methods             ////////
///////////////////////////////////////////////////////


LeafNode::LeafNode(FT_HSE* h, BV_Type bvtype)
    : hse{h}
{
    setBV(BVH_Node::createLeafBV(h,bvtype));
}

