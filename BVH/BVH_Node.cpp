#include "BVH_Node.h"

/////////////////////////////////////////////////////////
////////           BVH_Node Methods             ////////
///////////////////////////////////////////////////////


/////////////////////////////////////////////////////////
////////            LeafNode Methods             ////////
///////////////////////////////////////////////////////

/*
LeafNode::LeafNode(Hse* h)
    : hse{h} 
{
    setBV(BoundingVolume(h));
}
*/

/////////////////////////////////////////////////////////
////////           InternalNode Methods             ////////
///////////////////////////////////////////////////////

/*
InternalNode::InternalNode(BVH_Node* lc, BVH_Node* rc)
{
    assert(lc != nullptr && rc != nullptr);

    this->setLeft(lc);
    this->setRight(rc);
    rc->setParent(this);
    lc->setParent(this);
    //bv = BVH_Node::createParentBV(lc,rc);
}*/

/*
const InternalNode* const InternalNode::getSibling() const
{
    auto p = this->getParent();
    if( p != nullptr )
    {
        return p->getLeft() == this ?
            p->getRight() : p->getLeft();
    }
}
*/




