#include "BV_Node.h"


/////////////////////////////////////////////////////////
////////            BV_Node Methods             ////////
///////////////////////////////////////////////////////

BV_Node::BV_Node(BV_Node* lc, BV_Node* rc)
{
    assert(lc != nullptr && rc != nullptr);
    setLeft(lc);    lc->setParent(this);
    setRight(rc);   rc->setParent(this);
    bv = AABB(lc->getBV(),rc->getBV());
}

const BV_Node* const BV_Node::getSibling() const
{
    auto p = this->getParent();
    if( p != nullptr )
    {
        return p->getLeft() == this ?
            p->getRight() : p->getLeft();

    }
}


/////////////////////////////////////////////////////////
////////            BV_Leaf Methods             ////////
///////////////////////////////////////////////////////


BV_Leaf::BV_Leaf(FT_HSE* h)
    : hse{h}
{
    setBV(AABB(h));
}

