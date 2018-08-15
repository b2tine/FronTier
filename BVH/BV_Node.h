#include "BV.h"

#ifndef BV_NODE_H
#define BV_NODE_H

template<typename BV>
class BV_Node
{
    public:

        BV* bv{nullptr};
        //TODO: need parent and children nodes

        BV_Node() = default;
        virtual ~BV_Node() {delete bv;}

        //delete copy and move ops until there
        //is a good reason not to.
        BV_Node(const BV_Node&) = delete;
        BV_Node& operator=(const BV_Node&) = delete;
        BV_Node(BV_Node&&) = delete;
        BV_Node& operator=(BV_Node&&) = delete;

        //TODO: Need constructor with children as args.
        //      Requires a merge method for the AABB BV type.

        //returns a const pointer to a const BV
        const BV* const getBV() const {return bv;}

        virtual const bool isLeaf() const {return false;}
        virtual const bool isRoot() const
        {
            //TODO: implement this
        }

};

template<typename BV>
class BV_Leaf : public BV_Node<BV>
{
    public:

        FT_HSE* hse{nullptr};

        BV_Leaf() = default;
        ~BV_Leaf() = default;

        //delete copy and move ops until there
        //is a good reason not to.
        BV_Leaf(const BV_Leaf&) = delete;
        BV_Leaf& operator=(const BV_Leaf&) = delete;
        BV_Leaf(BV_Leaf&&) = delete;
        BV_Leaf& operator=(BV_Leaf&&) = delete;

        //NOTE: hse is a shallow copy of h
        explicit BV_Leaf(FT_HSE* h)
            : hse{h}
        {
            this->bv = new BV(h);
        }

        //returns a const pointer to a const FT_HSE
        const FT_HSE* const getHse() const {return hse;}
        const bool isRoot() const override {return false;}
        const bool isLeaf() const override {return true;}


};


#endif




