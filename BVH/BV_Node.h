#include "BV.h"

#ifndef BV_NODE_H
#define BV_NODE_H

template<typename BV>
class BV_Node
{
    public:

        BV_Node() = default;
        virtual ~BV_Node() = default;

        //delete copy and move ops until there
        //is a good reason not to.
        BV_Node(const BV_Node&) = delete;
        BV_Node& operator=(const BV_Node&) = delete;
        BV_Node(BV_Node&&) = delete;
        BV_Node& operator=(BV_Node&&) = delete;

        virtual const BV* const getBV() const = 0;
        virtual bool isLeaf() const = 0;
        virtual bool isRoot() const = 0;

};

template<typename BV>
class BV_Leaf : public BV_Node<BV>
{
    public:

        BV* bv{nullptr};
        FT_HSE* hse{nullptr};

        BV_Leaf() = default;
        
        ~BV_Leaf()
        {
            delete bv;
        }

        //delete copy and move ops until there
        //is a good reason not to.
        BV_Leaf(const BV_Leaf&) = delete;
        BV_Leaf& operator=(const BV_Leaf&) = delete;
        BV_Leaf(BV_Leaf&&) = delete;
        BV_Leaf& operator=(BV_Leaf&&) = delete;

        //NOTE: hse is a shallow copy of h
        explicit BV_Leaf(FT_HSE* h)
            : bv{new BV(h)}, hse{h}
        {}

        //returns a const pointer to a const FT_HSE
        const FT_HSE* const getHse() const
        {
            return hse;
        }

        //returns a const pointer to a const BV
        const BV* const getBV() const override
        {
            return bv;
        }

        bool isRoot() const override {return false;}
        bool isLeaf() const override {return true;}


};


//TODO:BV_iNode class











#endif




