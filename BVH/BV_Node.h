#include "BV.h"

#ifndef BV_NODE_H
#define BV_NODE_H

/* Using AABB as the Bounding Volume for now.
 * will make BV_Node and derived classes templates,
 * or create a class heirarchy of BV types ... */


//NOTE: Scott says concrete classes should not inherit from other concrete classes (item 33 of "More Effective C++") ... first design may be better.


//Does not have FT_HSE*.
class BV_Node
{
    private:

        BV_Node* left{nullptr};
        BV_Node* right{nullptr};
        BV_Node* parent{nullptr};
        //may want to use std::shared_ptr<BV_Node> but it
        //might not matter since there won't be any tree balancing
        //(at least I don't think there will be ...)

        void setLeft(BV_Node* lc) {left = lc;}
        void setRight(BV_Node* rc) {right = rc;}

    protected:

        AABB bv;

    public:

        BV_Node() = default;
        ~BV_Node() = default;

        //delete copy and move ops until there
        //is a good reason not to.
        BV_Node(const BV_Node&) = delete;
        BV_Node& operator=(const BV_Node&) = delete;
        BV_Node(BV_Node&&) = delete;
        BV_Node& operator=(BV_Node&&) = delete;

        BV_Node(BV_Node* lc, BV_Node* rc);

        void setParent(BV_Node* p) {parent = p;}

        const AABB& getBV() const {return bv;}
        const BV_Node* const getLeft() const {return left;}
        const BV_Node* const getRight() const {return right;}
        const BV_Node* const getParent() const {return parent;}
        const BV_Node* const getSibling() const;
        
        virtual const bool isLeaf() const {return false;}
        virtual const bool isRoot() const {return false;}
        //TODO: Implement isRoot(); is dummy right now.
        //
        //      Note: This should only be called after
        //      BVH construction has been completed.
};

//has a FT_HSE* and does not have children
class BV_Leaf : public BV_Node
{
    private:

        FT_HSE* hse{nullptr};
        void setBV(const AABB& BB) {this->bv = BB;}

    public:

        BV_Leaf() = default;
        ~BV_Leaf() = default;

        //delete copy and move ops until there
        //is a good reason not to.
        BV_Leaf(const BV_Leaf&) = delete;
        BV_Leaf& operator=(const BV_Leaf&) = delete;
        BV_Leaf(BV_Leaf&&) = delete;
        BV_Leaf& operator=(BV_Leaf&&) = delete;

        //NOTE: hse is a shallow copy of h
        explicit BV_Leaf(FT_HSE* h);

        const FT_HSE* const getHse() const {return hse;}
        const bool isLeaf() const override {return true;}
        const bool isRoot() const override {return false;}
};


#endif
