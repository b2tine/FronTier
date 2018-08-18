#include "BV.h"

#ifndef BV_NODE_H
#define BV_NODE_H

/* Using AABB as the Bounding Volume for now.
 * will make BV_Node and derived classes templates,
 * or create a class heirarchy of BV types ... */

class BV_Node
{
    private:

        //Moved AABB bv into subclasses to avoid having to
        //write copy constructor/assignment. Uncopyable BV
        //seems appropriate for now.
        
        BV_Node* parent{nullptr};
        //may want to use a shared_ptr, but also may not
        //matter since there won't be any tree balancing
        //(at least I don't think there will be ...)

    public:

        BV_Node() = default;
        virtual ~BV_Node() = default;//{delete bv;}

        //delete copy and move ops until there
        //is a good reason not to.
        BV_Node(const BV_Node&) = delete;
        BV_Node& operator=(const BV_Node&) = delete;
        BV_Node(BV_Node&&) = delete;
        BV_Node& operator=(BV_Node&&) = delete;

        void setParent(BV_Node* p) {parent = p;}
        const BV_Node* const getParent() const {return parent;}

        virtual const AABB& getBV() const = 0;
        virtual const bool isLeaf() const = 0;
        virtual const bool isRoot() const = 0;
};

//has a FT_HSE* and does not have children
class BV_Leaf : public BV_Node
{
    private:

        AABB bv;
        FT_HSE* hse{nullptr};

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
        explicit BV_Leaf(FT_HSE* h)
            : bv{h}, hse{h}
        {}

        const AABB& getBV() const override {return bv;}
        const FT_HSE* const getHse() const {return hse;}
        const bool isLeaf() const override {return true;}
        const bool isRoot() const override {return false;}
};

//has children, does not have FT_HSE*.
class BV_iNode : public BV_Node
{
    private:

        AABB bv;
        BV_Node* left{nullptr};
        BV_Node* right{nullptr};

    public:

        BV_iNode() = default;
        ~BV_iNode() = default;

        //delete copy and move ops until there
        //is a good reason not to.
        BV_iNode(const BV_iNode&) = delete;
        BV_iNode& operator=(const BV_iNode&) = delete;
        BV_iNode(BV_iNode&&) = delete;
        BV_iNode& operator=(BV_iNode&&) = delete;

        BV_iNode(BV_Node* lc, BV_Node* rc)
            : bv{lc->getBV(),rc->getBV()},
            left{lc}, right{rc}
        {
            lc->setParent(this);
            rc->setParent(this);
        }

        const AABB& getBV() const override {return bv;}
        const BV_Node* const getLeft() const {return left;}
        const BV_Node* const getRight() const {return right;}
        const bool isLeaf() const override {return false;}
        const bool isRoot() const override {return false;}
        //TODO: Implement isRoot(); is dummy right now.
        //
        //      Note: This should only be called after
        //      BVH construction has been completed.
};


#endif
