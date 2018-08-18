#include "BV.h"

#ifndef BV_NODE_H
#define BV_NODE_H

/* Using AABB as the Bounding Volume for now.
 * will make BV_NodeBase and derived classes templates,
 * or create a class heirarchy of BV types ... */

class BV_NodeBase
{
    private:

        //Moved AABB bv into subclasses to avoid having to
        //write copy constructor/assignment. Uncopyable BV
        //seems appropriate for now.
        
        BV_NodeBase* parent{nullptr};
        //may want to use a shared_ptr, but also may not
        //matter since there won't be any tree balancing
        //(at least I don't think there will be ...)

    public:

        BV_NodeBase() = default;
        virtual ~BV_NodeBase() = default;//{delete bv;}

        //delete copy and move ops until there
        //is a good reason not to.
        BV_NodeBase(const BV_NodeBase&) = delete;
        BV_NodeBase& operator=(const BV_NodeBase&) = delete;
        BV_NodeBase(BV_NodeBase&&) = delete;
        BV_NodeBase& operator=(BV_NodeBase&&) = delete;

        void setParent(BV_NodeBase* p) {parent = p;}
        const BV_NodeBase* const getParent() const {return parent;}

        virtual const AABB& getBV() const = 0;
        virtual const bool isLeaf() const = 0;
        virtual const bool isRoot() const = 0;
};

//has a FT_HSE* and does not have children
class BV_Leaf : public BV_NodeBase
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
class BV_Node : public BV_NodeBase
{
    private:

        AABB bv;
        BV_NodeBase* left{nullptr};
        BV_NodeBase* right{nullptr};

    public:

        BV_Node() = default;
        ~BV_Node() = default;

        //delete copy and move ops until there
        //is a good reason not to.
        BV_Node(const BV_Node&) = delete;
        BV_Node& operator=(const BV_Node&) = delete;
        BV_Node(BV_Node&&) = delete;
        BV_Node& operator=(BV_Node&&) = delete;

        BV_Node(BV_NodeBase* lc, BV_NodeBase* rc)
            : bv{lc->getBV(),rc->getBV()},
            left{lc}, right{rc}
        {
            lc->setParent(this);
            rc->setParent(this);
        }

        const AABB& getBV() const override {return bv;}
        const BV_NodeBase* const getLeft() const {return left;}
        const BV_NodeBase* const getRight() const {return right;}
        const bool isLeaf() const override {return false;}
        const bool isRoot() const override {return false;}
        //TODO: Implement isRoot(); is dummy right now.
        //
        //      Note: This should only be called after
        //      BVH construction has been completed.
};


#endif
