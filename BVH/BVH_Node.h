#include "BoundingVolume.h"

#ifndef BVH_NODE_H
#define BVH_NODE_H

class InternalNode;
using BoundingVolume = AABB;


class BVH_Node
{
    private:

        BoundingVolume bv;
        InternalNode* parent{nullptr};
        
    public:

        BVH_Node() = default;
        virtual ~BVH_Node() = default;

        BVH_Node(const BVH_Node&) = delete;
        BVH_Node& operator=(const BVH_Node&) = delete;
        BVH_Node(BVH_Node&&) = delete;
        BVH_Node& operator=(BVH_Node&&) = delete;

        virtual const bool isLeaf() const = 0;
        virtual const bool isRoot() const = 0;

        BoundingVolume& getBV() {return bv;}
        const InternalNode* const getParent() const {return parent;}
        
        void setBV(const BoundingVolume& BV) {bv = BV;}
        void setParent(InternalNode* p) {parent = p;}

        //static void createLeafBV(Hse*);
        //static void createParentBV(const BoundingVolume&,const BoundingVolume&);
    
};


//has a Hse* and does not have children
class LeafNode : public BVH_Node
{
    private:

        Hse* hse{nullptr};

    public:

        LeafNode(Hse*);
        LeafNode() = default;
        ~LeafNode() = default;

        LeafNode(const LeafNode&) = delete;
        LeafNode& operator=(const LeafNode&) = delete;
        LeafNode(LeafNode&&) = delete;
        LeafNode& operator=(LeafNode&&) = delete;

        const bool isLeaf() const override {return true;}
        const bool isRoot() const override {return false;}

        const Hse* const getHse() const {return hse;}
};


class InternalNode : public BVH_Node
{
    private:

        BVH_Node* left{nullptr};
        BVH_Node* right{nullptr};

        void setLeft(BVH_Node* lc) {left = lc;}
        void setRight(BVH_Node* rc) {right = rc;}

    public:

        InternalNode() = default;
        virtual ~InternalNode() = default;

        InternalNode(const InternalNode&) = delete;
        InternalNode& operator=(const InternalNode&) = delete;
        InternalNode(InternalNode&&) = delete;
        InternalNode& operator=(InternalNode&&) = delete;

        InternalNode(InternalNode* lc, InternalNode* rc);

        const BVH_Node* const getLeft() const {return left;}
        const BVH_Node* const getRight() const {return right;}
        
        const bool isLeaf() const override {return false;}
        const bool isRoot() const override {return false;}
};


class RootNode : public InternalNode
{

    public:

        const bool isRoot() const override {return true;}

};


#endif


