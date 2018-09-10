#include "BoundingVolume.h"

#ifndef BVH_NODE_H
#define BVH_NODE_H


class BVH_Node
{
    protected:
        BVH_Node* parent{nullptr};
        BoundingVolume* bv{nullptr};

        static void createLeafBV(FT_HSE*, const BV_Type&);
        static void createParentBV(BoundingVolume*, BoundingVolume*);
    
    public:
        const BoundingVolume* const getBV() const {return bv;}

        virtual const bool isLeaf() const;
        virtual const bool isRoot() const;
        virtual ~BVH_Node() {delete bv;}
};


//An internal node.
//Decided against the name BVHH_iNode since we have
//BVHH_Leaf and BVHH_Root
class InternalNode : public BVH_Node
{
    private:

        BVH_Node* left{nullptr};
        BVH_Node* right{nullptr};
        //may want to use std::shared_ptr<BVH_Node> but it
        //might not matter since there won't be any tree balancing
        //(at least I don't think there will be ...)

        void setLeft(BVH_Node* lc) {left = lc;}
        void setRight(BVH_Node* rc) {right = rc;}


    public:

        InternalNode() = default;
        ~InternalNode() = default;

        //delete copy and move ops until there
        //is a good reason not to.
        InternalNode(const InternalNode&) = delete;
        InternalNode& operator=(const InternalNode&) = delete;
        InternalNode(InternalNode&&) = delete;
        InternalNode& operator=(InternalNode&&) = delete;

        InternalNode(InternalNode* lc, InternalNode* rc);

        void setParent(BVH_Node* p) {parent = p;}

        const InternalNode* const getLeft() const {return left;}
        const InternalNode* const getRight() const {return right;}
        const InternalNode* const getParent() const {return parent;}
        const InternalNode* const getSibling() const;
        
        const bool isLeaf() const override {return false;}
        const bool isRoot() const override {return false;}
};

//TODO: implement BVH_Root subclass


//has a FT_HSE* and does not have children
class LeafNode : public BVH_Node
{
    private:

        FT_HSE* hse{nullptr};
        InternalNode* parent{nullptr};

        void setBV(BoundingVolume* BV) {this->bv = BV;}

    public:

        LeafNode() = default;
        ~LeafNode() = default;

        //delete copy and move ops until there
        //is a good reason not to.
        LeafNode(const LeafNode&) = delete;
        LeafNode& operator=(const LeafNode&) = delete;
        LeafNode(LeafNode&&) = delete;
        LeafNode& operator=(LeafNode&&) = delete;

        LeafNode(FT_HSE*, const BVH_Type&);

        const FT_HSE* const getHse() const {return hse;}
        const bool isLeaf() const override {return true;}
        const bool isRoot() const override {return false;}
};


#endif
