#include "BoundingVolume.h"

#ifndef BVH_NODE_H
#define BVH_NODE_H

#include <memory>

class InternalNode;
using BoundingVolume = AABB;


class BVH_Node : public std::enable_shared_from_this<BVH_Node>
{
    private:

        BoundingVolume bv;
        std::weak_ptr<InternalNode> parent;
        
    public:

        BVH_Node() = default;
        virtual ~BVH_Node() = default;

        BVH_Node(const BVH_Node&) = delete;
        BVH_Node& operator=(const BVH_Node&) = delete;
        BVH_Node(BVH_Node&&) = delete;
        BVH_Node& operator=(BVH_Node&&) = delete;

        virtual const bool isLeaf() const = 0;
        virtual const bool isRoot() const = 0;

        void setBV(const BoundingVolume& BV) {bv = BV;}
        const BoundingVolume& getBV() const {return bv;}
        
        void setParent(std::shared_ptr<InternalNode> P)
        {
            parent = P;
        }
        
        const std::weak_ptr<const InternalNode> getParent() const
        {
            return parent;
        }
};


//has a Hse* and does not have children
class LeafNode : public BVH_Node
{
    private:

        Hse* hse{nullptr};

    public:

        LeafNode(Hse* h)
            : hse{h}
        {
            setBV(BoundingVolume(h));
        }

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

        std::shared_ptr<BVH_Node> left;
        std::shared_ptr<BVH_Node> right;
        void setLeft(std::shared_ptr<BVH_Node> lc) {left = lc;}
        void setRight(std::shared_ptr<BVH_Node> rc) {right = rc;}

    public:

        InternalNode(InternalNode* lc, InternalNode* rc);

        InternalNode() = default;
        virtual ~InternalNode() = default;

        InternalNode(const InternalNode&) = delete;
        InternalNode& operator=(const InternalNode&) = delete;
        InternalNode(InternalNode&&) = delete;
        InternalNode& operator=(InternalNode&&) = delete;


        const std::weak_ptr<const BVH_Node> getLeft() const
        {
            return std::weak_ptr<BVH_Node>(left);
        }

        const std::weak_ptr<const BVH_Node> getRight() const
        {
            return std::weak_ptr<BVH_Node>(right);
        }
        
        const bool isLeaf() const override {return false;}
        const bool isRoot() const override {return false;}
        //TODO: need to check if parent is nullptr for root check
};



#endif


