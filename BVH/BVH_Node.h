#include "BoundingVolume.h"

#ifndef BVH_NODE_H
#define BVH_NODE_H

#include <memory>

//Design Considerations:
//
//  1. Do we really need smart_ptrs for Node linkage?
//     No deletion or rotations are going to be performed,
//     and the objects will persist for the entirety of
//     the program.
//  2. How much of an impact will the use of smart_ptrs
//     have on performance?

class InternalNode;
using BoundingVolume = AABB;


class BVH_Node
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

        void setBV(const BoundingVolume& BV);
        const BoundingVolume& getBV() const;
        
        void setParent(std::shared_ptr<InternalNode> P);
        const std::weak_ptr<const InternalNode> getParent() const;
};


class InternalNode : public BVH_Node,
    public std::enable_shared_from_this<InternalNode>
{
    private:

        std::shared_ptr<BVH_Node> left;
        std::shared_ptr<BVH_Node> right;
        
        void setLeftChild(std::shared_ptr<BVH_Node> lc);
        void setRightChild(std::shared_ptr<BVH_Node> rc);

    public:

        InternalNode(std::shared_ptr<BVH_Node> lc,
                std::shared_ptr<BVH_Node> rc);

        InternalNode() = default;
        virtual ~InternalNode() = default;

        InternalNode(const InternalNode&) = delete;
        InternalNode& operator=(const InternalNode&) = delete;
        InternalNode(InternalNode&&) = delete;
        InternalNode& operator=(InternalNode&&) = delete;
        
        const bool isLeaf() const override;

        //Would be better if this could be made private,
        //or coupled to the constructor. However, this does
        //not seem possible given the use of smart_ptrs for
        //node linkage. See InternalNode constructor in
        //BVH_Node.cpp for details.
        void setChildren(std::shared_ptr<BVH_Node> lc,
                std::shared_ptr<BVH_Node> rc);
       
        const std::weak_ptr<const BVH_Node> getLeftChild() const;
        const std::weak_ptr<const BVH_Node> getRightChild() const;
        
};


class LeafNode : public BVH_Node
{
    private:

        Hse* hse{nullptr};

    public:

        LeafNode(Hse* h);
        LeafNode() = default;
        ~LeafNode() = default;

        LeafNode(const LeafNode&) = delete;
        LeafNode& operator=(const LeafNode&) = delete;
        LeafNode(LeafNode&&) = delete;
        LeafNode& operator=(LeafNode&&) = delete;

        const bool isLeaf() const override;

        const Hse* const getHse() const;
};

#endif


