#include "gmock/gmock.h"
#include "../BVH_Node.h"

class BVH_NodeTests : public testing::Test
{
    protected:

    static POINT *a, *b, *c, *d, *e, *f, *g;
    static TRI *t1, *t2, *t3, *t4,*t5;
    static HsTri *T1, *T2, *T3, *T4, *T5;

    std::shared_ptr<LeafNode>  l1, l2, l3, l4, l5;

    static void SetUpTestCase()
    {
        a = new POINT;         b = new POINT;
        Coords(a)[0] = 0.0;    Coords(b)[0] = 1.0;
        Coords(a)[1] = 0.0;    Coords(b)[1] = 0.0;
        Coords(a)[2] = 0.0;    Coords(b)[2] = 0.0;

        c = new POINT;          d = new POINT;
        Coords(c)[0] = 0.0;     Coords(d)[0] = 0.0;
        Coords(c)[1] = 1.0;     Coords(d)[1] = 0.0;
        Coords(c)[2] = 0.0;     Coords(d)[2] = 1.0;

        e = new POINT;          f = new POINT;
        Coords(e)[0] = -1.0;     Coords(f)[0] = 0.0;
        Coords(e)[1] = 0.0;     Coords(f)[1] = -1.0;
        Coords(e)[2] = 0.0;     Coords(f)[2] = 0.0;

        g = new POINT;          //f = new POINT;
        Coords(g)[0] = 0.0;     //Coords(f)[0] = 1.0;
        Coords(g)[1] = 0.0;     //Coords(f)[1] = 1.0;
        Coords(g)[2] = -1.0;    //Coords(f)[2] = 0.0;

        t1 = new TRI;                t2 = new TRI;
        Point_of_tri(t1)[0] = a;     Point_of_tri(t2)[0] = b;
        Point_of_tri(t1)[1] = b;     Point_of_tri(t2)[1] = c;
        Point_of_tri(t1)[2] = c;     Point_of_tri(t2)[2] = d;

        t3 = new TRI;                t4 = new TRI;
        Point_of_tri(t3)[0] = c;     Point_of_tri(t4)[0] = d;
        Point_of_tri(t3)[1] = d;     Point_of_tri(t4)[1] = e;
        Point_of_tri(t3)[2] = e;     Point_of_tri(t4)[2] = f;

        t5 = new TRI;                //t6 = new TRI;
        Point_of_tri(t5)[0] = e;     //Point_of_tri(t6)[0] = f;
        Point_of_tri(t5)[1] = f;     //Point_of_tri(t6)[1] = g;
        Point_of_tri(t5)[2] = g;     //Point_of_tri(t6)[2] = h;

        T1 = new HsTri(t1);    T2 = new HsTri(t2);
        T3 = new HsTri(t3);    T4 = new HsTri(t4);
        T5 = new HsTri(t5);    //T6 = new HsTri(t6);
    }

    static void TearDownTestCase()
    {
        delete t1; delete t2; delete t3;
        delete t4; delete t5;
        delete T1; delete T2; delete T3;
        delete T4; delete T5;
        delete a; delete b; delete c;
        delete d; delete e; delete f;
        delete g;
    }

    void SetUp() override
    {
        l1 = std::make_shared<LeafNode>(T1);
        l2 = std::make_shared<LeafNode>(T2);
        l3 = std::make_shared<LeafNode>(T3);
        l4 = std::make_shared<LeafNode>(T4);
        l5 = std::make_shared<LeafNode>(T5);
    }

    void TearDown() override
    {

    }

    ~BVH_NodeTests() = default;
};

TRI* BVH_NodeTests::t1 = nullptr;
TRI* BVH_NodeTests::t2 = nullptr;
TRI* BVH_NodeTests::t3 = nullptr;
TRI* BVH_NodeTests::t4 = nullptr;
TRI* BVH_NodeTests::t5 = nullptr;
POINT* BVH_NodeTests::a = nullptr;
POINT* BVH_NodeTests::b = nullptr;
POINT* BVH_NodeTests::c = nullptr;
POINT* BVH_NodeTests::d = nullptr;
POINT* BVH_NodeTests::e = nullptr;
POINT* BVH_NodeTests::f = nullptr;
POINT* BVH_NodeTests::g = nullptr;
HsTri* BVH_NodeTests::T1 = nullptr;
HsTri* BVH_NodeTests::T2 = nullptr;
HsTri* BVH_NodeTests::T3 = nullptr;
HsTri* BVH_NodeTests::T4 = nullptr;
HsTri* BVH_NodeTests::T5 = nullptr;



using DISABLED_BVH_NodeTests = BVH_NodeTests;



TEST_F(DISABLED_BVH_NodeTests, GetSibling)
{
    //auto s = l1->getSibling();
    //ASSERT_EQ(s,l2);
}

TEST_F(BVH_NodeTests, ConstructorInternalNodeDeathTest)
{
    std::shared_ptr<LeafNode> l6;
    std::shared_ptr<InternalNode> p;
    ASSERT_DEATH(p = std::make_shared<InternalNode>(l1,l6),"");
}
      
TEST_F(BVH_NodeTests, InternalNodePrototypeFactoryTest)
{
    //See ../BVH_Node.cpp for why construction of a
    //shared_ptr<InternalNode> must be decoupled from
    //the linking of parent and children.
    std::shared_ptr<InternalNode> p1 =
        std::make_shared<InternalNode>(l1,l2);
    p1->setChildren(l1,l2);
    
    std::shared_ptr<InternalNode> p2 =
        std::make_shared<InternalNode>(l3,l4);
    p2->setChildren(l3,l4);

    std::shared_ptr<InternalNode> gp =
        std::make_shared<InternalNode>(p1,p2);
    gp->setChildren(p1,p2);
    
    ASSERT_FALSE(p1->isLeaf());
    ASSERT_EQ(p1->getLeftChild().lock(),l1);
    ASSERT_EQ(p1->getRightChild().lock(),l2);
    ASSERT_EQ(l1->getParent().lock(),p1);
    ASSERT_EQ(l2->getParent().lock(),p1);

    BoundingVolume bvp1 = p1->getBV();
    ASSERT_DOUBLE_EQ(bvp1.lower[0],0.0);
    ASSERT_DOUBLE_EQ(bvp1.upper[0],1.0);

    BoundingVolume bvp2 = p2->getBV();
    ASSERT_DOUBLE_EQ(bvp2.lower[1],-1.0);
    ASSERT_DOUBLE_EQ(bvp2.upper[1],1.0);

    ASSERT_EQ(gp->getLeftChild().lock(),p1);
    ASSERT_EQ(gp->getRightChild().lock(),p2);
    ASSERT_EQ(p1->getParent().lock(),gp);
    ASSERT_EQ(p2->getParent().lock(),gp);

    BoundingVolume bvgp = gp->getBV();
    ASSERT_DOUBLE_EQ(bvgp.lower[0],-1.0);
    ASSERT_DOUBLE_EQ(bvgp.upper[2],1.0);
}

TEST_F(BVH_NodeTests, ConstructorLeafNode)
{
    ASSERT_TRUE(l5->isLeaf());
    ASSERT_NE(l5->getHse(),nullptr);
    ASSERT_EQ(l5->getParent().lock(),nullptr);

    BoundingVolume bv5 = l5->getBV();
    ASSERT_DOUBLE_EQ(bv5.upper[0],0.0);
    ASSERT_DOUBLE_EQ(bv5.lower[1],-1.0);
}



