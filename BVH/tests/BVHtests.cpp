#include "gmock/gmock.h"
#include "../BVH.h"


class BVH_Tests : public testing::Test
{
    protected:

    static TRI *t1, *t2, *t3, *t4,*t5;
    static POINT *a, *b, *c, *d, *e, *f, *g;

    static HsTri *T1, *T2, *T3, *T4, *T5;

    //BVH bvh;
    std::shared_ptr<LeafNode>  l1, l2, l3, l4, l5;
    std::shared_ptr<InternalNode> p1, p2, gp;

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
        Coords(e)[0] = -1.0;    Coords(f)[0] = 0.0;
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
        l1 = BVH::createLeafNode(T1);
        l2 = BVH::createLeafNode(T2);
        l3 = BVH::createLeafNode(T3);
        l4 = BVH::createLeafNode(T4);
        l5 = BVH::createLeafNode(T5);
        p1 = BVH::createInternalNode(l1,l2);
        p2 = BVH::createInternalNode(l3,l4);
        gp = BVH::createInternalNode(p1,p2);
    }

    void TearDown() override
    {
    }

    ~BVH_Tests() = default;
};

TRI* BVH_Tests::t1 = nullptr;
TRI* BVH_Tests::t2 = nullptr;
TRI* BVH_Tests::t3 = nullptr;
TRI* BVH_Tests::t4 = nullptr;
TRI* BVH_Tests::t5 = nullptr;
POINT* BVH_Tests::a = nullptr;
POINT* BVH_Tests::b = nullptr;
POINT* BVH_Tests::c = nullptr;
POINT* BVH_Tests::d = nullptr;
POINT* BVH_Tests::e = nullptr;
POINT* BVH_Tests::f = nullptr;
POINT* BVH_Tests::g = nullptr;
HsTri* BVH_Tests::T1 = nullptr;
HsTri* BVH_Tests::T2 = nullptr;
HsTri* BVH_Tests::T3 = nullptr;
HsTri* BVH_Tests::T4 = nullptr;
HsTri* BVH_Tests::T5 = nullptr;


using DISABLED_BVH_Tests = BVH_Tests;



TEST_F(BVH_Tests, FactoryCreateInternalNode)
{
    ASSERT_NE(p1,nullptr);
    ASSERT_NE(p2,nullptr);
    ASSERT_NE(gp,nullptr);
}

TEST_F(BVH_Tests, FactoryCreateLeafNode)
{
    ASSERT_NE(l1,nullptr);
    ASSERT_NE(l2,nullptr);
    ASSERT_NE(l3,nullptr);
    ASSERT_NE(l4,nullptr);
}




