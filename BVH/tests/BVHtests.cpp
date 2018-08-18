#include "gmock/gmock.h"
#include "../BVH.h"


class BVH_TestData : public ::testing::Test
{
    protected:

    TRI *t1, *t2, *t3, *t4;
    POINT *a, *b, *c, *d, *e, *f;

    FT_TRI *T1, *T2, *T3, *T4;

    //can probably make this a static StartUp/TearDown method
    //in the derived class fixture to share the resources.
    //Do it when compilation/run time slowdown observed..
    BVH_TestData()
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

        t1 = new TRI;                t2 = new TRI;
        Point_of_tri(t1)[0] = a;     Point_of_tri(t2)[0] = b;
        Point_of_tri(t1)[1] = b;     Point_of_tri(t2)[1] = c;
        Point_of_tri(t1)[2] = c;     Point_of_tri(t2)[2] = d;

        t3 = new TRI;                t4 = new TRI;
        Point_of_tri(t3)[0] = c;     Point_of_tri(t4)[0] = d;
        Point_of_tri(t3)[1] = d;     Point_of_tri(t4)[1] = e;
        Point_of_tri(t3)[2] = e;     Point_of_tri(t4)[2] = f;

        T1 = new FT_TRI(t1);    T2 = new FT_TRI(t2);
        T3 = new FT_TRI(t3);    T4 = new FT_TRI(t4);
    }

    virtual void TearDown()
    {
        delete t1; delete t2; delete t3; delete t4;
        delete T1; delete T2; delete T3; delete T4;
        delete a; delete b; delete c;
        delete d; delete e; delete f;
    }

    virtual ~BVH_TestData() = default;

};

class BVH_Tests : public BVH_TestData
{
    public:

    BVH bvh;
    BV_Leaf  *l1, *l2, *l3, *l4;
    BV_Node *p1, *p2, *gp;

    BVH_Tests()
    {
        l1 = BVH::createLeaf(T1);
        l2 = BVH::createLeaf(T2);
        l3 = BVH::createLeaf(T3);
        l4 = BVH::createLeaf(T4);
        p1 = BVH::createNode(l1,l2);
        p2 = BVH::createNode(l3,l4);
        gp = BVH::createNode(p1,p2);
    }
    
    void TearDown() override
    {
        delete l1; delete l2;
        delete l3; delete l4;
        delete p1; delete p2; delete gp;
        BVH_TestData::TearDown();
    }

    ~BVH_Tests() = default;

};


using DISABLED_BVH_Tests = BVH_Tests;


//TODO: Add factory death tests


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

TEST_F(BVH_Tests, RootNullByDefault)
{
    ASSERT_EQ(bvh.getRoot(),nullptr);
}



