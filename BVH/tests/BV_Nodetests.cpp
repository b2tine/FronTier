#include "gmock/gmock.h"
#include "../BV_Node.h"

class BV_NodeTestData : public ::testing::Test
{
    public:

    TRI *t1, *t2;
    POINT *a, *b, *c, *d, *e, *f;

    FT_TRI *T1, *T2;

    //can probably make this a static StartUp/TearDown method
    //in the derived class fixture so share the resources.
    BV_NodeTestData()
    {
        a = new POINT;         b = new POINT;
        Coords(a)[0] = 0.0;    Coords(b)[0] = 1.0;
        Coords(a)[1] = 0.0;    Coords(b)[1] = 0.0;
        Coords(a)[2] = 0.0;    Coords(b)[2] = 0.0;

        c = new POINT;          d = new POINT;
        Coords(c)[0] = 0.0;     Coords(d)[0] = 0.0;
        Coords(c)[1] = 1.0;     Coords(d)[1] = 0.0;
        Coords(c)[2] = 0.0;     Coords(d)[2] = 1.0;

        t1 = new TRI;                t2 = new TRI;
        Point_of_tri(t1)[0] = a;     Point_of_tri(t2)[0] = b;
        Point_of_tri(t1)[1] = b;     Point_of_tri(t2)[1] = c;
        Point_of_tri(t1)[2] = c;     Point_of_tri(t2)[2] = d;

        T1 = new FT_TRI(t1);    T2 = new FT_TRI(t2);
    }

    virtual void TearDown()
    {
        delete t1; delete t2;
        delete T1; delete T2;
        delete a; delete b;
        delete c; delete d;
    }

    virtual ~BV_NodeTestData() = default;

};

class BV_NodeTests : public BV_NodeTestData
{
    public:
    
    BV_Leaf<AABB>* leaf;

    BV_NodeTests()
        : BV_NodeTestData{}
    {
        leaf = new BV_Leaf<AABB>(T1);
    }

    void TearDown() override
    {
        delete leaf;
        BV_NodeTestData::TearDown();
    }

    ~BV_NodeTests() = default;

};

//Parameterized by AABB
TEST_F(BV_NodeTests, Constructor_BV_Leaf)
{
    ASSERT_NE(leaf->getHse(),nullptr);
    ASSERT_NE(leaf->getBV(),nullptr);
    ASSERT_TRUE(leaf->isLeaf());
    ASSERT_FALSE(leaf->isRoot());
}



