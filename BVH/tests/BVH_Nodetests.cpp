#include "gmock/gmock.h"
#include "../BVH_Node.h"

class BVH_NodeTestData : public ::testing::Test
{
    public:

    TRI *t1, *t2;
    POINT *a, *b, *c, *d, *e, *f;

    FT_TRI *T1, *T2;

    BVH_NodeTestData()
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

    virtual ~BVH_NodeTestData() = default;

};

class BVH_NodeTests : public BVH_NodeTestData
{
    public:
    
    BVH_Node* node;

    BVH_NodeTests()
        : BVH_NodeTestData{}
    {
        node = new BVH_Node(T1);
    }

    void TearDown() override
    {
        delete node;
        BVH_NodeTestData::TearDown();
    }

    ~BVH_NodeTests() = default;

};

TEST_F(BVH_NodeTests, Constructor_FT_HSE)
{
    ASSERT_NE(node->getHse(),nullptr);
    ASSERT_NE(node->getBV(),nullptr);//not passed
}



