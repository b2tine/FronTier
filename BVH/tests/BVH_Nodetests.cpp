#include "gmock/gmock.h"
#include "../BVH_Node.h"

class TestData : public ::testing::Test
{
    public:

    TRI *t1, *t2;
    POINT *a, *b, *c, *d, *e, *f;

    FT_TRI *T1, *T2;

    TestData()
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

    void TearDown() override
    {
        delete t1; delete t2;
        delete T1; delete T2;
        delete a; delete b;
        delete c; delete d;
    }

};

/*
class BVH_NodeTests : public TestData
{
    public:
    
    AABB *bbP, *bbB1, *bbB2, *bbT1, *bbT2;

    BVH_NodeTests()
        : TestData{}
    {
        bbP = new AABB(P);
        bbB1 = new AABB(B1);    bbB2 = new AABB(B2);
        bbT1 = new AABB(T1);    bbT2 = new AABB(T2);
    }

    void TearDown() override
    {
        delete bbP;
        delete bbB1; delete bbB2;
        delete bbT1; delete bbT2;
        TestData::TearDown();
    }

};
*/


TEST(BVH_NodeTest, Constructor_FT_HSE)
{
    BVH_Node node;
} 




