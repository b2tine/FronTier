#include "gmock/gmock.h"
#include "../BoundingVolume.h"


class AABBTests : public testing::Test
{
    protected:

    static TRI *t1, *t2, *t3;
    static BOND *s1, *s2;
    static POINT *a, *b, *c, *d, *e, *f, *g;

    static HsTri *T1, *T2, *T3;
    static HsBond *B1, *B2;
    static HsPoint* P;

    static BV_Point L;
    static BV_Point U;

    AABB bbT1, bbT2, bbT3, bbB1, bbB2, bbBVPts;

    static void SetUpTestCase()
    {
        a = new POINT;         b = new POINT;
        Coords(a)[0] = 2.0;    Coords(b)[0] = 1.0;
        Coords(a)[1] = 2.0;    Coords(b)[1] = 1.0;
        Coords(a)[2] = 3.0;    Coords(b)[2] = 1.0;

        c = new POINT;          d = new POINT;
        Coords(c)[0] = 9.0;     Coords(d)[0] = 10.0;
        Coords(c)[1] = 0.5;     Coords(d)[1] = 0.0;
        Coords(c)[2] = 2.0;     Coords(d)[2] = 0.0;

        e = new POINT;          f = new POINT;
        Coords(e)[0] = 0.0;     Coords(f)[0] = 0.0;
        Coords(e)[1] = 10.0;    Coords(f)[1] = 0.0;
        Coords(e)[2] = 0.0;     Coords(f)[2] = 10.0;

        g = new POINT;          //h = new POINT;
        Coords(g)[0] = 0.0;     //Coords(h)[0] = 0.0;
        Coords(g)[1] = 0.0;     //Coords(h)[1] = 0.0;
        Coords(g)[2] = 0.0;     //Coords(h)[2] = 10.0;

        t1 = new TRI;                t2 = new TRI;
        Point_of_tri(t1)[0] = a;     Point_of_tri(t2)[0] = d;
        Point_of_tri(t1)[1] = b;     Point_of_tri(t2)[1] = e;
        Point_of_tri(t1)[2] = c;     Point_of_tri(t2)[2] = f;

        t3 = new TRI;               //t2 = new TRI;
        Point_of_tri(t3)[0] = e;    //Point_of_tri(t2)[0] = d;
        Point_of_tri(t3)[1] = f;    //Point_of_tri(t2)[1] = e;
        Point_of_tri(t3)[2] = g;    //Point_of_tri(t2)[2] = f;

        s1 = new BOND;      s2 = new BOND;
        s1->start = a;      s2->start = c;
        s1->end = b;        s2->end = f;

        T1 = new HsTri(t1);    T2 = new HsTri(t2);
        T3 = new HsTri(t3);

        B1 = new HsBond(s1);   B2 = new HsBond(s2);
        P = new HsPoint(a);
    }

    static void TearDownTestCase()
    {
        delete a; delete b; delete c;
        delete d; delete e; delete f;
        delete g; delete s1; delete s2;
        delete t1; delete t2; delete t3;
        delete P; delete B1; delete B2;
        delete T1; delete T2; delete T3;
    }

    void SetUp() override
    {
        bbT1 = AABB(T1);
        bbT2 = AABB(T2);
        bbT3 = AABB(T3);
        bbB1 = AABB(B1);
        bbB2 = AABB(B2);
        bbBVPts = AABB(L,U);
    }
  
    ~AABBTests() = default;

};

TRI* AABBTests::t1 = nullptr;
TRI* AABBTests::t2 = nullptr;
TRI* AABBTests::t3 = nullptr;
BOND* AABBTests::s1 = nullptr;
BOND* AABBTests::s2 = nullptr;
POINT* AABBTests::a = nullptr;
POINT* AABBTests::b = nullptr;
POINT* AABBTests::c = nullptr;
POINT* AABBTests::d = nullptr;
POINT* AABBTests::e = nullptr;
POINT* AABBTests::f = nullptr;
POINT* AABBTests::g = nullptr;
HsTri* AABBTests::T1 = nullptr;
HsTri* AABBTests::T2 = nullptr;
HsTri* AABBTests::T3 = nullptr;
HsBond* AABBTests::B1 = nullptr;
HsBond* AABBTests::B2 = nullptr;
HsPoint* AABBTests::P = nullptr;
BV_Point AABBTests::L = {-1,-1,-1};
BV_Point AABBTests::U = {1,1,1};


using DISABLED_AABBTests = AABBTests;


TEST_F(AABBTests, BoxesOverlapVsContain)
{
    //contains means strictly contained
    EXPECT_TRUE(bbT2.contains(bbT1));
    EXPECT_FALSE(bbT1.contains(bbB1));

    //shared surface is considered overlap
    EXPECT_TRUE(bbT1.overlaps(bbB1));
    EXPECT_TRUE(bbT1.overlaps(bbB2));
    //TODO: Is this the behavior we really want?
    //      May need to distinguish between:
    //          1. No overlap but share surface.
    //          2. Contained but share a surface
    //
    //      Too early to make a call.
}

TEST_F(AABBTests, ConstructorTwoAABBs)
{
    AABB parentbox(bbT2,bbT3);

    BV_Point lower = parentbox.lower;
    BV_Point upper = parentbox.upper;

    ASSERT_DOUBLE_EQ(lower[2],0.0);
    ASSERT_DOUBLE_EQ(upper[0],10.0);
}

TEST_F(AABBTests, ConstructorTwoBV_Points)
{
    BV_Point centroid = bbBVPts.getCentroid();
    ASSERT_DOUBLE_EQ(centroid[0],0.0);
    ASSERT_DOUBLE_EQ(centroid[1],0.0);
    ASSERT_DOUBLE_EQ(centroid[2],0.0);
}

TEST_F(AABBTests, ConstructorOneHse)
{
    ASSERT_DOUBLE_EQ(bbT1.upper[0],9.0);
    ASSERT_DOUBLE_EQ(bbT1.lower[2],1.0);
}



