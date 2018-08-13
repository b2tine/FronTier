#include "gmock/gmock.h"
#include "../FT_HSE.h"

/*
class FT_POINT_TEST : public ::testing::Test
{
    public:

        POINT* a;
};
*/

TEST(FT_POINT, OutOfRangeDeathTest)
{
    POINT* a = new POINT;
    Coords(a)[0] = 1.0;
    Coords(a)[1] = 0.0;
    Coords(a)[2] = 4.0;

    FT_POINT p(a);
    ASSERT_DEATH(p.Point_of_hse(1),".*Assertion.*");
    ASSERT_DEATH(p.Point_of_hse(-1),".*Assertion.*");
    delete a;
}

TEST(FT_POINT, PointOfHse)
{
    POINT* a = new POINT;
    Coords(a)[0] = 1.0;
    Coords(a)[1] = 0.0;
    Coords(a)[2] = 4.0;

    FT_POINT p(a);
    POINT* b = p.Point_of_hse(0);

    ASSERT_DOUBLE_EQ(Coords(b)[0],1.0);
    ASSERT_DOUBLE_EQ(Coords(b)[1],0.0);
    ASSERT_DOUBLE_EQ(Coords(b)[2],4.0);
    delete a;
}

TEST(FT_POINT, MaxCoordEqualMinCoord)
{
    POINT* a = new POINT;
    Coords(a)[0] = 1.0;
    Coords(a)[1] = 0.0;
    Coords(a)[2] = 4.0;

    FT_POINT p(a);

    ASSERT_DOUBLE_EQ(p.min_coord(0),p.max_coord(0));
    ASSERT_DOUBLE_EQ(p.min_coord(2),4.0);
    delete a;
}

//////////////////////////////////////////////////////

TEST(FT_BOND, OutOfRangeDeathTest)
{
    POINT* a = new POINT;   POINT* b = new POINT;
    Coords(a)[0] = 1.0;     Coords(b)[0] = 1.0;
    Coords(a)[1] = 0.0;     Coords(b)[1] = 1.0;
    Coords(a)[2] = 4.0;     Coords(b)[2] = 1.0;

    BOND* A = new BOND;
    A->start = a;   A->end = b;

    FT_BOND B(A);
    ASSERT_DEATH(B.Point_of_hse(2),".*Assertion.*");
    ASSERT_DEATH(B.Point_of_hse(-1),".*Assertion.*");
    delete a; delete b; delete A;
}

//////////////////////////////////////////////////////

TEST(FT_TRI, OutOfRangeDeathTest)
{
    POINT* a = new POINT;   POINT* b = new POINT;
    Coords(a)[0] = 1.0;     Coords(b)[0] = 1.0;
    Coords(a)[1] = 0.0;     Coords(b)[1] = 1.0;
    Coords(a)[2] = 4.0;     Coords(b)[2] = 1.0;

    POINT* c = new POINT;   TRI* t = new TRI;
    Coords(c)[0] = -1.0;    Point_of_tri(t)[0] = a;
    Coords(c)[1] = 0.0;     Point_of_tri(t)[1] = b;
    Coords(c)[2] = 2.0;     Point_of_tri(t)[2] = c;

    FT_TRI T(t);
    ASSERT_DEATH(T.Point_of_hse(3),".*Assertion.*");
    ASSERT_DEATH(T.Point_of_hse(-1),".*Assertion.*");
    delete a; delete b; delete c; delete t;
}

//////////////////////////////////////////////////////






