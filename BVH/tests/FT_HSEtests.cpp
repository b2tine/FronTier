#include "gmock/gmock.h"
#include "../FT_HSE.h"

/*
class FT_POINT_TEST : public ::testing::Test
{
    public:

        POINT* a;


};
*/

TEST(FT_POINT, DeathTest)
{
    POINT* a = new POINT;
    Coords(a)[0] = 1.0;
    Coords(a)[1] = 0.0;
    Coords(a)[2] = 4.0;

    FT_POINT p(a);
    ASSERT_DEATH(p.Point_of_hse(1),".*Assertion.*");
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

