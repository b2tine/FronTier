#include "gmock/gmock.h"
#include "../BV_Point.h"

class BV_Point_Tests :  public ::testing::Test
{
    protected:

        BV_Point defaultPoint;
        BV_Point xyzTriplePoint;

        BV_Point p1, p2, p3, p4;

    public:

        BV_Point_Tests()
            : defaultPoint(),
            xyzTriplePoint(1,2,3),
            p1(0,0,0), p2(1,2,3),
            p3(0,2,3), p4(1,2,0)
        {}

        ~BV_Point_Tests() = default;

};

using DISABLED_BV_Point_Tests = BV_Point_Tests;


TEST_F(BV_Point_Tests, ZeroVectorByDefault)
{
    ASSERT_DOUBLE_EQ(defaultPoint[0],0.0);
    ASSERT_DOUBLE_EQ(defaultPoint[1],0.0);
    ASSERT_DOUBLE_EQ(defaultPoint[2],0.0);
}

TEST_F(BV_Point_Tests, xyzTripleConstructor)
{
    ASSERT_DOUBLE_EQ(xyzTriplePoint[0],1.0);
    ASSERT_DOUBLE_EQ(xyzTriplePoint[1],2.0);
    ASSERT_DOUBLE_EQ(xyzTriplePoint[2],3.0);
}

TEST_F(BV_Point_Tests, OutOfRangeDeathTest)
{
    ASSERT_DEATH(xyzTriplePoint[-1],"");
    ASSERT_DEATH(xyzTriplePoint[3],"");
}

TEST_F(BV_Point_Tests, xCoordsNotEqualLessThan)
{
    ASSERT_TRUE(p1 < p2);
}

TEST_F(BV_Point_Tests, xCoordsEqualLessThan)
{
    ASSERT_TRUE(p1 < p3);
}

TEST_F(BV_Point_Tests, xyCoordsBothEqualLessThan)
{
    ASSERT_TRUE(p4 < p2);
}

TEST_F(BV_Point_Tests, PointsAreEqualLessThan)
{
    ASSERT_FALSE(defaultPoint < p1);
}

TEST_F(BV_Point_Tests, IsValidKeyForStdMap)
{
    std::map<BV_Point,int> imap;
    ASSERT_TRUE(imap.empty());

    imap[p1] = 1;
    ASSERT_FALSE(imap.empty());
    ASSERT_EQ(imap[p1],1);
}

//TODO: Start test driving HilbertSortingTraits


