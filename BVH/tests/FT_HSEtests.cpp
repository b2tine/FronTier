#include "gmock/gmock.h"
#include "../FT_HSE.h"

class FT_POINT_TESTS : public ::testing::Test
{
    protected:

    POINT* a;
    FT_POINT* p;

    FT_POINT_TESTS()
        : a{new POINT}
    {
        Coords(a)[0] = 1.0;
        Coords(a)[1] = 0.0;
        Coords(a)[2] = 4.0;
        p = new FT_POINT(a);
    }

    void TearDown() override
    {
        delete a;
        delete p;
    }

};


TEST_F(FT_POINT_TESTS, PointOfHse)
{
    POINT* b = p->Point_of_hse(0);

    ASSERT_DOUBLE_EQ(Coords(b)[0],1.0);
    ASSERT_DOUBLE_EQ(Coords(b)[1],0.0);
    ASSERT_DOUBLE_EQ(Coords(b)[2],4.0);
}

TEST_F(FT_POINT_TESTS, MaxCoordEqualMinCoord)
{
    ASSERT_DOUBLE_EQ(p->min_coord(0),p->max_coord(0));
    ASSERT_DOUBLE_EQ(p->min_coord(2),4.0);
}

TEST_F(FT_POINT_TESTS, OutOfRangeDeathTest)
{
    ASSERT_DEATH(p->Point_of_hse(1),".*Assertion.*");
    ASSERT_DEATH(p->Point_of_hse(-1),".*Assertion.*");
}

//////////////////////////////////////////////////////



class FT_BOND_TESTS : public ::testing::Test
{
    protected:
    
    BOND* s;
    POINT *a, *b;
    FT_BOND* B;

    FT_BOND_TESTS() 
        : a{new POINT}, b{new POINT}, s{new BOND}
    {
        Coords(a)[0] = 1.0;     Coords(b)[0] = 1.0;
        Coords(a)[1] = 0.0;     Coords(b)[1] = 1.0;
        Coords(a)[2] = 4.0;     Coords(b)[2] = 1.0;
        
        s->start = a;   s->end = b;

        B = new FT_BOND(s);
    }

    void TearDown() override
    {
        delete a;   delete b;
        delete s;   delete B;
    }

};


TEST_F(FT_BOND_TESTS, OutOfRangeDeathTest)
{
    ASSERT_DEATH(B->Point_of_hse(2),".*Assertion.*");
    ASSERT_DEATH(B->Point_of_hse(-1),".*Assertion.*");
}



//////////////////////////////////////////////////////


class FT_TRI_TESTS : public ::testing::Test
{
    protected:
    
    TRI* t;
    POINT *a, *b, *c;
    FT_TRI* T;

    FT_TRI_TESTS() 
        : a{new POINT}, b{new POINT},
        c{new POINT}, t{new TRI}
    {
        Coords(a)[0] = 1.0;     Coords(b)[0] = 1.0;
        Coords(a)[1] = 0.0;     Coords(b)[1] = 1.0;
        Coords(a)[2] = 4.0;     Coords(b)[2] = 1.0;

        Coords(c)[0] = -1.0;    Point_of_tri(t)[0] = a;
        Coords(c)[1] = 0.0;     Point_of_tri(t)[1] = b;
        Coords(c)[2] = 2.0;     Point_of_tri(t)[2] = c;

        T = new FT_TRI(t);
    }

    void TearDown() override
    {
        delete a; delete b;
        delete c; delete t;
    }

};

TEST_F(FT_TRI_TESTS, OutOfRangeDeathTest)
{
    ASSERT_DEATH(T->Point_of_hse(3),".*Assertion.*");
    ASSERT_DEATH(T->Point_of_hse(-1),".*Assertion.*");
}







