#include "gmock/gmock.h"
#include "../HyperSurfElement.h"



class HsPoint_Tests : public testing::Test
{
    protected:

    POINT* a;
    HsPoint* A;

    HsPoint_Tests()
        : a{new POINT}
    {
        Coords(a)[0] = 1.0;     //Coords(b)[0] = 1.0;
        Coords(a)[1] = 0.0;     //Coords(b)[0] = 1.0;
        Coords(a)[2] = 4.0;     //Coords(b)[0] = 1.0;
        A = new HsPoint(a);    //B = new HsPoint(b);
    }

    void TearDown() override
    {
        delete a;   delete A;
        //delete b;   delete B;
    }

};

using DISABLED_HsPoint_Tests = HsPoint_Tests;

TEST_F(HsPoint_Tests, MaxCoordIsEqualToMinCoord)
{
    ASSERT_DOUBLE_EQ(A->min_coord(0),A->max_coord(0));
    ASSERT_DOUBLE_EQ(A->min_coord(2),A->max_coord(2));
}

TEST_F(HsPoint_Tests, PointOfHseDefaultArgIsZero)
{
    POINT* b = A->Point_of_hse();
    ASSERT_DOUBLE_EQ(Coords(b)[0],1.0);
    ASSERT_DOUBLE_EQ(Coords(b)[1],0.0);
    ASSERT_DOUBLE_EQ(Coords(b)[2],4.0);
}

TEST_F(HsPoint_Tests, OutOfRangeDeathTest)
{
    ASSERT_DEATH(A->Point_of_hse(1),"");
    ASSERT_DEATH(A->Point_of_hse(-1),"");
}

//////////////////////////////////////////////////////



class HsBond_Tests : public testing::Test
{
    protected:
    
    BOND* s;
    POINT *a, *b;
    HsBond* B;

    HsBond_Tests() 
        : a{new POINT}, b{new POINT}, s{new BOND}
    {
        Coords(a)[0] = 1.0;     Coords(b)[0] = 1.0;
        Coords(a)[1] = 0.0;     Coords(b)[1] = 1.0;
        Coords(a)[2] = 4.0;     Coords(b)[2] = 1.0;
        
        s->start = a;   s->end = b;

        B = new HsBond(s);
    }

    void TearDown() override
    {
        delete a;   delete b;
        delete s;   delete B;
    }

};


TEST_F(HsBond_Tests, OutOfRangeDeathTest)
{
    ASSERT_DEATH(B->Point_of_hse(2),"");
    ASSERT_DEATH(B->Point_of_hse(-1),"");
}



//////////////////////////////////////////////////////


class HsTri_Tests : public testing::Test
{
    protected:
    
    TRI* t;
    POINT *a, *b, *c;
    HsTri* T;

    HsTri_Tests() 
        : a{new POINT}, b{new POINT},
        c{new POINT}, t{new TRI}
    {
        Coords(a)[0] = 1.0;     Coords(b)[0] = 1.0;
        Coords(a)[1] = 0.0;     Coords(b)[1] = 1.0;
        Coords(a)[2] = 4.0;     Coords(b)[2] = 1.0;

        Coords(c)[0] = -1.0;    Point_of_tri(t)[0] = a;
        Coords(c)[1] = 0.0;     Point_of_tri(t)[1] = b;
        Coords(c)[2] = 2.0;     Point_of_tri(t)[2] = c;

        T = new HsTri(t);
    }

    void TearDown() override
    {
        delete a; delete b;
        delete c; delete t;
    }

};

TEST_F(HsTri_Tests, OutOfRangeDeathTest)
{
    ASSERT_DEATH(T->Point_of_hse(3),"");
    ASSERT_DEATH(T->Point_of_hse(-1),"");
}







