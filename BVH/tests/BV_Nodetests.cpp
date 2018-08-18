#include "gmock/gmock.h"
#include "../BV_Node.h"

class BV_NodeTestData : public ::testing::Test
{
    protected:

    TRI *t1, *t2, *t3, *t4,*t5;
    POINT *a, *b, *c, *d, *e, *f, *g;

    FT_TRI *T1, *T2, *T3, *T4, *T5;

    //can probably make this a static StartUp/TearDown method
    //in the derived class fixture to share the resources.
    //Do it when compilation/run time slowdown observed..
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

        e = new POINT;          f = new POINT;
        Coords(e)[0] = -1.0;     Coords(f)[0] = 0.0;
        Coords(e)[1] = 0.0;     Coords(f)[1] = -1.0;
        Coords(e)[2] = 0.0;     Coords(f)[2] = 0.0;

        g = new POINT;          //f = new POINT;
        Coords(g)[0] = 0.0;     //Coords(f)[0] = 1.0;
        Coords(g)[1] = 0.0;     //Coords(f)[1] = 1.0;
        Coords(g)[2] = -1.0;    //Coords(f)[2] = 0.0;

        t1 = new TRI;                t2 = new TRI;
        Point_of_tri(t1)[0] = a;     Point_of_tri(t2)[0] = b;
        Point_of_tri(t1)[1] = b;     Point_of_tri(t2)[1] = c;
        Point_of_tri(t1)[2] = c;     Point_of_tri(t2)[2] = d;

        t3 = new TRI;                t4 = new TRI;
        Point_of_tri(t3)[0] = c;     Point_of_tri(t4)[0] = d;
        Point_of_tri(t3)[1] = d;     Point_of_tri(t4)[1] = e;
        Point_of_tri(t3)[2] = e;     Point_of_tri(t4)[2] = f;

        t5 = new TRI;                //t6 = new TRI;
        Point_of_tri(t5)[0] = e;     //Point_of_tri(t6)[0] = f;
        Point_of_tri(t5)[1] = f;     //Point_of_tri(t6)[1] = g;
        Point_of_tri(t5)[2] = g;     //Point_of_tri(t6)[2] = h;

        T1 = new FT_TRI(t1);    T2 = new FT_TRI(t2);
        T3 = new FT_TRI(t3);    T4 = new FT_TRI(t4);
        T5 = new FT_TRI(t5);    //T6 = new FT_TRI(t6);
    }

    virtual void TearDown()
    {
        delete t1; delete t2; delete t3;
        delete t4; delete t5;
        delete T1; delete T2; delete T3;
        delete T4; delete T5;
        delete a; delete b; delete c;
        delete d; delete e; delete f;
        delete g;
    }

    virtual ~BV_NodeTestData() = default;

};

//AABB is currently the only BV type
class AABB_NodeTests : public BV_NodeTestData
{
    public:
    
    BV_Leaf  *l1, *l2, *l3, *l4, *l5;
    BV_iNode *p1, *p2, *gp;

    AABB_NodeTests()
        : BV_NodeTestData{}
    {
        l1 = new BV_Leaf(T1);
        l2 = new BV_Leaf(T2);
        l3 = new BV_Leaf(T3);
        l4 = new BV_Leaf(T4);
        l5 = new BV_Leaf(T5);
        p1 = new BV_iNode(l1,l2);
        p2 = new BV_iNode(l3,l4);
        gp = new BV_iNode(p1,p2);
    }

    void TearDown() override
    {
        delete l1; delete l2;
        delete l3; delete l4; delete l5;
        delete p1; delete p2; delete gp;
        BV_NodeTestData::TearDown();
    }

    ~AABB_NodeTests() = default;

};

using DISABLED_AABB_NodeTests = AABB_NodeTests;


TEST_F(DISABLED_AABB_NodeTests, GetSibling)
{
    //BV_Node* s = l1->getSibling();
}

TEST_F(AABB_NodeTests, ConstructorBV_iNodeDeathTest)
{
    BV_Node* achild;
    ASSERT_DEATH(BV_iNode* node = new BV_iNode(l1,achild),"");
}

TEST_F(AABB_NodeTests, ConstructorBV_Node)
{
    BV_iNode* p1 = new BV_iNode(l1,l2);
    ASSERT_NE(p1->getLeft(),nullptr);
    ASSERT_EQ(l1->getParent(),p1);
    ASSERT_EQ(l2->getParent(),p1);
    ASSERT_FALSE(p1->isLeaf());
    delete p1;
}

TEST_F(AABB_NodeTests, ConstructorBV_Leaf)
{
    ASSERT_NE(l5->getHse(),nullptr);
    ASSERT_EQ(l5->getParent(),nullptr);
    ASSERT_TRUE(l5->isLeaf());
    ASSERT_FALSE(l5->isRoot());
}



