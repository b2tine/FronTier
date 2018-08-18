#include "gmock/gmock.h"
#include "../BV_Node.h"

class BV_NodeTestData : public ::testing::Test
{
    protected:

    TRI *t1, *t2;
    POINT *a, *b, *c, *d, *e, *f;

    FT_TRI *T1, *T2;

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

//AABB is currently the only BV type
class AABB_NodeTests : public BV_NodeTestData
{
    public:
    
    BV_Leaf *lchild, *rchild;

    AABB_NodeTests()
        : BV_NodeTestData{}
    {
        lchild = new BV_Leaf(T1);
        rchild = new BV_Leaf(T2);
    }

    void TearDown() override
    {
        delete lchild;
        delete rchild;
        BV_NodeTestData::TearDown();
    }

    ~AABB_NodeTests() = default;

};

using DISABLED_AABB_NodeTests = AABB_NodeTests;

TEST_F(AABB_NodeTests, ConstructorBV_NodeDeathTest)
{
    BV_Node* achild;
    ASSERT_DEATH(BV_Node* node = new BV_Node(lchild,achild),"");
}

TEST_F(AABB_NodeTests, ConstructorBV_Node)
{
    BV_Node* parentnode = new BV_Node(lchild,rchild);
    ASSERT_NE(parentnode->getLeft(),nullptr);
    ASSERT_EQ(lchild->getParent(),parentnode);
    ASSERT_EQ(rchild->getParent(),parentnode);
    ASSERT_FALSE(parentnode->isLeaf());
    delete parentnode;
}

TEST_F(AABB_NodeTests, ConstructorBV_Leaf)
{
    ASSERT_NE(lchild->getHse(),nullptr);
    ASSERT_EQ(lchild->getParent(),nullptr);
    ASSERT_TRUE(lchild->isLeaf());
    ASSERT_FALSE(lchild->isRoot());
}



