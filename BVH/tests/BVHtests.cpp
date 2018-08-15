#include "gmock/gmock.h"
#include "../BVH.h"

TEST(BVH, BVH_EmptyByDefault)
{
    BVH bvh;
    ASSERT_EQ(bvh.getRoot(),nullptr);
}
