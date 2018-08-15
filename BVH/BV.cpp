#include "BV.h"

///////////////////////////////////
//////     AABB methods     //////
/////////////////////////////////

AABB::AABB(FT_HSE* h)
{
    for( int i = 0; i < 3; i++ )
    {
        lower.push_back(h->min_coord(i));
        upper.push_back(h->max_coord(i));
        centroid.push_back(0.5*(lower[i]+upper[i]));
    }
}

AABB::AABB(const BV_Point& L, const BV_Point& U)
    : lower{L}, upper{U}
{
    for( int i = 0; i < 3; i++ )
        centroid.push_back(0.5*(lower[i]+upper[i]));
}

BV_Point AABB::getCentroid() const
{
    return centroid;
}

bool AABB::contains(AABB* BB) const
{
    for( int i = 0; i < 3; i++ )
    {
        if( lower[i] >= BB->lower[i] ) return false;
        if( upper[i] <= BB->upper[i] ) return false;
    }
    return true;
}

bool AABB::overlaps(AABB* BB) const
{
    for( int i = 0; i < 3; i++ )
    {
        if( lower[i] > BB->upper[i] ) return false;
        if( upper[i] < BB->lower[i] ) return false;
    }
    return true;
}

void AABB::print() const
{
    printf("\n   upper: (%3g,%3g,%3g) \n", upper[0], upper[1], upper[2]);
    printf("centroid: (%3g,%3g,%3g) \n", centroid[0], centroid[1], centroid[2]);
    printf("   lower: (%3g,%3g,%3g) \n\n", lower[0], lower[1], lower[2]);
}

BV_Type AABB::getTypeBV() const
{
    return BV_Type::AABB;
}

//static factory like function
AABB* AABB::merge(AABB* A, AABB* B)
{
    BV_Point L;
    BV_Point U;

    for( int i = 0; i < 3; i++ )
    {
        L.push_back(std::min(A->lower[i],B->lower[i]));
        U.push_back(std::max(A->upper[i],B->upper[i]));
    }

    return new AABB(L,U);
}

