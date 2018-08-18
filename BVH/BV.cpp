#include "BV.h"

///////////////////////////////////
//////     AABB methods     //////
/////////////////////////////////

AABB::AABB(FT_HSE* h)
{
    for( int i = 0; i < 3; i++ )
    {
        lower[i] = h->min_coord(i);
        upper[i] = h->max_coord(i);
        //lower.push_back(h->min_coord(i));
        //upper.push_back(h->max_coord(i));
        //centroid.push_back(0.5*(lower[i]+upper[i]));
    }
}

AABB::AABB(const BV_Point& L, const BV_Point& U)
    : lower{L}, upper{U}
{
    /*for( int i = 0; i < 3; i++ )
        centroid.push_back(0.5*(lower[i]+upper[i]));*/
}

AABB::AABB(const AABB& A, const AABB& B)
{
    for( int i = 0; i < 3; i++ )
    {
        lower[i] = std::min(A.lower[i],B.lower[i]);
        upper[i] = std::max(A.upper[i],B.upper[i]);
        //lower.push_back(std::min(A.lower[i],B.lower[i]));
        //upper.push_back(std::max(A.upper[i],B.upper[i]));
        //centroid.push_back(0.5*(lower[i]+upper[i]));
    }
}

const BV_Point AABB::centroid() const
{
    BV_Point centroid;
    for( int i = 0; i < 3; i++ )
        centroid[i] = 0.5*(lower[i]+upper[i]);
    return centroid;
}

//bool AABB::contains(AABB* BB) const
bool AABB::contains(const AABB& BB) const
{
    for( int i = 0; i < 3; i++ )
    {
        if( lower[i] >= BB.lower[i] ) return false;
        if( upper[i] <= BB.upper[i] ) return false;
    }
    return true;
}

//bool AABB::overlaps(AABB* BB) const
bool AABB::overlaps(const AABB& BB) const
{
    for( int i = 0; i < 3; i++ )
    {
        if( lower[i] > BB.upper[i] ) return false;
        if( upper[i] < BB.lower[i] ) return false;
    }
    return true;
}

void AABB::print() const
{
    BV_Point centroid = this->centroid();
    printf("\n   upper: (%3g,%3g,%3g) \n", upper[0], upper[1], upper[2]);
    printf("centroid: (%3g,%3g,%3g) \n", centroid[0], centroid[1], centroid[2]);
    printf("   lower: (%3g,%3g,%3g) \n\n", lower[0], lower[1], lower[2]);
}

BV_Type AABB::getTypeBV() const
{
    return BV_Type::AABB;
}


