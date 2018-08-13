#include "FT_HSE.h"

#ifndef BV_H
#define BV_H

#include <vector>

using BV_Point = std::vector<double>;


class AABB
{
    public:

        BV_Point lower;
        BV_Point upper;
        BV_Point centroid;

        explicit AABB(FT_HSE*);

        bool contains(AABB*) const;
        bool overlaps(AABB*) const;
        void print() const;
};



#endif


//Every BV type must be constructable from a pointer to
//a FT_HSE and have a centroid of type BV_Point.
//More to come...  


//Potential Methods:
//
//  1. AABB of 2 AABBs (constructor)
//  2. volume
//  3. inflation
