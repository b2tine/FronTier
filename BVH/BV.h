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

        bool contains(AABB*);
        bool overlaps(AABB*);
        void print();
};



#endif


//Every BV type must be constructable from a pointer to
//a FT_HSE and have a centroid of type BV_Point.
//More to come...  
