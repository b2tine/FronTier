#include "HyperSurfElement.h"

#ifndef BOUNDING_VOLUME_H
#define BOUNDING_VOLUME_H

#include <vector>


enum class BV_Type {AABB};//, OBB, KDOP, SPHERE};

using BV_Point = std::vector<double>;


//Axis Aligned Bounding Box (AABB)
class AABB
{
    public:

        BV_Point lower;
        BV_Point upper;

        AABB();
        explicit AABB(Hse*);
        AABB(const AABB&,const AABB&);
        AABB(const BV_Point&,const BV_Point&);

        AABB(const AABB&) = default;
        AABB& operator=(const AABB&) = default;
        AABB(AABB&&) = default;
        AABB& operator=(AABB&&) = default;
        ~AABB() = default;

        const BV_Type getBvType() const;
        const BV_Point getCentroid() const;

        bool contains(const AABB&) const;
        bool overlaps(const AABB&) const;
        //void inflate();

        void print() const;
};


//TODO: may need notion of containment with shared surfaces


#endif
