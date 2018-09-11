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

        AABB() = default;
        AABB(const AABB&) = default;
        AABB& operator=(const AABB&) = default;
        AABB(AABB&&) = default;
        AABB& operator=(AABB&&) = default;
        ~AABB() = default;

        explicit AABB(Hse*);
        AABB(const AABB&,const AABB&);
        AABB(const BV_Point&,const BV_Point&);

        const BV_Point getLower() const {return lower;}
        const BV_Point getUpper() const {return upper;}
        
        const BV_Point centroid() const;
        bool contains(const AABB&) const;
        bool overlaps(const AABB&) const;
        //void inflate();
        BV_Type getBvType() const;

        void print() const;
};


//TODO: may need notion of containment with shared surfaces


#endif
