/*
 *                  BV.h
 *
 *
 * Every BV type must be constructable from a FT_HSE* and
 * have the following member functions.
 * 
 *      1. BV_Type typeBV()
 *
 *      2. BV_Point centroid()
 *
 *      3. bool contains()
 *      
 *      4. bool overlaps()
 *
 *      5. void inflate()
 *
 *
 *      More to come ...
 *      
*/

#include "FT_HSE.h"

#ifndef BV_H
#define BV_H

#include <vector>


enum class BV_Type {AABB, OBB, KDOP, SPHERE};
using BV_Point = std::vector<double>;


//Axis Aligned Bounding Box (AABB)
class AABB
{
    public:

        BV_Point lower;
        BV_Point upper;
        BV_Point centroid;

        AABB() = default;
        ~AABB() = default;

        //Delete copy and move ops until there
        //is a good reason not to.
        AABB(const AABB&) = delete;
        AABB& operator=(const AABB&) = delete;
        AABB(AABB&&) = delete;
        AABB& operator=(AABB&&) = delete;

        explicit AABB(FT_HSE*);
        AABB(const BV_Point&,const BV_Point&);

        BV_Type getTypeBV() const;
        BV_Point getCentroid() const;
        bool contains(AABB*) const;
        bool overlaps(AABB*) const;
        //void inflate() override;
        void print() const;

        //create new AABB containing the two given as args
        static AABB* merge(AABB*,AABB*);
};

//TODO: may need notion of containment with shared surfaces


#endif
