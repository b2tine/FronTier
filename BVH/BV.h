/*
 *
 *                  BV.h
 *
 *  
 *  Bounding Volume (BV) types used to parameterize the
 *  Bounding Volume Heirarchy template class, BVH, and
 *  the associated BVH_Node template class.
 *  
 *
 *  Current BV types: AABB
 *
 * */


#include "FT_HSE.h"

#ifndef BV_H
#define BV_H

#include <vector>


enum class BV_Type {AABB, OBB, KDOP, SPHERE};

using BV_Point = std::vector<double>;



class AABB
{
    public:

        BV_Point lower;
        BV_Point upper;
        BV_Point centroid;

        AABB() = default;
        ~AABB() = default;

        //delete copy and move ops until there
        //is a good reason not to.
        AABB(const AABB&) = delete;
        AABB& operator=(const AABB&) = delete;
        AABB(AABB&&) = delete;
        AABB& operator=(AABB&&) = delete;

        explicit AABB(FT_HSE*);

        bool contains(AABB*) const;
        bool overlaps(AABB*) const;

        void print() const;
        BV_Type getType() const;

};

//Potential AABB Methods:
//
//  1. AABB of 2 AABBs (constructor)
//  2. volume
//  3. inflation


//Every BV type must be constructable from a pointer to
//a FT_HSE and have a centroid of type BV_Point.
//
//  - contains()
//  - overlaps()
//
//  methods that operate on other BVs of the same type.  
//
//  More to come...


#endif
