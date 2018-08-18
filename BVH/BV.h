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

#include "BV_Point.h"
#include "FT_HSE.h"

#ifndef BV_H
#define BV_H

//#include <vector>


enum class BV_Type {AABB, OBB, KDOP, SPHERE};

//using BV_Point = std::vector<double>;
/*
 * Note: BV_Point is now its own class that
 *      is compatible with as a key type in
 *      std::map<key,val>.
 *
 *      TODO: Test functors needed for
 *      compatability with the CGAL function
 *      hilbert_sort().
 *
 *      If we end up implementing our own
 *      hilbert_sort() we will probably go
 *      back to using a std::vector<double>.
 *
 */

//Axis Aligned Bounding Box (AABB)
class AABB
{
    public:

        BV_Point lower;
        BV_Point upper;
        //BV_Point centroid;

        AABB() = default;
        AABB(const AABB&) = default;
        AABB& operator=(const AABB&) = default;
        AABB(AABB&&) = default;
        AABB& operator=(AABB&&) = default;
        ~AABB() = default;

        explicit AABB(FT_HSE*);
        AABB(const AABB&,const AABB&);
        AABB(const BV_Point&,const BV_Point&);

        const BV_Point centroid() const;
        bool contains(const AABB&) const;
        bool overlaps(const AABB&) const;
        //void inflate();
        BV_Type getTypeBV() const;
        void print() const;
};


//TODO: may need notion of containment with shared surfaces


#endif
