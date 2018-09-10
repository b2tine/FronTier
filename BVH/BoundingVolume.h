#include "BV_Point.h"
#include "FT_HSE.h"

#ifndef BOUNDING_VOLUME_H
#define BOUNDING_VOLUME_H


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
 *      NOTE: If we end up implementing our own
 *      hilbert_sort() we will probably go
 *      back to using a std::vector<double>.
 *
 */



enum class BV_Type {AABB, OBB, KDOP, SPHERE};




class BoundingVolume
{
    public:
        virtual const BV_Point centroid() const;
        virtual bool contains(const AABB&) const;
        virtual bool overlaps(const AABB&) const;
        //virtual void inflate();
        virtual BV_Type TypeBV() const;
        virtual ~BoundingVolume();
};


//Axis Aligned Bounding Box (AABB)
class AABB : public BoundingVolume
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

        explicit AABB(FT_HSE*);
        AABB(const AABB&,const AABB&);
        AABB(const AABB* const ,const AABB* const);
        AABB(const BV_Point&,const BV_Point&);

        const BV_Point centroid() const override;
        bool contains(const AABB&) const override;
        bool overlaps(const AABB&) const override;
        //void inflate() override;
        BV_Type getBvType() const override;

        void print() const;
};


//TODO: may need notion of containment with shared surfaces


#endif
