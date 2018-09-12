#ifndef HYPER_SURF_ELEMENT_H
#define HYPER_SURF_ELEMENT_H

#include <algorithm>
#include <iostream>
#include <exception>
#include <cassert>

#include "FronTier.h"

//Abstract base class for FronTier hypersurface element wrappers
class Hse
{
    public:

        Hse() = default;
        virtual ~Hse() = default;

        Hse(const Hse&) = delete;
        Hse& operator=(const Hse&) = delete;
        Hse(Hse&&) = delete;
        Hse& operator=(Hse&&) = delete;

        virtual double max_coord(int) const = 0;
        virtual double min_coord(int) const = 0;
        virtual POINT* Point_of_hse(int) const = 0;
        virtual int num_pts() const = 0;
};

//Wrapper for FronTier POINT
class HsPoint : public Hse
{
    private:

        POINT* point{nullptr};

    public:

        //TODO: throw exception if p is nullptr
        explicit HsPoint(POINT* p)
            : point{p}
        {}

        ~HsPoint() = default;

        HsPoint() = delete;
        HsPoint(const HsPoint&) = delete;
        HsPoint& operator=(const HsPoint&) = delete;
        HsPoint(HsPoint&&) = delete;
        HsPoint& operator=(HsPoint&&) = delete;

        POINT* Point_of_hse(int i = 0) const override;
        double min_coord(int) const override;
        double max_coord(int) const override;
        int num_pts() const override { return 1; }

        //TODO: Add method to check if a Hse (HsPoint, HsBond, or HsTri)
        //      is incident to an instance of HsPoint.
        //      May not need this one since we don't associate
        //      a bounding volume to points.
};

//Wrapper for FronTier BOND
class HsBond : public Hse
{
    private:

        BOND* bond{nullptr};

    public:

        //TODO: throw exception if b is nullptr
        explicit HsBond(BOND* b)
            : bond{b}
        {}

        ~HsBond() = default;

        HsBond() = delete;
        HsBond(const HsBond&) = delete;
        HsBond& operator=(const HsBond&) = delete;
        HsBond(HsBond&&) = delete;
        HsBond& operator=(HsBond&&) = delete;
        
        POINT* Point_of_hse(int) const override;
        double min_coord(int) const override;
        double max_coord(int) const override;
        int num_pts() const override { return 2; }
        
        //TODO: Add method to check if a Hse (HsPoint, HsBond, or HsTri)
        //      is incident to an instance of HsBond.
};


//Wrapper for FronTier TRI
class HsTri : public Hse
{
    private:

        TRI* tri{nullptr};

    public:

        //TODO: throw exception if t is nullptr
        explicit HsTri(TRI* t)
            : tri{t}
        {}

        ~HsTri() = default;

        HsTri() = delete;
        HsTri(const HsTri&) = delete;
        HsTri& operator=(const HsTri&) = delete;
        HsTri(HsTri&&) = delete;
        HsTri& operator=(HsTri&&) = delete;
        
        POINT* Point_of_hse(int) const override;
        double min_coord(int) const override;
        double max_coord(int) const override;
        int num_pts() const override { return 3; }

        //TODO: Add method to check if a Hse (HsPoint, HsBond, or HsTri)
        //      is incident to an instance of HsTri.
};




#endif
