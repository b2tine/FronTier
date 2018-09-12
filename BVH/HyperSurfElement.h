#ifndef HYPER_SURF_ELEMENT_H
#define HYPER_SURF_ELEMENT_H

#include <algorithm>
#include <iostream>
#include <exception>
#include <cassert>

#include "FronTier.h"


//FronTier hypersurface elements wrapper abstract base class.
class Hse
{
    public:

        Hse() = default;
        virtual ~Hse() = default;

        //delete copy and move ops since no member data
        Hse(const Hse&) = delete;
        Hse& operator=(const Hse&) = delete;
        Hse(Hse&&) = delete;
        Hse& operator=(Hse&&) = delete;

        virtual double max_coord(int) const = 0;
        virtual double min_coord(int) const = 0;
        virtual POINT* Point_of_hse(int) const = 0;
        virtual int num_pts() const = 0;
};

//Wrapper class for FronTier POINT
class HsPoint : public Hse
{
    private:

        POINT* point{nullptr};

    public:

        //TODO: throw exception if p is nullptr
        explicit HsPoint(POINT* p)
            : point{p}
        {}

        HsPoint() = default;
        ~HsPoint() = default;

        HsPoint(const HsPoint&) = delete;
        HsPoint& operator=(const HsPoint&) = delete;
        HsPoint(HsPoint&&) = delete;
        HsPoint& operator=(HsPoint&&) = delete;

        POINT* Point_of_hse(int i = 0) const override;
        double min_coord(int) const override;
        double max_coord(int) const override;
        int num_pts() const override { return 1; }
};

//Wrapper class for FronTier BOND
class HsBond : public Hse
{
    private:

        BOND* bond{nullptr};

    public:

        //TODO: throw exception if b is nullptr
        explicit HsBond(BOND* b)
            : bond{b}
        {}

        HsBond() = default;
        ~HsBond() = default;

        HsBond(const HsBond&) = delete;
        HsBond& operator=(const HsBond&) = delete;
        HsBond(HsBond&&) = delete;
        HsBond& operator=(HsBond&&) = delete;
        
        POINT* Point_of_hse(int) const override;
        double min_coord(int) const override;
        double max_coord(int) const override;
        int num_pts() const override { return 2; }
};

//TODO: Add method to check if an Hse (HsBond or HsTri)
//      is HsBond's neighbor.

//Wrapper class for FronTier TRI
class HsTri : public Hse
{
    private:

        TRI* tri{nullptr};

    public:

        //TODO: throw exception if t is nullptr
        explicit HsTri(TRI* t)
            : tri{t}
        {}

        HsTri() = default;
        ~HsTri() = default;

        HsTri(const HsTri&) = delete;
        HsTri& operator=(const HsTri&) = delete;
        HsTri(HsTri&&) = delete;
        HsTri& operator=(HsTri&&) = delete;
        
        POINT* Point_of_hse(int) const override;
        double min_coord(int) const override;
        double max_coord(int) const override;
        int num_pts() const override { return 3; }
};

//TODO: Add method to check if an Hse (HsBond or HsTri)
//      is HsTri's neighbor.



#endif
