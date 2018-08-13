#ifndef FT_HSE_H
#define FT_HSE_H

#include <iostream>
#include <cassert>
#include <algorithm>
#include "FronTier.h"


//FronTier hypersurface elements wrapper base class.
class FT_HSE
{
    public:

        FT_HSE() = default;
        FT_HSE(const FT_HSE&) = delete;
        FT_HSE& operator=(const FT_HSE&) = delete;
        FT_HSE(FT_HSE&&) = delete;
        FT_HSE& operator=(FT_HSE&&) = delete;
        virtual ~FT_HSE() = default;

        virtual double max_coord(int) const = 0;
        virtual double min_coord(int) const = 0;
        virtual POINT* Point_of_hse(int) const = 0;
        virtual int num_pts() const = 0;

};

//Wrapper class for FronTier POINT
class FT_POINT : public FT_HSE
{
    private:

        POINT* point;

    public:

        FT_POINT() = default;
        FT_POINT(const FT_POINT&) = delete;
        FT_POINT& operator=(const FT_POINT&) = delete;
        FT_POINT(FT_POINT&&) = delete;
        FT_POINT& operator=(FT_POINT&&) = delete;
        ~FT_POINT() = default;

        explicit FT_POINT(POINT* p)
            : point{p}
        {}

        int num_pts() const override { return 1; }
        POINT* Point_of_hse(int) const override;
        double min_coord(int) const override;
        double max_coord(int) const override;
};

//Wrapper class for FronTier BOND
class FT_BOND : public FT_HSE
{
    private:

        BOND* bond;

    public:

        FT_BOND() = default;
        FT_BOND(const FT_BOND&) = delete;
        FT_BOND& operator=(const FT_BOND&) = delete;
        FT_BOND(FT_BOND&&) = delete;
        FT_BOND& operator=(FT_BOND&&) = delete;
        ~FT_BOND() = default;
        
        explicit FT_BOND(BOND* b)
            : bond{b}
        {}

        int num_pts() const override { return 2; }
        POINT* Point_of_hse(int) const override;
        double min_coord(int) const override;
        double max_coord(int) const override;
};

//Wrapper class for FronTier TRI
class FT_TRI : public FT_HSE
{
    private:

        TRI* tri;

    public:

        FT_TRI() = default;
        FT_TRI(const FT_TRI&) = delete;
        FT_TRI& operator=(const FT_TRI&) = delete;
        FT_TRI(FT_TRI&&) = delete;
        FT_TRI& operator=(FT_TRI&&) = delete;
        ~FT_TRI() = default;

        explicit FT_TRI(TRI* t)
            : tri{t}
        {}

        int num_pts() const override { return 3; }
        POINT* Point_of_hse(int) const override;
        double min_coord(int) const override;
        double max_coord(int) const override;
};





#endif
