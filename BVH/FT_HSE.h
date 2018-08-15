#ifndef FT_HSE_H
#define FT_HSE_H

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/hilbert_sort.h>
#include <algorithm>
#include <iostream>
#include <cassert>
#include "FronTier.h"


//FronTier hypersurface elements wrapper base class.
class FT_HSE
{
    public:

        FT_HSE() = default;
        virtual ~FT_HSE() = default;

        //delete copy and move ops since no member data
        FT_HSE(const FT_HSE&) = delete;
        FT_HSE& operator=(const FT_HSE&) = delete;
        FT_HSE(FT_HSE&&) = delete;
        FT_HSE& operator=(FT_HSE&&) = delete;

        virtual int num_pts() const = 0;
        virtual double max_coord(int) const = 0;
        virtual double min_coord(int) const = 0;
        virtual POINT* Point_of_hse(int) const = 0;

};

//Wrapper class for FronTier POINT
class FT_POINT : public FT_HSE
{
    private:

        POINT* point{nullptr};

    public:

        FT_POINT() = default;
        ~FT_POINT() = default;

        //delete copy and move ops until there
        //is a good reason not to.
        FT_POINT(const FT_POINT&) = delete;
        FT_POINT& operator=(const FT_POINT&) = delete;
        FT_POINT(FT_POINT&&) = delete;
        FT_POINT& operator=(FT_POINT&&) = delete;

        //NOTE: point is a shallow copy of p
        explicit FT_POINT(POINT* p)
            : point{p} {}

        POINT* Point_of_hse(int) const override;
        double min_coord(int) const override;
        double max_coord(int) const override;
        int num_pts() const override { return 1; }
};

//Wrapper class for FronTier BOND
class FT_BOND : public FT_HSE
{
    private:

        BOND* bond{nullptr};

    public:

        FT_BOND() = default;
        ~FT_BOND() = default;

        //delete copy and move ops until there
        //is a good reason not to.
        FT_BOND(const FT_BOND&) = delete;
        FT_BOND& operator=(const FT_BOND&) = delete;
        FT_BOND(FT_BOND&&) = delete;
        FT_BOND& operator=(FT_BOND&&) = delete;
        
        //NOTE: bond is a shallow copy of p
        explicit FT_BOND(BOND* b)
            : bond{b} {}

        POINT* Point_of_hse(int) const override;
        double min_coord(int) const override;
        double max_coord(int) const override;
        int num_pts() const override { return 2; }
};

//Wrapper class for FronTier TRI
class FT_TRI : public FT_HSE
{
    private:

        TRI* tri{nullptr};

    public:

        FT_TRI() = default;
        ~FT_TRI() = default;

        //delete copy and move ops until there
        //is a good reason not to.
        FT_TRI(const FT_TRI&) = delete;
        FT_TRI& operator=(const FT_TRI&) = delete;
        FT_TRI(FT_TRI&&) = delete;
        FT_TRI& operator=(FT_TRI&&) = delete;
        
        //NOTE: tri is a shallow copy of t
        explicit FT_TRI(TRI* t)
            : tri{t} {}

        POINT* Point_of_hse(int) const override;
        double min_coord(int) const override;
        double max_coord(int) const override;
        int num_pts() const override { return 3; }
};





#endif
