#include "FT_HSE.h"

POINT* FT_POINT::Point_of_hse(int i) const
{
    assert( i >= 0 && i < this->num_pts() );
    assert( point != nullptr );
    return point;
}

double FT_POINT::min_coord(int dim) const
{
    assert( dim >= 0 && dim < 3 );
    return Coords(this->Point_of_hse(0))[dim];
}

double FT_POINT::max_coord(int dim) const
{
    assert( dim >= 0 && dim < 3 );
    return Coords(this->Point_of_hse(0))[dim];
}

POINT* FT_BOND::Point_of_hse(int i) const
{
    assert( i >= 0 && i < this->num_pts() );
    assert( bond != nullptr );
    return i == 0 ? bond->start : bond->end;
}

double FT_BOND::min_coord(int dim) const
{
    assert( dim >= 0 && dim < 3 );
    return std::min(Coords(this->Point_of_hse(0))[dim],
            Coords(this->Point_of_hse(1))[dim]);
}

double FT_BOND::max_coord(int dim) const
{
    assert( dim >= 0 && dim < 3 );
    return std::max(Coords(this->Point_of_hse(0))[dim],
            Coords(this->Point_of_hse(1))[dim]);
}

POINT* FT_TRI::Point_of_hse(int i) const
{
    assert( i >= 0 && i < this->num_pts() );
    assert( tri != nullptr );
    return Point_of_tri(tri)[i];
}

double FT_TRI::min_coord(int dim) const
{
    assert( dim >= 0 && dim < 3 );
    double val = HUGE;
    for( int i = 0; i < 3; i++ )
    {
        val = std::min(val,
                Coords(this->Point_of_hse(i))[dim]);
    }
    return val;
}

double FT_TRI::max_coord(int dim) const
{
    assert( dim >= 0 && dim < 3 );
    double val = -HUGE;
    for( int i = 0; i < 3; i++ )
    {
        val = std::max(val,
                Coords(this->Point_of_hse(i))[dim]);
    }
    return val;
}

