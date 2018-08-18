#ifndef BV_POINT_H
#define BV_POINT_H

#include <CGAL/hilbert_sort.h>

//TODO: figure out good way to test the hilbert_sort feature



//Point type that can be used in std::map 
class BV_Point
{
    private:

        double x{0}, y{0}, z{0};
    
    public:

        //need all default ctors for hilbert_sort()
        BV_Point() = default;
        BV_Point(const BV_Point&) = default;
        BV_Point& operator=(const BV_Point&) = default;
        BV_Point(BV_Point&&) = default;
        BV_Point& operator=(BV_Point&&) = default;
        ~BV_Point() = default;

        BV_Point(double X, double Y, double Z)
            : x{X}, y{Y}, z{Z}
        {}

        const double& operator[](const std::size_t i) const
        {
            assert(i >=0 && i <= 2);
            switch(i)
            {
                case 0: return x; break;
                case 1: return y; break;
                case 2: return z; break;
            }
        }

        //calls const version through typecasts
        double& operator[](const std::size_t i)
        {
            const_cast<double&>(
                    static_cast<const BV_Point&>(*this)[i]);
        }

        //needed for use as key in std::map<key,val>
        bool operator < (const BV_Point& p)
        {
            if( x == p[0] )
            {
                if( y == p[1] )
                {
                    return z < p[2];
                }
                return y < p[1];
            }
            return x < p[0];
        }
};


//additional structures allowing BV_Point to
//be used in the CGAL function hilbert_sort();

struct BV_LessX
{
    bool operator()(const BV_Point& p, const BV_Point& q) const
    {
        return p[0] < q[0];
    }
};

struct BV_LessY
{
    bool operator()(const BV_Point& p, const BV_Point& q) const
    {
        return p[1] < q[1];
    }
};

struct BV_LessZ
{
    bool operator()(const BV_Point& p, const BV_Point& q) const
    {
        return p[2] < q[2];
    }
};

struct BV_HilberSortingTraits
{
    using Point_3 = BV_Point;
    using Less_x_3 = BV_LessX;
    using Less_y_3 = BV_LessY;
    using Less_z_3 = BV_LessZ;

    Less_x_3 less_x_3_object() const
    {
        return Less_x_3();
    }

    Less_y_3 less_y_3_object() const
    {
        return Less_y_3();
    }

    Less_z_3 less_z_3_object() const
    {
        return Less_z_3();
    }
};


#endif
