#ifndef BV_POINT_H
#define BV_POINT_H

#include <CGAL/hilbert_sort.h>
#include <vector>


class BV_Point
{
    private:

        std::vector<double> point;
    
    public:

        BV_Point()
            : point(3,0.0)
        {}

        BV_Point(double x, double y, double z)
            : point{x,y,z}
        {}

        ~BV_Point() = default;
        BV_Point(const BV_Point&) = default;

        BV_Point& operator=(const BV_Point&) = delete;
        BV_Point(BV_Point&&) = delete;
        BV_Point& operator=(BV_Point&&) = delete;

        const double& operator[](const std::size_t i) const
        {
            assert(i >=0 && i <= 2);
            return point[i];
        }

        double& operator[](const std::size_t i)
        {
            const_cast<double&>(
                    static_cast<const BV_Point&>(*this)[i]);
        }

        bool operator < (const BV_Point& rhs) const
        {
            if( point[0] == rhs[0] )
            {
                if( point[1] == rhs[1] )
                {
                    return point[2] < rhs[2];
                }
                return point[1] < rhs[1];
            }
            return point[0] < rhs[0];
        }

};

//additional structures allowing BV_Point to
//be used in the CGAL function hilbert_sort();

/*
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

struct BV_HilbertSortingTraits
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
*/

#endif
