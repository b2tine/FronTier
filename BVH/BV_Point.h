#ifndef BV_POINT_H
#define BV_POINT_H

#include <CGAL/hilbert_sort.h>

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

//TODO: figure out good way to test the hilber_sort feature

#endif
