#include "gmock/gmock.h"
#include "FronTier.h"

//Base class for hypersurface elements (HSE).
//Can refer to a point, bond, or triangle
//depending on the dimension of the hypersurface.
class FT_HSE
{
    public:
        ~FT_HSE() = default;
};

class FT_POINT
    : public FT_HSE
{
    public:
};



TEST(FT_HSE, ReturnNumberOfPointsInHSE)
{
    FT_POINT ft_point;
}
