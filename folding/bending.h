#include <FronTier.h>
#include "spring_solver.h"

class BendingForce: public SpringSolver::ExtForceInterface
{
    void calculateBendingForce3d2003(POINT* p1, TRI* t1, TRI* t2);
    void calculateBendingForce3d2006(POINT* p1, TRI* t1, TRI* t2);
    void calculateBendingForce3dparti(POINT*, TRI*, TRI*);
    INTERFACE* intfc;
    // bend stiffness
    double bends;
    // bend damping
    double bendd;
    // specify which method to be used
    // 0: 2003
    // 1: 2006
    // 2: particle
    // default particle
    int index;
    void clear_surf_point_force(SURFACE*);
    static const int num = 3;
    // pointer to member function
    void (BendingForce::*method[num])(POINT*,TRI*, TRI*);
public:
    double* getExternalForce(SpringVertex* sv);
    double calOriLeng(int, int, TRI*, TRI*);
    void computeExternalForce();
    void getParaFromFile(const char*);
    double& getBendStiff() { return bends; }
    double& getBendDamp() { return bendd; }
    int methodIndex() { return index; }
    BendingForce(INTERFACE* _intfc, double s = 0.01, double d = 0.0);

};

