#ifndef SPRING_SOLVER_H
#define SPRING_SOLVER_H
#include <vector>
#include <algorithm>
#include <string>
#include "drag.h"

#define SPRING_EPS 1e-10
struct SpringVertex 
{
    double* x; //pointing to Point->x
    double* v; //pointing to left_state(Point)->vel
    void* org_vtx;//pointing to external point
    double accel[3];
    double ext_accel[3];
    std::vector<size_t> index_nb;
    std::vector<double> length0; 
    bool is_registered;
    //this could be used to identify
    //point type at runtime, should be 
    //set up before applying drag
    int point_type;
    
    SpringVertex():is_registered(false){ 
	for (int i = 0; i < 3; i++) { 
	     accel[i] = 0.0; 
	     ext_accel[i] = 0.0; 
	}
    }
    SpringVertex(SpringVertex&);
    virtual ~SpringVertex(){}
    double* getVel() {return v;}
    double* getCoords() {return x;}
    double* getExternalAccel() {return ext_accel;}
    void addNeighbor(size_t,double);
    bool isRegistered(){return is_registered;}
    void setRegistered(){is_registered = true;}
    void unsetRegistered(){is_registered = false;}
};

class Drag;
class SpringSolver 
{
private:
    double m_dt;  //time intrement
    Drag* m_drag; //control the prescribed points
    void checkVertexNeighbors();
public:
    struct SpringParameter {
	double k;
	double lambda;
	double m;
	SpringParameter(){}
	SpringParameter(double stf, double fri, double mas) :
		k(stf), lambda(fri), m(mas) {}
	//copy constructor
	SpringParameter(SpringParameter& params) {
	    k = params.k;
	    lambda = params.lambda;
	    m = params.m;
	}
    };
    enum ODE_SCHEME {EXPLICIT, IMPLICIT};
    void setParameters(double,double,double);
    static SpringSolver* createSpringSolver(ODE_SCHEME);
    void resetVelocity();
    void printAdjacencyList(std::string);
    void printPointList(std::string);
    virtual double getTimeStepSize() {return m_dt;}
    virtual void setTimeStepSize(double dt) {m_dt = dt;}
    void setParameters(SpringParameter&); 
    void setDrag(Drag*);
    virtual ~SpringSolver(){}
    void solve(double);
    //this function retrieves spring vertex 
    //you can fill in the spring vertex from outside
    std::vector<SpringVertex*>& getSpringMesh() {return pts;}

    //function for computing acceleration
    void computeAccel(SpringVertex*);
    void computeJacobian();

    //provide interface for computing acceleration externally
    class ExtForceInterface 
    {
	public:
	virtual void computeExternalForce() = 0;
	virtual double* getExternalForce(SpringVertex*) = 0;	
    };
    std::vector<ExtForceInterface*> ext_forces;

protected:
    virtual void doSolve(double) = 0;
    std::vector<SpringVertex*> pts;
    void setPresetVelocity(SpringVertex*);
    SpringParameter springParameter;
    //singleton pattern
    //using virtual constructor: createSpringSolver()
    SpringSolver(){}
};

//solve spring-mass ode system using 
//4th-order explicit Runge-Kutta method
class EX_SPRING_SOLVER: public SpringSolver {
public:
    void doSolve(double);
private:
    double getTimeStepSize();
    std::vector<std::vector<double> > x_old;
    std::vector<std::vector<double> > x_new;
    std::vector<std::vector<double> > v_old;
    std::vector<std::vector<double> > v_new;
};

//solve spring-mass ode system using
//2nd-order implicit BDF method
class IM_SPRING_SOLVER: public SpringSolver {
public:
   void doSolve(double);
};

#endif
