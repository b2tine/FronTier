#include <CGAL/Simple_cartesian.h>
#include <CGAL/box_intersection_d.h>
#include <FronTier.h>
#include "../iFluid/ifluid_state.h"
#include <functional>
#include <map>
#if defined(isnan)
#undef isnan
#endif

#define DEBUGGING false
const double ROUND_EPS = 1e-10;
const double EPS = 1e-6;
const double DT = 0.001;
/*
user-defined state should include the following
struct UF{
	POINT* next_pt;
	POINT* root;
	POINT* tail;
	int num_pts;
};

struct STATE{
	double vel[3];
	double collsnImpulse[3];
	double friction[3];
	double avgVel[3];
	double x_old[3];
	int    collsn_num;
	bool   has_collsn;
	bool   is_fixed;
	UF     impZone;
};*/

//abstract base class for hypersurface element(HSE)
//can be a point or a bond or a triangle
class CD_HSE{
public:
	virtual double max_static_coord(int) = 0;
	virtual double min_static_coord(int) = 0;
	virtual double max_moving_coord(int,double) = 0;
	virtual double min_moving_coord(int,double) = 0;
	virtual POINT* Point_of_hse(int) const  = 0;
	virtual int num_pts() const= 0;
	virtual ~CD_HSE(){};
};

//wrap class for triangle
class CD_TRI: public CD_HSE{
public:
	TRI* m_tri;
	CD_TRI(TRI* tri):m_tri(tri){}
	double max_static_coord(int);
	double min_static_coord(int);
	double max_moving_coord(int,double);
	double min_moving_coord(int,double);
	POINT* Point_of_hse(int) const;
	int num_pts() const {return 3;}
};

//wrap class for bond
class CD_BOND: public CD_HSE{
public:
	BOND* m_bond;
	int m_dim;
	CD_BOND(BOND* bond, int dim):m_bond(bond), m_dim(dim){}
	double max_static_coord(int);
	double min_static_coord(int);
	double max_moving_coord(int,double);
	double min_moving_coord(int,double);
	POINT* Point_of_hse(int) const;
	int num_pts()const{return 2;}
};

//wrap class for point
class CD_POINT: public CD_HSE{
public:
	POINT* m_point;
	CD_POINT(POINT* point):m_point(point){}
	double max_static_coord(int);
	double min_static_coord(int);
	double max_moving_coord(int,double);
	double min_moving_coord(int,double);
	POINT* Point_of_hse(int) const;
	int num_pts()const {return 1;}
};

//box traits structure for proximity detection 
struct traitsForProximity{
	typedef double 		NT;
	typedef CD_HSE*		Box_parameter;

	static int m_dim;
	static double s_eps;
	static int dimension(){ return m_dim;}
	static double min_coord(Box_parameter b, int d)
		{return b->min_static_coord(d)-s_eps;}
	static double max_coord(Box_parameter b, int d)
		{return b->max_static_coord(d)+s_eps;}
	static std::ptrdiff_t id(Box_parameter b) { return (std::ptrdiff_t)(b);}
};

//box traits structure for collision detection
struct traitsForCollision{
	typedef double          		NT;
        typedef CD_HSE*      Box_parameter;

	static int m_dim;
	static double s_eps;
	static double s_dt;
        static int dimension(){ return m_dim;}
        static double min_coord(Box_parameter b, int d)
		{return b->min_moving_coord(d,s_dt)-s_eps;}
        static double max_coord(Box_parameter b, int d)
		{return b->max_moving_coord(d,s_dt)+s_eps;}
	static std::ptrdiff_t id(Box_parameter b) { return (std::ptrdiff_t)(b);}
};

//abstract base class for collision detection and handling
class CollisionSolver {
private:
	//global parameters
	static double s_eps;
	static double s_thickness;
	static double s_dt;
	static double s_m;
	static double s_k;
	static double s_lambda;
	static double s_cr;
	bool has_collision;
	double Boundary[3][2]; //domain boundary[dir][side]
	static void turnOffImpZone();
	static void turnOnImpZone();
	bool reduceSuperelastOnce(int&);
	void computeAverageVelocity();
	void updateFinalPosition();
	void reduceSuperelast();
	void updateFinalVelocity();
	void updateAverageVelocity();
	void computeImpactZone();
	void updateImpactZoneVelocity(int&);
	void updateImpactZoneVelocityForRG();
	void setTraitsDimension();
	void detectProximity();
	void detectCollision();
	void detectDomainBoundaryCollision();
	void updateFinalForRG();
	void setHasCollision(bool judge) {has_collision = judge;}

	virtual void updateImpactListVelocity(POINT*) = 0;
	virtual bool BondToBond(const BOND*,const BOND*,double) = 0;
	virtual bool TriToTri(const TRI*,const TRI*,double) = 0;
	virtual bool TriToBond(const TRI*,const BOND*,double)=0;
	virtual bool MovingBondToBond(const BOND*,const BOND*,double) = 0;
	virtual bool MovingTriToTri(const TRI*,const TRI*,double) = 0;
	virtual bool MovingTriToBond(const TRI*,const BOND*,double)=0;

protected:
	int m_dim;
	std::vector<CD_HSE*> hseList;
	std::map< int, std::vector<double> > mrg_com;
	static bool s_detImpZone;
	void clearHseList();
public:
	CollisionSolver(int dim):m_dim(dim){}
	CollisionSolver(){}
	static void setRoundingTolerance(double);
	static double getRoundingTolerance();
	static void setFabricThickness(double);
	static double getFabricThickness();
	static void setTimeStepSize(double);
	static double getTimeStepSize();
	static void setSpringConstant(double);
	static double getSpringConstant();
	static void setFrictionConstant(double);
	static double getFrictionConstant();
	static void setPointMass(double);
	static double getPointMass();
	static void setRestitutionCoef(double);
	static double getRestitutionCoef();
	static bool getImpZoneStatus();	
	virtual ~CollisionSolver(){} //virtual destructor
	//pure virtual functions
	virtual void assembleFromInterface(const INTERFACE*,double dt) = 0;
	virtual void createImpZoneForRG(const INTERFACE*) = 0;
	bool isProximity(const CD_HSE*,const CD_HSE*);	
	bool isCollision(const CD_HSE*,const CD_HSE*);
	void resolveCollision();
	void recordOriginPosition();	
	void setDomainBoundary(double* L,double *U);
	double getDomainBoundary(int dir,int side) {return Boundary[dir][side];}
	bool hasCollision() {return has_collision;}

	//for debugging
	static void printDebugVariable();
	static int moving_edg_to_edg;
	static int moving_pt_to_tri;
	static int is_coplanar;
	static int edg_to_edg;
	static int pt_to_tri;
};

//derived 2D-class for collision detection and handling
class CollisionSolver2d : public CollisionSolver {
private:
	void updateImpactListVelocity(POINT*);
	bool BondToBond(const BOND*,const BOND*,double);
	bool TriToTri(const TRI*,const TRI*,double);
	bool TriToBond(const TRI*,const BOND*,double);
public:
	CollisionSolver2d():CollisionSolver(2){}
	void assembleFromInterface(const INTERFACE*,double dt);
	void createImpZoneForRG(const INTERFACE*);
};

//derived 3D-class for collision detection and handling
class CollisionSolver3d : public CollisionSolver {
private:
	//input tris
	void updateImpactListVelocity(POINT*);
	bool BondToBond(const BOND*,const BOND*,double);
	bool TriToTri(const TRI*,const TRI*,double);
	bool TriToBond(const TRI*,const BOND*,double);
	bool MovingBondToBond(const BOND*,const BOND*,double);
	bool MovingTriToTri(const TRI*,const TRI*,double);
	bool MovingTriToBond(const TRI*,const BOND*,double);
public:
	CollisionSolver3d():CollisionSolver(3){}
	void assembleFromInterface(const INTERFACE*,double dt);
	void createImpZoneForRG(const INTERFACE*);
};

//callback functor to identify real collision
struct reportProximity{
    int& num_pairs;
    CollisionSolver* collsn_solver;
    reportProximity(int &npair,CollisionSolver* solver): 
			 	 num_pairs(npair = 0),
				 collsn_solver(solver){}
    void operator()( const CD_HSE* a, const CD_HSE* b) {
	if(collsn_solver->isProximity(a,b)){
	    num_pairs++;
	}
    }
};

struct reportCollision{
    bool& is_collision;
    int&  num_pairs;
    CollisionSolver* collsn_solver;
    reportCollision(bool &status, int &npairs,CollisionSolver* solver): 
		     is_collision(status), 
		     num_pairs(npairs = 0), 
		     collsn_solver(solver){}
    void operator()( const CD_HSE* a, const CD_HSE* b) {
	if (collsn_solver->isCollision(a,b)){
	    num_pairs ++;
	    is_collision = true;
	}
    }
};

void initSurfaceState(SURFACE*,const double*);
void initCurveState(CURVE*,const double*);
void initTestModule(Front&, char*);
void Pts2Vec(const POINT*, const POINT*, double*); 
void scalarMult(double a,double* v, double* ans); 
void addVec(double* v1, double* v2, double* ans); 
void minusVec(double* v1, double* v2, double* ans); 
double myDet3d(double[][3]);
double distBetweenCoords(double* v1, double* v2);
extern void printPointList(POINT**, const int);
extern void createImpZone(POINT*[],int num = 4,bool first = NO);
extern void makeSet(std::vector<CD_HSE*>&);
void unsortHseList(std::vector<CD_HSE*>&);
POINT*& next_pt(POINT*);
int& weight(POINT*);
bool isStaticRigidBody(const POINT*);
bool isStaticRigidBody(const CD_HSE*);
bool isMovableRigidBody(const POINT*);
bool isMovableRigidBody(const CD_HSE*);
bool isRigidBody(const POINT*);
bool isRigidBody(const CD_HSE*);
extern void SpreadImpactZoneImpulse(POINT*, double, double*);

void vtkplotVectorSurface(std::vector<CD_HSE*>&,const char*);
