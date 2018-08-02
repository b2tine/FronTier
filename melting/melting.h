/**********************************************************************
 * 		liquid.h
 * the code is a direct modification of the code by Leveque[1].
 *
 * References:
 * [1] R.J. Leveque, High-resolution conservative algorithms for
 *     advection in incompressible flow.
 *
 **********************************************************************/

#include <FronTier.h>
#include <vector>
#include <petscksp.h>
#include <assert.h>

#define         LIQUID_COMP            3
#define		alternate_comp(comp) 					\
		(comp) == LIQUID_COMP ? SOLID_COMP : LIQUID_COMP

enum _NUM_SCHEME {
	UNSPLIT_EXPLICIT = 1,
	UNSPLIT_EXPLICIT_CIM,
	UNSPLIT_IMPLICIT,
	UNSPLIT_IMPLICIT_CIM,
	CRANK_NICOLSON
};
typedef enum _NUM_SCHEME NUM_SCHEME;

struct _PHASE_FIELD {
        double *temperature;
        double **vel;
};
typedef struct _PHASE_FIELD PHASE_FIELD;

struct _PARAMS {
        int dim;
	NUM_SCHEME num_scheme;
	PHASE_FIELD *field;
	int pde_order;
	int num_phases;
        double *Ti;    /* melting temperature at the interface */
	double *L;       /* latent heat at the interface */
	double *T0;	/* Ambient temperature of the phase */
        double *rho;    /* density of the phase */
        double *Cp;     /* specific heat of the phase */
        double *k;   /* thermal conductivity of the phase */
	double D;
	double min_temperature;
	double max_temperature;
	boolean no_fluid;
};
typedef struct _PARAMS PARAMS;

typedef class CARTESIAN CARTESIAN_EB;

/**************************************************************************
 *		vector/matrix operation functions 
 **************************************************************************/

void VectorIntPrint(int my_rank, char *name, int n, int *vector);
void VectorPrint(int my_rank, char *name, int n, double *vector);
void VectorZero(int n, double *vector);
void VectorZero(int n, int *vector);
void VectorCopy(int n, double *vector1, double *vector2);
void MatrixPrint(int my_rank, char *name, int m, int n, double **matrix);
void MatrixZero(int m, int n, double **matrix);
void MatrixIdentity(int n, double **matrix);
void MatrixCopy(int m, int n, double **U2, double **U1);
void ArrayZero(int m, int n, int l, double ***matrix);
void Vector2Matrix(int m, int n, double **matrix, double *vector);
void ArrayPrint(int my_rank, char*name, int m, int n, int l, double ***array);
void MatrixMultiply(int n, double **C, double **A, double **B);	// C = AB
void MatrixMultiply(int n, double **C, double **A);		// C = A'A 
void MatrixMultiply(double C[3][3], double A[3][3], double B[3][3]);  // C = AB
void MatrixVectorMultiply(int n, double *C, double **A, double *b);   // c = Ab
void SymmetricMatrixInverse(double B[3][3], double A[3][3]); // B = inverse(A); 
double MatrixMax(int m, int n, double **matrix);
double MatrixMin(int m, int n, double **matrix);

/******************************************************************************
 * 		lcartsn.h
 * A simple incompressible flow solver using the ghost fluid method and the
 * projection method.
 *
 * the main function is 
 * 	CARTESIAN::solve().
 *
 * References:
 ******************************************************************************/

class SOLVER;
class CARTESIAN;

class RECTANGLE {
public:
	int index;			// rectangle index
	int comp;			 
	double area;
	double coords[MAXD];	
	int icoords[MAXD];

	RECTANGLE();

	void setCoords(double*,int);
};


class CARTESIAN{
	Front *front;
public:
	CARTESIAN(Front &front);

	// member data: RECT_GRID
	int dim;

	// On topological grid
	RECT_GRID *top_grid;
	double *array;		// for scatter states;
	double *source;		// for source of parabolic solver;
	double *top_L,*top_U,*top_h,hmin;
	int *top_gmax;
	COMPONENT *top_comp;
	PARAMS *eqn_params;
	PHASE_FIELD *field;
	int comp_size;

	int *lbuf,*ubuf,*gmax;
	int *i_to_I,*I_to_i;		// Index mapping for 1D
	int **ij_to_I,**I_to_ij;	// Index mapping for 2D
	int ***ijk_to_I,**I_to_ijk;	// Index mapping for 3D

	// Sweeping limites
	int imin,jmin,kmin;
	int imax,jmax,kmax;

	enum BC_TYPE { 									// used by ADVECTION
		BC_PERIODIC = 1,
		BC_Extrapolation0 = 2,
		BC_Extrapolation1 = 3,
		BC_InflowOutflow = 4};	
	BC_TYPE m_bc[4];								// down, right, up, left 		

	// member data: mesh storage
	std::vector<RECTANGLE> 	cell_center;

	double m_t;                     // time
	double m_dt;			// time increment
	double min_dt;			// Minimum dt to use non-explicit

	// constructor
	~CARTESIAN();

	// for parallel partition
	int             NLblocks,ilower,iupper;
        int             *n_dist;

	// mesh: full cells mesh
	void initMesh(void);		// setup the cartesian grid
	void setDomain(void);
	void setComponent(void);	// init components	
	void setAdvectionDt(void);	// compute time step 	
	void setBoundary(void);		// set up boundary conditions 	
	void readFrontInteriorState(char*);
	void printFrontInteriorState(char*);

	void computeAdvection();

	void computeAdvectionExplicitCim(COMPONENT);
	void computeAdvectionExplicit(COMPONENT);
	void computeAdvectionImplicit(COMPONENT);
	void computeAdvectionCN(COMPONENT);
	void computeAdvectionCim();

	void pointExplicitCimSolver(int*,COMPONENT);

	// interface functions
	void makeGridIntfc();
	void deleteGridIntfc();

	// Extra plot functions
	void oneDimPlot(char*);
	void xgraphOneDimPlot(char*);
	void initMovieVariables(void);
	void vtk_plot_temperature2d(char*);
        void vtk_plot_temperature3d(char*);
	void initSampleTemperature(char *in_name);
        void sampleTemperature();
	void sampleTemperature2d();
	void sampleTemperature3d();

	// Extra movie functions
	void temperatureMovie(char*);

	void checkStates();

	// parallelization related functions
	//
	void setGlobalIndex(COMPONENT);

	// physics calculation
	void setInitialCondition(void);

	void setIndexMap(COMPONENT);
		// for compProjWithSmoothProperty(), 
		// should be changed to use setIndexMap() only.

	// -------------------------------------------------------
	// 		incompressible solver functions

	// main step function
	void solve(double dt);		

	void getVelocity(double *p, double *U);

	// ----------------------------------------------------------
	// 		utility functions
	// ----------------------------------------------------------

	void getRectangleIndex(int indexRectangle, int &i, int &j);
	void getRectangleIndex(int indexRectangle, int &i, int &j, int &k);
	void getRectangleCenter(int index, double *coords);
	void getRectangleCenter(int index0, int index1, double *coords);
	int getRectangleComponent(int index);	// the center component
	
	double getDistance(double *coords0, double *coords1);
	
			// incompletely implemented
	void getNearestInterfacePoint(double *q,double *p); 
		
	int  getComponent(int *icoords);	
	int  getComponent(double *coords);	
	void save(char *filename);
};

extern void readPhaseParams(Front*);
extern void initPhaseIntfc(char*,int,LEVEL_FUNC_PACK*,PARAMS*);
extern void ifluid_point_propagate(Front*,POINTER,POINT*,POINT*,
			HYPER_SURF_ELEMENT*,HYPER_SURF*,double,double*);
extern void melting_point_propagate(Front*,POINTER,POINT*,POINT*,
                    HYPER_SURF_ELEMENT*,HYPER_SURF*,double,double*);
extern void read_crt_dirichlet_bdry_data(char*,Front*,F_BASIC_DATA);
extern void read_melt_dirichlet_bdry_data(char*,Front*,F_BASIC_DATA);
extern void melt_flowThroughBoundaryState(double*,HYPER_SURF*,Front*,
			POINTER,POINTER);
extern void read_fluid_params(Front*);
extern void init_fluid_state_func(Front*,Incompress_Solver_Smooth_Basis*);		
extern void assignStateTemperature(double,POINTER);
extern double getStateTemperature(POINTER);
extern double jumpT(POINTER,int,double*);
extern double jumpEpsGradDotNorm(POINTER,int,double*,double*);
extern double jumpGradDotTan(POINTER,int,int,double*,double*);
