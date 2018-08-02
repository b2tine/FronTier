#include "spring_solver.h"
#include <iostream>
#include <stdexcept>

#include <cvode/cvode.h>             /* main integrator header file */
#include <cvode/cvode_spgmr.h>       /* prototypes & constants for CVSPGMR solver */
#include <cvode/cvode_spbcgs.h>      /* prototypes & constants for CVSPBCG solver */
#include <cvode/cvode_sptfqmr.h>     /* prototypes & constants for CVSPTFQMR solver */
#include <nvector/nvector_serial.h>  /* serial N_Vector types, fct. and macros */
#include <sundials/sundials_dense.h> /* use generic DENSE solver in preconditioning */
#include <sundials/sundials_types.h> /* definition of realtype */
#include <sundials/sundials_math.h>  /* contains the macros ABS, SUNSQR, and EXP */


static void updateSpringVertex(N_Vector, std::vector<SpringVertex*>&);
static void setInitialCondition(N_Vector, std::vector<SpringVertex*>&);
static int f(realtype, N_Vector, N_Vector, void*);
/*static int Jac(long int N, long int mu, long int ml,
               realtype t, N_Vector u, N_Vector fu,
               DlsMat J, void *user_data,
               N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);*/
static int check_flag(void *flagvalue, const char *funcname, int opt);

void IM_SPRING_SOLVER::doSolve(double t) {
    int NEQ = pts.size() * 6; //number of equations
    //const int MAX_NB = 7; //number of neighbours
    //create vector
    N_Vector u = N_VNew_Serial(NEQ);
    
    //set tolerances
    realtype reltol = RCONST(0.0);
    realtype abstol = RCONST(1.0e-5);

    //set initial condition
    setInitialCondition(u, pts);

    //Create ode solver
    void *cvode_mem = CVodeCreate(CV_BDF, CV_NEWTON);

    //Initialize
    int flag = CVodeInit(cvode_mem, f, 0, u);

    //set tolerance
    flag = CVodeSStolerances(cvode_mem, reltol, abstol);

    //sIMet pointer to user-defined data
    flag = CVodeSetUserData(cvode_mem, (void*)this);

    //call CVSpgmr to specify the KSP linear solver
    flag = CVSpbcg(cvode_mem, PREC_NONE, 0);

    //try to add preconditioner if iteration too many times 
    //flag = CVSpilsSetPreconditioner(cvode_mem, Precond, PSolve);
    realtype tret;
    long nst;
    //realtype umax = N_VMaxNorm(u);

    for (size_t i = 0; i < ext_forces.size(); ++i)	
	ext_forces[i]->computeExternalForce();
    //solve ode
    flag = CVode(cvode_mem, t, u, &tret, CV_NORMAL);

    flag = CVodeGetNumSteps(cvode_mem, &nst);

    if (check_flag(&flag, "CVodeGetNumSteps", 1)) 
	throw std::runtime_error("CVodeGetNumSteps");
    updateSpringVertex(u, pts);

    N_VDestroy_Serial(u);
    CVodeFree(&cvode_mem);
}

static void updateSpringVertex(N_Vector u, std::vector<SpringVertex*>& pts) {
    realtype *udata = N_VGetArrayPointer_Serial(u);
    for (size_t i = 0; i < pts.size(); ++i) {
        pts[i]->x[0] = udata[i*6];
        pts[i]->x[1] = udata[i*6+1];
        pts[i]->x[2] = udata[i*6+2];
        pts[i]->v[0] = udata[i*6+3];
        pts[i]->v[1] = udata[i*6+4];
        pts[i]->v[2] = udata[i*6+5];
    }
}

static void setInitialCondition(N_Vector u, std::vector<SpringVertex*>& pts) {
    realtype *udata = N_VGetArrayPointer_Serial(u);
    for (size_t i = 0; i < pts.size(); ++i) {
	udata[i*6]   = pts[i]->x[0];
        udata[i*6+1] = pts[i]->x[1];
        udata[i*6+2] = pts[i]->x[2];
        udata[i*6+3] = pts[i]->v[0];
        udata[i*6+4] = pts[i]->v[1];
        udata[i*6+5] = pts[i]->v[2];
    }
}

static int f(realtype t, N_Vector u, N_Vector udot, void* user_data) {
    realtype *dudata;
    IM_SPRING_SOLVER* sp_solver = (IM_SPRING_SOLVER*)user_data;
    std::vector<SpringVertex*>& pts = sp_solver->getSpringMesh();

    dudata = N_VGetArrayPointer_Serial(udot);

    //update the coords and velocity in spring vertex
    //and compute acceleration
    updateSpringVertex(u, pts);
    for (size_t i = 0; i < pts.size(); ++i)
	sp_solver->computeAccel(pts[i]);

    for (size_t i = 0; i < pts.size(); ++i) {
	//dx = v
	dudata[i*6]   = pts[i]->v[0];
	dudata[i*6+1] = pts[i]->v[1];
	dudata[i*6+2] = pts[i]->v[2];
	//dv = f
	dudata[i*6+3] = pts[i]->accel[0];
	dudata[i*6+4] = pts[i]->accel[1];
	dudata[i*6+5] = pts[i]->accel[2];
    }
    return 0;
}

/*static int Jac(long int N, long int mu, long int ml,
               realtype t, N_Vector u, N_Vector fu,
               DlsMat J, void *user_data,
               N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
    return 0;
}*/

static int check_flag(void *flagvalue, const char *funcname, int opt)
{
  int *errflag;

  /* Check if SUNDIALS function returned NULL pointer - no memory allocated */
  if (opt == 0 && flagvalue == NULL) {
    fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed - returned NULL pointer\n\n",
            funcname);
    return(1); }

  /* Check if flag < 0 */
  else if (opt == 1) {
    errflag = (int *) flagvalue;
    if (*errflag < 0) {
      fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed with flag = %d\n\n",
              funcname, *errflag);
      return(1); }}

  /* Check if function returned NULL pointer - no memory allocated */
  else if (opt == 2 && flagvalue == NULL) {
    fprintf(stderr, "\nMEMORY_ERROR: %s() failed - returned NULL pointer\n\n",
            funcname);
    return(1); }

  return(0);
}

