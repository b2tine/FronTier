#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/box_intersection_d.h>
#include <CGAL/intersections.h>
#include <vector>
#include <fstream>
#include <math.h>
#include <algorithm>
#include <FronTier.h>
#include "collid.h"
#include "../iFluid/ifluid_state.h"
#include <omp.h>

/*****declaration of static functions starts here********/
//static void makeSet(std::vector<CD_HSE*>&);
static POINT* findSet(POINT*);
static void mergePoint(POINT*,POINT*);
inline POINT*& root(POINT*);
inline POINT*& tail(POINT*);
/*******************end of declaration*******************/

typedef CGAL::Exact_predicates_inexact_constructions_kernel   Kernel;
typedef Kernel::Point_3                                       Point_3;
typedef Kernel::Triangle_3                                    Triangle_3;

//define default parameters for collision detection
bool   CollisionSolver::s_detImpZone = false;
double CollisionSolver::s_eps = EPS;
double CollisionSolver::s_thickness = 0.0001;
double CollisionSolver::s_dt = DT;
double CollisionSolver::s_k = 1000;
double CollisionSolver::s_m = 0.01;
double CollisionSolver::s_lambda = 0.02;
double CollisionSolver::s_cr = 0.0;
int traitsForProximity::m_dim = 3;
int traitsForCollision::m_dim = 3;
double traitsForProximity::s_eps = EPS;	
double traitsForCollision::s_eps = EPS;
double traitsForCollision::s_dt = DT;

//debugging variables
int CollisionSolver::moving_edg_to_edg = 0;
int CollisionSolver::moving_pt_to_tri = 0;
int CollisionSolver::is_coplanar = 0;
int CollisionSolver::edg_to_edg = 0;
int CollisionSolver::pt_to_tri = 0;

//functions in the abstract base class
void CollisionSolver::clearHseList(){
	for (unsigned i = 0; i < hseList.size(); ++i){
		delete hseList[i];
	}
	hseList.clear();
}
//set rounding tolerance
void CollisionSolver::setRoundingTolerance(double neweps){
	s_eps = neweps;
	traitsForProximity::s_eps = neweps;	
	traitsForCollision::s_eps = neweps;
}
double CollisionSolver::getRoundingTolerance(){return s_eps;}

//set fabric thickness
void CollisionSolver::setFabricThickness(double h){s_thickness = h;}
double CollisionSolver::getFabricThickness(){return s_thickness;}

//this function should be called at every time step
void CollisionSolver::setTimeStepSize(double new_dt){	
	s_dt = new_dt;
	traitsForCollision::s_dt = new_dt;
}
double CollisionSolver::getTimeStepSize(){return s_dt;}

//set spring constant
void   CollisionSolver::setSpringConstant(double new_k){s_k = new_k;}
double CollisionSolver::getSpringConstant(){return s_k;}

//set spring friction 
void   CollisionSolver::setFrictionConstant(double new_la){s_lambda = new_la;}
double CollisionSolver::getFrictionConstant(){return s_lambda;}

//set mass of fabric point
void   CollisionSolver::setPointMass(double new_m){s_m = new_m;}
double CollisionSolver::getPointMass(){return s_m;}

//set restitution coefficient between rigid bodies
void   CollisionSolver::setRestitutionCoef(double new_cr){s_cr = new_cr;}
double CollisionSolver::getRestitutionCoef(){return s_cr;}

void CollisionSolver::recordOriginPosition(){
	POINT* pt;
	STATE* sl;
	//#pragma omp parallel for private(pt,sl)
	for (std::vector<CD_HSE*>::iterator it = hseList.begin();
	     it < hseList.end(); ++it){
	    for (int i = 0; i < (*it)->num_pts(); ++i){
		pt = (*it)->Point_of_hse(i);
		sl = (STATE*)left_state(pt); 
		sl->has_collsn = false;
		if (isMovableRigidBody(pt)) continue;
		for (int j = 0; j < 3; ++j)
		    sl->x_old[j] = Coords(pt)[j];
		if (std::isnan(sl->x_old[0])) std::cout<<"nan_x_old"<<std::endl;
	    }
	}
}

void CollisionSolver::setDomainBoundary(double* L, double* U) {
	for (int i = 0; i < m_dim; ++i) {
	    Boundary[i][0] = L[i];
	    Boundary[i][1] = U[i];
	}
}

void CollisionSolver::detectDomainBoundaryCollision() {
	double dt = getTimeStepSize();
	double mu = getFrictionConstant();
	for (std::vector<CD_HSE*>::iterator it = hseList.begin();
                it < hseList.end(); ++it) {
	    for (int i = 0; i < (*it)->num_pts(); ++i) {
		POINT* pt = (*it)->Point_of_hse(i);
		if (isMovableRigidBody(pt)) continue;
                STATE* sl = (STATE*)left_state(pt);
		double cand_coords[3]; //candidate position
		//try to modify the average velocity 
		//according to the new candidate position
		double dv = 0;
		for (int j = 0; j < m_dim; ++j) {
		    cand_coords[j] = sl->x_old[j] + dt*sl->avgVel[j];
		    double L, U;
		    L = getDomainBoundary(j,0);
		    U = getDomainBoundary(j,1);
		    if (cand_coords[j] <= L) 
		    {
			sl->has_collsn = true;
			cand_coords[j] = L + s_thickness;
			dv = fabs(sl->avgVel[j]);
		    	sl->avgVel[j] = 0.0;
		        Coords(pt)[j] = cand_coords[j];
		    }
		    else if (cand_coords[j] >= U)
		    {
			sl->has_collsn = true;
			cand_coords[j] = U - s_thickness;
			dv = fabs(sl->avgVel[j]);
		    	sl->avgVel[j] = 0.0;
		        Coords(pt)[j] = cand_coords[j];
		    }
		}
		//reduce tangential velocity with friction
		double preVt = Mag3d(sl->avgVel);
		if (preVt > MACH_EPS)
		for (int j = 0; j < m_dim; ++j) 
		    sl->avgVel[j] *= std::max(1.0-mu*dv/preVt,0.0); 
	    }
	}
}

void CollisionSolver::computeAverageVelocity(){
	POINT* pt;
        STATE* sl; 
	double dt = getTimeStepSize();
	double max_speed = 0, *max_vel = NULL;
	POINT* max_pt;
	//#pragma omp parallel for private(pt,sl)
        for (std::vector<CD_HSE*>::iterator it = hseList.begin();
                it < hseList.end(); ++it){
            for (int i = 0; i < (*it)->num_pts(); ++i){
                pt = (*it)->Point_of_hse(i);
                sl = (STATE*)left_state(pt); 
                for (int j = 0; j < 3; ++j)
		{
		    if (dt > ROUND_EPS)
                        sl->avgVel[j] = (Coords(pt)[j] - sl->x_old[j])/dt;
		    else
		        sl->avgVel[j] = 0.0;
		    if (std::isnan(sl->avgVel[j]) || std::isinf(sl->avgVel[j]))
		    {
			std::cout<<"nan avgVel" << std::endl;
			printf("dt = %e, x_old = %f, x_new = %f\n",
			dt,sl->x_old[j],Coords(pt)[j]);
			clean_up(ERROR);
		    }
		}
		if (debugging("collision"))
		if (Mag3d(sl->avgVel) >= max_speed){
		    max_speed = Mag3d(sl->avgVel);
		    max_vel = sl->avgVel;
		    max_pt = pt;
		}
            }
        }
	if (debugging("collision"))
	{
	    std::cout << "Maximum average velocity is " 
		      << max_vel[0] << " "
		      << max_vel[1] << " "
		      << max_vel[2] << std::endl; 
	    sl = (STATE*)left_state(max_pt);
	    printf("x_old = [%f %f %f]\n",
		   sl->x_old[0],sl->x_old[1],sl->x_old[2]);
	    printf("x_new = [%f %f %f]\n",
		   Coords(max_pt)[0],Coords(max_pt)[1],Coords(max_pt)[2]);
	    printf("dt = %f\n",dt);
	}
	//restore coords of points to old coords !!!
	//x_old is the only valid coords for each point 
	//Coords(point) is for temporary judgement
	//#pragma omp parallel private(pt,sl)
	for (std::vector<CD_HSE*>::iterator it = hseList.begin();
                it < hseList.end(); ++it){
            for (int i = 0; i < (*it)->num_pts(); ++i){
                pt = (*it)->Point_of_hse(i);
                sl = (STATE*)left_state(pt);
                for (int j = 0; j < m_dim; ++j)
                    Coords(pt)[j] =  sl->x_old[j];
            }
        }
}

void CollisionSolver::turnOffImpZone(){s_detImpZone = false;}
void CollisionSolver::turnOnImpZone(){s_detImpZone = true;}
bool CollisionSolver::getImpZoneStatus(){ return s_detImpZone;}
//this function is needed if collision still happens
//after several iterations;
void CollisionSolver::computeImpactZone()
{
	bool is_collision = true;
        int numZones = 0;
	int niter = 0;
	int cd_pair = 0;

        std::cout<<"Starting compute Impact Zone: "<<std::endl;
	turnOnImpZone();
	//makeSet(hseList);
        while(is_collision){
            is_collision = false;

	    //start UF alogrithm
	    //merge four pts if collision happens

	    start_clock("cgal_impactzone");
	    CGAL::box_self_intersection_d(hseList.begin(),
                  hseList.end(),reportCollision(is_collision,cd_pair,this),
                  traitsForCollision());
	    stop_clock("cgal_impactzone");
 
            updateAverageVelocity();

	    updateImpactZoneVelocity(numZones);
            std::cout <<"    #"<<niter++ << ": " << cd_pair 
                      << " pair of collision tris" << std::endl;
	    std::cout <<"     "<< numZones
		      <<" zones of impact" << std::endl;
        }
	turnOffImpZone();
	return;
}

void CollisionSolver::updateImpactZoneVelocityForRG()
{
	POINT* pt;
	unsortHseList(hseList);

	for (std::vector<CD_HSE*>::iterator it = hseList.begin();
	     it < hseList.end(); ++it){
	    for (int i = 0; i < (*it)->num_pts(); ++i){
		pt = (*it)->Point_of_hse(i);
		//skip traversed or isolated pts
		if (sorted(pt) ||
		    weight(findSet(pt)) == 1) continue;
		else if (!isMovableRigidBody(pt))
		{
		    sorted(pt) = YES;
		    continue;
		}
		else
		    updateImpactListVelocity(findSet(pt));
	    }
	}	
}

void CollisionSolver::updateImpactZoneVelocity(int &nZones)
{
	POINT* pt;
	int numZones = 0;
	unsortHseList(hseList);

	for (std::vector<CD_HSE*>::iterator it = hseList.begin();
	     it < hseList.end(); ++it){
	    for (int i = 0; i < (*it)->num_pts(); ++i){
		pt = (*it)->Point_of_hse(i);
		//skip traversed or isolated pts
		if (sorted(pt) ||
		    weight(findSet(pt)) == 1) continue;
		else{
		    updateImpactListVelocity(findSet(pt));
		    numZones++;
		}
	    }
	}	
	nZones = numZones;
}

void CollisionSolver::setTraitsDimension(){
	traitsForProximity::m_dim = m_dim;
	traitsForCollision::m_dim = m_dim;
}
//resolve collision in the input tris list
void CollisionSolver::resolveCollision()
{
	//catch floating point exception: nan/inf
	feenableexcept(FE_INVALID | FE_OVERFLOW);

	setTraitsDimension();

	start_clock("computeAverageVelocity");
	//compute average velocity
	computeAverageVelocity();
	stop_clock("computeAverageVelocity");

	start_clock("detectProximity");
	//test proximity for tris on surfaces
	//or bond on curves
	detectProximity();
	stop_clock("detectProximity");

	if (debugging("collision"))
	    printDebugVariable();

	start_clock("detectCollision");
	//test collision for tri-tri
	//or bond-tri or bond-bond
	detectCollision();
	stop_clock("detectCollision");

	if (debugging("collision"))
	    printDebugVariable();
	
	start_clock("detectDomainBoundaryCollision");
	detectDomainBoundaryCollision();
	stop_clock("detectDomainBoundaryCollision");

	start_clock("updateFinalPosition");
	//update position using average velocity
	updateFinalPosition();
	stop_clock("updateFinalPosition");

	reduceSuperelast();
	
	start_clock("updateFinalVelocity");
	//update velocity using average velocity
	updateFinalVelocity();
	stop_clock("updateFinalVelocity");
}

void CollisionSolver::detectProximity()
{
	int num_pairs = 0;
	CGAL::box_self_intersection_d(hseList.begin(),hseList.end(),
                                     reportProximity(num_pairs,this),
				     traitsForProximity());
	updateAverageVelocity();
	if (debugging("collision"))
	std::cout << num_pairs << " pair of proximity" << std::endl;
}

void CollisionSolver::detectCollision()
{
	bool is_collision = true; 
	const int MAX_ITER = 5;
	int niter = 1;
	int cd_pair = 0;

	std::cout<<"Starting collision handling: "<<std::endl;
	//record if has an actual collision
	//this is useful for adpative dt
	int cd_count = 0;
	setHasCollision(false);
	//
	while(is_collision){
	    is_collision = false;
	    start_clock("cgal_collision");
	    CGAL::box_self_intersection_d(hseList.begin(),
		  hseList.end(),reportCollision(is_collision,cd_pair,this),
		  traitsForCollision());
	    stop_clock("cgal_collision");
	
	    if (cd_count++ == 0 && is_collision) 
		setHasCollision(true);

	    updateAverageVelocity();
	    std::cout<<"    #"<<niter << ": " << cd_pair 
		     << " pair of collision tris" << std::endl;
	    if (++niter > MAX_ITER) break;
	}
	start_clock("computeImpactZone");
	if (is_collision) 
	    computeImpactZone();
	stop_clock("computeImpactZone");
}

//helper function to detect a collision between 
//a moving point and a moving triangle
//or between two moving edges 
extern void createImpZone(POINT* pts[], int num, bool first){
	for (int i = 0; i < num; ++i)
	{
	    for (int j = 0; j < i; ++j)
	    {
	        if ((!first) && (isMovableRigidBody(pts[i]) || 
				 isMovableRigidBody(pts[j])))
		    continue;
		mergePoint(pts[i],pts[j]); 
	    }
	}
}

bool CollisionSolver::reduceSuperelastOnce(int& num_edges)
{
	double dt = getTimeStepSize();
	const double superelasTol = 0.10;
	bool has_superelas = false;
	num_edges = 0;
	for (unsigned i = 0; i < hseList.size(); ++i){
	    CD_HSE* hse = hseList[i];
	    int np = hse->num_pts();
	    if (isRigidBody(hse)) continue;
	    for (int j = 0; j < ((np == 2) ? 1 : np); ++j){
		POINT* p[2];
		STATE* sl[2];
		p[0] = hse->Point_of_hse(j%np);	
		p[1] = hse->Point_of_hse((j+1)%np);
		sl[0]= (STATE*)left_state(p[0]);
		sl[1]= (STATE*)left_state(p[1]);

		double x_cand[2][3];
		for (int k = 0; k < 2; ++k){
		    double tmp[3];
		    scalarMult(dt,sl[k]->avgVel,tmp);
		    addVec(sl[k]->x_old,tmp,x_cand[k]);
		}	
		double len_new = distance_between_positions(x_cand[0],x_cand[1],3);
		double len_old = distance_between_positions(sl[0]->x_old,sl[1]->x_old,3);
		double len0;
		if (CD_TRI* cd_tri = dynamic_cast<CD_TRI*>(hse))
		    len0 = cd_tri->m_tri->side_length0[j];
		else if (CD_BOND* cd_bond = dynamic_cast<CD_BOND*>(hse))
		    len0 = cd_bond->m_bond->length0;
		else{
		    std::cout<<"Unknown type"<<std::endl;
		    clean_up(ERROR);
		}
		double vec[3], v_rel[3];

		minusVec(sl[0]->x_old,sl[1]->x_old,vec);
	        minusVec(sl[0]->avgVel,sl[1]->avgVel,v_rel);
		if (len_old > ROUND_EPS && len_new > ROUND_EPS)
		{
		    scalarMult(1/len_old,vec,vec); //normalize
		    double strain_rate = (len_new-len_old)/len_old;
		    double strain = (len_new-len0)/len0;
		    if (fabs(strain) > superelasTol || fabs(strain_rate) > superelasTol){
			double v_tmp[3];
			addVec(sl[0]->avgVel,sl[1]->avgVel,v_tmp);
                        scalarMult(0.5,v_tmp,v_tmp);
                        memcpy((void*)sl[0]->avgVel,(void*)v_tmp,3*sizeof(double)); 
                        memcpy((void*)sl[1]->avgVel,(void*)v_tmp,3*sizeof(double));
			num_edges ++;
		        has_superelas = true;
		    }
		}
		else
		{
		    	printf("Warning: len0 = %e, len_new = %e, len_old = %e\n",
			    len0,len_new,len_old);
			printf("p0 = %p, p1 = %p\n",(void*)p[0],(void*)p[1]);
			printf("x_old[0] = [%f %f %f]\n",sl[0]->x_old[0],sl[0]->x_old[1],sl[0]->x_old[2]);
			printf("avgVel[0] = [%f %f %f]\n",sl[0]->avgVel[0],sl[0]->avgVel[1],sl[0]->avgVel[2]);
			printf("x_old[1] = [%f %f %f]\n",sl[1]->x_old[0],sl[1]->x_old[1],sl[1]->x_old[2]);
			printf("avgVel[1] = [%f %f %f]\n",sl[1]->avgVel[0],sl[1]->avgVel[1],sl[1]->avgVel[2]);
		        double v_tmp[3];
                        addVec(sl[0]->avgVel,sl[1]->avgVel,v_tmp);
                        scalarMult(0.5,v_tmp,v_tmp);
                        memcpy((void*)sl[0]->avgVel,(void*)v_tmp,3*sizeof(double)); 
                        memcpy((void*)sl[1]->avgVel,(void*)v_tmp,3*sizeof(double));
			num_edges ++;
		        has_superelas = true;
		}	
	    }
	}
	return has_superelas;
}

void CollisionSolver::updateFinalPosition()
{
	POINT* pt;
	STATE* sl;
	double dt = getTimeStepSize();

	//#pragma omp parallel for private(sl,pt)
	for (std::vector<CD_HSE*>::iterator it = hseList.begin();
	     it < hseList.end(); ++it)
	{
	    for (int i = 0; i < (*it)->num_pts(); ++i){
		pt = (*it)->Point_of_hse(i);
		sl = (STATE*)left_state(pt);
	    	for (int j = 0; j < 3; ++j)
		{
		    Coords(pt)[j] = sl->x_old[j]+sl->avgVel[j]*dt;
		    if (std::isnan(Coords(pt)[j]))
			printf("nan coords, x_old = %f, avgVel = %f\n",
				sl->x_old[j],sl->avgVel[j]);
		}
	    }
	}
}

void CollisionSolver::reduceSuperelast()
{
	bool has_superelas = true;
	int niter = 0, num_edges;
	const int max_iter = 10;
	while(has_superelas && niter++ < max_iter){
	    has_superelas = reduceSuperelastOnce(num_edges);
	}
	if (debugging("collision"))
	printf("    %d edges are over strain limit after %d iterations\n",num_edges,niter);
}

void CollisionSolver::updateFinalVelocity()
{
	//TODO:avgVel is actually the velocity at t(n+1/2)
	//need to call spring solver to get velocity at t(n+1)
	//for simplicity now set v(n+1) = v(n+1/2)
	POINT* pt;
	STATE* sl;
	double dt = getTimeStepSize();
	//#pragma omp parallel for private(pt,sl)
	for (std::vector<CD_HSE*>::iterator it = hseList.begin();
             it < hseList.end(); ++it)
        {
            for (int i = 0; i < (*it)->num_pts(); ++i){
                pt = (*it)->Point_of_hse(i);
                sl = (STATE*)left_state(pt);
		if (!sl->has_collsn) 
		    continue;
                for (int j = 0; j < 3; ++j){
                    pt->vel[j] = sl->avgVel[j];
                    sl->vel[j] = sl->avgVel[j];
		    if (std::isnan(pt->vel[j]))
			printf("nan vel and avgVel\n");
		}
            }
        }
	updateFinalForRG();
}

void CollisionSolver::updateFinalForRG()
{
	POINT* pt;
        STATE* sl;
        double dt = getTimeStepSize();
	std::vector<int> mrg;
	std::map<int, bool> visited;

	for (std::vector<CD_HSE*>::iterator it = hseList.begin();
             it < hseList.end(); ++it)
        {
            for (int i = 0; i < (*it)->num_pts(); ++i){
                pt = (*it)->Point_of_hse(i);
                sl = (STATE*)left_state(pt);
                if (!isMovableRigidBody(pt)) continue;
		int rg_index = body_index(pt->hs);
                if (sl->has_collsn && 
		    std::find(mrg.begin(), mrg.end(), rg_index) == mrg.end())
                {
                    mrg.push_back(rg_index); 
                    for (int j = 0; j < 3; ++j)
                    {
                        center_of_mass_velo(pt->hs)[j] = sl->avgVel[j];
                        center_of_mass(pt->hs)[j] = sl->avgVel[j] * dt + 
                                                (mrg_com[rg_index])[j];
                    }
		    visited[rg_index] = false;
		    if (debugging("rigid_body"))
		    {
			printf("After collision handling: \n");
			printf("Body Index: %d\n", rg_index);
			printf("center_of_mass = %f %f %f\n", 
				center_of_mass(pt->hs)[0], 
				center_of_mass(pt->hs)[1], 
				center_of_mass(pt->hs)[2]);
			printf("center_of_mass_velo = %f %f %f\n", 
				center_of_mass_velo(pt->hs)[0], 
				center_of_mass_velo(pt->hs)[1], 
				center_of_mass_velo(pt->hs)[2]);
		    }
                }
		if ((visited.count(rg_index) == 0) || (!visited[rg_index]))
		{
		    double* com = center_of_mass(pt->hs);
		    mrg_com[rg_index] = std::vector<double>(com, com+3);
		    visited[rg_index] = true;
		}
            }
        }
}

void CollisionSolver::updateAverageVelocity()
{
	POINT *p;
	STATE *sl;
	double maxSpeed = 0;
	double* maxVel = NULL;

#ifdef __VTK__
	if (debugging("CollisionImpulse")){
       	  char fname[200] = "vtk_test";
       	  static int count = 0;
	  updateFinalPosition();
       	  if (create_directory(fname,NO)){
       	  	sprintf(fname,"%s/surf-%03d.vtp",fname,count++);
       		vtkplotVectorSurface(hseList,fname);	
       	  }
	}
#endif
	unsortHseList(hseList);
	for (unsigned i = 0; i < hseList.size(); ++i)
	{
	    CD_HSE* hse = hseList[i];
	    int np = hse->num_pts(); 

	    for (int j = 0; j < np; ++j)
	    {
		p = hse->Point_of_hse(j);
		if (isStaticRigidBody(p)) continue;
		if (sorted(p)) continue;
		sl = (STATE*)left_state(p);
		if (sl->collsn_num > 0)
		{
		    sl->has_collsn = true;
		    for (int k = 0; k < 3; ++k)
		    {
			sl->avgVel[k] += (sl->collsnImpulse[k] + sl->friction[k])/sl->collsn_num;
			if (std::isinf(sl->avgVel[k]) || std::isnan(sl->avgVel[k])) 
			{
			    printf("inf/nan vel[%d]: impulse = %f, friction = %f, collsn_num = %d\n",
				k,sl->collsnImpulse[k],sl->friction[k],sl->collsn_num);
			    clean_up(ERROR);
			}
			//reset impulse and fricition to 0
			//collision handling will recalculate the impulse
			sl->collsnImpulse[k] = sl->friction[k] = 0.0;
		    }
		    sl->collsn_num = 0;
		}
		/* test for RG */
		if (sl->collsn_num_RG > 0)
		{
		    sl->has_collsn = true;
		    for (int k = 0; k < 3; ++k)
			sl->avgVel[k] += sl->collsnImpulse_RG[k]/
							sl->collsn_num_RG;
		    sl->collsn_num_RG = 0;
		}

		if (debugging("collision")){
		    //debugging: print largest speed
		    double speed = Mag3d(sl->avgVel);
		    if (speed > maxSpeed)
			maxVel = sl->avgVel;
		}
		sorted(p) = YES;
	    }
	}
	if (getTimeStepSize() > 0.0)
	    updateImpactZoneVelocityForRG(); // test for moving objects
	if (debugging("collision"))
	if (maxVel != NULL)
	    printf("    max velocity = [%f %f %f]\n",maxVel[0],maxVel[1],maxVel[2]);
	if (debugging("collision"))
	    printDebugVariable();
}

bool CollisionSolver::isCollision(const CD_HSE* a, const CD_HSE* b){
	const CD_BOND *cd_b1, *cd_b2;
	const CD_TRI  *cd_t1, *cd_t2;
	double h = CollisionSolver3d::getRoundingTolerance();
	if ((cd_t1 = dynamic_cast<const CD_TRI*>(a)) && 
	    (cd_t2 = dynamic_cast<const CD_TRI*>(b)))
	{
	    TRI* t1 = cd_t1->m_tri;
	    TRI* t2 = cd_t2->m_tri;
	    if ((t1->surf == t2->surf) && isRigidBody(a))
		return false;
	    return MovingTriToTri(t1,t2,h);
	}
	else if ((cd_b1 = dynamic_cast<const CD_BOND*>(a)) && 
	         (cd_b2 = dynamic_cast<const CD_BOND*>(b)))
	{
	    BOND* b1 = cd_b1->m_bond;
	    BOND* b2 = cd_b2->m_bond;
	    return MovingBondToBond(b1,b2,h);
	}
	else if ((cd_b1 = dynamic_cast<const CD_BOND*>(a)) &&
		 (cd_t1 = dynamic_cast<const CD_TRI*>(b)))
	{
	    BOND* b1 = cd_b1->m_bond;
	    TRI* t1  = cd_t1->m_tri;
	    return MovingTriToBond(t1,b1,h);
	}
	else if ((cd_t1 = dynamic_cast<const CD_TRI*>(a)) &&
                 (cd_b1 = dynamic_cast<const CD_BOND*>(b)))
	{
	    BOND* b1 = cd_b1->m_bond;
	    TRI* t1  = cd_t1->m_tri;
	    return MovingTriToBond(t1,b1,h);
	}
	else
	{
	    std::cout<<"This case has not been implemented"<<std::endl;
	    clean_up(ERROR);
	}
	return false;
}

bool CollisionSolver::isProximity(const CD_HSE* a, const CD_HSE* b){
	const CD_BOND *cd_b1, *cd_b2;
	const CD_TRI  *cd_t1, *cd_t2;
	double h = CollisionSolver3d::getFabricThickness();

	if ((cd_t1 = dynamic_cast<const CD_TRI*>(a)) && 
	    (cd_t2 = dynamic_cast<const CD_TRI*>(b)))
	{
	    TRI* t1 = cd_t1->m_tri;
	    TRI* t2 = cd_t2->m_tri;
	    if ((t1->surf == t2->surf) && isRigidBody(a))
		return false;
	    return TriToTri(t1,t2,h);
	}
	else if ((cd_b1 = dynamic_cast<const CD_BOND*>(a)) && 
	         (cd_b2 = dynamic_cast<const CD_BOND*>(b)))
	{
	    BOND* b1 = cd_b1->m_bond;
	    BOND* b2 = cd_b2->m_bond;
	    return BondToBond(b1,b2,h);
	}
	else if ((cd_b1 = dynamic_cast<const CD_BOND*>(a)) &&
		 (cd_t1 = dynamic_cast<const CD_TRI*>(b)))
	{
	    BOND* b1 = cd_b1->m_bond;
	    TRI* t1  = cd_t1->m_tri;
	    return TriToBond(t1,b1,h);
	}
	else if ((cd_t1 = dynamic_cast<const CD_TRI*>(a)) &&
                 (cd_b1 = dynamic_cast<const CD_BOND*>(b)))
	{
	    BOND* b1 = cd_b1->m_bond;
	    TRI* t1  = cd_t1->m_tri;
	    return TriToBond(t1,b1,h);
	}
	else
	{
	    std::cout<<"This case has not been implemented"<<std::endl;
	    clean_up(ERROR);
	}
	return false;
}

void CollisionSolver::printDebugVariable(){
	std::cout << "Enter EdgeToEdge " << edg_to_edg 
		  << " times"<< std::endl;
	std::cout << "Enter PointToTri " << pt_to_tri 
		  << " times"<< std::endl;
	std::cout << "Enter isCoplanar " << is_coplanar
		  << " times"<< std::endl;
	moving_edg_to_edg = moving_pt_to_tri = is_coplanar = 0;
	edg_to_edg = pt_to_tri = 0;
}

/********************************
* implementation for CD_HSE     *
*********************************/
double CD_BOND::max_static_coord(int dim){
    return std::max(Coords(m_bond->start)[dim],
		    Coords(m_bond->end)[dim]);
}

double CD_BOND::min_static_coord(int dim){
    return std::min(Coords(m_bond->start)[dim],
		    Coords(m_bond->end)[dim]);
}

double CD_BOND::max_moving_coord(int dim,double dt){
    double ans = -HUGE;
    for (int i = 0; i < 2; ++i){
	POINT* pt = (i == 0)? m_bond->start : m_bond->end;
	STATE* sl = (STATE*)left_state(pt);
	ans = std::max(ans,sl->x_old[dim]);
	ans = std::max(ans,sl->x_old[dim]+sl->avgVel[dim]*dt); 
    }    
    return ans;
}

double CD_BOND::min_moving_coord(int dim,double dt){
    double ans = HUGE;
    for (int i = 0; i < 2; ++i){
	POINT* pt = (i == 0)? m_bond->start : m_bond->end;
	STATE* sl = (STATE*)left_state(pt);
	ans = std::min(ans,sl->x_old[dim]);
	ans = std::min(ans,sl->x_old[dim]+sl->avgVel[dim]*dt); 
    }    
    return ans;
}

POINT* CD_BOND::Point_of_hse(int i) const{
    if (i >= num_pts())
	return NULL;
    else
        return (i == 0) ? m_bond->start : 
			  m_bond->end;
}

double CD_TRI::max_static_coord(int dim){
    double ans = -HUGE;
    for (int i = 0; i < 3; ++i){
	POINT* pt = Point_of_tri(m_tri)[i];
	STATE* sl = (STATE*)left_state(pt);
	ans = std::max(sl->x_old[dim],ans);
    }
    return ans;
}

double CD_TRI::min_static_coord(int dim){
    double ans = HUGE;
    for (int i = 0; i < 3; ++i){
	POINT* pt = Point_of_tri(m_tri)[i];
	STATE* sl = (STATE*)left_state(pt);
	ans = std::min(sl->x_old[dim],ans);
    }
    return ans;
}

double CD_TRI::max_moving_coord(int dim,double dt){
    double ans = -HUGE;
    for (int i = 0; i < 3; ++i){
	POINT* pt = Point_of_tri(m_tri)[i];
	STATE* sl = (STATE*)left_state(pt);
	ans = std::max(ans,sl->x_old[dim]);
	ans = std::max(ans,sl->x_old[dim]+sl->avgVel[dim]*dt);
    }
    return ans;
}

double CD_TRI::min_moving_coord(int dim,double dt){
    double ans = HUGE;
    for (int i = 0; i < 3; ++i){
	POINT* pt = Point_of_tri(m_tri)[i];
	STATE* sl = (STATE*)left_state(pt);
	ans = std::min(ans,sl->x_old[dim]);
	ans = std::min(ans,sl->x_old[dim]+sl->avgVel[dim]*dt);
    }
    return ans;
}

POINT* CD_TRI::Point_of_hse(int i) const{
    if (i >= num_pts())
	return NULL;
    else
        return Point_of_tri(m_tri)[i];
}

/*******************************
* utility functions start here *
*******************************/
/* The followings are helper functions for vector operations. */
void Pts2Vec(const POINT* p1, const POINT* p2, double* v){
	for (int i = 0; i < 3; ++i)	
	    v[i] = Coords(p1)[i] - Coords(p2)[i];
}

double distBetweenCoords(double* v1, double* v2)
{
	double dist = 0.0;
	for (int i = 0; i < 3; ++i){
		dist += sqr(v1[i]-v2[i]);
	}
	return std::sqrt(dist);
}

void addVec(double* v1, double* v2, double* ans)
{
	for (int i = 0; i < 3; ++i)
	    ans[i] = v1[i]+v2[i];
}

void minusVec(double* v1, double* v2, double* ans)
{
	for (int i = 0; i < 3; ++i)
	    ans[i] = v1[i]-v2[i];
}

void scalarMult(double a,double* v, double* ans)
{
	for (int i = 0; i < 3; ++i)
            ans[i] = a*v[i];	
}

extern double myDet3d(double a[3][3]){
    return  a[0][0]*(a[1][1]*a[2][2] - a[2][1]*a[1][2]) 
	  - a[0][1]*(a[1][0]*a[2][2] - a[2][0]*a[1][2]) 
	  + a[0][2]*(a[1][0]*a[2][1] - a[2][0]*a[1][1]);
}

void unsortHseList(std::vector<CD_HSE*>& hseList){
	for (unsigned j = 0; j < hseList.size(); ++j)
	{
	    CD_HSE* hse = hseList[j];
	    int np = hse->num_pts();
	    for (int i = 0; i < np; ++i){
		sorted(hse->Point_of_hse(i)) = NO;
	    }
	}
}

//functions for UF alogrithm
int& weight(POINT* p){
	STATE* sl = (STATE*)left_state(p);
	return sl->impZone.num_pts;
}

inline POINT*& root(POINT* p){
	STATE* sl = (STATE*)left_state(p);
	return sl->impZone.root;
}

POINT*& next_pt(POINT* p){
	STATE* sl = (STATE*)left_state(p);
        return sl->impZone.next_pt;
}

inline POINT*& tail(POINT* p){
	STATE* sl = (STATE*)left_state(p);
	return sl->impZone.tail;
}

extern void makeSet(std::vector<CD_HSE*>& hseList){
	STATE* sl;
	POINT* pt;
	//#pragma omp parallel for private(sl,pt)
        for (std::vector<CD_HSE*>::iterator it = hseList.begin();
                it < hseList.end(); ++it){
            for (int i = 0; i < (*it)->num_pts(); ++i){
		pt = (*it)->Point_of_hse(i);
                sorted(pt) = NO;
		sl = (STATE*)left_state(pt);
		sl->impZone.next_pt = NULL;
		sl->impZone.tail = pt;
		sl->impZone.root = pt;
		sl->impZone.num_pts = 1;
            }
        }
}

POINT* findSet(POINT* p){
	if (root(p) != p)
		root(p) = findSet(root(p));
	return root(p);
}

void mergePoint(POINT* X, POINT* Y){
	POINT* PX = findSet(X);
	POINT* PY = findSet(Y);
	if (PX == PY) return;
	if (weight(PX) > weight(PY)){
	    //update root after merge
	    weight(PX) += weight(PY);
	    root(PY) = PX;
	    //link two list, update tail
	    next_pt(tail(PX)) = PY;
	    tail(PX) = tail(PY); 
	}
	else{
	    //update root after merge
	    weight(PY) += weight(PX);
	    root(PX) = PY;
	    //link two list, update tail
	    next_pt(tail(PY)) = PX;
	    tail(PY) = tail(PX); 
	}
}
//end of UF functions

void printPointList(POINT** plist,const int n){
	for (int i = 0; i < n; ++i){
	    printf("pt[%d] = [%f %f %f]\n",i,Coords(plist[i])[0],
		Coords(plist[i])[1],Coords(plist[i])[2]);
	}
}

bool isStaticRigidBody(const POINT* p){
    STATE* sl = (STATE*)left_state(p);
    return sl->is_fixed;
}

bool isStaticRigidBody(const CD_HSE* hse){
    for (int i = 0; i < hse->num_pts(); ++i)
   	if (isStaticRigidBody(hse->Point_of_hse(i)))
	    return true;
    return false;
}

bool isMovableRigidBody(const POINT* p){
    STATE* sl = (STATE*)left_state(p);
    return sl->is_movableRG;
}

bool isMovableRigidBody(const CD_HSE* hse){
    for (int i = 0; i < hse->num_pts(); ++i)
        if (isMovableRigidBody(hse->Point_of_hse(i)))
            return true;
    return false;
}

bool isRigidBody(const POINT* p){
    return isStaticRigidBody(p) || isMovableRigidBody(p);
}

bool isRigidBody(const CD_HSE* hse){
    return isStaticRigidBody(hse) || isMovableRigidBody(hse);
}

extern void SpreadImpactZoneImpulse(
        POINT* p,
        double impulse,
        double* nor)
{
        POINT* root = findSet(p);
        while (root)
        {
            STATE *sl = (STATE*)left_state(root);
            for (int i = 0; i < 3; ++i)
                sl->collsnImpulse_RG[i] += impulse * nor[i];
            sl->collsn_num_RG += 1;
            root = next_pt(root);
        }
}
