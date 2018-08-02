#include "collid.h"
static void unsort_surface_point(SURFACE *surf);
static bool MovingPointToTri(POINT**,double);
static bool MovingEdgeToEdge(POINT**,double);
static bool PointToTri(POINT**,double,double root = 0.0);
static bool EdgeToEdge(POINT**,double,double root = 0.0);
static bool isCoplanar(POINT**,double,double*);
static void EdgeToEdgeImpulse(POINT**, double*, double, double, double,double);
static void PointToTriImpulse(POINT**, double*, double*, double,double);

//functions in CollisionSolver3d
void CollisionSolver3d::assembleFromInterface(
	const INTERFACE* intfc,const double dt)
{
	//assemble tris list from input intfc
	//this function should be called before
	//spring interior dynamics computed
	SURFACE** s;
	CURVE** c;
	TRI *tri;
	BOND *b;
	int n_tri = 0, n_bond = 0;
	setTimeStepSize(dt);
	clearHseList();
	intfc_surface_loop(intfc,s)
	{
	    if (is_bdry(*s)) continue;
	    unsort_surface_point(*s);
	    surf_tri_loop(*s,tri)
	    {
	        hseList.push_back(new CD_TRI(tri));
		n_tri++;
	    }
	}
	intfc_curve_loop(intfc,c)
	{
	    if (hsbdry_type(*c) != STRING_HSBDRY) continue; 
	    curve_bond_loop(*c,b)
	    {
		hseList.push_back(new CD_BOND(b,m_dim));
		n_bond++;
	    }
	}
	makeSet(hseList);
	createImpZoneForRG(intfc);
	setDomainBoundary(intfc->table->rect_grid.L,
			  intfc->table->rect_grid.U);
	if (debugging("collision")){
	    printf("%d num of tris, %d num of bonds\n",n_tri,n_bond);
	    printf("%lu number of elements is assembled\n",hseList.size());
	}
}

void CollisionSolver3d::createImpZoneForRG(const INTERFACE* intfc)
{	// test function for creating impact zone for each movable RG
	SURFACE** s;
	TRI* tri;

	intfc_surface_loop(intfc, s)
	{
	    if (is_bdry(*s)) continue;
	    if (!isMovableRigidBody(Point_of_tri(first_tri(*s))[0])) continue;
	    surf_tri_loop(*s, tri)
	    {
		createImpZone(Point_of_tri(tri), 3, YES);
	    }
	}
}

void CollisionSolver3d::updateImpactListVelocity(POINT* head){
	STATE* sl = NULL;
	POINT* p = head;
	double m = getPointMass();
	double x_cm[3] = {0.0}, v_cm[3] = {0.0};
	double L[3] = {0.0}; //angular momentum
	double I[3][3] = {0.0}; //inertia tensor
	double tmp[3][3];
	int num_pts = 0;

	while(p){
		num_pts++; //debug
		sorted(p) = YES;
		sl = (STATE*)left_state(p);
		for (int i = 0; i < m_dim; ++i){
		    x_cm[i] += sl->x_old[i]; 
		    v_cm[i] += sl->avgVel[i];
		}
                p = next_pt(p);
        }
	if (debugging("collision"))
	    printf("%d number of points in this zone\n",num_pts);
	//compute center and veclocity of impact Zone
	for (int i = 0; i < m_dim; ++i){
	    x_cm[i] /= num_pts;
	    v_cm[i] /= num_pts;
	}

	//compute angular momentum
	p = head;
	while(p){
	    double dx[3], dv[3], Li[3];
	    sl = (STATE*)left_state(p);
	    minusVec(sl->x_old,x_cm,dx);
	    minusVec(sl->avgVel,v_cm,dv); 	
	    Cross3d(dx,dv,Li);
	    scalarMult(m,Li,Li);
	    addVec(Li,L,L);    
	    p = next_pt(p);
	}
	//compute Inertia tensor
	p = head;
	while(p){
	    double dx[3], mag_dx = 0.0;
	    sl = (STATE*)left_state(p);
	    minusVec(sl->x_old,x_cm,dx);
	    mag_dx = Mag3d(dx);
	    for (int i = 0; i < 3; ++i)
	    for (int j = 0; j < 3; ++j){
		tmp[i][j] = -dx[i]*dx[j];
		if (i == j)
		    tmp[i][j] += mag_dx*mag_dx; 
	 	I[i][j] += tmp[i][j]*m;
	    } 
	    p = next_pt(p);
	}

	//compute angular velocity w: I*w = L;
	double w[3], mag_w = 0;
	for (int i = 0; i < 3; ++i){
	    memcpy(tmp,I,9*sizeof(double));
	    for (int j = 0; j < 3; j++)
		tmp[j][i] = L[j];
	    if (myDet3d(I) < ROUND_EPS)
		w[i] = 0.0;
	    else
	        w[i] = myDet3d(tmp)/myDet3d(I);
	}
	mag_w = Mag3d(w);
	
	//compute average velocity for each point
	double dt = getTimeStepSize();
	p = head;
        while(p){
	    if (isStaticRigidBody(p)) {
		p = next_pt(p);
		continue;
	    }
	    double x_new[3],dx[3];
	    double xF[3], xR[3];
	    sl = (STATE*)left_state(p);
	    minusVec(sl->x_old,x_cm,dx);
	    double wxR[3],tmpV[3];
	    if (mag_w < ROUND_EPS){
	        for (int i = 0; i < 3; ++i){
		    xF[i] = dx[i];
		    wxR[i] = 0.0;
		}
		minusVec(dx,xF,xR);
	    }
	    else{
	        scalarMult(Dot3d(dx,w)/Dot3d(w,w),w,xF);
	        minusVec(dx,xF,xR);
	        scalarMult(sin(dt*mag_w)/mag_w,w,tmpV);
	        Cross3d(tmpV,xR,wxR);
	    }
	    for (int i = 0; i < 3; ++i)
	    {
	    	x_new[i] = x_cm[i] + dt*v_cm[i]
			 + xF[i] + cos(dt*mag_w)*xR[i]
			 + wxR[i];
		sl = (STATE*)left_state(p);
		sl->avgVel[i] = (x_new[i] - sl->x_old[i])/dt;
	    	if (std::isnan(sl->avgVel[i]))
		{ 
			printf("coords[3], vel[3]\n");
			p = head;
			while(p){
			    sl = (STATE*)left_state(p);
			    printf("%f %f %f %f %f %f;\n",
				sl->x_old[0],sl->x_old[1],sl->x_old[2],
				sl->avgVel[0],sl->avgVel[1],sl->avgVel[2]);
			    p = next_pt(p);
			}
			printf("num_pts = %d, weight = %d\n",
			num_pts,weight(head));
			printf("nan vel, w = %f, mag_w = %f\n",
			w[i],mag_w);
			printf("L = [%f %f %f]\n",L[0],L[1],L[2]);
			printf("I = [%f %f %f;  %f %f %f; %f %f %f]\n",
			I[0][0],I[0][1],I[0][2],I[1][0],I[1][1],I[1][2],
			I[2][0],I[2][1],I[2][2]);
			printf("xF = %f %f %f, xR = %f %f %f\n",
			xF[0],xF[1],xF[2],xR[0],xR[1],xR[2]);
			clean_up(ERROR);
		}
	    }
	    p = next_pt(p);
	}
	//done!!!
}

//helper function to detect collision between elements 
bool CollisionSolver3d::MovingTriToBond(const TRI* tri,const BOND* bd, double h){
	bool is_detImpZone = CollisionSolver3d::getImpZoneStatus();
	POINT* pts[4];
	bool status = false;

	/* do not consider bond point that is a tri vertex */
	for (int i = 0; i < 3; ++i)
	{
	    if (Point_of_tri(tri)[i] == bd->start ||
		Point_of_tri(tri)[i] == bd->end)
		return false;
	}

	for (int i = 0; i < 3; ++i)
	    pts[i] = Point_of_tri(tri)[i];

	/* detect collision of start point of bond w.r.t to tri */
	pts[3] = bd->start;
	if (MovingPointToTri(pts, h)) status = true;
	if (status && is_detImpZone)
	    createImpZone(pts,4);

	/* detect collision of end point of bond to w.r.t. tri */
	pts[3] = bd->end;
	if (MovingPointToTri(pts, h)) status = true;
	if (status && is_detImpZone)
	    createImpZone(pts,4);

	/* detect collision of each of tri edge w.r.t to bond */
	pts[2] = bd->start;
	pts[3] = bd->end;
	for (int i = 0; i < 3; ++i)
	{
	    pts[0] = Point_of_tri(tri)[i];
	    pts[1] = Point_of_tri(tri)[(i+1)%3];
	    if (MovingEdgeToEdge(pts,h)) status = true;
	    if (status && is_detImpZone)
		createImpZone(pts,4);
	}

	return status;
}

bool CollisionSolver3d::MovingBondToBond(const BOND* b1, const BOND* b2, double h){
	bool is_detImpZone = CollisionSolver3d::getImpZoneStatus();
	POINT* pts[4];
	bool status = false;

	pts[0] = b1->start;
	pts[1] = b1->end;
	pts[2] = b2->start;
	pts[3] = b2->end;

	/* do not consider two bonds that share a common point */
	for (int i = 0; i < 4; ++i)
	{
	    for (int j = i + 1; j < 4; ++j)
	    {
		if (pts[i] == pts[j])
		    return false;
	    }
	}

	/* detect collision between two bonds */	
	if(MovingEdgeToEdge(pts,h)) status = true;
	if (status && is_detImpZone)
	    createImpZone(pts,4);
	
	return status;
}

bool CollisionSolver3d::MovingTriToTri(const TRI* a,const TRI* b, double h)
{
	bool is_detImpZone = CollisionSolver3d::getImpZoneStatus();
	POINT* pts[4];
	bool status = false;
	for (int i = 0; i < 3; ++i)
	for (int j = 0; j < 3; ++j)
	{
	    if (Point_of_tri(a)[i] == Point_of_tri(b)[j])
		return false;
	}
	//detect point to tri collision
	for (int k = 0; k < 2; ++k)
	for (int i = 0; i < 3; ++i){
	    const TRI* tmp_tri1 = (k == 0) ? a : b;
	    const TRI* tmp_tri2 = (k == 0) ? b : a;
	    for (int j = 0; j < 3; ++j)
	    	pts[j] = Point_of_tri(tmp_tri1)[j];
	    pts[3] = Point_of_tri(tmp_tri2)[i];

	    //we do not consider point against 
	    //a triangle that cotains it
	    if (pts[3] == pts[0] || pts[3] == pts[1] ||
		pts[3] == pts[2])
		continue; 

	    if(MovingPointToTri(pts,h)) status = true;
	    if (status && is_detImpZone)
		createImpZone(pts,4);
	}

	//detect edge to edge collision
	for (int i = 0; i < 3; ++i)
        {
            pts[0] = Point_of_tri(a)[i];
            pts[1] = Point_of_tri(a)[(i+1)%3];
            for (int j = 0; j < 3; ++j){
                pts[2] = Point_of_tri(b)[j];
                pts[3] = Point_of_tri(b)[(j+1)%3];
		
		//we do not consider edges share a point
		if (pts[0] == pts[2] || pts[0] == pts[3] ||
		    pts[1] == pts[2] || pts[1] == pts[3])
		    continue;

                if(MovingEdgeToEdge(pts,h)) status = true;
	        if (status && is_detImpZone)
		    createImpZone(pts,4);
	    }
        }
	return status;
}

static bool MovingPointToTri(POINT* pts[],const double h){

	double dt = CollisionSolver3d::getTimeStepSize();
	double roots[4] = {-1,-1,-1,dt};
	STATE* sl;
	if (isCoplanar(pts,dt,roots)){
	    for (int i = 0; i < 4; ++i){
		if (roots[i] < 0) continue;
		for (int j = 0; j < 4; ++j){
		    sl = (STATE*)left_state(pts[j]);
		    for (int k = 0; k < 3; ++k)
		        Coords(pts[j])[k] = sl->x_old[k]+roots[i]*sl->avgVel[k];
		}
		if (PointToTri(pts,h,roots[i])) 
		    return true;
	    }
	    return false;
	}
	else
	    return false;
}

static bool MovingEdgeToEdge(POINT* pts[],const double h){
	double dt = CollisionSolver3d::getTimeStepSize();
	double roots[4] = {-1,-1,-1,dt};
	STATE* sl;

	if (isCoplanar(pts,dt,roots)){
            for (int i = 0; i < 4; ++i){
                if (roots[i] < 0) continue;
                for (int j = 0; j < 4; ++j){
                    sl = (STATE*)left_state(pts[j]);
                    for (int k = 0; k < 3; ++k)
                        Coords(pts[j])[k] = sl->x_old[k]+roots[i]*sl->avgVel[k];
                }
                if (EdgeToEdge(pts,h,roots[i]))
                        return true;
            }
	    return false;
        }
	else
	    return false;
}

static void isCoplanarHelper(double* s[], double v[][3]) {
    	v[0][0] = s[0][0];         v[0][1] = s[0][1];         v[0][2] = s[0][2];
	v[1][0] = s[1][0]-s[0][0]; v[1][1] = s[1][1]-s[0][1]; v[1][2] = s[1][2]-s[0][2];
	v[2][0] = s[2][0]-s[0][0]; v[2][1] = s[2][1]-s[0][1]; v[2][2] = s[2][2]-s[0][2];
	v[3][0] = s[3][0]-s[0][0]; v[3][1] = s[3][1]-s[0][1]; v[3][2] = s[3][2]-s[0][2];
}
static bool isCoplanar(POINT* pts[], const double dt, double roots[])
{
	if (debugging("collision"))
	    CollisionSolver::is_coplanar++;

	double v[4][3] = {0}, x[4][3] = {0};
	double* tmp[4] = {0};
	//for performance, unrolling the loop
	tmp[0] = ((STATE*)left_state(pts[0]))->avgVel;
	tmp[1] = ((STATE*)left_state(pts[1]))->avgVel;
	tmp[2] = ((STATE*)left_state(pts[2]))->avgVel;
	tmp[3] = ((STATE*)left_state(pts[3]))->avgVel;
	isCoplanarHelper(tmp, v);

	tmp[0] = ((STATE*)left_state(pts[0]))->x_old;
	tmp[1] = ((STATE*)left_state(pts[1]))->x_old;
	tmp[2] = ((STATE*)left_state(pts[2]))->x_old;
	tmp[3] = ((STATE*)left_state(pts[3]))->x_old;
	isCoplanarHelper(tmp, x);

	//get roots "t" of a cubic equation
	//(x1+tv1)x(x2+tv2)*(x3+tv3) = 0
	//transform to at^3+bt^2+ct+d = 0
	double a, b, c, d;
	double vv[3], vx[3], xx[3];
	vv[0] = v[1][1]*v[2][2]-v[1][2]*v[2][1];
	vv[1] = v[1][0]*v[2][2]-v[1][2]*v[2][0];
	vv[2] = v[1][0]*v[2][1]-v[1][1]*v[2][0];
	
	vx[0] = v[1][1]*x[2][2]-v[1][2]*x[2][1]-v[2][1]*x[1][2]+v[2][2]*x[1][1];
	vx[1] = v[1][0]*x[2][2]-v[1][2]*x[2][0]-v[2][0]*x[1][2]+v[2][2]*x[1][0];
	vx[2] = v[1][0]*x[2][1]-v[1][1]*x[2][0]-v[2][0]*x[1][1]+v[2][1]*x[1][0];

	xx[0] = x[1][1]*x[2][2]-x[1][2]*x[2][1];
	xx[1] = x[1][0]*x[2][2]-x[1][2]*x[2][0];
	xx[2] = x[1][0]*x[2][1]-x[1][1]*x[2][0];

	a = v[3][0]*vv[0] - v[3][1]*vv[1] + v[3][2]*vv[2];

	b = x[3][0]*vv[0] - x[3][1]*vv[1] + x[3][2]*vv[2] + 
	    v[3][0]*vx[0] - v[3][1]*vx[1] + v[3][2]*vx[2];

	c = x[3][0]*vx[0] - x[3][1]*vx[1] + x[3][2]*vx[2] +
            v[3][0]*xx[0] - v[3][1]*xx[1] + v[3][2]*xx[2];

	d = x[3][0]*xx[0] - x[3][1]*xx[1] + x[3][2]*xx[2]; 
	//solve equation using method from "Art of Scientific Computing"
	//transform equation to t^3+at^2+bt+c = 0
	if (fabs(a) > MACH_EPS){
	    b /= a; c /= a; d /= a;
	    a = b; b = c; c = d;
	    double Q, R, theta;
	    double Q3, R2;
	    Q = (a*a-3*b)/9;
	    R = (2*a*a*a-9*a*b+27*c)/54;
	    Q3 = Q*Q*Q;
	    R2 = R*R;
	    if (R2 < Q3){
	        double Qsqrt = sqrt(Q);
		theta = acos(R/sqrt(Q3));
		roots[0] = -2*Qsqrt*cos(theta/3)-a/3;
		roots[1] = -2*Qsqrt*cos((theta+2*M_PI)/3)-a/3;
		roots[2] = -2*Qsqrt*cos((theta-2*M_PI)/3)-a/3;	
	    }
	    else{
		double A, B;
		double sgn = (R > 0) ? 1.0 : -1.0;
		A = -sgn*pow(fabs(R)+sqrt(R2-Q3),1.0/3.0);
		B = (fabs(A) < ROUND_EPS) ? 0.0 : Q/A;
		roots[0] = (A+B)-a/3.0;
		if (fabs(A-B) < ROUND_EPS)
		    roots[1] = roots[2] = -0.5*(A+B)-a/3.0; //multiple roots
	    }
	}
	else{
		a = b; b = c; c = d;
	   	double delta = b*b-4.0*a*c;
	   	if (fabs(a) > ROUND_EPS && delta > 0){
		    double delta_sqrt = sqrt(delta);
		    roots[0] = (-b+delta_sqrt)/(2.0*a);
	    	    roots[1] = (-b-delta_sqrt)/(2.0*a);
	   	}
		else if (fabs(a) < ROUND_EPS && fabs(b) > ROUND_EPS)
		{
		    roots[0] = -c/b;
	        }
	}
	//elimiate invalid roots;
	for (int i = 0; i < 3; ++i){
	        roots[i] = roots[i]-MACH_EPS;
	    	if (roots[i] < 0 || roots[i] > dt) 
		    roots[i] = -1;
	}
	//sort the roots
	if (roots[0] > roots[1])
	    std::swap(roots[0], roots[1]);
	if (roots[0] > roots[2])
	    std::swap(roots[0], roots[2]);
	if (roots[1] > roots[2])
	    std::swap(roots[1], roots[2]);

	if (roots[0] > MACH_EPS || roots[1] > MACH_EPS || roots[2] > MACH_EPS)
	    return true;
	else
	    return false;
}

//helper function to detect proximity between elements 
bool CollisionSolver3d::TriToBond(const TRI* tri,const BOND* bd, double h){
	POINT* pts[4];
	STATE* sl;
	bool status = false;

	/* do not consider bond point that is a tri vertex */
	for (int i = 0; i < 3; ++i)
	{
	    if (Point_of_tri(tri)[i] == bd->start ||
		Point_of_tri(tri)[i] == bd->end)
		return false;
	}

	/* make sure the coords are old coords */
	for (int i = 0; i < 3; ++i)
	    pts[i] = Point_of_tri(tri)[i];
	for (int i = 0; i < 3; ++i)
	{
	    sl = (STATE*)left_state(pts[i]);
	    for (int j = 0; j < 3; ++j)
		Coords(pts[i])[j] = sl->x_old[j];
	}

	/* detect proximity of start point of bond w.r.t. tri */
	pts[3] = bd->start;
	sl = (STATE*)left_state(pts[3]);
	for (int j = 0; j < 3; ++j)
	    Coords(pts[3])[j] = sl->x_old[j];
	if (PointToTri(pts,h)) status = true;

	/* detect proximity of end point of bond w.r.t. tri */
	pts[3] = bd->end;
	sl = (STATE*)left_state(pts[3]);
	for (int j = 0; j < 3; ++j)
	    Coords(pts[3])[j] = sl->x_old[j];
	if (PointToTri(pts,h)) status = true;
	
	/* detect proximity each edge of tri w.r.t. bond */
	pts[2] = bd->start;
	pts[3] = bd->end;
	for (int i = 0; i < 3; ++i)
	{
	    pts[0] = Point_of_tri(tri)[i];
	    pts[1] = Point_of_tri(tri)[(i+1)%3];
	    if (EdgeToEdge(pts, h)) status = true;
	}

	return status;
}

bool CollisionSolver3d::BondToBond(const BOND* b1, const BOND* b2, double h){
	POINT* pts[4];
	STATE* sl;
	bool status = false;

	pts[0] = b1->start;
	pts[1] = b1->end;
	pts[2] = b2->start;
	pts[3] = b2->end;

	/* do not consider two bonds that share a common point */
	for (int i = 0; i < 4; ++i)
	{
	    for (int j = i + 1; j < 4; ++j)
	    {
		if (pts[i] == pts[j])
		    return false;
	    }
	}

	/* make sure the coords are old coords */
	for (int i = 0; i < 4; ++i)
	{
	    sl = (STATE*)left_state(pts[i]);
	    for (int j = 0; j < 3; ++j)
		Coords(pts[i])[j] = sl->x_old[j];
	}

	/* detect proximity between two bonds */
	if (EdgeToEdge(pts, h))
		status = true;

	return status;
}

bool CollisionSolver3d::TriToTri(const TRI* tri1, const TRI* tri2, double h){
	POINT* pts[4];
	STATE* sl[2];
	bool status = false;
	for (int i = 0; i < 3; ++i)
	for (int j = 0; j < 3; ++j)
	{
	    if (Point_of_tri(tri1)[i] == Point_of_tri(tri2)[j])
		return false;
	}

	//make sure the coords are old coords;
	for (int i = 0; i < 3; ++i){
	    pts[0] = Point_of_tri(tri1)[i];
	    sl[0] = (STATE*)left_state(pts[0]);
	    pts[1] = Point_of_tri(tri2)[i];
	    sl[1] = (STATE*)left_state(pts[1]);
	    for (int j = 0; j < 2; ++j)
	    for (int k = 0; k < 3; ++k)
	        Coords(pts[j])[k] = sl[j]->x_old[k];
	}

	for (int k = 0; k < 2; ++k)
	for (int i = 0; i < 3; ++i)
	{
	    const TRI* tmp_tri1 = (k == 0) ? tri1 : tri2;
	    const TRI* tmp_tri2 = (k == 0) ? tri2 : tri1;
	    for (int j = 0; j < 3; ++j)
		pts[j] = Point_of_tri(tmp_tri2)[j];
	    pts[3] = Point_of_tri(tmp_tri1)[i];
	
	    //we do not consider point against 
            //a triangle that cotains it
            if (pts[3] == pts[0] || pts[3] == pts[1] ||
                pts[3] == pts[2])
                continue;
	    if (PointToTri(pts,h))
		status = true;
	}
	for (int i = 0; i < 3; ++i)
	{
	    pts[0] = Point_of_tri(tri1)[i];
	    pts[1] = Point_of_tri(tri1)[(i+1)%3];
	    for (int j = 0; j < 3; ++j){
		pts[2] = Point_of_tri(tri2)[j];
		pts[3] = Point_of_tri(tri2)[(j+1)%3];
		
		//we do not consider edges share a point
                if (pts[0] == pts[2] || pts[0] == pts[3] ||
                    pts[1] == pts[2] || pts[1] == pts[3])
                    continue;

	  	if (EdgeToEdge(pts, h))
		    status = true;
	    }  
	}
	return status;
}

static void PointToLine(POINT* pts[],double &a)
{
/*
*	x1 -----projP---- x2
*		  |
*		  |  dist
*		  *x3
*/
	double x12[3], x13[3];
	Pts2Vec(pts[0],pts[1],x12);
        Pts2Vec(pts[0],pts[2],x13);
        a = Dot3d(x13,x12)/Dot3d(x12,x12);
}

static bool EdgeToEdge(POINT** pts, double h, double root)
{
/*	x1	x3
 *	/	 \
 *     /	  \
 * x2 /		   \ x4
 * solve equation
 * x21*x21*a - x21*x43*b = x21*x31
 * -x21*x43*a + x43*x43*b = -x43*x31
 */
	double x21[3], x43[3], x31[3];
	double a, b;
	double tmp[3];
	double v1[3],v2[3];
	double nor[3], nor_mag, dist;

	Pts2Vec(pts[1],pts[0],x21);    
	Pts2Vec(pts[3],pts[2],x43);
	Pts2Vec(pts[2],pts[0],x31);
	Cross3d(x21,x43,tmp);
	if (Mag3d(tmp) < ROUND_EPS)
	{
	    return false; //ignore the case where two edges are parallel
	    //degenerate cases to parallel line segments
	    if (Mag3d(x21) > ROUND_EPS || Mag3d(x43) > ROUND_EPS){
	    	POINT* plist[3];
		double tmp_min_dist = HUGE;
		for (int i = 0; i < 4; ++i){
		    // p0-> p2--p3
		    // p1-> p2--p3
		    // p2-> p0--p1
		    // p3-> p0--p1
		    plist[0] = (i%2 == 0) ? pts[(i+2)%4] : pts[(i+1)%4];
		    plist[1] = (i%2 == 0) ? pts[(i+3)%4] : pts[(i+2)%4];
		    plist[2] = pts[i];
		    double tmp_vec[3], tmp_nor[3], tmp_dist, tmp_a;
		    minusVec(Coords(plist[1]),Coords(plist[0]),tmp_vec);
		    if (Mag3d(tmp_vec) < ROUND_EPS) continue;

		    PointToLine(plist,tmp_a);
		    tmp_a = std::max(std::min(tmp_a,1.0),0.0);
		    scalarMult(tmp_a,tmp_vec,tmp_vec);
		    addVec(Coords(plist[0]),tmp_vec,tmp_vec);
		    if (i/2 == 0)
		        minusVec(Coords(plist[2]),tmp_vec,tmp_nor);
		    else 
		        minusVec(tmp_vec,Coords(plist[2]),tmp_nor);
		    tmp_dist = distance_between_positions(
				tmp_vec,Coords(plist[2]),3);
		    if (tmp_dist < tmp_min_dist){
			memcpy((void*)nor,(void*)tmp_nor,3*sizeof(double));
			dist = tmp_min_dist = tmp_dist;
			if	(i == 0){a = 0.0; b = tmp_a;}
			else if (i == 1){a = 1.0; b = tmp_a;}
			else if (i == 2){a = tmp_a; b = 0.0;}
			else if (i == 3){a = tmp_a; b = 1.0;}
		    }
		}		
		if (dist < ROUND_EPS){
		    memcpy((void*)nor,
			   (Mag3d(x21) < ROUND_EPS) ? (void*)x43 : (void*)x21,
			    3*sizeof(double));
		}
	    }
	    else{
	   	//both x21 and x43 degenerate to points
		minusVec(Coords(pts[0]),Coords(pts[2]),nor);
		dist = Mag3d(nor);
		if (dist < ROUND_EPS){
		    nor[0] = nor[1] = nor[2] = 1.0;
		}
		a = 0.0; b = 0.0;
	    }
	}
	else
	{
	    a = (Dot3d(x43,x43)*Dot3d(x21,x31)-Dot3d(x21,x43)*Dot3d(x43,x31))/
	        (Dot3d(x21,x21)*Dot3d(x43,x43)-Dot3d(x21,x43)*Dot3d(x21,x43)); 
	    b = (Dot3d(x21,x43)*Dot3d(x21,x31)-Dot3d(x21,x21)*Dot3d(x43,x31))/
                (Dot3d(x21,x21)*Dot3d(x43,x43)-Dot3d(x21,x43)*Dot3d(x21,x43));
	    a = std::max(std::min(a,1.0),0.0);	
	    b = std::max(std::min(b,1.0),0.0);	
	    scalarMult(a,x21,v1);
	    scalarMult(b,x43,v2);
	    addVec(Coords(pts[0]),v1,v1);
	    addVec(Coords(pts[2]),v2,v2);
	    minusVec(v2,v1,nor);
	    nor_mag = Mag3d(nor);
	    if (nor_mag < 1000 * MACH_EPS)
	    {
		//v1 == v2;
                //two edges intersect with each other
                //normal direction is calculated with old position
		STATE* sl[4];
		for (int i = 0; i < 4; ++i)
		    sl[i] = (STATE*)left_state(pts[i]);
		for (int j = 0; j < 3; ++j)
		{
		    nor[j]  = (1.0-b) * sl[2]->x_old[j] + b * sl[3]->x_old[j];
		    nor[j] -= (1.0-a) * sl[0]->x_old[j] + a * sl[1]->x_old[j];
		}
	    }
	    dist = distBetweenCoords(v1,v2);
	}
	if (dist > h) 
	    return false;
	nor_mag = Mag3d(nor);
	if (nor_mag < MACH_EPS)
	{
            printf("Normal vector is nan");
            printf("a = %f, b = %f\n",a,b);
            printf("x_old:\n");
            for (int i = 0; i < 4; ++i){
                STATE* sl1 = (STATE*)left_state(pts[i]);
                printf("%f %f %f\n",sl1->x_old[0],sl1->x_old[1],sl1->x_old[2]);
            }
            printf("x_new:\n");
            for (int i = 0; i < 4; ++i){
                printf("%f %f %f\n",Coords(pts[i])[0],Coords(pts[i])[1],Coords(pts[i])[2]);
            }
            printf("avgVel:\n");
            for (int i = 0; i < 4; ++i){
                STATE* sl1 = (STATE*)left_state(pts[i]);
                printf("%f %f %f\n",sl1->avgVel[0],sl1->avgVel[1],sl1->avgVel[2]);
            }
            clean_up(ERROR);
	}
	else
	for (int i = 0; i < 3; ++i)
            nor[i] /= nor_mag;

	EdgeToEdgeImpulse(pts, nor, a, b, dist, root);
	return true;
}

static bool PointToTri(POINT** pts, double h, double root)
{
/*	x1
 *  	/\     x4 *
 *     /  \
 * x2 /____\ x3
 *
 * solve equation
 * x13*x13*w1 + x13*x23*w2 = x13*x43
 * x13*x23*w1 + x23*x23*w2 = x23*x43
 */
	double w[3] = {0.0};
	double x13[3], x23[3], x43[3];
	double nor[3] = {0.0}, nor_mag = 0.0, dist, det;

	Pts2Vec(pts[0],pts[2],x13);
	Pts2Vec(pts[1],pts[2],x23);
	Pts2Vec(pts[3],pts[2],x43);
	
	det = Dot3d(x13,x13)*Dot3d(x23,x23)-Dot3d(x13,x23)*Dot3d(x13,x23);
	if (fabs(det) < 1000 * MACH_EPS){ // change ROUND_EPS to 1000*MACH_EPS
	    return false; // ignore cases where tri reduces to a line or point
	    /*consider the case when det = 0*/
	    /*x13 and x23 are collinear*/
	    POINT* tmp_pts[3]; 
	    double max_len = 0;
	    //find longest side
	    for (int i = 0; i < 3; ++i){
		double tmp_dist = distance_between_positions(Coords(pts[i]),
				Coords(pts[(i+1)%3]),3);
		if (tmp_dist > max_len){
		    tmp_pts[0] = pts[i];
		    tmp_pts[1] = pts[(i+1)%3];
		    max_len = tmp_dist;
		}
	    }
	    if (max_len < ROUND_EPS){
		//triangle degenerate to a point
		minusVec(Coords(pts[0]),Coords(pts[3]),nor);
		dist = Mag3d(nor);
		if (dist < ROUND_EPS)
		for (int i = 0; i < 3; ++i){
		    w[i] = 1.0/3.0;
		    nor[i] = 1.0; 
		}
	    }
	    else{
		//triangle degnerate to a line between tmp_pts
		tmp_pts[2] = pts[3];
		double a;
		PointToLine(tmp_pts,a);
		for (int i = 0; i < 3; ++i)
		for (int j = 0; j < 2; ++j){
		    if (tmp_pts[j] == pts[i]){
			w[i] = (j == 0) ? a : 1.0-a;
		    }
		}
		double v[3];
		minusVec(Coords(tmp_pts[1]),Coords(tmp_pts[0]),v);
		scalarMult(a,v,v);
		addVec(Coords(tmp_pts[0]),v,v);
		minusVec(Coords(pts[3]),v,nor);
		dist = Mag3d(nor);
		//define nor vec if dist == 0
		if (Mag3d(nor) < ROUND_EPS)
		    nor[0] = nor[1] = nor[2] = 1.0;
	    }
	}
	else{
	    /*det != 0*/
	    /*x13 and x23 are non-collinear*/
	    Cross3d(x13, x23, nor);
	    nor_mag = Mag3d(nor);

	    /*compute the old direction*/
	    double x43_old[3];
	    STATE* tmp_sl[2];
	    tmp_sl[0] = (STATE*)left_state(pts[3]);
	    tmp_sl[1] = (STATE*)left_state(pts[2]);
	    minusVec(tmp_sl[0]->x_old,tmp_sl[1]->x_old,x43_old);

	    /*correct the normal direction*/
	    /*always pointing from triangle to p4*/
	    dist = Dot3d(x43_old, nor);
	    for (int i = 0; i < 3; ++i)
	        nor[i] /= nor_mag * ((dist >= 0)? 1.0:-1.0);
	    dist = fabs(Dot3d(x43, nor));

	    w[0] = (Dot3d(x13,x43)*Dot3d(x23,x23)-Dot3d(x23,x43)*Dot3d(x13,x23))/det;
	    w[1] = (Dot3d(x13,x13)*Dot3d(x23,x43)-Dot3d(x13,x23)*Dot3d(x13,x43))/det;
	    w[2] = 1 - w[0] - w[1];
	    
	}
	/* test for corner cases */
	if (fabs(w[0]) < ROUND_EPS || fabs(w[1]) < ROUND_EPS 
	    || fabs(w[2]) < ROUND_EPS)
	{
	    double vec[3] = {0, 0, 0};
	    STATE* tmp_sl = (STATE*)left_state(pts[3]);
	    for (int j = 0; j < 3; ++j)
		vec[j] = tmp_sl->x_old[j];
	    for (int i = 0; i < 3; ++i)
	    {
		tmp_sl = (STATE*)left_state(pts[i]);
		for (int j = 0; j < 3; ++j)
		    vec[j] -= w[i] * tmp_sl->x_old[j];
	    }
	    double vec_mag = Mag3d(vec);
	    if (vec_mag > ROUND_EPS)
	    {
		for (int j = 0; j < 3; ++j)
		    nor[j] = vec[j];
	    }
	}
	/* end of the test */
	
	nor_mag = Mag3d(nor);
	if (nor_mag > ROUND_EPS)
	    for (int i = 0; i < 3; ++i)
	        nor[i] /= nor_mag;
	else{
	    std::cout << "nan nor vec" << std::endl;
	    printPointList(pts,4);
	    clean_up(ERROR);
	}

	double c_len = 0;	
	for (int i = 0; i < 3; ++i){
	    double tmp_dist = distance_between_positions(Coords(pts[i]),
                                Coords(pts[(i+1)%3]),3);
	    if (tmp_dist > c_len) 
		c_len = tmp_dist;
	}
	if (dist > h)
	    return false;
	for (int i = 0; i < 3; ++i)
	{
	    double eps = CollisionSolver3d::getRoundingTolerance(); 
	    //test, use eps instead of h/c_len
	    if (w[i] > 1+eps || w[i] < -eps) 
		return false;
	}
	PointToTriImpulse(pts, nor, w, dist,root);
	return true;
}

/* repulsion and friction functions, update velocity functions */
static void PointToTriImpulse(POINT** pts, double* nor, double* w, double dist, double root)
{
	if (debugging("collision"))
	    CollisionSolver::pt_to_tri++;
	STATE *sl[4];
	for (int i = 0; i < 4; ++i)
	    sl[i] = (STATE*)left_state(pts[i]);

	double v_rel[3] = {0.0}, vn = 0.0, vt = 0.0;
	double impulse = 0.0, m_impulse = 0.0;
	double k, m, lambda, dt, h, cr, sum_w = 0.0;
	k      = CollisionSolver::getSpringConstant();
	m      = CollisionSolver::getPointMass();
	dt     = CollisionSolver::getTimeStepSize();
	lambda = CollisionSolver::getFrictionConstant(); 
	h      = CollisionSolver::getFabricThickness();
	cr     = CollisionSolver::getRestitutionCoef();
	dist   = h - dist;
	double rigid_impulse[2] = {0.0};

	/* it is supposed to use the average velocity*/
	for (int i = 0; i < 3; ++i)
	{
	    v_rel[i] += sl[3]->avgVel[i];
	    for (int j = 0; j < 3; ++j)
		v_rel[i] -= w[j] * sl[j]->avgVel[i];
	}
	vn = Dot3d(v_rel, nor);
	if (Dot3d(v_rel, v_rel) > sqr(vn))
	    vt = sqrt(Dot3d(v_rel, v_rel) - sqr(vn));
	else
	    vt = 0.0;
	if (vn < 0)
	{
	    if (isStaticRigidBody(pts[3]) ||
	       (isStaticRigidBody(pts[0]) && isStaticRigidBody(pts[1]) &&
		isStaticRigidBody(pts[2]))){
		impulse = vn;
		rigid_impulse[0] = vn;
		rigid_impulse[1] = vn;
	    }
	    else if (isMovableRigidBody(pts[0]) && isMovableRigidBody(pts[1]) 
		&& isMovableRigidBody(pts[2]) && isMovableRigidBody(pts[3]))
	    {
		double m1 = total_mass(pts[0]->hs);
		double m2 = total_mass(pts[3]->hs);
		rigid_impulse[0] = vn * m2 / (m1 + m2);
		rigid_impulse[1] = vn * m1 / (m1 + m2);
	    }
	    else if (isMovableRigidBody(pts[0]) && isMovableRigidBody(pts[1])
		&& isMovableRigidBody(pts[2]))
	    {
		rigid_impulse[0] = 0.5 * vn;
		impulse = 0.5 * vn;
	    }
	    else if (isMovableRigidBody(pts[3]))
	    {
		impulse = 0.5 * vn;
		rigid_impulse[1] = 0.5 * vn;
	    }
	    else
	        impulse = vn * 0.5;
	    for (int i = 0; i < 3; ++i)
	    {
		if (isStaticRigidBody(pts[i]))
		        w[i] = 0.0;
		sum_w += w[i];
	    }
	    if (fabs(sum_w) > MACH_EPS)
	        scalarMult(1.0/sum_w,w,w);
	}
	if (vn * dt < 0.1 * dist)
	{
	    if (isRigidBody(pts[0]) && isRigidBody(pts[1]) &&
		isRigidBody(pts[2]) && isRigidBody(pts[3]))
	    {
		rigid_impulse[0] *= 1.0 + cr;
		rigid_impulse[1] *= 1.0 + cr;
	    }
	    else
	    {
		double tmp = - std::min(dt*k*dist/m, (0.1*dist/dt - vn));
		impulse += tmp;
		rigid_impulse[0] += tmp;
		rigid_impulse[1] += tmp;
	    }
	}
	if (fabs(sum_w) < MACH_EPS)
	    m_impulse = impulse;
	else
	    m_impulse = 2.0 * impulse / (1.0 + Dot3d(w, w));

//uncomment the following the debugging purpose
if (debugging("CollisionImpulse"))
if (fabs(m_impulse) > 0.0){
	printf("real PointToTri collision, dist = %e\n",dist);
	printf("vt = %f, vn = %f, dist = %f\n",vt,vn,dist);
	printf("v_rel = %f %f %f\n",v_rel[0],v_rel[1],v_rel[2]);
	printf("nor = %f %f %f\n",nor[0],nor[1],nor[2]);
	printf("m_impuse = %f, impulse = %f, w = [%f %f %f]\n",
		m_impulse,impulse,w[0],w[1],w[2]);
	printf("dt = %f, root = %f\n",dt,root);
	printf("k = %f, m = %f\n",k,m);
	printf("x_old:\n");
	for (int i = 0; i < 4; ++i){
	    STATE* sl1 = (STATE*)left_state(pts[i]);
	    printf("%f %f %f\n",sl1->x_old[0],sl1->x_old[1],sl1->x_old[2]);
	}
	printf("x_new:\n");
	for (int i = 0; i < 4; ++i){
	    printf("%f %f %f\n",Coords(pts[i])[0],Coords(pts[i])[1],Coords(pts[i])[2]);
	}
	printf("avgVel:\n");
	for (int i = 0; i < 4; ++i){
	    STATE* sl1 = (STATE*)left_state(pts[i]);
	    printf("%f %f %f\n",sl1->avgVel[0],sl1->avgVel[1],sl1->avgVel[2]);
	}
}

	if (isRigidBody(pts[0]) && isRigidBody(pts[1]) && 
	    isRigidBody(pts[2]) && isRigidBody(pts[3]))
	{
	    if (isMovableRigidBody(pts[0]))
		SpreadImpactZoneImpulse(pts[0], rigid_impulse[0], nor);
	    if (isMovableRigidBody(pts[3]))
		SpreadImpactZoneImpulse(pts[3], -1.0 * rigid_impulse[1], nor);
	    return;
	}
	for (int i = 0; i < 3; ++i)
	{
	    if (isStaticRigidBody(pts[i])) continue;
	    double t_impulse = m_impulse;
	    if (isMovableRigidBody(pts[i])) t_impulse = rigid_impulse[0];
	    for(int j = 0; j < 3; ++j)
	    {
	        sl[i]->collsnImpulse[j] += w[i] * t_impulse * nor[j];
	        if (fabs(vt) > ROUND_EPS)
	            sl[i]->friction[j] += std::max(-fabs(lambda * w[i] * 
	    	    t_impulse/vt), -1.0) * (v_rel[j] - vn * nor[j]);
	    }
	    sl[i]->collsn_num += 1;
	}
	if (!isStaticRigidBody(pts[3]))
	{
	    double t_impulse = m_impulse;
	    if (isMovableRigidBody(pts[3])) t_impulse = rigid_impulse[1];
	    for (int j = 0; j < 3; ++j)
	    {
	        sl[3]->collsnImpulse[j] -= t_impulse * nor[j];
	        if (fabs(vt) > ROUND_EPS)
	            sl[3]->friction[j] += std::max(-fabs(lambda * 
			    t_impulse/vt), -1.0) * (v_rel[j] - vn * nor[j]);
	    }
	    sl[3]->collsn_num += 1;
	}
	for (int j = 0; j < 4; ++j)
	     if(isStaticRigidBody(pts[j])) 
		memset((void*)sl[j]->collsnImpulse,0,3*sizeof(double));

	if (debugging("CollisionImpulse"))
	for (int i = 0; i < 4; ++i){
	    printf("pt[%d], collsnImp = [%f %f %f], friction = [%f %f %f]\n",
	i, sl[i]->collsnImpulse[0],sl[i]->collsnImpulse[1],sl[i]->collsnImpulse[2],
	sl[i]->friction[0],sl[i]->friction[1],sl[i]->friction[2]);
	}

	for (int kk = 0; kk < 4; kk++)
	for (int j = 0; j < 3; ++j){
	    if (std::isnan(sl[kk]->collsnImpulse[j]) ||
		std::isinf(sl[kk]->collsnImpulse[j])){
		printf("PointToTri: sl[%d]->impl[%d] = nan\n",kk,j);
		for (int i = 0; i < 4; ++i){
		printf("points[%d] = %p\n",i,(void*)pts[i]);
		printf("coords = [%f %f %f]\n",Coords(pts[i])[0],
			Coords(pts[i])[1],Coords(pts[i])[2]);
		}
		printf("w = [%f %f %f]\nnor = [%f %f %f]\ndist = %f\n",
		w[0],w[1],w[2],nor[0],nor[1],nor[2],dist);
		printf("v_rel = [%f %f %f]\n",v_rel[0],v_rel[1],v_rel[2]);
	        clean_up(ERROR);
	    }
	}
}

static void EdgeToEdgeImpulse(POINT** pts, double* nor, double a, double b, double dist, double root)
{
	if (debugging("collision"))
	    CollisionSolver::edg_to_edg++;

	STATE *sl[4];
	for (int i = 0; i < 4; ++i)
	    sl[i] = (STATE*)left_state(pts[i]);

	double v_rel[3] = {0.0, 0.0, 0.0}, vn = 0.0, vt = 0.0;
	double impulse = 0.0, m_impulse = 0.0;
	double k, m, lambda, dt, h, cr;
	k      = CollisionSolver::getSpringConstant();
	m      = CollisionSolver::getPointMass();
	dt     = CollisionSolver::getTimeStepSize();
	lambda = CollisionSolver::getFrictionConstant(); 
	h      = CollisionSolver::getFabricThickness();
	cr     = CollisionSolver::getRestitutionCoef();
	dist   = h - dist;
	double rigid_impulse[2] = {0.0};
	double wa[2] = {1.0 - a, a}, wb[2] = {1.0 - b, b};

	/* it is supposed to use the average velocity*/
	for (int j = 0; j < 3; ++j)
	{
	    v_rel[j]  = (1.0-b) * sl[2]->avgVel[j] + b * sl[3]->avgVel[j];
	    v_rel[j] -= (1.0-a) * sl[0]->avgVel[j] + a * sl[1]->avgVel[j];
	}
	vn = Dot3d(v_rel, nor);
	if (Dot3d(v_rel, v_rel) > sqr(vn))
	    vt = sqrt(Dot3d(v_rel, v_rel) - sqr(vn));
	else
	    vt = 0.0;

	if (vn < 0.0)
	{
	    if ((isStaticRigidBody(pts[0]) && isStaticRigidBody(pts[1])) ||
	    	(isStaticRigidBody(pts[2]) && isStaticRigidBody(pts[3])))
	    {
		impulse = vn;
		rigid_impulse[0] = vn;
		rigid_impulse[1] = vn;
	    }
	    else if (isMovableRigidBody(pts[0]) && isMovableRigidBody(pts[1]) 
		&& isMovableRigidBody(pts[2]) && isMovableRigidBody(pts[3]))
	    {
		double m1 = total_mass(pts[0]->hs);
		double m2 = total_mass(pts[2]->hs);
		rigid_impulse[0] = vn * m2 / (m1 + m2);
		rigid_impulse[1] = vn * m1 / (m1 + m2);
	    }
	    else if (isMovableRigidBody(pts[0]) && isMovableRigidBody(pts[1]))
	    {
		rigid_impulse[0] = 0.5 * vn;
		impulse = 0.5 * vn; 
	    }
	    else if (isMovableRigidBody(pts[2]) && isMovableRigidBody(pts[3]))
	    {
		impulse = 0.5 * vn;
		rigid_impulse[1] = 0.5 * vn;
	    }
	    else
		impulse = vn * 0.5;
    	    if (isStaticRigidBody(pts[0])) wa[0] = 0.0;
	    if (isStaticRigidBody(pts[1])) wa[1] = 0.0;
	    if (isStaticRigidBody(pts[2])) wb[0] = 0.0;
	    if (isStaticRigidBody(pts[3])) wb[1] = 0.0;
	}
	if (vn * dt < 0.1 * dist)
	{
	    if (isRigidBody(pts[0]) && isRigidBody(pts[1]) &&
		isRigidBody(pts[2]) && isRigidBody(pts[3]))
	    {
		rigid_impulse[0] *= 1.0 + cr;
		rigid_impulse[1] *= 1.0 + cr;
	    }
	    else
	    {
		double tmp = - std::min(dt*k*dist/m, (0.1*dist/dt - vn));
		impulse += tmp;
		rigid_impulse[0] += tmp;
		rigid_impulse[1] += tmp;
	    }
	}
	if (wa[0] + wa[1] < MACH_EPS || wb[0] + wb[1] < MACH_EPS)
	    m_impulse = impulse;
	else
	    m_impulse = 2.0 * impulse / (wa[0]*wa[0] + wa[1]*wa[1] 
					+ wb[0]*wb[0] + wb[1]*wb[1]);

//uncomment the following for the debugging purpose
if (debugging("CollisionImpulse"))
if (fabs(m_impulse) > 0){
	printf("real EdgeToEdge collision\n");
	printf("vt = %f, vn = %f, dist = %f\n",vt,vn,dist);
	printf("v_rel = %f %f %f\n",v_rel[0],v_rel[1],v_rel[2]);
	printf("nor = %f %f %f\n",nor[0],nor[1],nor[2]);
	printf("m_impuse = %f, impulse = %f, a = %f, b = %f\n",
		m_impulse,impulse,a,b);
	printf("root = %e,h = %e, dt = %e\n",root,h,dt);
	printf("x_old:\n");
	for (int i = 0; i < 4; ++i){
	    STATE* sl1 = (STATE*)left_state(pts[i]);
	    printf("%f %f %f\n",sl1->x_old[0],sl1->x_old[1],sl1->x_old[2]);
	}
	printf("x_new:\n");
	for (int i = 0; i < 4; ++i){
	    printf("%f %f %f\n",Coords(pts[i])[0],Coords(pts[i])[1],Coords(pts[i])[2]);
	}
	printf("avgVel:\n");
	for (int i = 0; i < 4; ++i){
	    STATE* sl1 = (STATE*)left_state(pts[i]);
	    printf("%f %f %f\n",sl1->avgVel[0],sl1->avgVel[1],sl1->avgVel[2]);
	}
	printf("\n");
}

	if (isRigidBody(pts[0]) && isRigidBody(pts[1]) && 
	    isRigidBody(pts[2]) && isRigidBody(pts[3]))
	{
	    if (isMovableRigidBody(pts[0]))
		SpreadImpactZoneImpulse(pts[0], rigid_impulse[0], nor);
	    if (isMovableRigidBody(pts[2]))
		SpreadImpactZoneImpulse(pts[2], -1.0 * rigid_impulse[1], nor);
	    return;
	}
	for (int j = 0; j < 3; ++j)
	{
	    //p[0]
	    if (!isStaticRigidBody(pts[0]))
	    {
		double t_impulse = m_impulse;
		if (isMovableRigidBody(pts[0])) t_impulse = rigid_impulse[0];
	        sl[0]->collsnImpulse[j] += wa[0] * t_impulse * nor[j];
	        if (fabs(vt) > ROUND_EPS)
	            sl[0]->friction[j] += std::max(-fabs(lambda * wa[0] *
			    t_impulse/vt), -1.0) * (v_rel[j] - vn * nor[j]);
		if (j == 0) sl[0]->collsn_num += 1;
	    }

	    //p[1]
	    if (!isStaticRigidBody(pts[1]))
	    {
		double t_impulse = m_impulse;
		if (isMovableRigidBody(pts[1])) t_impulse = rigid_impulse[0];
	        sl[1]->collsnImpulse[j] += wa[1] * t_impulse * nor[j];
	        if (fabs(vt) > ROUND_EPS)
	            sl[1]->friction[j] += std::max(-fabs(lambda * wa[1] * 
			    t_impulse/vt), -1.0) * (v_rel[j] - vn * nor[j]);
		if (j == 0) sl[1]->collsn_num += 1;
	    }

	    //p[2]
	    if (!isStaticRigidBody(pts[2]))
	    {
		double t_impulse = m_impulse;
		if (isMovableRigidBody(pts[2])) t_impulse = rigid_impulse[1];
	        sl[2]->collsnImpulse[j] -= wb[0] * t_impulse * nor[j];
	        if (fabs(vt) > ROUND_EPS)
	            sl[2]->friction[j] += std::max(-fabs(lambda * wb[0] *
			    t_impulse/vt), -1.0)*(v_rel[j] - vn * nor[j]);
		if (j == 0) sl[2]->collsn_num += 1;
	    }

	    //p[3]
	    if (!isStaticRigidBody(pts[3]))
	    {
		double t_impulse = m_impulse;
		if (isMovableRigidBody(pts[3])) t_impulse = rigid_impulse[1];
	        sl[3]->collsnImpulse[j] -= wb[1] * t_impulse * nor[j];
	        if (fabs(vt) > ROUND_EPS)
	            sl[3]->friction[j] += std::max(-fabs(lambda * wb[1] * 
			    t_impulse/vt), -1.0) * (v_rel[j] - vn * nor[j]);
		if (j == 0) sl[3]->collsn_num += 1;
	    }
	}
	for (int j = 0; j < 4; ++j)
	     if(isStaticRigidBody(pts[j])) 
		memset((void*)sl[j]->collsnImpulse,0,3*sizeof(double));

	for (int i = 0; i < 4; i++)
	for (int j = 0; j < 3; ++j){
	    if (std::isnan(sl[i]->collsnImpulse[j]) ||
		std::isinf(sl[i]->collsnImpulse[j])){
		printf("EdgeToEdge: sl[%d]->impl[%d] = nan\n",i,j);
		printf("a b = %f %f, nor = [%f %f %f], dist = %f\n",
		a,b,nor[0],nor[1],nor[2],dist);
	        clean_up(ERROR);
	    }
	}

}

static void unsort_surface_point(SURFACE *surf)
{
        TRI *tri;
        POINT *p;
        int i;

        for (tri = first_tri(surf); !at_end_of_tri_list(tri,surf);
                        tri = tri->next)
        {
            for (i = 0; i < 3; ++i)
            {
                p = Point_of_tri(tri)[i];
                sorted(p) = NO;
            }
        }
}       /* end unsort_surface_point */

