#include "drag.h"
#include <iostream>
#include <cmath>
#include <limits>

const double EPS = 1e-15;
static bool isPointInBall(double p[], double c[], double r);
static double Dot3d(const double*,const double*);
static void Cross3d(const double*, const double*, double*);
static double Mag3d(const double*);
double Drag::m_thickness = 0.001;
double Drag::m_tol = 0.001;

bool PLANE::isFront(double *p)
{
    double v[3] = {0};
    for (int i = 0; i < 3; ++i)
	v[i] = p[i] - center[i];
    if (Dot3d(v, nor) > 0)
	return true;
    else
	return false;
}

bool PLANE::isBack(double *p)
{
    return !isFront(p);
}

double PLANE::distance(double *p)
{
    double v[3] = {0};
    for (int i = 0; i < 3; ++i)
    {
	v[i] = p[i] - center[i];
    }
    return fabs(Dot3d(v, nor));
}

PLANE::PLANE(const double*p, const double* v)
{
    for (int i = 0; i < 3; ++i)
    {
	center[i] = p[i];
	dir[i] = v[i];
    }
    double v_m = Mag3d(v);
    if (v_m < EPS)
    {
	std::cout << "v_m is too small" << std::endl;
	return;
    }
    for (int i = 0; i < 3; ++i)
    {
	nor[i] = dir[i]/v_m;
    }
}

LINE::LINE(const double* p1, const double* p2)
{
    std::copy(p1, p1 + 3, pOnLine1); 
    std::copy(p2, p2 + 3, pOnLine2); 
}

// calculate the shortest distance from a point to a line
double LINE::distance(double* p)
{
    double p21[3], p10[3];

    std::transform(pOnLine2, pOnLine2 + 3, pOnLine1, p21, 
		std::minus<double>());
    std::transform(pOnLine1, pOnLine1 + 3, p, p10, std::minus<double>());
	
    double dtemp[3]; 
    
    Cross3d(p21, p10, dtemp); 

    return Mag3d(dtemp) / Mag3d(p21);     
}

static void computeRotationMatrix(double* dir,
                           double theta,
                           double R[][3])
{
    double u,v,w;
    //dir should be normalized
    double M = std::sqrt(dir[0]*dir[0]
        + dir[1]*dir[1] + dir[2]*dir[2]);
    u = dir[0]/M; v = dir[1]/M; w = dir[2]/M;

    R[0][0] = u*u + (1-u*u)*cos(theta);
    R[0][1] = u*v*(1-cos(theta)) - w*sin(theta);
    R[0][2] = u*w*(1-cos(theta))+v*sin(theta);

    R[1][0] = u*v*(1-cos(theta)) + w*sin(theta);
    R[1][1] = v*v + (1 - v*v)*cos(theta);
    R[1][2] = v*w*(1-cos(theta)) - u*sin(theta);

    R[2][0] = u*w*(1-cos(theta)) - v*sin(theta);
    R[2][1] = v*w*(1-cos(theta)) + u*sin(theta);
    R[2][2] = w*w+(1-w*w)*cos(theta);
}

static void Normalize(double* v)
{
    double mv = Mag3d(v);
    for (int i = 0; i < 3; ++i)
        v[i] /= mv;
}

static double distance_between(const double p1[], 
		const double p2[], int dim) {
    double dist = 0;
    for (int i = 0; i < dim; ++i) {
	dist += (p1[i]-p2[i])*(p1[i]-p2[i]);
    }
    return sqrt(dist);
}

static double angleBetween(const double* v1, const double *v2)
{
    return acos(Dot3d(v1, v2)/(Mag3d(v1)*Mag3d(v2)));
}

static void Cross3d(const double B[], const double C[], double ans[])
{
    (ans)[0] = ((B)[1])*((C)[2]) - ((B)[2])*((C)[1]);       
    (ans)[1] = ((B)[2])*((C)[0]) - ((B)[0])*((C)[2]);      
    (ans)[2] = ((B)[0])*((C)[1]) - ((B)[1])*((C)[0]);     
}


static double Dot3d(const double* v1, const double* v2) 
{ 
    return v1[0]*v2[0]+v1[1]*v2[1]+v1[2]*v2[2];
}

static double Mag3d(const double* v) 
{ 
    return sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);
}

static void M3xV3(const double M[][3], const double* V, double* ans)
{
    for (int i = 0; i < 3; ++i)
    {
        ans[i] = 0;
        for (int j = 0; j < 3; ++j)
            ans[i] += M[i][j] * V[j]; 
    }
}

static double pointToLine(
	const double* p, 
	const double *o, 
	const double *dir, 
	double* ans_v = NULL)
{
    double v[3];
    for (int i = 0; i < 3; ++i)
    {
	v[i] = p[i] - o[i];
    }
    double nor_mag = Mag3d(dir);
    double nor[3];
    for (int i = 0; i < 3; ++i)
    {
	nor[i] = dir[i]/nor_mag;
    }
    double v1[3];
    for (int i = 0; i < 3; ++i)
    {
	v1[i] = v[i] - nor[i] * Dot3d(v, nor);
    }
    if (ans_v != NULL)
	std::copy(v1, v1+3, ans_v);
    return Mag3d(v1);
}

//Drag
bool Drag::validateData(const Drag::Info& info) {
    if (info.data().size() != this->dataSize()) 
    {
	std::cout << "Warning: "<< info.id() << " needs " 
	<< this->dataSize() << " parameters, but " << info.data().size() 
	<< " are given!" << std::endl;
	return false;
    }
    return true;
}

Drag* Drag::dragFactory(const Drag::Info& info) {
    for (size_t i = 0; i < prototypes.size(); ++i) {
	if (prototypes[i]->id() == info.id()) 
	{
	    std::cout << "Adding drag: " 
		<< prototypes[i]->id() << std::endl;
	    for (size_t j = 0; j < info.data().size(); ++j)
		std::cout << info.data()[j] << " ";
	    std::cout << std::endl;
	    return prototypes[i]->clone(info);
	}
    }
    std::cout << "Warning: unknown drag type: " 
	      << info.id() << std::endl;
    return NULL;
}

void Drag::spinToAxis(double* cen, double* dir, double theta, double* p)
{
    double a,b,c,u,v,w,x,y,z;
    //normalize dir
    double len = Mag3d(dir);
    if (len < EPS)
    {
	std::cout << "Warning: len < EPS" << std::endl;
        std::cout << "dir = " << dir[0] << ' ' << dir[1] << ' ' << dir[2] 
		<< std::endl;
	return;
    }
    else
    for (int i = 0; i < 3; ++i)
	dir[i] /= len;

    a = cen[0]; b = cen[1]; c = cen[2];
    u = dir[0]; v = dir[1]; w = dir[2];
    x = p[0]; y = p[1]; z = p[2];
    p[0] = (a*(v*v+w*w)-u*(b*v+c*w-u*x-v*y-w*z))*(1-cos(theta)) + 
            x*cos(theta) + (-c*v+b*w-w*y+v*z)*sin(theta);
    p[1] = (b*(u*u+w*w)-v*(a*u+c*w-u*x-v*y-w*z))*(1-cos(theta)) +
            y*cos(theta) + (c*u-a*w+w*x-u*z)*sin(theta);
    p[2] = (c*(u*u+v*v)-w*(a*u+b*v-u*x-v*y-w*z))*(1-cos(theta)) +
            z*cos(theta) + (-b*u+a*v-v*x+u*y)*sin(theta);
}

// Point Drag
Drag* PointDrag::clone(const Drag::Info& info) {
    if (!validateData(info))
	return NULL;
    else {
	const std::vector<double>& v = info.data();
	const double *it = &(v.front());
	return new PointDrag(it,*(it+3),it+4,it+7,*(it+10));
    }
}

void PointDrag::preprocess(std::vector<SpringVertex*>& pts) {
    if (!first) return;
    first = false;
    for (size_t i = 0; i < pts.size(); ++i)
    {
        double *p = pts[i]->getCoords();
	if (rad*rad > (p[0]-cent[0])*(p[0]-cent[0]) +
                     (p[1]-cent[1])*(p[1]-cent[1]) +
                     (p[2]-cent[2])*(p[2]-cent[2]))
	    pts[i]->setRegistered();
	else
	    pts[i]->unsetRegistered();
    }
}

void PointDrag::setVel(SpringVertex* sv) {
    double *v = sv->getVel();
    if (sv->isRegistered())
        std::copy(m_v,m_v+3,v);
}

void PointDrag::setAccel(SpringVertex* sv) {
    double *a = sv->getExternalAccel();
    if (sv->isRegistered())
        std::copy(m_a,m_a+3,a);
}    

PointDrag::PointDrag(
	const double c[], 
	double r, 
	const double v[], 
	const double a[],  
	double t) : Drag(t), rad(r) {
    if (v) 
	std::copy(v,v+3,m_v);
    else 
	std::fill(m_v,m_v+3,0);

    if (a) 
	std::copy(a,a+3,m_a);
    else 
	std::fill(m_a,m_a+3,0);

    if (c)
        std::copy(c,c+3,cent);
}

//Line drag
void LineDrag::preprocess(std::vector<SpringVertex*>& pts)
{
    if (!first) return ; 
    first = false; 
    for (size_t i = 0; i < pts.size(); i++)
    {
	 SpringVertex* sv = pts[i]; 
         double *p = sv->getCoords(); 
	 
	 if (dragLine->distance(p) < 1.1 * getTolerance())
         {
	     sv->point_type = TRANS_DPOINT; 
  	     sv->setRegistered(); 
         }
         else if (controlLine1->distance(p) < 1.1 * getTolerance())
	 {
	     sv->point_type = TRANS_CPOINT1; 
	     sv->setRegistered(); 
	 }
	 else if (controlLine2->distance(p) < 1.1 * getTolerance())
	 {
	     sv->point_type = TRANS_CPOINT2;
             sv->setRegistered();
	 }
	 else 
	 {
	     sv->point_type = FREE_POINT; 
	     sv->unsetRegistered(); 
	 }
    }
}

void LineDrag::setVel(SpringVertex* sv) 
{
    if (sv->point_type >= TRANS_DPOINT)
    {
	if (m_dt < EPS)
            std::cout << "Warning: m_dt = " << m_dt
                      << " is too small. "
                      << "Have you called setTimeStepSize() before using it?"
                      << std::endl;
    }
    if (sv->point_type == TRANS_DPOINT)
        std::copy(veld, veld + 3, sv->getVel());
    else if (sv->point_type == TRANS_CPOINT1)
	std::copy(velc1, velc1 + 3, sv->getVel());
    else if (sv->point_type == TRANS_CPOINT2)
	std::copy(velc2, velc2 + 3, sv->getVel());
    else if (sv->point_type == STATIC_POINT)
        std::fill(sv->getVel(), sv->getVel() + 3, 0);  
    else
	return; 
}

void LineDrag::setAccel(SpringVertex* sv) 
{
}

Drag* LineDrag::clone(const Drag::Info& info) 
{
    if (!validateData(info))
        return NULL;
    else
    {
	LineDrag* ld = new LineDrag(); 
	const double* it = &(info.data().front());
        ld->dragLine = new LINE(it, it + 3); 
        std::copy(it + 6, it + 9, ld->veld);
        ld->controlLine1 = new LINE(it + 9, it + 12);
        std::copy(it + 15, it + 18, ld->velc1);
        std::copy(it + 18, it + 21, ld->accelc1);
        ld->accelStartTime1 = *(it + 21);
        ld->controlLine2 = new LINE(it + 22, it + 25);
        std::copy(it + 28, it + 31, ld->velc2); 
	std::copy(it + 31, it + 34, ld->accelc2);
        ld->accelStartTime2 = *(it + 34);
	ld->m_t = *(it + 35); 
	return ld; 
    }
}

void LineDrag::postprocess(std::vector<SpringVertex*>& pts)
{
    accumCurTime(); 
    if (m_ct > accelStartTime1)
    {
	for (size_t i = 0; i < pts.size(); i++)
    	{
             double deltav[3] = {0.0};

             if (pts[i]->point_type == TRANS_CPOINT1)
             {
                 std::transform(accelc1, accelc1 + 3, deltav,
                        bind1st(std::multiplies<double>(), m_dt));
                 std::transform(velc1, velc1 + 3, deltav, velc1,
                        std::plus<double>());
             }
        }
    }
    if (m_ct > accelStartTime2)
    {
	double deltat = m_ct - accelStartTime2; 

        for (size_t i = 0; i < pts.size(); i++)
        {
             double deltav[3];

	     if (pts[i]->point_type == TRANS_CPOINT2)
             {
                 std::transform(accelc2, accelc2 + 3, deltav,
                        bind1st(std::multiplies<double>(), m_dt));
                 std::transform(velc2, velc2 + 3, deltav, velc2,
                        std::plus<double>());
             }
        }
    }
}

//Dravity drag
void GravityDrag::setAccel(SpringVertex* sv) {
    double *a = sv->getExternalAccel();
    if (sv->isRegistered())
        std::fill(a,a+3,0);
    else
        std::copy(m_a,m_a+3,a);
}

GravityDrag::GravityDrag(const double c[], const double r, const double a[], double t)
             : PointDrag(c,r,NULL,a,t) {}

Drag* GravityDrag::clone(const Drag::Info& info) {
    if (!validateData(info)) 
	return NULL;
    else {
	const std::vector<double>& v = info.data();
	const double *it = &(v.front());
	return new GravityDrag(it,*(it+3),it+4,*(it+7));
    }
}

//MultiplePointDrag
void MultiplePointDrag::preprocess(std::vector<SpringVertex*>& pts) {
    if (!first) return;
    first = false;
    for (size_t i = 0; i < pts.size(); ++i)
    {
	int iid = getDragPointNumber(pts[i]->getCoords());
    	if (iid != -1)
	{
	    pts[i]->point_type = iid;
	    pts[i]->setRegistered();
	}
	else
	{
	    pts[i]->unsetRegistered();
	}
    }
}

int MultiplePointDrag::getDragPointNumber(double p[]) {
    for (size_t i = 0; i < dragPoints.size(); ++i) 
        if (distance_between(p,dragPoints[i]->cent,3) 
		< dragPoints[i]->rad)
	    return i;
    return -1;
}

void MultiplePointDrag::setVel(SpringVertex* sv) {
    double *vel = sv->getVel();

    if (sv->isRegistered()) {
	size_t i = static_cast<size_t>(sv->point_type);
	std::copy(dragPoints[i]->m_v,dragPoints[i]->m_v+3,vel);
    }
} 

void MultiplePointDrag::setAccel(SpringVertex* sv) {
    double *accel = sv->getExternalAccel();
    if (sv->isRegistered()) {
	size_t i = static_cast<size_t>(sv->point_type);
        std::copy(dragPoints[i]->m_a,dragPoints[i]->m_a+3,accel);
    }
}

void MultiplePointDrag::addDragPoint(DragPoint* dragPoint) {
    dragPoints.push_back(dragPoint);
}

MultiplePointDrag::MultiplePointDrag(double t): Drag(t) {}

MultiplePointDrag::DragPoint::DragPoint(double c[], double r,
		double v[], double a[]) {
    std::copy(v,v+3,m_v);
    std::copy(a,a+3,m_a);
    std::copy(c,c+3,cent);
    rad = r;
}

Drag* MultiplePointDrag::clone(const Drag::Info& info) {
    if ((info.data().size() - 1)%10 != 0) {
	std::cout << "Multiple point drag should have "
	          << "10x + 1 parameters, " 
		  << "but " << info.data().size() 
		  << " parameters are given" 
		  << std::endl; 
	return NULL;
    }
    else {
	MultiplePointDrag* mpd = new MultiplePointDrag(info.data().back());
	int num_of_point = (info.data().size()-1)/10;
	std::vector<double> v = info.data();
	for (int i = 0; i < num_of_point; ++i) {
	    double *it = &(v.front()) + i*10;
	    mpd->addDragPoint(new DragPoint(it,*(it+3),it+4,it+7));
	}
	return mpd;
    }
}

//FoldDrag
void FoldDrag::preprocess(std::vector<SpringVertex*>& pts) {
    if (!first) return;
    first = false;
    PLANE rotate_plane(spinOrig, static_plane->dir);

    for (size_t i = 0; i < pts.size(); ++i)
    {
	SpringVertex* sv = pts[i];
        double *p =sv->getCoords();
    	if (static_plane->isBack(p)) 
	{
	    sv->point_type = STATIC_POINT;
	    sv->setRegistered();
    	}
    	else if (rotate_plane.isFront(p)) 
	{
	    sv->point_type = ROTATE_POINT;
            sv->setRegistered();
    	}
        else
        {
	    sv->point_type = FREE_POINT;
	    sv->unsetRegistered();
     	}
    }
}

void FoldDrag::setVel(SpringVertex* sv) {
    double *p0 = sv->getCoords();
    if (sv->point_type == ROTATE_POINT) {
	double p1[3];
	std::copy(p0,p0+3,p1);
	if (m_dt < EPS)
	    std::cout << "Warning: m_dt = " << m_dt
		      << " is too small. "
		      << "Have you called setTimeStepSize() before using it?" 
		      << std::endl;
	double theta = angVel * m_dt;
	spinToAxis(spinOrig,spinDir,theta,p1);
	double *v = sv->getVel();
	for (size_t i = 0; i < 3; ++i) {
	    v[i] = (p1[i] - p0[i])/m_dt;
	}
    }
    else if (sv->point_type == STATIC_POINT) {
	std::fill(sv->getVel(),sv->getVel()+3,0);
    }
    else 
	return;	
}

void FoldDrag::setAccel(SpringVertex* sv) {}

Drag* FoldDrag::clone(const Drag::Info& info) {
    if (!validateData(info)) 
        return NULL;
    else {
        FoldDrag* td = new FoldDrag();
        const double* it = &(info.data().front());
        std::copy(it,it+3,td->spinOrig);
        std::copy(it+3,it+6,td->spinDir);
        td->angVel = *(it+6);
	td->static_plane = new PLANE(it+7, it+10);
        td->m_t = *(it+13);
        return td;
    }
}

//ZFoldDrag
void ZFoldDrag::preprocess(std::vector<SpringVertex*>& pts)
{
    if (!first) return; 
    first = false; 

    for (size_t i = 0; i < pts.size(); i++)
    {
	 SpringVertex* sv = pts[i]; 
	 double *p = sv->getCoords(); 

         if (static_plane1->isBack(p))
         {
	     sv->point_type = STATIC_POINT; 
	     sv->setRegistered(); 
	 }
	 else if (static_plane1->isFront(p) && 
		static_plane2->isBack(p))
	 {
	     sv->point_type = ROTATE_POINT; 
	     sv->setRegistered(); 
	 }
	 else 
 	 {
	     sv->point_type = FREE_POINT; 
	     sv->unsetRegistered(); 	
	 }
    }
}

void ZFoldDrag::setVel(SpringVertex* sv)
{
    if (sv->point_type == STATIC_POINT)
	std::fill(sv->getVel(), sv->getVel() + 3, 0.0);
    else if (sv->point_type == ROTATE_POINT)
    {
	if (m_dt < EPS)
            std::cout << "Warning: m_dt = " << m_dt
                      << " is too small. "
                      << "Have you called setTimeStepSize() before using it?"
                      << std::endl;
	double *p = sv->getCoords(); 
	double p1[3]; 
	double theta = angVel * m_dt; 

	std::copy(p, p + 3, p1); 
	spinToAxis(spinOrig, spinDir, theta, p1);
	std::transform(p1, p1 + 3, p, sv->getVel(), std::minus<double>()); 
	std::transform(sv->getVel(), sv->getVel() + 3, sv->getVel(), 
		bind2nd(std::divides<double>(), m_dt)); 
    }
    else return; 
}

void ZFoldDrag::setAccel(SpringVertex* sv) {}

Drag* ZFoldDrag::clone(const Drag::Info& info)
{
    if (!validateData(info)) 
        return NULL;
    else
    {
	ZFoldDrag *zfd = new ZFoldDrag(); 
	const double *it = &(info.data().front()); 
	std::copy(it, it + 3, zfd->spinOrig); 
	std::copy(it + 3, it + 6, zfd->spinDir);
	zfd->angVel = *(it + 6); 
	zfd->static_plane1 =  new PLANE (it + 7, it + 10); 
	zfd->static_plane2 = new PLANE (it + 13, it + 16); 
	zfd->m_t = *(it + 19); 
	return zfd; 
    }
}

static bool isPointInBall(double p[], double c[], double r) {
    double v[3];
    std::transform(p,p+3,c,v,std::minus<double>());
    return (Mag3d(v) < r) ? true : false;
}

//closeUmbrellaDrag
void CloseUmbrellaDrag::preprocess(std::vector<SpringVertex*>& pts) {
    if (!first) 
	return;
    first = false;
    for (size_t i = 0; i < pts.size(); ++i)
    {
        double d_theta = 2*M_PI/num_tuck_line;
    	double theta = 0;
	SpringVertex* sv = pts[i];
    	double *p = sv->getCoords();
    	double len = distance_between(p,spinOrig,3);

    	if (!isPointInBall(p,spinOrig,1.1*getTolerance())) 
	{
            theta = (p[1] > spinOrig[1]) ? acos((p[0]-spinOrig[0])/len) :
				       acos((p[0]-spinOrig[0])/len) + M_PI;
	    double theta_l = floor(theta/d_theta)*d_theta;
	    double theta_r = ceil(theta/d_theta)*d_theta;
	    double dist = std::min(len*fabs(sin(theta-theta_l)),len*fabs(sin(theta-theta_r)));
	    if (fabs(dist) < 1.1*getTolerance()) 
	    {
	        sv->point_type = ROTATE_POINT;
	        sv->setRegistered();
	    }
	    else 
	    {
	        sv->point_type = FREE_POINT;
		sv->unsetRegistered();
	    }
    	}
    	else 
	{
	    sv->point_type = STATIC_POINT;
	    sv->setRegistered();
    	}
    }
}

void CloseUmbrellaDrag::setVel(SpringVertex* sv) {
    double *p0 = sv->getCoords();
    double *v = sv->getVel();
    double p1[3];
    std::copy(p0,p0+3,p1);
    double rotateDir[3];
    double p_to_o[3];
    std::transform(p0,p0+3,spinOrig,p_to_o,std::minus<double>());

    if (sv->point_type == ROTATE_POINT) 
    {
        double theta = angVel * m_dt;
	Cross3d(p_to_o,spinDir,rotateDir);
	if (Mag3d(rotateDir) < EPS) return;
        spinToAxis(spinOrig,rotateDir,theta,p1);
        for (size_t i = 0; i < 3; ++i)
            v[i] = (p1[i] - p0[i])/m_dt;
    }
    else
        return;
}

Drag* CloseUmbrellaDrag::clone(const Drag::Info& info) {
    if (!validateData(info)) 
        return NULL;
    else {
        CloseUmbrellaDrag* cud = new CloseUmbrellaDrag();
        const double* it = &(info.data().front());
	cud->num_tuck_line = *it;
	std::copy(it+1,it+4,cud->spinOrig);	
	std::copy(it+4,it+7,cud->spinDir);	
	cud->angVel = *(it+7);
        cud->m_t = *(it+8);
        return cud;
    }
}

//RelaxDrag
Drag* RelaxDrag::clone(const Drag::Info& info) {
    RelaxDrag* rd = new RelaxDrag();
    const double* it = &(info.data().front());
    rd->m_t = *it;
    return rd;
}

//GravityBoxDrag
Drag* GravityBoxDrag::clone(const Drag::Info& info) {
    if (!validateData(info)) 
        return NULL;
    else {
        const std::vector<double>& v = info.data();
        const double *it = &(v.front());
	GravityBoxDrag* gbd = new GravityBoxDrag();
	std::copy(it, it+3, gbd->L);
	std::copy(it+3, it+6, gbd->U);
	std::copy(it+6, it+9, gbd->g);
	gbd->m_t = *(it+9);
        return gbd;
    }
}

void GravityBoxDrag::preprocess(std::vector<SpringVertex*>& sv) {
    if (!first) return;
    first = false;
    for (size_t i = 0; i < sv.size(); ++i)
    {
        double *p = sv[i]->getCoords();
	bool is_out_of_box = false;
        for (size_t j = 0; j < 3; ++j)
	    if (p[j] < L[j] || p[j] > U[j])
		is_out_of_box = true;;
    	if (is_out_of_box)
	    sv[i]->unsetRegistered();
	else
	    sv[i]->setRegistered();
    }
}

void GravityBoxDrag::setAccel(SpringVertex* sv) {
    double *a = sv->getExternalAccel();
    if (!sv->isRegistered())
        std::copy(g,g+3,a);
}

//CompressDrag
//add force to compress the surface from several direction
Drag* CompressDrag::clone(const Drag::Info& info)
{
    if (!validateData(info))
      	return NULL;
    else
    {
	const std::vector<double> & v = info.data();
	CompressDrag* cd = new CompressDrag();
	for (size_t i = 0; i < 3; ++i)
	{
	    cd->center[i] = v[i];
	    cd->accel[i] = v[i+3];
	}
	cd->thickness = v[6];
	cd->m_t = v.back();
	return cd;
    }
}

void CompressDrag::preprocess(std::vector<SpringVertex*>& sv)
{
    if (!first) return;
    first = false;
    for (size_t i = 0; i < sv.size(); ++i)
    {
        double* p = sv[i]->getCoords();
        double v[3] = {0};
        for (size_t j = 0; j < 3; ++j)
        {
	    v[j] = p[j] - center[j];
        }
        if (fabs(Dot3d(v, accel)/Mag3d(accel)) < thickness*0.5)
	{
	    sv[i]->setRegistered();
	}
        else
	    sv[i]->unsetRegistered();
    }
}

void CompressDrag::setAccel(SpringVertex* sv)
{
    double* a = sv->getExternalAccel();
    std::fill(a, a+3, 0);
    if (sv->isRegistered())
	return;

    double v[3] = {0};
    double* p = sv->getCoords();
    for (size_t j = 0; j < 3; ++j)
	v[j] = p[j] - center[j];
    double dist = Dot3d(v, accel)/Mag3d(accel);
    if (dist < -0.5*thickness) 
    {
	//back of the plane, push in the same direction
	std::copy(accel, accel+3, a);
    }
    else if (dist > 0.5*thickness)
    {
	//front of the plane, push in opposite direction
	for (int i = 0; i < 3; ++i)
	    a[i] = -accel[i];
    }
}

//RelaxDrag
void RelaxDrag::preprocess(std::vector<SpringVertex*>& pts)
{
    if (!first) return;
    first = false;
    double L = std::numeric_limits<double>::max();
    double S = std::numeric_limits<double>::min();
    double min_coords[3] = {L, L, L};
    double max_coords[3] = {S, S, S};
    SpringVertex* b_pts[6];

    for (size_t i = 0; i < pts.size(); i++)
    {
	SpringVertex* sv = pts[i]; 

	sv->unsetRegistered();
	for (int j = 0; j < 3; ++j)
	{
	    if (sv->getCoords()[j] < min_coords[j])
	    {
		b_pts[j] = sv;
		min_coords[j] = sv->getCoords()[j];
	    }
	    if (sv->getCoords()[j] > max_coords[j])
	    {
		b_pts[3+j] = sv;
		max_coords[j] = sv->getCoords()[j];
	    }
	}
    }
}

//SeparateDrag
//add force to separate the surface from several direction
Drag* SeparateDrag::clone(const Drag::Info& info)
{
    if (!validateData(info))
      	return NULL;
    else
    {
	const std::vector<double> & v = info.data();
	SeparateDrag* cd = new SeparateDrag();
	for (size_t i = 0; i < 3; ++i)
	{
	    cd->spin_center[i] = v[i];
	    cd->spin_dir[i] = v[i+3];
	    cd->nor[i] = v[i+6];
	}
	cd->angVel = v[9];
	cd->m_t = v.back();
	return cd;
    }
}

static void computeCenter(std::vector<SpringVertex*>& sv, double* c)
{
    if (c == NULL)
    {
	std::cout << "memory is not allocated" << std::endl;
    }
    std::fill(c, c+3, 0.0);
    for (size_t i = 0; i < sv.size(); i++)
    {
	SpringVertex* p = sv[i];
	double* coords = p->getCoords();

	for (int j = 0; j < 3; ++j)
	{
	    c[j] += coords[j];
	}
    }
    for (int j = 0; j < 3; ++j)
    {
	c[j] /= (double)sv.size();
    }
}

void SeparateDrag::preprocess(std::vector<SpringVertex*>& sv)
{
    if (!first) return;
    first = false;
    double max_radius = 0.0;
    for (size_t i = 0; i < sv.size(); ++i)
    {
	double r = pointToLine(sv[i]->getCoords(), spin_center, spin_dir);
	max_radius = std::max(max_radius, r);
    }
    for (size_t i = 0; i < sv.size(); ++i)
    {
	if (pointToLine(sv[i]->getCoords(), spin_center, spin_dir) < 0.8*max_radius)
       	    sv[i]->unsetRegistered();
	else
	    sv[i]->setRegistered();
    }
    computeCenter(sv, old_body_center);
}

void SeparateDrag::postprocess(std::vector<SpringVertex*>& sv)
{
    double body_center[3] = {0};
    computeCenter(sv, body_center);
    for (size_t i = 0; i < sv.size(); i++)
    {
	SpringVertex* p = sv[i];
	double* coords = p->getCoords();
	for (int j = 0; j < 3; ++j)
	{
	    coords[j] = coords[j] - body_center[j] + old_body_center[j];
	}
    }
}

void SeparateDrag::setVel(SpringVertex* sv)
{
    if (!sv->isRegistered())
	return;

    double* vel = sv->getVel();
    double* p = sv->getCoords();
    double v1[3] = {0};
    pointToLine(p, spin_center, spin_dir, v1);
    double d_theta = m_dt * angVel;
    
    double cross_v[3]; 
    Cross3d(nor, spin_dir, cross_v);

    double p1[3];
    std::copy(p, p+3, p1);
    double theta = angleBetween(v1, nor);
    double eps = M_PI/36.0;
    if (theta < M_PI/2.0 && theta > eps)
    {
	//spin towards 0
	if (angleBetween(cross_v, v1) > M_PI/2.0)
	{
	    d_theta *= -1;
	}
   	spinToAxis(spin_center, spin_dir, d_theta, p1);
    }
    else if (theta >= M_PI/2.0 && theta < M_PI-eps)
    {
	//spin towards PI
	if (angleBetween(cross_v, v1) < M_PI/2.0)
	{
            d_theta *= -1;
	}
	spinToAxis(spin_center, spin_dir, d_theta, p1);	
    }
    
    for (int i = 0; i < 3; ++i)
	vel[i] = (p1[i] - p[i])/m_dt;    
}

//Roll Drag
Drag* RollDrag::clone(const Drag::Info& info)
{
    const std::vector<double> & v = info.data();
    if (!this->validateData(info))
	return NULL;

    RollDrag* drag = new RollDrag();
    for (size_t i = 0; i < 3; ++i)
    {   
        drag->spin_center[i] = v[i];
        drag->spin_dir[i] = v[i+3];
        drag->mov_center[i] = v[i+6];
        drag->mov_dir[i] = v[i+9];
    }
    drag->ang_vel = v[12];
    drag->m_t = v.back();

    double dist = pointToLine(drag->mov_center, 
			     drag->spin_center, drag->spin_dir);
    drag->num_layers = 3;
    drag->spacing = dist/(drag->num_layers*M_PI*2.0);
    return drag;
}

void RollDrag::preprocess(std::vector<SpringVertex*>& pts)
{
    double bound = spacing * 2.0 * M_PI * num_layers + getTolerance();
    for(size_t i = 0; i < pts.size(); i++)
    {
        SpringVertex* sv = pts[i];
	if (pointToLine(sv->getCoords(), spin_center, spin_dir) < bound)
	{
	     sv->setRegistered();
	     sv->point_type = ROTATE_POINT;
	}
	else
	{
	    sv->unsetRegistered();
	    sv->point_type = FREE_POINT;
	}
    }
}

void RollDrag::setVel(SpringVertex* sv)
{
    if (sv->point_type != ROTATE_POINT)
    {
	return;
    }
    double* p0 = sv->getCoords();
    double p1[3];
    std::copy(p0,p0+3,p1);
    double d_theta = ang_vel * m_dt;
    double v[3] = {0};
    spinToAxis(spin_center,spin_dir,d_theta,p1);
    double r = pointToLine(p1, spin_center, spin_dir, v);
    r = r - spacing * fabs(d_theta);
    
    double v_mag = Mag3d(v);
    for (size_t i = 0; i < 3; ++i)
    {
	p1[i] = p1[i] - v[i] + r * v[i]/v_mag;
    }
    double *vel = sv->getVel();
    for (size_t i = 0; i < 3; ++i) {
        vel[i] = (p1[i] - p0[i])/m_dt;
    }
}

//AlignDrag
//Align the geometric center to required direction
Drag* AlignDrag::clone(const Drag::Info& info)
{
    if (!validateData(info))
        return NULL;
    else {
        const std::vector<double>& v = info.data();
        const double *it = &(v.front());
        return new AlignDrag(it, it+3);
    }
}

AlignDrag::AlignDrag(const double* rotate_center,
                     const double* dir)
{
    rotate_angle = 0; 
    std::copy(rotate_center, rotate_center+3, this->rotate_center);
    std::copy(dir, dir+3, this->dir);
    m_t = 1;
}

void AlignDrag::preprocess(std::vector<SpringVertex*>& pts)
{
    //calculate geomeric center
    std::fill(gravity_center, gravity_center+3, 0);
    for (int i = 0; i < (int)pts.size(); i++)
    {
        SpringVertex* sv = pts[i];  
	sv->setRegistered();
        double* crds = sv->getCoords();
        for (int i = 0; i < 3; ++i)
            gravity_center[i] += crds[i];
    }
    for (int i = 0; i < 3; ++i)
            gravity_center[i] /= pts.size();

    double v1[3] = {0};
    double v2[3] = {0};
    for (int i = 0; i < 3; ++i)
    {
        v1[i] = gravity_center[i] - rotate_center[i];
        v2[i] = dir[i];
    }
    Cross3d(v1, v2, rotate_axis);
    if (Mag3d(rotate_axis) < 1e-15)
    {
        rotate_axis[0] = v2[1];
        rotate_axis[1] = -v2[0];
        rotate_axis[2] = 0;
    }
    Normalize(rotate_axis);
    rotate_angle = std::acos(std::max(
                   std::min(Dot3d(v1, v2)/(Mag3d(v1)*Mag3d(v2)), 1.0), -1.0));
    computeRotationMatrix(rotate_axis, rotate_angle, R);
}

void AlignDrag::postprocess(std::vector<SpringVertex*>& pts)
{
    for (int i = 0; i < (int)pts.size(); i++)
    {
	SpringVertex* sv = pts[i]; 
        double crds[3], new_crds[3];
        std::copy(sv->getCoords(), sv->getCoords()+3, crds);
        std::copy(sv->getCoords(), sv->getCoords()+3, new_crds);
        for (int i = 0; i < 3; ++i)
            crds[i] -= rotate_center[i];
        M3xV3(R, crds, new_crds);
        for (int i = 0; i < 3; ++i)
        sv->getCoords()[i] = new_crds[i] + rotate_center[i];
    }
    m_t = 0; // terminate drag
}

