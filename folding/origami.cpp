#include "origami.h"
#include <iostream>
#include <iterator>
#include <map>
#include <fstream>
static std::ostream& operator << (std::ostream& os, const std::vector<double>&);
static std::ostream& operator << (std::ostream& os, const std_matrix&);

struct Compare {
    Compare(std::vector<double> xs, Vertex* v) : x_axis(xs), vtx(v) {}
    bool operator()(const Crease* c1, const Crease* c2) const {
	double a1 = Math::angleBetweenWithDir(x_axis, c1->dir, vtx->nor); 
	double a2 = Math::angleBetweenWithDir(x_axis, c2->dir, vtx->nor); 
	return a1 < a2; 
    }
private : 
    std::vector<double> x_axis;
    Vertex* vtx; 
}; 
//class OgmPoint 
OgmPoint::OgmPoint(SpringVertex& sv) : SpringVertex(sv)
{
    x0.resize(3,0);
    std::copy(sv.getCoords(), sv.getCoords()+3, x0.begin());
}

//class OrigamiFold
void OrigamiFold :: preprocess(std::vector<SpringVertex*>& pts)
{
    //perform only once
    if (first)
    {
	first = false;
    	//set all points registered
	//replace all spring vertices with new OgmPoints
	for (size_t i = 0; i < pts.size(); ++i)
	{
	    SpringVertex* sv = new OgmPoint(*pts[i]);
	    delete pts[i];
	    pts[i] = sv;
	    pts[i]->setRegistered();
	}
    	//assign id to vertex, each vertex uniquely relates to a crease
    	assignVertexId(pts);
    }
    //recompute folding matrix for each crease
    findNextFoldingAngle();
    int start = 0;
    for (size_t i = 0; i < vertices.size(); ++i)
    {
	Vertex* v = vertices[i];
	std::vector<double> myrho(rho_delta.begin() + start,
				  rho_delta.begin() + start + v->creases.size());
	v->updateCreaseFoldingMatrix(myrho); 
	start += v->creases.size();
    }
    
}

void OrigamiFold::assignVertexId(std::vector<SpringVertex*>& pts)
{
    sortCreases(pts[rand()%pts.size()]->getCoords());
    //assign crease to each spring point
    for (int i = 0; i < (int)pts.size(); i++)
    {
	SpringVertex* sv = pts[i]; 
	for (int k = 0; k < (int)vertices.size(); k++)
	{
	    Vertex* vtx = vertices[k]; 
	    if (vtx->creases.size() == 0)
		continue;
	    if (vtx->creases.size() == 1)
	    {
		OgmPoint* og_pt = static_cast<OgmPoint*>(sv);
                og_pt->crs.push_back(vtx->creases.back());
		continue;
	    }

	    std::vector<double> v0(sv->getCoords(), sv->getCoords()+3);
	    Math::V3mV3(v0, vtx->crds, v0);
	    double alpha0 = Math::angleBetweenWithDir(
					vtx->creases[0]->dir,
					v0, vtx->nor);
	    if (Math::Mag(v0) < 1e-12)
	    {
		OgmPoint* og_pt = static_cast<OgmPoint*>(sv);
                og_pt->crs.push_back(vtx->creases.back());
	    }
	    else if (alpha0 > Math::angleBetweenWithDir(
					vtx->creases[0]->dir,
					vtx->creases.back()->dir, vtx->nor))
	    {
		OgmPoint* og_pt = static_cast<OgmPoint*>(sv);
                og_pt->crs.push_back(vtx->creases.back());
	    }
	    else
	    for (int j = 1; j < vtx->creases.size(); ++j)
	    {
		double alpha1 = Math::angleBetweenWithDir(
					vtx->creases[0]->dir,	
					vtx->creases[j-1]->dir, vtx->nor);
		double alpha2 = Math::angleBetweenWithDir(
					vtx->creases[0]->dir,
					vtx->creases[j]->dir, vtx->nor);
		if (alpha0 >= alpha1 && alpha0 <= alpha2)
		{
		    OgmPoint* og_pt = static_cast<OgmPoint*>(sv);
		    og_pt->crs.push_back(vtx->creases[j-1]);
		}
	    }
	    if (static_cast<OgmPoint*>(sv)->crs.size() == 0)
	    {
		std::cout << "alpha0 = " << alpha0 << std::endl;
	    	std::cout << "ERROR: crease not assigned" << std::endl;
	    	exit(-1);
	    }
	}
    }
}

std::vector<double> OrigamiFold::findFoldable(std::vector<double>& input, double& err)
{
    if (m_opt == NULL)
    {
	size_t N = rho.size();
	m_opt = new nlopt::opt(nlopt::LN_COBYLA, N);
	this->m_opt->set_stopval(1e-6);
        this->m_opt->set_ftol_abs(1e-6);
        this->m_opt->set_xtol_rel(1e-6);
        this->m_opt->set_maxeval(3000);  //maxiteration
	nlopt::srand(1234);
	std::vector<double> lowerbounds(N, -M_PI);
        std::vector<double> upperbounds(N, M_PI);

	for (size_t i = 0; i < N; ++i)
	{
	    if (rho_T[i] > 0)
		lowerbounds[i] = 0;
	    if (rho_T[i] < 0)
		upperbounds[i] = 0;
	}
	m_opt->set_lower_bounds(lowerbounds); 
	m_opt->set_upper_bounds(upperbounds);
	this->m_opt->set_min_objective(staticTargetFunction, this);
    }    
    std::vector<double> ans = m_opt->optimize(input);
    err = targetFunction(ans);
    return ans;
}

double OrigamiFold::staticTargetFunction(
		const std::vector<double>& x,
                std::vector<double> &grad,
                void *my_func_data)
{
    OrigamiFold* obj = static_cast<OrigamiFold*>(my_func_data);
    return obj->targetFunction(x);
}

double OrigamiFold::targetFunction(
		const std::vector<double>& x)
{
    double ans = 0; 
    int start = 0;
    for (size_t i = 0; i < vertices.size(); ++i)
    {
	Vertex* v = vertices[i];
	std::vector<double> myrho(x.begin() + start, 
				  x.begin() + start + v->creases.size());
	v->updateCreaseFoldingMatrix(myrho);
	std_matrix M = Math::Mat(3, 3);
	std_matrix I = Math::Eye(3);
        Math::M3mM3(v->creases.back()->fd_matrix, I, M);
	ans += Math::Norm(M, 1);
	start += v->creases.size();
    }
    return ans;
}

std::ostream& operator << (std::ostream& os, const std::vector<double>& v)
{
    std::copy(v.begin(), v.end(), std::ostream_iterator<double>(os, " "));
    return os;
}

std::ostream& operator << (std::ostream& os, const std_matrix& m)
{
    os << "[" << std::endl;
    for (int i = 0; i < (int)m.size(); i++)
    {
	std::vector<double> v = m[i]; 
    	std::copy(v.begin(), v.end(), std::ostream_iterator<double>(os, " "));
	os << std::endl;
    }
    os << "]" << std::endl;
    return os;
}

void OrigamiFold::findNextFoldingAngle()
{
    const int maxIter = 2000;
    int iter = 0;
    bool success = false;
    double w0, w1, w2, w, D;
    w0 = 0.8; w1 = 0.2; w2 = 0.01; w = w0; D = 0.015;
    while (!success && iter++ < maxIter)
    {
	int N = rho.size();
	std::vector<double> rho_rand(N, 0);
	std::vector<double> dir(N, 0);
	for (int i = 0; i < N; ++i)
	{
	    rho_rand[i] = (rand()/(double)RAND_MAX * 2 - 1) * M_PI;
	    dir[i] = (1 - w) * rho_rand[i] + w * rho_T[i];
	    rho_tau[i] = rho_delta[i] + D * dir[i];
	    //hard limit for input
	    rho_tau[i] = std::min(M_PI, 
			 std::max(rho_tau[i], -M_PI));
	}

	double err = 0;
	rho = findFoldable(rho_tau, err);
	if (isValid(rho), Math::dist3d(rho_T, rho) < Math::dist3d(rho_T, rho_delta))
	{
	    rho_delta = rho;
	    w = w + w1;	
	    success = true;
	    std::cout << "rho = " << rho << std::endl;
	}
	else
	{
	    w = w - w2;
	    success = false;
	}
	w = std::max(std::min(w, 1.0), 0.0);
    }    
    if (iter >= maxIter)
    {
	std::cout << "Maxium iteration reached\n" 
	<< "folding angle = [" << rho_delta << "]" 
	<< std::endl;
	std::cout << "err = " << targetFunction(rho_delta) << std::endl;
	std::cout << "Folding process is terminated" 
	<< std::endl;
	m_t = 0;
    }
}

bool OrigamiFold::isValid(const std::vector<double>& rho)
{
    for (int i = 0; i < (int)rho.size(); i++)
    {
	double a = rho[i]; 
	if (fabs(a) > M_PI)
	    return false;
    }
    return true;
}

void OrigamiFold::postprocess(std::vector<SpringVertex*>& pts) {}

void OrigamiFold:: setVel(SpringVertex* sv)
{  
    std::vector<double> new_crds(3, 0); 
    ogmComputeNewPosition(sv, new_crds);
    double* vel = sv->getVel();
    for (int i = 0; i < 3; ++i)
	vel[i] = (new_crds[i] - sv->getCoords()[i])/getTimeStepSize();
}

void OrigamiFold::setAccel(SpringVertex*){}

size_t OrigamiFold::dataSize()
{
    return totalDataSize;
}

Drag* OrigamiFold::clone(const Info & info)
{
    if (info.data().size() == 0)
    {
	std::cout << "Warning: OrigamiFold is not created, no data provided"
		  << std::endl;
	return NULL;
    }
    if (info.data().size() > 0 && info.data().size() < 2)
    {
	std::cout << "Warning: OrigamiFold is not created, " 
		  << "insufficient data is given"
		  << std::endl;
    }
    //parse info to create drag object
    std::vector<double> v = info.data();
    std::vector< std::vector<double> > points;
    std::vector< std::vector<int> > creases; 
    std::vector<double> angles;
    size_t n_pt = (size_t)v[0];
    size_t n_crs = (size_t)v[1];
    totalDataSize = 2 + n_pt * 3 + n_crs * 2 + n_crs * 1;
    
    if(!validateData(info))
	return NULL;

    double *it = &(v.front()) + 2;
    //create points
    for (int i = 0; i < n_pt; ++i)
    {
	points.push_back(std::vector<double> (it, it + 3));
	it += 3;
    }
    //create creases
    for (int i = 0; i < n_crs; ++i)
    {
        creases.push_back(std::vector<int> (it, it + 2));
        it += 2;
    }

    //create angles
    for (int i = 0; i < n_crs; ++i)
    {
	angles.push_back(*it);
  	it ++;
    }
    return new OrigamiFold(points, creases, angles);
}

OrigamiFold::OrigamiFold(): m_opt(NULL), first(true) {}

OrigamiFold::OrigamiFold(const std::vector< std::vector<double> >& points, 
			 const std::vector< std::vector<int> >& creases,
			 const std::vector<double>& angles)
{
    first = true; 
    m_t = 2.0;
    m_opt = NULL;
    std::map<int, Vertex*> vertex_map;
    if (creases.size() != angles.size())
    {
	std::cerr << "Number of creases does not equal to number of angles" 
		  << std::endl; 
    }

    size_t N = creases.size();
    //allocate memory for rho
    rho.resize(N, 0);
    rho_delta.resize(N, 0);
    rho_tau.resize(N, 0);

    //create creases and vertices
    for (size_t i = 0; i < creases.size(); ++i)
    {
	int vindex = creases[i][0];
	int pindex = creases[i][1];
	Vertex* vtx = NULL;
	if (vertex_map.find(vindex) == vertex_map.end())
	{
	    vtx = new Vertex(points[vindex]);
	    vertices.push_back(vtx);
	    vertex_map[vindex] = vtx;
	}
	else
	{
	    vtx = vertex_map[vindex];
	}
	Crease* crs = new Crease(points[vindex], 
				     points[pindex], angles[i]);
	crs->vtx = vtx;
	vtx->creases.push_back(crs);
    }

}

void OrigamiFold::sortCreases(double* rand_crds)
{
    //compute normal vector
    for (int i = 0; i < (int)vertices.size(); i++)
    {
	Vertex* vtx = vertices[i]; 
	for (int j = 0; j < 3; j++)
    	     vtx->nor.push_back(0);
	if (vtx->creases.size() == 0)
	{
	
	    vtx->nor[2] = 1;
	}
	else if (vtx->creases.size() > 1)
	{
	    Math::Cross3d(vtx->creases[0]->dir, vtx->creases[1]->dir, vtx->nor);
	}
	if (Math::Mag(vtx->nor) < 1e-15)
	{
	    //one crease line or two parallel crease lines
	    std::vector<double> rand_v(rand_crds, rand_crds+3);
	    Math::V3mV3(rand_v, vtx->crds, rand_v);
	    Math::Cross3d(vtx->creases[0]->dir, rand_v, vtx->nor);
	}
	Math::Normalize(vtx->nor);
    }


    //sort creases
    for (int i = 0; i < (int)vertices.size(); i++)
    {
	Vertex* vtx = vertices[i]; 
        std::vector<double> x_axis = vtx->creases[0]->dir;
	struct Compare customLess(x_axis, vtx); 

	sort(vtx->creases.begin(), vtx->creases.end(), customLess); 
    }
    //assemble rho_T according to the order of creases
    for (int i = 0; i < (int)vertices.size(); i++)
    {
	Vertex* vtx = vertices[i]; 
	for (int j = 0; j < (int)vtx->creases.size(); j++)
        {
	    Crease* crs = vtx->creases[j]; 
	    rho_T.push_back(crs->rho_T);
     	}
    } 

}

void OrigamiFold::ogmComputeNewPosition(SpringVertex* sv, std::vector<double>& new_crds)
{
    OgmPoint* ov = static_cast<OgmPoint*>(sv);
    
    for (int i = 0; i < ov->crs.size(); i++)
    {
	Crease* crs = ov->crs[i]; 
        std::copy(ov->x0.begin(), ov->x0.end(), new_crds.begin());
        crs->crsFoldCrds(new_crds);
    }
}

//destructor
OrigamiFold::~OrigamiFold() {
    if (this->m_opt)
    {
	delete this->m_opt;
	this->m_opt = NULL;
    }
    for (int i = 0; i < vertices.size(); i++)
    {
	Vertex* v = vertices[i]; 
	for (int j = 0; j < v->creases.size(); j++)
	{
	    Crease* c = v->creases[j]; 
	    delete c;
	}
	delete v;
    }
}

//Class Crease
Crease::Crease(const std::vector<double>& p1, 
	       const std::vector<double>& p2, 
	       const double angle)
{
    dir.resize(3, 0);
    for (int i = 0; i < 3; ++i)
    {
	dir[i] = p2[i] - p1[i];
    }
    Math::Normalize(dir);
    rho_T = angle;
    fd_matrix = Math::Mat(3,3);
}

void Crease::crsFoldCrds(std::vector<double>& crds)
{
     std::vector<double> ans(3, 0);
     const std::vector<double>& center = vtx->crds;
     Math::V3mV3(crds, center, crds);
     Math::M3xV3(fd_matrix, crds, ans);
     Math::V3pV3(ans, center, crds);
}

void Crease::crsRotationMatrix(std::vector< std::vector<double> >& R, 
				double theta)
{
    double u,v,w;
    //dir should be normalized
    u = dir[0]; v = dir[1]; w = dir[2];

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

void Crease::crsRotationMatrixDerivative(
		std::vector< std::vector<double> >& dR, double theta)
{
    double u,v,w;
    //dir should be normalized
    u = dir[0]; v = dir[1]; w = dir[2];

    dR[0][0] = -(1-u*u)*sin(theta);
    dR[0][1] = u*v*sin(theta) - w*cos(theta);
    dR[0][2] = u*w*sin(theta)+v*cos(theta);
             
    dR[1][0] = u*v*sin(theta) + w*cos(theta);
    dR[1][1] = -(1-v*v)*sin(theta);
    dR[1][2] = v*w*sin(theta) - u*cos(theta);
             
    dR[2][0] = u*w*sin(theta) - v*cos(theta);
    dR[2][1] = v*w*sin(theta) + u*cos(theta);
    dR[2][2] = -(1-w*w)*sin(theta);
}

void Crease::updateFoldingMatrix(const std::vector< std::vector<double> >& M)
{
    for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 3; ++j)
    {
	fd_matrix[i][j] = M[i][j];
    }
}

Vertex::Vertex(const std::vector<double>& coords) : crds(coords)
{}

void Vertex::updateCreaseFoldingMatrix(const std::vector<double>& rho)
{
    if (creases.size() == 0)
	return;
    std::vector< std::vector<double> > M(4,std::vector<double>(4,0)), ans = M;

    Crease* crs = creases[0];
    crs->crsRotationMatrix(M, rho[0]);
    crs->updateFoldingMatrix(M);

    for (int i = 1; i < creases.size(); ++i)
    {
	crs = creases[i];
	crs->crsRotationMatrix(M, rho[i]);
	Math::M3xM3(creases[i-1]->fd_matrix, M, ans);	
	crs->updateFoldingMatrix(ans);
    }
}

//Math class
/*
void Math::MatxMat(const std::vector<std::vector<double>>& M1, 
		const std::vector<std::vector<double>>& M2, int dim, 
		std::vector<std::vector<double>>& ans) {
    for (int i = 0; i < dim; i++)
	 for (int j = 0; j < dim; j++)
	      for (int k = 0; k < dim; k++)
	           ans[i][j] += M1[i][k] * M2[k][i];
}

*/
void Math::M3xM3(const std::vector< std::vector<double> > & M1, 
		   const std::vector< std::vector<double> > & M2,
		   std::vector< std::vector<double> > & ans)
{
    for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 3; ++j)
    {
	ans[i][j] = 0.0;
	for (int k = 0; k < 3; ++k)
	    ans[i][j] += M1[i][k] * M2[k][j];
    }
}

void Math::M3xV3(const std::vector< std::vector<double> >& M,
                   const std::vector<double>& v,
                   std::vector<double>& ans)
{
    for (int i = 0; i < 3; ++i)
    {
        ans[i] = 0;
        for (int j = 0; j < 3; ++j)
            ans[i] += M[i][j] * v[j];
    }
}

void Math::M3pM3(const std::vector< std::vector<double> >& M1,
                 const std::vector< std::vector<double> >& M2,
                 std::vector< std::vector<double> > & ans)
{
    for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 3; ++j)
    {
        ans[i][j] = M1[i][j] + M2[i][j];
    }
}

void Math::M3mM3(const std::vector< std::vector<double> >& M1,
                 const std::vector< std::vector<double> >& M2,
                 std::vector< std::vector<double> > & ans)
{
    for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 3; ++j)
    {
        ans[i][j] = M1[i][j] - M2[i][j];
    }
}

void Math::V3mV3(const std::vector<double>& v1,
                 const std::vector<double>& v2,
                       std::vector<double>& ans)
{
    for (size_t i = 0; i < v1.size(); ++i)
	ans[i] = v1[i] - v2[i];
}

void Math::V3pV3(const std::vector<double>& v1,
                 const std::vector<double>& v2,
                       std::vector<double>& ans)
{
    for (size_t i = 0; i < v1.size(); ++i)
        ans[i] = v1[i] + v2[i];
}

double Math::angleBetween(const std::vector<double>& v1, 
			  const std::vector<double>& v2)
{
    double dp = 0;
    for (int i = 0; i < 3; ++i)
    {
	dp += v1[i]*v2[i];
    }
    double m1 = Math::Mag(v1);
    double m2 = Math::Mag(v2);
    return acos(std::max(std::min(dp/(m1*m2), 1.0), -1.0));
}

double Math::angleBetweenWithDir(
		const std::vector<double>& v1,
                const std::vector<double>& v2,
		const std::vector<double>& nor)
{
    //angle increases from v1 to v2
    std::vector<double> n = nor;
    Normalize(n);
    double det = TriProd(n, v1, v2);
    double dot = Math::Dot3d(v1, v2);
    return det >= 0 ? atan2(det, dot) : 2*M_PI + atan2(det, dot);
}

void Math::Cross3d(const std::vector<double>& u,
                     const std::vector<double>& v,
		     std::vector<double>& ans)
{
    ans[0] = u[1]*v[2] - u[2]*v[1];
    ans[1] = u[2]*v[0] - u[0]*v[2];
    ans[2] = u[0]*v[1] - u[1]*v[0];
}

double Math::TriProd(const std::vector<double>& a,
                     const std::vector<double>& b,
                     const std::vector<double>& c)
{
    std::vector<double> crx(3, 0);
    Math::Cross3d(b,c,crx);
    return Math::Dot3d(a, crx);
}

double Math::Dot3d(const std::vector<double>& a,
                   const std::vector<double>& b)
{
    double ans = 0;
    for (size_t i = 0; i < a.size(); ++i)
    {
	ans += a[i] * b[i];
    }
    return ans;
}

double Math::Mag(const std::vector<double>& v)
{
    double ans = 0;
    for (int i = 0; i < (int)v.size(); i++)
	ans += v[i]*v[i];
    return sqrt(ans); 
}

std_matrix Math::Eye(int N)
{
    std_matrix res(N, std::vector<double>(N, 0));
    for (int i = 0; i < N; ++i)
	res[i][i] = 1.0;
    return res;
}

std_matrix Math::Mat(int m, int n)
{
   return std_matrix(m, std::vector<double>(n, 0));
}

double Math::Norm(const std_matrix& m, int p)
{
   if (p == 0)
	return 1.0;
   double res = 0;
   for (int i = 0; i < (int)m.size(); i++) {
	std::vector<double> v = m[i]; 
        for (int j = 0; j < (int)v.size(); j++)
        {
	     double a = v[j]; 
	     res += std::pow(fabs(a), p);
        }
   }
   return std::pow(res, 1/p);
}

double Math::Norm(const std::vector<double>& v, int p)
{
    if (p == 0)
	return 1.0;
    double res = 0;
    for (int i = 0; i < (int)v.size(); i++)
	res += std::pow(fabs(v[i]), p);
    return std::pow(res, 1/p);
}

double Math::dist3d(const std::vector<double>& p1, 
		    const std::vector<double>& p2)
{
    double ans = 0;
    for (size_t i = 0; i < p1.size(); ++i)
    {
	ans += (p1[i]-p2[i])*(p1[i]-p2[i]);
    }
    return sqrt(ans);
}

void Math::Normalize(std::vector<double>& v)
{
   double vm = Mag(v);
   if (vm == 0)
	return;
   for (size_t i = 0; i < v.size(); ++i)
	v[i] /= vm;	
}
