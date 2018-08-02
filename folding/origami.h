#ifndef __ORIGAMI_H__
#define __ORIGAMI_H__
#include "spring_solver.h"
#include "drag.h"
#include "nlopt.hpp"

struct Vertex;

struct Crease //folding line
{
    Vertex *vtx;
    double rho_T;  //final rotation angle    
    std::vector<double> dir; //axis direction
    std::vector< std::vector<double> > fd_matrix;
    void crsUpdateFoldingMatrix(double rho);

    //utility functions
    void crsRotateCrdsByAngle(double*, double);
    void crsRotationMatrix(std::vector< std::vector<double> >&, double);
    void crsRotationMatrixDerivative(std::vector< std::vector<double> >&, 
		double);
    void updateFoldingMatrix(const std::vector< std::vector<double> >&);
    void crsFoldCrds(std::vector<double>&);

    //constructor
    Crease(const std::vector<double>&, const std::vector<double>&,const double);
};

struct OgmPoint : public SpringVertex
{
    std::vector<Crease*> crs;
    std::vector<double> x0; //initial coords 
public:
    OgmPoint(SpringVertex&);
};

struct Vertex 
{
    std::vector<double> crds;
    std::vector<Crease*> creases;
    std::vector<double> nor;
    //constructor
    Vertex(const std::vector<double>&);
    void updateCreaseFoldingMatrix(const std::vector<double>& rho);
};

class OrigamiFold : public Drag
{
public:
    void preprocess(std::vector<SpringVertex*>&);
    void postprocess(std::vector<SpringVertex*>&);
    void setVel(SpringVertex*);
    void updateVel(std::vector<SpringVertex*>& pts, double t) {}
    void setAccel(SpringVertex*);
    std::string id() {return "OrigamiDrag";} 
    OrigamiFold(const std::vector< std::vector<double> >&,
		const std::vector< std::vector<int> >&, 
		const std::vector<double>&);
    OrigamiFold();
    ~OrigamiFold();
    size_t dataSize();
private:
    bool first;
    std::vector<Vertex*> vertices;
    std::vector<double> rho, rho_delta, rho_tau, rho_T;
    nlopt::opt *m_opt = NULL; //optimization tool
    void assignVertexId(std::vector<SpringVertex*>&);
    std::vector<double> findFoldable(std::vector<double>&, double&);
    static double staticTargetFunction(const std::vector<double>&, std::vector<double>&, void*);
    double targetFunction(const std::vector<double>&);
    void findNextFoldingAngle();
    void ogmComputeNewPosition(SpringVertex* sv, std::vector<double>&);
    void sortCreases(double*);
    bool isValid(const std::vector<double>&);
    Drag * clone(const Info&);
    size_t totalDataSize;
};

typedef std::vector< std::vector<double> > std_matrix;
class Math
{
public:

    static void M3xM3(const std_matrix &,
                   const std_matrix &,
                   std_matrix &);

    static void M3xV3(const std_matrix &,
                   const std::vector<double>&,
                   std::vector<double>&);

    static void M3pM3(const std_matrix &,
                   const std_matrix &,
                   std_matrix &);

    static void M3mM3(const std_matrix &,
                   const std_matrix &,
                   std_matrix &);
  
    static void V3mV3(const std::vector<double>&,
		      const std::vector<double>&,
		      std::vector<double>&);

    static void V3pV3(const std::vector<double>&,
		      const std::vector<double>&,
		      std::vector<double>&);

    static double angleBetween(const std::vector<double>&,
		               const std::vector<double>&);

    static double angleBetweenWithDir(
                const std::vector<double>& v1,
                const std::vector<double>& v2,
		const std::vector<double>& n);
 
    static double Mag(const std::vector<double>&);

    static std_matrix Eye(int);

    static std_matrix Mat(int, int);

    static double Norm(const std_matrix &, int);
 
    static double Norm(const std::vector<double>&, int);

    static double dist3d(const std::vector<double>& p1, 
			 const std::vector<double>& p2);

    static void Normalize(std::vector<double>& v);

    static double TriProd(const std::vector<double>& a,
                     const std::vector<double>& b,
                     const std::vector<double>& c);

    static double Dot3d(const std::vector<double>& a,
                   const std::vector<double>& b);

    static void Cross3d(const std::vector<double>& u,
                     const std::vector<double>& v,
                     std::vector<double>& ans);

private:
    Math();
};
#endif

