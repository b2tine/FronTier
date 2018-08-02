#ifndef CGAL_H_
#define CGAL_H_

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Delaunay_mesher_2.h>
#include <CGAL/Delaunay_mesh_face_base_2.h>
#include <CGAL/Delaunay_mesh_size_criteria_2.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>
#include <iostream>
#include <fstream>
#include <FronTier.h>
#include <vector>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Triangulation_vertex_base_with_info_2<unsigned, K> Vb;
typedef CGAL::Delaunay_mesh_face_base_2<K> Fb;
typedef CGAL::Triangulation_data_structure_2<Vb, Fb> Tds;
typedef CGAL::Exact_predicates_tag Itag;
typedef CGAL::Constrained_Delaunay_triangulation_2<K, Tds,Itag> CDT;
typedef CGAL::Delaunay_mesh_size_criteria_2<CDT> Criteria;

typedef CDT::Vertex_handle Vertex_handle;
typedef CDT::Point Cgal_Point;

bool findAndLocate(std::ifstream&, const char*);

class cgalSurf
{
    INTERFACE* c_intfc;
    SURFACE** c_surf; 
    CDT cdt; 
    int num_finite_face; 
    double k_dx; 
    // used to provide upper bound for minimum angle 
    double c_bound; 
    double height; 
    std::vector<double> extConPoint;
protected :
    void setHeight(double h) { height = h; } 
    double readHeight() const { return height; }
    void setCGALCoeffRestricSize(double k) { k_dx = k; }
    double readCGALCoeffRestricSize() const { return k_dx; }
    void setCGALMinAngleUb(double c) { c_bound = c; }
    double readCGALMinAngleUb() const { return c_bound; }
    Vertex_handle insertPointToCDT(Cgal_Point a) { return cdt.insert(a); }
    void insertConstraintToCDT(Vertex_handle v1, Vertex_handle v2) { 
	cdt.insert_constraint(v1, v2); }
    CDT readCDT() const { return cdt; }
    CDT& readCDT() { return cdt; }
    std::vector<double> readExtConPoint() const { return extConPoint; }
    void setExtConPoint(double p, int index) { extConPoint[index] = p; }
    INTERFACE* readIntface() { return c_intfc; }
    int readNumFinFace() const { return num_finite_face; }
    void setNumFinFace(int num) { num_finite_face = num; }
    SURFACE** readSurface() { return c_surf; }
public :
    // std::ifstream& can't be used in constructor
    // any constructor involving it been added =delete in the end
    cgalSurf(INTERFACE*, SURFACE**);   
    // add extra constraint
    void getExtraConstPoint(std::ifstream&); 
    void setSurfZeroMesh(); 
    void setCurveZeroLength(CURVE*, double); 
    void setMonoCompBdryZeroLength(); 
    // by generated tri mesh generates surface
    void cgalGenSurf(); 
    // add boundary constraint. different for differernt shape
    virtual void addCgalConst() = 0;
    // preprocess and generate tri mesh 
    virtual void cgalTriMesh(std::ifstream&) = 0; 
    virtual void getParaFromFile(std::ifstream&) = 0; 
    virtual ~cgalSurf() {}; 
};

class cgalRectangleSurf : public cgalSurf
{
    double lower[2]; 
    double upper[2]; 
public :
    cgalRectangleSurf(INTERFACE* intfc, SURFACE** surf) :
    		cgalSurf(intfc, surf) {}
    virtual void addCgalConst(); 
    virtual void cgalTriMesh(std::ifstream&); 
    virtual void getParaFromFile(std::ifstream&); 
    virtual ~cgalRectangleSurf() {}
};

class cgalCircleSurf : public cgalSurf
{
    double cen[2]; 
    double radius; 
    int num_reg_const;
    static const double eps;
protected : 
    double* getCenter() { return cen; }
    void setCenter(double c1, double c2) { cen[0] = c1; cen[1] = c2; }
    double getRadius() const { return radius; }
    void setRadius(double r) { radius = r; }
    int getNumRegConst() const { return num_reg_const; }
    void setNumRegConst(int num) { num_reg_const = num; } 
public : 
    cgalCircleSurf(INTERFACE* intfc, SURFACE** surf) : 
		cgalSurf(intfc, surf) { num_reg_const = 150; }
    virtual void addCgalConst();
    virtual void cgalTriMesh(std::ifstream&);
    virtual void getParaFromFile(std::ifstream&);
    double distance(double*, double*, int); 
    virtual ~cgalCircleSurf() {}
};

class cgalParaSurf : public cgalCircleSurf {
    int num_lines; 
public : 
    cgalParaSurf(INTERFACE* intfc, SURFACE** surf) : 
		cgalCircleSurf(intfc, surf) {}
    virtual void addCgalConst(); 
    virtual void getParaFromFile(std::ifstream&); 
    virtual ~cgalParaSurf() {}; 
};

#endif
