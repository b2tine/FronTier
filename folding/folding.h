#ifndef FOLDING_H
#define FOLDING_H

#include<vector>
#include <set>
#include "../Collision/collid.h"
#include "drag.h"
#include "spring_solver.h"
#include "bending.h"
#include <string>

class Folder {
public:
    virtual void doFolding() = 0;
    virtual ~Folder(){}
    void addDrag(Drag* drag) { drags.push_back(drag);}
    void addDragsFromFile(std::string);
    void setThickness(double h){this->m_thickness = h;}
    void setFrameStepSize(double dt) {max_dt = dt;}
    virtual void setupMovie(std::string, std::string, double) = 0;
    void setOdeScheme(SpringSolver::ODE_SCHEME scheme) {ode_scheme = scheme;}
    void setSpringParameters(double, double, double);
    SpringSolver::SpringParameter& getSpringParams(){return spring_params;}
    SpringSolver::ODE_SCHEME getOdeScheme() {return ode_scheme;}
    void setParaFromFile(const char*);
    void setInputFile(std::string& s) { inname = s; }
    std::string& getInputFile() { return inname; }
protected:
    Folder(){ SpringSolver::ODE_SCHEME ode_scheme = SpringSolver::EXPLICIT; }
    double getThickness(){return m_thickness;}
    double getFrameStepSize() {return max_dt;}
    std::vector<Drag*> drags;
    SpringSolver::SpringParameter spring_params;
private:
    static double max_dt;
    static double m_thickness;
    std::string inname; 
    SpringSolver::ODE_SCHEME ode_scheme;
};

struct Movie;
class Folder3d:public Folder {
public:
    Folder3d(INTERFACE*, SURFACE*);
    ~Folder3d();
    void doFolding();
    Folder3d(){}
    void setupMovie(std::string, std::string, double);
    void check_force(SpringSolver*);
private:
    void doFolding(Drag*,SpringSolver*,CollisionSolver*);
    void straightenStrings();
    double computePotentialEnergy();
    void recordData(double, std::string);
    void appendDataToFile(double, double, std::string);
    std::set<std::string> dataFileSet;
    INTERFACE* m_intfc;
    Movie* movie;
};

struct Movie {
    double mv_dt;
    int mv_count;
    std::string mv_gv_dir;
    std::string mv_vtk_dir; 
    std::string out_dir;
    INTERFACE* mv_intfc;
    bool doMakeMovie;
    void recordMovieFrame();
    bool isMovieTime(double t);
    Movie() : mv_dt(0),mv_count(1),doMakeMovie(false) {}
};

//helper functions
extern void FT_Intfc2SpringMesh(INTERFACE*, std::vector<SpringVertex*>&);
#endif
