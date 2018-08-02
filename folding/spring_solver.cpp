#include "spring_solver.h"
#include <assert.h>
#include <iterator>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <cmath>
//class SpringVeterx
void SpringVertex::addNeighbor(size_t i_nb, double len0) 
{
    index_nb.push_back(i_nb);
    length0.push_back(len0);
}

SpringVertex::SpringVertex(SpringVertex& spv)
{
    x = spv.x;
    v = spv.v;
    org_vtx = spv.org_vtx;
    std::copy(spv.accel, spv.accel+3, accel);
    std::copy(spv.ext_accel, spv.ext_accel+3, ext_accel);
    index_nb = spv.index_nb;
    length0 = spv.length0;
    is_registered = spv.is_registered;
    point_type = spv.point_type;
}

//class SpringSolver
void SpringSolver::printAdjacencyList(std::string fname) {
    std::ofstream ofs(fname.c_str());
    for (size_t i = 0; i < pts.size(); ++i) {
	ofs << i << " ";
	std::copy(pts[i]->index_nb.begin(),pts[i]->index_nb.end(),
		std::ostream_iterator<size_t>(ofs," "));
	ofs << "\n";
    }
}

void SpringSolver::printPointList(std::string fname) {
    printf("%lu number of points in list\n",pts.size());
    std::ofstream ofs(fname.c_str());
    ofs << std::setprecision(9);
    for (size_t i = 0; i < pts.size(); ++i) {
	double* crds = pts[i]->getCoords();
	std::copy(crds,crds+3,std::ostream_iterator<double>(ofs," "));
	ofs << pts[i]<< "\n";
    }
}

void SpringSolver::checkVertexNeighbors() {
    for (size_t i = 0; i < pts.size(); ++i) {
	SpringVertex* sv = pts[i];
	for (size_t j = 0; j < sv->index_nb.size(); ++j) {
	    SpringVertex* sv_nb = pts[static_cast<size_t>(sv->index_nb[j])];
	    sv_nb->addNeighbor(i,sv->length0[j]);
	}
    }
}
void SpringSolver::setParameters(double k, double lambda, double m) {
    if (k == 0 || m == 0)
	throw std::invalid_argument("tensile stiffness and mass cannot zero");
    springParameter.k = k;
    springParameter.lambda = lambda;
    springParameter.m = m;
}

void SpringSolver::setParameters(SpringSolver::SpringParameter& params) {
    if (params.k == 0 || params.m == 0)
	throw std::invalid_argument("tensile stiffness and mass cannot zero");
    springParameter = params;
}

SpringSolver* SpringSolver::createSpringSolver(ODE_SCHEME scheme) {
    switch (scheme) {
	case EXPLICIT:
	    return new EX_SPRING_SOLVER();
	case IMPLICIT:
	    return new IM_SPRING_SOLVER();
	default:
	    return new EX_SPRING_SOLVER();
    }
}

void SpringSolver::computeAccel(SpringVertex* sv) {
        double k = springParameter.k;
        double m = springParameter.m;
        double lambda = springParameter.lambda;
        const int dim = 3;
        double vec[3] = {0};
        double v_rel[3] = {0};

	//reset acceleration
	std::fill(sv->accel, sv->accel+3, 0.0);

	m_drag->setVel(sv);
	m_drag->setAccel(sv);
        if (sv->isRegistered())
            return;
	
 	//add external acceleration: bending
	for (size_t j = 0; j < ext_forces.size(); ++j)
	{
	    double* ext_f = ext_forces[j]->getExternalForce(sv);   
	    for (size_t i = 0; i < dim; ++i)
                sv->accel[i] += ext_f[i]/m; 
	}

	//add internal acceleration: tensile stiffness
        for (size_t i = 0; i < sv->index_nb.size(); ++i)
        {
            double len = 0.0;
            SpringVertex* sv_nb = pts[sv->index_nb[i]];
            for (int j = 0; j < dim; ++j)
            {
                vec[j] = sv_nb->x[j] - sv->x[j];
                len += vec[j]*vec[j];
                v_rel[j] = sv_nb->v[j] - sv->v[j];
            }
            len = sqrt(len);
            for (int j = 0; j < dim; ++j)
            {
                if (len > SPRING_EPS)
                    vec[j] /= len;
                else
                    vec[j] = 0.0;
                sv->accel[j] += k*((len - sv->length0[i])*vec[j])/m;
                sv->accel[j] += lambda*v_rel[j]/m;
            }
            if (sv->length0[i] < 0) exit(-1);
        }

        for (int j = 0; j < dim; ++j)
        {
            sv->accel[j] += (sv->getExternalAccel())[j];
        }
}

void SpringSolver::computeJacobian() {
}

void SpringSolver::solve(double dt)
{
    if (m_drag) 
    {
	m_drag->setTimeStepSize(dt);
	m_drag->preprocess(pts);
	int count = 0;
	for (int i = 0; i < (int)pts.size(); i++)
	{
	    SpringVertex* sv = pts[i]; 
	    if (sv->isRegistered())
	    {
		count++;
	    }
	}
    }
    doSolve(dt);
    if (m_drag)
    { 
	m_drag->postprocess(pts);
    }
}

void SpringSolver::setDrag(Drag* dg) {
    m_drag = dg;
}

void SpringSolver::resetVelocity() {
    for (size_t i = 0; i < pts.size(); ++i)
    {
        std::fill(pts[i]->v,pts[i]->v+3,0);
        std::fill(pts[i]->ext_accel,pts[i]->ext_accel+3,0);
    }
}
