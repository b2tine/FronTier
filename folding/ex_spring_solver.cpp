#include "spring_solver.h"
#include <cmath>
#include <iostream>

void EX_SPRING_SOLVER::doSolve(double t) {
    const size_t size = pts.size();
    if (x_old.size() < pts.size()) {
        x_old.resize(size,std::vector<double>(3,0));
        x_new.resize(size,std::vector<double>(3,0));
        v_old.resize(size,std::vector<double>(3,0));
        v_new.resize(size,std::vector<double>(3,0));
    }

    size_t i, j;
    const int dim = 3;

    double dt = getTimeStepSize();
    int n_loop = t/dt+1;
    dt = t/n_loop;

    std::cout << "Starting spring solver:\n";
    for (int n = 0; n < n_loop; ++n)
    {
	for (i = 0; i < ext_forces.size(); ++i)
	    ext_forces[i]->computeExternalForce();

        std::cout << "    #sub_step  = " << n+1 << ' ' << n_loop << std::endl;
        for (i = 0; i < size; ++i)
            computeAccel(pts[i]);

        for (i = 0; i < size; ++i)
        for (j = 0; j < dim; ++j)
        {
            x_old[i][j] = pts[i]->x[j];
            v_old[i][j] = pts[i]->v[j];
        }

        for (i = 0; i < size; ++i)
        for (j = 0; j < dim; ++j)
        {
            x_new[i][j] = x_old[i][j] + dt*pts[i]->v[j]/6.0;
            v_new[i][j] = v_old[i][j] + dt*pts[i]->accel[j]/6.0;
            pts[i]->x[j] = x_old[i][j] + 0.5*pts[i]->v[j]*dt;
            pts[i]->v[j] = v_old[i][j] + 0.5*pts[i]->accel[j]*dt;
        }

        for (i = 0; i < size; ++i)
            computeAccel(pts[i]);

        for (i = 0; i < size; ++i)
        for (j = 0; j < dim; ++j)
        {
            x_new[i][j] += dt*pts[i]->v[j]/3.0;
            v_new[i][j] += dt*pts[i]->accel[j]/3.0;
            pts[i]->x[j] = x_old[i][j] + 0.5*pts[i]->v[j]*dt;
            pts[i]->v[j] = v_old[i][j] + 0.5*pts[i]->accel[j]*dt;
        }

        for (i = 0; i < size; ++i)
            computeAccel(pts[i]);
	
	for (i = 0; i < size; ++i)
        for (j = 0; j < dim; ++j)
        {
            x_new[i][j] += dt*pts[i]->v[j]/3.0;
            v_new[i][j] += dt*pts[i]->accel[j]/3.0;
            pts[i]->x[j] = x_old[i][j] + pts[i]->v[j]*dt;
            pts[i]->v[j] = v_old[i][j] + pts[i]->accel[j]*dt;
        }

        for (i = 0; i < size; ++i)
            computeAccel(pts[i]);

        for (i = 0; i < size; ++i)
        for (j = 0; j < dim; ++j)
        {
            x_new[i][j] += dt*pts[i]->v[j]/6.0;
            v_new[i][j] += dt*pts[i]->accel[j]/6.0;
        }

        for (i = 0; i < size; ++i)
        for (j = 0; j < dim; ++j)
        {
            pts[i]->x[j] = x_new[i][j];
            pts[i]->v[j] = v_new[i][j];
        }
    }

}

double EX_SPRING_SOLVER::getTimeStepSize() {
    return 0.1*sqrt(springParameter.m
                /springParameter.k);
}

