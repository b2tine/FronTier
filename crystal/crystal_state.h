/**********************************************************************
 * 		crystal_basic.h					      *
 **********************************************************************/

#ifndef _FT_CRYSTAL_BASIC_H_
#define _FT_CRYSTAL_BASIC_H_


struct _STATE {
        double dens;                    /* Density */
        double pres;                    /* Pressure */
        double phi;                     /* Potential */
        double vel[MAXD];               /* Velocities */
        double vort;                    /* Vorticity in 2D */
        double solute;                  /* For subsurface problem */
        double mu;                      /* For eddy viscosity */
};
typedef struct _STATE STATE;

#endif
