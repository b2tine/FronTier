/***************************************************************
FronTier is a set of libraries that implements different types of 
Front Traking algorithms. Front Tracking is a numerical method for 
the solution of partial differential equations whose solutions have 
discontinuities.  

Copyright (C) 1999 by The University at Stony Brook. 

This library is free software; you can redistribute it and/or
modify it under the terms of the GNU Lesser General Public
lcense as published by the Free Software Foundation; either
version 2.1 of the License, or (at your option) any later version.

This library is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public
License along with this library; if not, write to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
****************************************************************/


/*
*				example14.c:
*
*		User initialization example for Front Package:
*
*	Copyright 1999 by The University at Stony Brook, All rights reserved.
*	
*	This example generate a grenade type of interface.
*
*/

#include <FronTier.h>

	/*  Function Declarations */
static double grenade_func(POINTER,double*);


char *in_name,*restart_state_name,*restart_name,*out_name;
boolean RestartRun;
int RestartStep;

/********************************************************************
 *	Level function parameters for the initial interface 	    *
 ********************************************************************/
//
typedef struct {
	/*parameters for generating two ellipsoids  */
        double cen[3];
	double a1;
	double a2;
	double c1;
	double c2;
	/*parameters of cylinder coordinates for generating 
	  certain part of the grenade  */
	double phi1; // left angle
	double phi2; // right angle
	double z1;   // upper bound
	double z2;   // lower bound
} GRENADE_PARAMS;

int main(int argc, char **argv)
{
	static Front front;
	static RECT_GRID comp_grid;
	static F_BASIC_DATA f_basic;
	GRENADE_PARAMS s_params;
	static LEVEL_FUNC_PACK level_func_pack;
	static VELO_FUNC_PACK velo_func_pack;
	Locstate  sl;

	f_basic.dim = 3;	
	FT_Init(argc,argv,&f_basic);

	/* Initialize basic computational data */

	f_basic.L[0] = 0.0;	f_basic.L[1] = 0.0; 	f_basic.L[2] = 0.0;
	f_basic.U[0] = 2.0;	f_basic.U[1] = 2.0; 	f_basic.U[2] = 2.0;
	f_basic.gmax[0] = 100;	f_basic.gmax[1] = 100; f_basic.gmax[2] = 100;
	f_basic.boundary[0][0] = f_basic.boundary[0][1] = PERIODIC_BOUNDARY;
	f_basic.boundary[1][0] = f_basic.boundary[1][1] = PERIODIC_BOUNDARY;
	f_basic.boundary[2][0] = f_basic.boundary[2][1] = PERIODIC_BOUNDARY;
	f_basic.size_of_intfc_state = 0;

        in_name                 = f_basic.in_name;
        restart_state_name      = f_basic.restart_state_name;
        out_name                = f_basic.out_name;
        restart_name            = f_basic.restart_name;
        RestartRun              = f_basic.RestartRun;
        RestartStep             = f_basic.RestartStep;

        RestartRun              = f_basic.RestartRun;
        RestartStep             = f_basic.RestartStep;

        sprintf(restart_name,"%s/intfc-ts%s",restart_name,
                        right_flush(RestartStep,7));
        if (pp_numnodes() > 1)
            sprintf(restart_name,"%s-nd%s",restart_name, 
                                right_flush(pp_mynode(),4));

	FT_StartUp(&front,&f_basic);

	if (!RestartRun)
	{
	    /* Initialize interface through level function */
	    s_params.cen[0] = 1.0;
	    s_params.cen[1] = 1.0;
	    s_params.cen[2] = 1.0;


	    s_params.a1=0.5;
	    s_params.a2=0.4;
            s_params.c1=0.8;
	    s_params.c2=0.64;
 	
	    s_params.z1=0.2;
            s_params.z2=0.4;
	    s_params.phi1=PI*(4.0/3.0);
	    s_params.phi2=(11.0/6.0)*PI;
   
	    level_func_pack.neg_component = 1;
	    level_func_pack.pos_component = 2;
	    level_func_pack.func_params = (POINTER)&s_params;
	    level_func_pack.func = grenade_func;
	    level_func_pack.wave_type = FIRST_PHYSICS_WAVE_TYPE;

	    FT_InitIntfc(&front,&level_func_pack);
	}

	FT_Draw(&front);	 

	clean_up(0);
}

static double grenade_func(
        POINTER func_params,
        double *coords)
{
	/* parameters for the ellipoid */
        GRENADE_PARAMS *e_params = (GRENADE_PARAMS*)func_params;
        const double *cen;                   
	double x0,y0,z0;
	double x,y,z;
	/* upper bound and lower bound of cylinder coordinates of the 
	   certain parts of the grenade */
	double a1,c1,a2,c2;
	/* left angle and right angle of the cylinder coordinates of 
	   the certain parts fothe grenade */
	double phi1, phi2;
	double z1,z2;
 
	x=coords[0];
	y=coords[1];
	z=coords[2];
	
	x0 = e_params->cen[0];
        y0 = e_params->cen[1];
        z0 = e_params->cen[2];
	a1=e_params->a1;
	c1=e_params->c1;
	a2=e_params->a2;
	c2=e_params->c2;
	phi1=e_params->phi1;
	phi2=e_params->phi2;
	z1=e_params->z1;
	z2=e_params->z2;

	
	double dist1;
	double dist2;

	double theta;	
	int    i;
	double arg;

        dist1=sqr(x-x0)/sqr(a1)+sqr(y-y0)/sqr(a1)+sqr(z-z0)/sqr(c1)-1.0;
	dist2=sqr(x-x0)/sqr(a2)+sqr(y-y0)/sqr(a2)+sqr(z-z0)/sqr(c2)-1.0;

	
	if (x-x0 >= 0.0 && y-y0 >= 0.0)
	{
	    theta=asin((y-y0)/sqrt(sqr(x-x0)+sqr(y-y0)));
	}
	else if (x-x0 >= 0.0 && y-y0 <= 0.0)
	{
	    theta = 2*PI+asin((y-y0)/sqrt(sqr(x-x0)+sqr(y-y0)));
	}
	else
	{
	    theta = PI-asin((y-y0)/sqrt(sqr(x-x0)+sqr(y-y0)));
        }	


	if (dist1<0 && dist2>0 && z<z0+z2 && z>z0+z1 && 
	    theta<phi2 && theta>phi1)
	{
	    arg=-1;
	}
	else
	{
	    arg=1;
	}
	
	return arg;	
	
}       /* end grenade_func */
