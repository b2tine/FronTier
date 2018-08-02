/***************************************************************
FronTier is a set of libraries that implements different types of 
Front Traking algorithms. Front Tracking is a numerical method for 
the solution of partial differential equations whose solutions have 
discontinuities.  

Copyright (C) 1999 by The University at Stony Brook. 

This library is free software; you can redistribute it and/or
modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation; either
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
*				example00.c:
*
*		User initialization example for Front Package:
*
*	Copyright 1999 by The University at Stony Brook, All rights reserved.
*	
*	This example shows a slotted disk in rotation. It demonstrates
*	the geometry preservation of the front tracking method.
*
*/

#include <vector>
#include <FronTier.h>

	/*  Function Declarations */
static void split_curve_at_horizontal_line(CURVE*,double);

char *in_name,*restart_state_name,*restart_name,*out_name;
boolean RestartRun;
int RestartStep;

/********************************************************************
 *	Velocity function parameters for the front	 	    *
 ********************************************************************/

int main(int argc, char **argv)
{
	static Front front;
	static F_BASIC_DATA f_basic;
	static LEVEL_FUNC_PACK level_func_pack;
	ELLIP2D_PARAMS ellipse_params;
	CURVE **c,*curve;

	FT_Init(argc,argv,&f_basic);
	f_basic.dim = 2;	

	/* Initialize basic computational data */

	f_basic.L[0] = 0.0;	f_basic.L[1] = 0.0;
	f_basic.U[0] = 1.0;	f_basic.U[1] = 1.0;
	f_basic.gmax[0] = 128;	f_basic.gmax[1] = 128;
	f_basic.boundary[0][0] = f_basic.boundary[0][1] = PERIODIC_BOUNDARY;
	f_basic.boundary[1][0] = f_basic.boundary[1][1] = PERIODIC_BOUNDARY;
	f_basic.size_of_intfc_state = 0;

        in_name                 = f_basic.in_name;
        restart_state_name      = f_basic.restart_state_name;
        out_name                = f_basic.out_name;
        restart_name            = f_basic.restart_name;
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

	    ellipse_params.x0 = 0.5;
	    ellipse_params.y0 = 0.5;
	    ellipse_params.a = 0.4;
	    ellipse_params.b = 0.3;

	    level_func_pack.neg_component = 1;
	    level_func_pack.pos_component = 2;
	    level_func_pack.func_params = (POINTER)&ellipse_params;
	    level_func_pack.func = ellipse_func;
	    level_func_pack.wave_type = FIRST_PHYSICS_WAVE_TYPE;
	    FT_InitIntfc(&front,&level_func_pack);
	    if (f_basic.dim < 3)
                FT_ClipIntfcToSubdomain(&front);
	}

	intfc_curve_loop(front.interf,c)
	{
	    curve = *c;
	    if (wave_type(curve) == FIRST_PHYSICS_WAVE_TYPE) break;
	}
	split_curve_at_horizontal_line(curve,0.5);

	gview_plot_color_interface("test",front.interf,YES);
	xgraph_2d_intfc("intfc.xg",front.interf);

	clean_up(0);
	return 0;
}

static void split_curve_at_horizontal_line(
	CURVE *curve,
	double H)
{
	CURVE **c;
	BOND *b;
	POINT *pts[4];
	int nx;
	double coords[MAXD];
	double *ps,*pe;

	nx = 0;
	curve_bond_loop(curve,b)
	{
	    if ((Coords(b->start)[1] < H && Coords(b->end)[1] > H) ||
		(Coords(b->start)[1] > H && Coords(b->end)[1] < H))
	    {
		/* Interpolating the cut point */
		ps = Coords(b->start);
		pe = Coords(b->end);
		coords[0] = ps[0] + (H - ps[1])/(pe[1] - ps[1])*(pe[0] - ps[0]);
		coords[1] = H;
		pts[nx] = Point(coords);
		insert_point_in_bond(pts[nx],b,curve);
		nx++;
	    }
	    else if (Coords(b->start)[1] == H)
	    {
		/* use existing point */
		pts[nx] = b->start;
		nx++;
	    }
	}
	if (nx != 2)
	{
	    printf("Something wrong!\n");
	    clean_up(ERROR);
	}

	I_MoveNodeToPoint(pts[0],curve); // move the original node to pts[0]
	c = I_SplitCurve(pts[1],curve);  // c is an array of the split curves
}
