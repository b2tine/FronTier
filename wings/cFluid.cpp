/***************************************************************
FronTier is a set of libraries that implements different types of 
Front Traking algorithms. Front Tracking is a numerical method for 
the solution of partial differential equations whose solutions 
have discontinuities.  

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
*				cFluid.c:
*
*	Copyright 1999 by The University at Stony Brook, All rights reserved.
*	
*/

#include "cFluid.h"

	/*  Function Declarations */
static void gas_driver(Front*,G_CARTESIAN&);
static int g_cartesian_vel(POINTER,Front*,POINT*,HYPER_SURF_ELEMENT*,
                        HYPER_SURF*,double*);
static boolean compare_with_base_data(Front *front);

char *in_name,*restart_state_name,*restart_name,*out_name;
boolean RestartRun;
boolean ReadFromInput;
int RestartStep;

int main(int argc, char **argv)
{
	static Front front;
	static F_BASIC_DATA f_basic;
	static LEVEL_FUNC_PACK level_func_pack;
	static VELO_FUNC_PACK velo_func_pack;
	static EQN_PARAMS eqn_params;
	int i;

	G_CARTESIAN	g_cartesian(front);


	/* Initialize basic computational data */

	FT_Init(argc,argv,&f_basic);
	f_basic.size_of_intfc_state = sizeof(STATE);
	
        in_name                 = f_basic.in_name;
        restart_state_name      = f_basic.restart_state_name;
        out_name                = f_basic.out_name;
        restart_name            = f_basic.restart_name;
        RestartRun              = f_basic.RestartRun;
        RestartStep             = f_basic.RestartStep;

        sprintf(restart_state_name,"%s/state.ts%s",restart_name,
                        right_flush(RestartStep,7));
        sprintf(restart_name,"%s/intfc-ts%s",restart_name,
			right_flush(RestartStep,7));
	if (pp_numnodes() > 1)
	{
            sprintf(restart_name,"%s-nd%s",restart_name,
			right_flush(pp_mynode(),4));
            sprintf(restart_state_name,"%s-nd%s",restart_state_name,
                        right_flush(pp_mynode(),4));
	}

	FT_ReadSpaceDomain(in_name,&f_basic);
	FT_StartUp(&front,&f_basic);
	FT_InitDebug(in_name);
	if (debugging("sample_velocity"))
            g_cartesian.initSampleVelocity(in_name);

	if (debugging("trace")) printf("Passed FT_StartUp()\n");

	eqn_params.dim = f_basic.dim;
	read_cFluid_params(in_name,&eqn_params);

	if (eqn_params.use_base_soln == YES)
	{
            for (i = 0; i < f_basic.dim; ++i)
                eqn_params.f_basic->subdomains[i] = f_basic.subdomains[i];
	}

	front.extra1 = (POINTER)&eqn_params;
	read_movie_options(&front);
	if (debugging("trace")) printf("Passed read_cFluid_params()\n");

	/* Initialize interface through level function */

	level_func_pack.pos_component = GAS_COMP1;
	if (!RestartRun)
	{
	    FT_InitIntfc(&front,&level_func_pack);
	    initModules(&front);
	    initStateParams(&front);
	    read_dirichlet_bdry_data(in_name,&front);
	    if (f_basic.dim < 3)
	    	FT_ClipIntfcToSubdomain(&front);
	}
	else
	{
	    initStateParams(&front);
	    restart_set_dirichlet_bdry_function(&front);
	}
	initMotionParams(&front);

	/* Initialize velocity field function */

	front._compute_force_and_torque = preset_wing_motion;
	velo_func_pack.func_params = (POINTER)&g_cartesian;
	velo_func_pack.func = g_cartesian_vel;
	velo_func_pack.point_propagate = cFluid_point_propagate;
	FT_InitFrontVeloFunc(&front,&velo_func_pack);
	if (debugging("trace"))
	    printf("Passed FT_InitFrontVeloFunc()\n");

	g_cartesian.initMesh();
	if (RestartRun)
	{
	    readFrontStates(&front,restart_state_name);
	    g_cartesian.readInteriorStates(restart_state_name);
	}
	else
	{
	    g_cartesian.initFrontInteriorStates();
	}

	/* Propagate the front */

	gas_driver(&front, g_cartesian);

	//PetscFinalize();
	clean_up(0);
}

static  void gas_driver(
        Front *front,
	G_CARTESIAN &g_cartesian)
{
        double CFL;

	Curve_redistribution_function(front) = full_redistribute;

	FT_ReadTimeControl(in_name,front);
	CFL = Time_step_factor(front);

	if (!RestartRun || ReSetTime(front))
	{
	    FT_ResetTime(front);

	    FrontPreAdvance(front);
	    FT_Propagate(front);
	    g_cartesian.solve(front->dt);

	    FT_Save(front);
            g_cartesian.printFrontInteriorStates(out_name);
            if (compare_with_base_data(front))
            {
                g_cartesian.compareWithBaseData(out_name);
                g_cartesian.freeBaseFront();
            }
            g_cartesian.initMovieVariables();
            FT_Draw(front);
	    recordWingVar(front);

	    FT_SetTimeStep(front);
	    front->dt = std::min(front->dt,CFL*g_cartesian.max_dt);
	    FT_SetOutputCounter(front);
        }
        else
	    FT_SetOutputCounter(front);

	FT_TimeControlFilter(front);
	FT_PrintTimeStamp(front);

	if (debugging("trace"))
	{
	    printf("CFL = %f\n",CFL);
	    printf("Frequency_of_redistribution(front,GENERAL_WAVE) = %d\n",
			Frequency_of_redistribution(front,GENERAL_WAVE));
	}

	if (debugging("trace")) printf("Before time loop\n");
        for (;;)
        {
            /* Propagating interface for time step dt */

	    start_clock("time_loop");
	    print_storage("Storage at start of time step","trace");
	    if (debugging("trace")) printf("Begin a time step\n");
	    FrontPreAdvance(front);
	    FT_Propagate(front);

	    g_cartesian.solve(front->dt);
	    if (debugging("trace")) 
	    {
		print_storage("Storage after time step","trace");
	    }

	    FT_AddTimeStepToCounter(front);
				
            //Next time step determined by maximum speed of previous
            //step, assuming the propagation is hyperbolic and
            //is not dependent on second order derivatives of
            //the interface such as curvature, and etc.

	    FT_SetTimeStep(front);
	    if (debugging("step_size"))
	    {
		(void) printf("Step size from front:    %20.14f\n",front->dt);
		(void) printf("Step size from interior: %20.14f\n",
					CFL*g_cartesian.max_dt);
	    }
            front->dt = std::min(front->dt,CFL*g_cartesian.max_dt);
	
            /* Output section */

            if (FT_IsSaveTime(front))
	    {
            	FT_Save(front);
		g_cartesian.printFrontInteriorStates(out_name);
		if (compare_with_base_data(front))
		{
		    g_cartesian.compareWithBaseData(out_name);
		    g_cartesian.freeBaseFront();
		}
	    }
            if (FT_IsDrawTime(front))
	    {
	        g_cartesian.initMovieVariables();
            	FT_Draw(front);
	    	recordWingVar(front);
	    }

            if (FT_TimeLimitReached(front))
	    {
		if (!FT_IsSaveTime(front))
		{
            	    FT_Save(front);
		    g_cartesian.printFrontInteriorStates(out_name);
		}
		if (!FT_IsDrawTime(front))
		{
                    g_cartesian.initMovieVariables();
                    FT_Draw(front);
	    	    recordWingVar(front);
		}
	    	FT_PrintTimeStamp(front);
	    	stop_clock("time_loop");
                break;
	    }
	    FT_TimeControlFilter(front);
	    FT_PrintTimeStamp(front);
	    stop_clock("time_loop");
        }
	if (debugging("trace")) printf("After time loop\n");
}       /* end gas_driver */

static int g_cartesian_vel(
	POINTER params,
	Front *front,
	POINT *p,
	HYPER_SURF_ELEMENT *hse,
	HYPER_SURF *hs,
	double *vel)
{
	double *coords = Coords(p);
	((G_CARTESIAN_EB*)params)->getVelocity(coords, vel);
	return YES;
}	/* end g_cartesian_vel */

static boolean compare_with_base_data(Front *front)
{
	EQN_PARAMS *eqn_params = (EQN_PARAMS*)front->extra1;
	return eqn_params->use_base_soln;
}	/* end compare_with_base_data */
