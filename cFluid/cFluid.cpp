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


static void gas_driver(Front*,G_CARTESIAN*);
static int g_cartesian_vel(POINTER,Front*,
        POINT*,HYPER_SURF_ELEMENT*,HYPER_SURF*,double*);
static boolean compare_with_base_data(Front *front);
static void rgb_init(Front*,RG_PARAMS);

char *in_name,*restart_state_name,*restart_name,*out_name;
int RestartStep;
boolean RestartRun;
boolean ReadFromInput;

int main(int argc, char **argv)
{
	static Front front;
	static F_BASIC_DATA f_basic;
	static LEVEL_FUNC_PACK level_func_pack;
	static VELO_FUNC_PACK velo_func_pack;
	static EQN_PARAMS eqn_params;
	static RG_PARAMS rgb_params;
    char test_name[100];


	FT_Init(argc,argv,&f_basic);
	f_basic.size_of_intfc_state = sizeof(STATE);

    in_name = f_basic.in_name;
    restart_state_name = f_basic.restart_state_name;
    out_name = f_basic.out_name;
    restart_name = f_basic.restart_name;
    RestartRun = f_basic.RestartRun;
    RestartStep = f_basic.RestartStep;

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

	if (debugging("trace"))
        printf("Passed FT_StartUp()\n");

	
    eqn_params.dim = f_basic.dim;
	read_cFluid_params(in_name,&eqn_params);
	
    if (debugging("trace"))
        printf("Passed read_cFluid_params()\n");

	if (eqn_params.use_base_soln == YES)
	{
        for(int i = 0; i < f_basic.dim; ++i)
            eqn_params.f_basic->subdomains[i] = f_basic.subdomains[i];
	}

    

	front.extra1 = (POINTER)&eqn_params;

	/* Initialize interface through level function */

	if (!RestartRun)
	{
        cFluid_InitIntfc(&front,&level_func_pack);
	    
        if (f_basic.dim == 3)
            level_func_pack.set_3d_bdry = YES;
	    
        FT_InitIntfc(&front,&level_func_pack);
	    insert_objects(&front);
	    rgb_init(&front,rgb_params);
	    FT_PromptSetMixedTypeBoundary2d(in_name,&front);
 
        read_dirichlet_bdry_data(in_name,&front);
	    read_open_end_bdry_data(in_name,&front);
	    
        if (f_basic.dim < 3)
	    	FT_ClipIntfcToSubdomain(&front);
	    
        FT_RedistMesh(&front);
	    
        if (debugging("trace"))
                (void) printf("Passed FT_InitIntfc()\n");
	    
        if (debugging("init_intfc"))
        {
            switch (f_basic.dim)
            {
                case 2:
                    sprintf(test_name,"init_intfc-%d.xg",pp_mynode());
                    xgraph_2d_intfc(test_name,front.interf);
                    break;
                case 3:
                    gview_plot_interface("gv-init",front.interf);
                    break;
            }
        }

	}
	else
    {
	    restart_set_dirichlet_bdry_function(&front);
    }



    //Compressible fluid solver
    G_CARTESIAN	g_cartesian(&front);


    if (RestartRun)
	{
	    readFrontStates(&front,restart_state_name);
	    g_cartesian.readInteriorStates(restart_state_name);
	}
	else
	{
	    g_cartesian.setInitialStates();
	}

	if (debugging("trace"))
	    printf("Passed state initialization()\n");


	/* Initialize velocity field function */
	front._compute_force_and_torque = cfluid_compute_force_and_torque;
	velo_func_pack.func_params = (POINTER)&g_cartesian;
	velo_func_pack.func = g_cartesian_vel;
	velo_func_pack.point_propagate = cFluid_point_propagate;
	FT_InitFrontVeloFunc(&front,&velo_func_pack);

	if (debugging("trace"))
	    printf("Passed FT_InitFrontVeloFunc()\n");

    if (debugging("sample_velocity"))
        initSampleVelocity(&front,in_name);

	/* Propagate the front */

	gas_driver(&front,&g_cartesian);

	clean_up(0);
}

void gas_driver(Front *front,
        G_CARTESIAN* g_cartesian)
{
    double CFL;

	Curve_redistribution_function(front) = full_redistribute;

	FT_ReadTimeControl(in_name,front);
	CFL = Time_step_factor(front);

	if (!RestartRun)
	{
	    FT_ResetTime(front);

	    FrontPreAdvance(front);
	    FT_Propagate(front);

        g_cartesian->solve();

	    FT_Save(front);
        g_cartesian->printFrontInteriorStates(out_name);
        
        if (compare_with_base_data(front))
        {
            g_cartesian->compareWithBaseData(out_name);
            g_cartesian->freeBaseFront();
        }

        FT_Draw(front);

	    FT_SetTimeStep(front);
	    front->dt = std::min(front->dt,CFL*g_cartesian->max_dt);
	    FT_SetOutputCounter(front);
        
    }
    else
    {
        FT_SetOutputCounter(front);
    }

	FT_TimeControlFilter(front);
	FT_PrintTimeStamp(front);

	if (debugging("trace"))
	{
	    printf("CFL = %f\n",CFL);
	    printf("Frequency_of_redistribution(front,GENERAL_WAVE) = %d\n",
                Frequency_of_redistribution(front,GENERAL_WAVE));
	}

    /* Propagating interface for time step dt */
	if (debugging("trace"))
        printf("Entering time loop\n");

    for (;;)
    {

        start_clock("time_loop");
        print_storage("Storage at start of time step","trace");
        
        if (debugging("trace"))
            printf("Begin a time step\n");
        
        FrontPreAdvance(front);
        FT_Propagate(front);

        g_cartesian->solve();
        
        if (debugging("trace")) 
            print_storage("Storage after time step","trace");
        
        FT_AddTimeStepToCounter(front);
        FT_SetTimeStep(front);

        if (debugging("step_size"))
        {
            (void) printf("Step size from front:    %20.14f\n",front->dt);
            (void) printf("Step size from interior: %20.14f\n",
                    CFL*g_cartesian->max_dt);
        }

        front->dt = std::min(front->dt,CFL*g_cartesian->max_dt);

        /* Output section */

        start_clock("output");
        if (FT_IsSaveTime(front))
        {
            FT_Save(front);
            g_cartesian->printFrontInteriorStates(out_name);
        
            if (compare_with_base_data(front))
            {
                //TODO: remedy static bool memory allocation
                //      in compareWithBaseData()
                g_cartesian->compareWithBaseData(out_name);
                g_cartesian->freeBaseFront();
            }
        }
        
        if (FT_IsDrawTime(front))
            FT_Draw(front);

        stop_clock("output");

        if (FT_TimeLimitReached(front))
        {
            start_clock("exit-output");
            if (!FT_IsSaveTime(front))
            {
                FT_Save(front);
                g_cartesian->printFrontInteriorStates(out_name);
            }
    
            if (!FT_IsDrawTime(front))
                FT_Draw(front);

            FT_PrintTimeStamp(front);
            stop_clock("exit-output");
            stop_clock("time_loop");
            break;
        }

        FT_TimeControlFilter(front);
        FT_PrintTimeStamp(front);
        stop_clock("time_loop");

    }

	if (FT_Dimension() == 1)
    {
        g_cartesian->errFunction();
    }

	if (debugging("trace"))
        printf("After time loop\n");

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
	((G_CARTESIAN*)params)->getVelocity(coords, vel);
	return YES;
}	/* end g_cartesian_vel */

static boolean compare_with_base_data(Front *front)
{
	EQN_PARAMS *eqn_params = (EQN_PARAMS*)front->extra1;
	return eqn_params->use_base_soln;
}	/* end compare_with_base_data */

static void rgb_init(Front *front,
	RG_PARAMS rgb_params)
{
	CURVE **c;
	SURFACE **s;

	if (FT_Dimension() == 1) return;
	else if (FT_Dimension() == 2)
	{
	    for (c = front->interf->curves; c && *c; ++c)
	    {
		if (wave_type(*c) == MOVABLE_BODY_BOUNDARY)
		{
		    prompt_for_rigid_body_params(front->f_basic->dim,
				front->f_basic->in_name,&rgb_params);
		    body_index(*c) = 0;
		    set_rgbody_params(rgb_params,Hyper_surf(*c));
		}
	    }
	}
	else
	{
	    for (s = front->interf->surfaces; s && *s; ++s)
	    {
		if (wave_type(*s) == MOVABLE_BODY_BOUNDARY)
		{
		    prompt_for_rigid_body_params(front->f_basic->dim,
				front->f_basic->in_name,&rgb_params);
		    body_index(*s) = 0;
		    set_rgbody_params(rgb_params,Hyper_surf(*s));
		}
	    }
	}
} /* end rgb_init */
