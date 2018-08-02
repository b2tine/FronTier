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

#include "inverse.h"
#include "solver.h"

static boolean cim_driver(Front*,C_CARTESIAN&);
static void initCimTestParams(char*,Front*);

char *in_name,*restart_state_name,*restart_name,*out_name;
boolean RestartRun;
boolean ReSetTime;
int RestartStep;

static void CIM_flowThroughBoundaryState(double*,HYPER_SURF*,Front*,
			POINTER,POINTER);
static void CIM_timeDependBoundaryState(double*,HYPER_SURF*,Front*,
                        POINTER,POINTER);
static void read_dirichlet_bdry_data(char*,Front*,F_BASIC_DATA);
static void inverse_point_propagate(Front*,POINTER,POINT*,POINT*,	
			HYPER_SURF_ELEMENT*,HYPER_SURF*,double,double*);
static void neumann_point_propagate(Front*,POINTER,POINT*,POINT*,	
			HYPER_SURF_ELEMENT*,HYPER_SURF*,double,double*);
static void dirichlet_point_propagate(Front*,POINTER,POINT*,POINT*,	
			HYPER_SURF_ELEMENT*,HYPER_SURF*,double,double*);
static void interior_point_propagate(Front*,POINTER,POINT*,POINT*,	
			HYPER_SURF_ELEMENT*,HYPER_SURF*,double,double*);
static void computeDirichletGradient(char*,Front*);
static void integrateGradient(INTERFACE*,BDRY_INTEGRAL*);
static void integrateGradient2d(INTERFACE*,BDRY_INTEGRAL*);
static void integrateGradient3d(INTERFACE*,BDRY_INTEGRAL*);
static double dirichletPointGradient(Front*,POINT*,HYPER_SURF_ELEMENT*,
				HYPER_SURF*);


int main(int argc, char **argv)
{
	static Front front;
        static F_BASIC_DATA f_basic;
        static LEVEL_FUNC_PACK level_func_pack;
	CIRCLE_PARAMS circle_params;
	C_CARTESIAN       c_cartesian(front);
	static CIM_PARAMS cim_params;
	static RADIAL_MOTION_PARAMS rv_params;
	static VELO_FUNC_PACK velo_func_pack;

	FT_Init(argc,argv,&f_basic);
	f_basic.size_of_intfc_state = sizeof(STATE);

        PetscInitialize(&argc,&argv,PETSC_NULL,PETSC_NULL);

	/* Initialize basic computational data */

        in_name                 = f_basic.in_name;
        restart_state_name      = f_basic.restart_state_name;
        out_name                = f_basic.out_name;
        restart_name            = f_basic.restart_name;
        RestartRun              = f_basic.RestartRun;
        RestartStep             = f_basic.RestartStep;
        ReSetTime               = f_basic.ReSetTime;

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
	front.extra1 = (POINTER)&cim_params;
	initCimTestParams(in_name,&front);
	if (debugging("trace"))
            (void) printf("Passed FT_StartUp()\n");

	initCimIntfcParams(in_name,&front,&level_func_pack);

	c_cartesian.w_type = level_func_pack.wave_type;
	c_cartesian.neg_comp = level_func_pack.neg_component;
	c_cartesian.pos_comp = level_func_pack.pos_component;

	if (f_basic.dim == 3) level_func_pack.set_3d_bdry = YES;
	FT_InitIntfc(&front,&level_func_pack);
	read_dirichlet_bdry_data(in_name,&front,f_basic);
        if (f_basic.dim == 2)
	{
	    char xg_name[100];
            FT_ClipIntfcToSubdomain(&front);
	    sprintf(xg_name,"%s/init_intfc-%d",out_name,pp_mynode());
            xgraph_2d_intfc(xg_name,front.interf);
	}
	else if (f_basic.dim == 3)
	{
	    char gv_name[100];
	    sprintf(gv_name,"%s/init_intfc-%d",out_name,pp_mynode());
	    gview_plot_interface(gv_name,front.interf);
	}
	c_cartesian.initMesh();	
	rv_params.dim = f_basic.dim;
	rv_params.cen[0] = 0.0;
	rv_params.cen[1] = 0.0;
	rv_params.speed = 1.0;
	velo_func_pack.func_params = (POINTER)&rv_params;
        velo_func_pack.func = radial_motion_vel;
	velo_func_pack.point_propagate = inverse_point_propagate;
	FT_InitFrontVeloFunc(&front,&velo_func_pack);

	cim_driver(&front,c_cartesian);
	clean_up(0);
}

static boolean cim_driver(
	Front *front,
	C_CARTESIAN &c_cartesian)
{
	CIM_PARAMS *params;
	FIELD *field;
	int count;

	Curve_redistribution_function(front) = full_redistribute;
	FT_RedistMesh(front);

	FT_ResetTime(front);

        FT_Propagate(front);
	c_cartesian.solve();
	computeDirichletGradient(out_name,front);

	params = (CIM_PARAMS*)front->extra1;
	field = params->field;

	if (params->compare_method == CAUCHY_COMPARE)
	    c_cartesian.compareWithBaseSoln();
	else if (params->compare_method == EXACT_COMPARE)
	    c_cartesian.compareWithExacySoln();

	c_cartesian.initMovieVariables();
	FT_Draw(front);
	FT_Save(front);
	c_cartesian.printFrontInteriorStates(out_name);
	viewTopVariable(front,field->u,NO,0.0,0.0,
			out_name,(char*)"solution");

	if (params->prop_intfc == NO)
	    return YES;

	count = 0;
	for (;;)
        {
	    front->dt = 0.02;
            FT_Propagate(front);
	    c_cartesian.solve();
	    FT_AddTimeStepToCounter(front);
	    computeDirichletGradient(out_name,front);
	    c_cartesian.initMovieVariables();
	    FT_Draw(front);
	    if (count++ == 10) break;
	}

	return YES;
}	/* end cim_driver */

static void read_dirichlet_bdry_data(
	char *inname,
	Front *front,
	F_BASIC_DATA f_basic)
{
	char msg[100],s[100];
	int i,dim = front->rect_grid->dim;
	FILE *infile = fopen(inname,"r");
	static STATE *state;
	POINTER func_params;
	HYPER_SURF *hs;
	int i_surf;

	for (i = 0; i < dim; ++i)
	{
	    if (f_basic.boundary[i][0] == DIRICHLET_BOUNDARY)
	    {
		hs = NULL;
	        if (rect_boundary_type(front->interf,i,0) == DIRICHLET_BOUNDARY)
		    hs = FT_RectBoundaryHypSurf(front->interf,DIRICHLET_BOUNDARY,
					i,0);
		sprintf(msg,"For lower boundary in %d-th dimension",i);
		CursorAfterString(infile,msg);
		(void) printf("\n");
		CursorAfterString(infile,"Enter type of Dirichlet boundary:");
		fscanf(infile,"%s",s);
		(void) printf("%s\n",s);
		switch (s[0])
		{
		case 'c':			// Constant state
		case 'C':
		    FT_ScalarMemoryAlloc((POINTER*)&state,sizeof(STATE));
		    CursorAfterString(infile,"Enter boundary u:");
		    fscanf(infile,"%lf",&state->u);
		    (void) printf("%f\n",state->u);
		    FT_InsertDirichletBoundary(front,NULL,NULL,
					NULL,(POINTER)state,hs,i_surf);
		    break;
		case 'f':			// Flow through state
		case 'F':
		    FT_InsertDirichletBoundary(front,CIM_flowThroughBoundaryState,
					"flowThroughBoundaryState",NULL,
					NULL,hs,i_surf);
		    break;
		case 't':			// Flow through state
		case 'T':
		    //get_time_dependent_params(dim,infile,&func_params);
		    FT_InsertDirichletBoundary(front,CIM_timeDependBoundaryState,
					"CIM_timeDependBoundaryState",
					func_params,NULL,hs,i_surf);
		    break;
		default:
		    (void) printf("Unknown Dirichlet boundary!\n");
		    clean_up(ERROR);
		}
	    }
            if (f_basic.boundary[i][1] == DIRICHLET_BOUNDARY)
	    {
		hs = NULL;
                if (rect_boundary_type(front->interf,i,1) == DIRICHLET_BOUNDARY)                    hs = FT_RectBoundaryHypSurf(front->interf,DIRICHLET_BOUNDARY,
                                                i,1);
		sprintf(msg,"For upper boundary in %d-th dimension",i);
		CursorAfterString(infile,msg);
		(void) printf("\n");
		CursorAfterString(infile,"Enter type of Dirichlet boundary:");
		fscanf(infile,"%s",s);
		(void) printf("%s\n",s);
		switch (s[0])
		{
		case 'c':			// Constant state
		case 'C':
		    FT_ScalarMemoryAlloc((POINTER*)&state,sizeof(STATE));
		    CursorAfterString(infile,"Enter boundary u:");
		    fscanf(infile,"%lf",&state->u);
		    (void) printf("%f\n",state->u);
		    FT_InsertDirichletBoundary(front,NULL,NULL,
					NULL,(POINTER)state,hs,i_surf);
		    break;
		case 'f':			// Flow through state
		case 'F':
		    FT_InsertDirichletBoundary(front,CIM_flowThroughBoundaryState,
					"flowThroughBoundaryState",NULL,
					NULL,hs,i_surf);
		    break;
		case 't':			// Flow through state
		case 'T':
		    //get_time_dependent_params(dim,infile,&func_params);
		    FT_InsertDirichletBoundary(front,CIM_timeDependBoundaryState,
					"CIM_timeDependBoundaryState",
					func_params,NULL,hs,i_surf);
		    break;
		default:
		    (void) printf("Unknown Dirichlet boundary!\n");
		    clean_up(ERROR);
		}
	    }
	}
	fclose(infile);
}	/* end read_dirichlet_bdry_data */


static void CIM_flowThroughBoundaryState(
        double          *p0,
        HYPER_SURF      *hs,
        Front           *front,
        POINTER         params,
        POINTER         state)
{
	printf("Entering CIM_flowThroughBoundaryState()\n");
}	/* CIM_flowThroughBoundaryState */

static void CIM_timeDependBoundaryState(
        double          *p0,
        HYPER_SURF      *hs,
        Front           *front,
        POINTER         params,
        POINTER         state)
{
}	/* end CIM_timeDependBoundaryState */

static void initCimTestParams(
	char *inname,
	Front *front)
{
	CIM_PARAMS *cim_params = (CIM_PARAMS*)front->extra1;
	FILE *infile = fopen(inname,"r");
	int i,dim = front->rect_grid->dim;
	char string[100];

	cim_params->dim = dim;
	for (i = 0; i < dim; ++i)
	    cim_params->h[i] = front->rect_grid->h[i];

	CursorAfterString(infile,"Enter jump type:");
	fscanf(infile,"%s",string);
	(void) printf("%s\n",string);
	switch (string[0])
	{
	case 'E':
	case 'e':
	    cim_params->jump_type = EXACT_JUMP;
	    break;
	case 'C':
	case 'c':
	    cim_params->jump_type = CONSTANT_JUMP;
	    break;
	case 'F':
	case 'f':
	    cim_params->jump_type = FUNCTIONAL_JUMP;
	    break;
	}

	if (cim_params->jump_type == EXACT_JUMP)
	{
	    CursorAfterString(infile,"Enter run case number:");
	    fscanf(infile,"%d",&cim_params->Run_case);
	    (void) printf("%d\n",cim_params->Run_case);
	    if (dim == 2)
	    {
	    	if (cim_params->Run_case == 8)
	    	{
		    (void) printf("Run case %d not implemented in 2D!\n",
				cim_params->Run_case);
		    clean_up(ERROR);
	    	}
	    }
	    else if (dim == 3)
	    {
	    	if (cim_params->Run_case == 3 ||
		    cim_params->Run_case == 4 ||
		    cim_params->Run_case == 5 ||
		    cim_params->Run_case == 6 ||
		    cim_params->Run_case == 7) 
	    	{
		    (void) printf("Run case %d not implemented in 3D!\n",
				cim_params->Run_case);
		    clean_up(ERROR);
	    	}
	    }
	}
	else if (cim_params->jump_type == CONSTANT_JUMP)
	{
	    CursorAfterString(infile,"Enter jump of u:");
	    fscanf(infile,"%lf",&cim_params->jump_u);
	    (void) printf("%f\n",cim_params->jump_u);
	    CursorAfterString(infile,"Enter jump of eps gradient u dot n:");
	    fscanf(infile,"%lf",&cim_params->jump_eps_grad_u_dot_n);
	    (void) printf("%f\n",cim_params->jump_eps_grad_u_dot_n);
	    CursorAfterString(infile,"Enter jump of gradient u dot t:");
	    fscanf(infile,"%lf",&cim_params->jump_grad_u_dot_t);
	    (void) printf("%f\n",cim_params->jump_grad_u_dot_t);
	}

	CursorAfterString(infile,"Enter interface number:");
	fscanf(infile,"%d",&cim_params->intfc_num);
	(void) printf("%d\n",cim_params->intfc_num);

	CursorAfterString(infile,"Enter comparison method:");
	fscanf(infile,"%s",string);
	(void) printf("%s\n",string);
	switch (string[0])
	{
	case 'N':
	case 'n':
	    cim_params->compare_method = NO_COMPARE;
	    break;
	case 'E':
	case 'e':
	    cim_params->compare_method = EXACT_COMPARE;
	    if (cim_params->jump_type != EXACT_JUMP)
	    {
		(void) printf("Can't do exact comparison without"
			      " exact jump!\n");
		clean_up(ERROR);
	    }
	    break;
	case 'C':
	case 'c':
	    cim_params->compare_method = CAUCHY_COMPARE;
	    break;
	}

	if (cim_params->compare_method == CAUCHY_COMPARE)
	{
	    CursorAfterString(infile,"Enter base directory name:");
	    fscanf(infile,"%s",cim_params->base_dir_name);
	    (void) printf("%s\n",cim_params->base_dir_name);
	    CursorAfterString(infile,"Enter number of comparing steps:");
	    fscanf(infile,"%d",&cim_params->num_step);
	    (void) printf("%d\n",cim_params->num_step);
	    FT_VectorMemoryAlloc((POINTER*)&cim_params->steps,
				cim_params->num_step,sizeof(int));
	    for (i = 0; i < cim_params->num_step; ++i)
	    {
	    	sprintf(string,"Enter index of step %d:",i+1);
	    	CursorAfterString(infile,string);
	    	fscanf(infile,"%d",&cim_params->steps[i]);
	    	(void) printf("%d\n",cim_params->steps[i]);
	    }
	    FT_ScalarMemoryAlloc((POINTER*)&cim_params->f_basic,
				sizeof(F_BASIC_DATA));
	    cim_params->f_basic->dim = dim;
	    for (i = 0; i < dim; ++i)
		cim_params->f_basic->subdomains[i] = front->pp_grid->gmax[i];
	}
	cim_params->prop_intfc = NO;
    	CursorAfterStringOpt(infile,"Enter yes to propagate interface:");
	fscanf(infile,"%s",string);
	(void) printf("%s\n",string);
	if (string[0] == 'y' || string[0] == 'Y')
	{
	    cim_params->prop_intfc = YES;
	}
	fclose(infile);

	if (cim_params->compare_method == CAUCHY_COMPARE)
	    FT_ReadComparisonDomain(in_name,cim_params->f_basic);
}	/* end initCimTestParams */

static  void inverse_point_propagate(
        Front *front,
        POINTER wave,
        POINT *oldp,
        POINT *newp,
        HYPER_SURF_ELEMENT *oldhse,
        HYPER_SURF         *oldhs,
        double              dt,
        double              *V)
{
	switch(wave_type(oldhs))
	{
	case SUBDOMAIN_BOUNDARY:
	case NEUMANN_BOUNDARY:
	    return;
	case DIRICHLET_BOUNDARY:
	    dirichlet_point_propagate(front,wave,oldp,newp,oldhse,
                                        oldhs,dt,V);
	    return;
	case FIRST_PHYSICS_WAVE_TYPE:
	    interior_point_propagate(front,wave,oldp,newp,oldhse,
                                        oldhs,dt,V);
	    return;
	}
}	/* end inverse_point_propagate */

static  void dirichlet_point_propagate(
        Front *front,
        POINTER wave,
        POINT *oldp,
        POINT *newp,
        HYPER_SURF_ELEMENT *oldhse,
        HYPER_SURF         *oldhs,
        double              dt,
        double              *V)
{
        int i, dim = front->rect_grid->dim;
        STATE *sl,*sr,*state,*bstate;
	CIM_PARAMS *cim_params = (CIM_PARAMS*)front->extra1;

        for (i = 0; i < dim; ++i)
            Coords(newp)[i] = Coords(oldp)[i];
	if ((bstate = (STATE*)boundary_state(oldhs)) == NULL)
	    return;

        FT_GetStatesAtPoint(oldp,oldhse,oldhs,(POINTER*)&sl,(POINTER*)&sr);

	state =  (STATE*)left_state(newp);
	state->u = bstate->u;
	state =  (STATE*)right_state(newp);
	state->u = bstate->u;

        return;
}       /* dirichlet_point_propagate */


static  void interior_point_propagate(
        Front *front,
        POINTER wave,
        POINT *oldp,
        POINT *newp,
        HYPER_SURF_ELEMENT *oldhse,
        HYPER_SURF         *oldhs,
        double              dt,
        double              *V)
{
	fourth_order_point_propagate(front,wave,oldp,newp,oldhse,oldhs,dt,V);
}       /* interior_point_propagate */

static void computeDirichletGradient(
	char *outname,
	Front *front)
{
	POINT *p;
	HYPER_SURF_ELEMENT *hse;
        HYPER_SURF         *hs;	
	INTERFACE *intfc = front->interf;
	STATE *sl,*sr;
	BDRY_INTEGRAL bdry_integral;
	char fname[256];
	static FILE *bdry_file;

	sprintf(fname,"%s/intGradU.xg",outname);
	if (bdry_file == NULL)
	    bdry_file = fopen(fname,"w");

	next_point(intfc,NULL,NULL,NULL);
	while (next_point(intfc,&p,&hse,&hs))
	{
	    if (wave_type(hs) != DIRICHLET_BOUNDARY) continue;
            FT_GetStatesAtPoint(p,hse,hs,(POINTER*)&sl,(POINTER*)&sr);
	    sl->gradU = sr->gradU = dirichletPointGradient(front,p,hse,hs);
	}
	integrateGradient(intfc,&bdry_integral);	
	fprintf(bdry_file,"%f %f\n",front->time,bdry_integral.gradU);

}	/* end computeDirichletGradient */

static double dirichletPointGradient(
	Front *front,
	POINT *p,
	HYPER_SURF_ELEMENT *hse,
	HYPER_SURF *hs)
{
	CIM_PARAMS *cim_params = (CIM_PARAMS*)front->extra1;
	FIELD *field = cim_params->field;
	double *u = field->u;
	double u0,u1;
	STATE *bstate;
	double u_in;
	int comp;
	double p1[MAXD],nor[MAXD];
	int i,dim = front->rect_grid->dim;
	double dn,*h = front->rect_grid->h;
	double gradU;

	if ((bstate = (STATE*)boundary_state(hs)) == NULL)
	    return 0.0;
	comp = (positive_component(hs) != exterior_component(front->interf)) ?
		positive_component(hs) : negative_component(hs);
	u0 = bstate->u;

	FT_NormalAtPoint(p,front,nor,comp);
	dn = grid_size_in_direction(nor,h,dim);
	for (i = 0; i < dim; ++i)
            p1[i] = Coords(p)[i] + nor[i]*dn;
        FT_IntrpStateVarAtCoords(front,comp,p1,u,getStateU,&u1,&u0);
	gradU = (u1 - u0)/dn;
	return gradU;
}	/* end dirichletPointGradient */

static void integrateGradient(
	INTERFACE *intfc,
	BDRY_INTEGRAL *bdry_integral)
{
	int dim = Dimension(intfc);
	switch (dim)
	{
	case 2:
	    integrateGradient2d(intfc,bdry_integral);
	    return;
	case 3:
	    integrateGradient3d(intfc,bdry_integral);
	    return;
	}
}	/* end integrateGradient */

static void integrateGradient2d(
	INTERFACE *intfc,
	BDRY_INTEGRAL *bdry_integral)
{
	CURVE **c,*curve;
	BOND *b;
	POINT *ps,*pe;
	HYPER_SURF_ELEMENT *hse;
	HYPER_SURF *hs;
	STATE *sl,*sr;
	double gu_s,gu_e;

	bdry_integral->gradU = 0.0;
	for (c = intfc->curves; c && *c; ++c)
	{
	    if (wave_type(*c) != DIRICHLET_BOUNDARY) continue;
	    curve = *c;
	    hs = Hyper_surf(curve);
	    for (b = curve->first; b != NULL; b = b->next)
	    {
		ps = b->start;
		pe = b->end;
	 	hse = Hyper_surf_element(b);
            	FT_GetStatesAtPoint(ps,hse,hs,(POINTER*)&sl,(POINTER*)&sr);
		gu_s = sl->gradU;
            	FT_GetStatesAtPoint(pe,hse,hs,(POINTER*)&sl,(POINTER*)&sr);
		gu_e = sl->gradU;
		bdry_integral->gradU += 0.5*bond_length(b)*(gu_s + gu_e);
	    }
	}
}	/* end integrateGradient2d */

static void integrateGradient3d(
	INTERFACE *intfc,
	BDRY_INTEGRAL *bdry_integral)
{
}	/* end integrateGradient3d */
