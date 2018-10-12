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

#include <cFluid.h>

static void getRTState(STATE*,EQN_PARAMS*,double*,COMPONENT);
static void getRMState(STATE*,EQN_PARAMS*,double*,COMPONENT);
static void getBubbleState(STATE*,EQN_PARAMS*,double*,COMPONENT);
static void getAmbientState(STATE*,EQN_PARAMS*,double*,COMPONENT);
static void getBlastState(STATE*,EQN_PARAMS*,double*,COMPONENT);
static void getShockSineWaveState(STATE*,EQN_PARAMS*,double*,COMPONENT);
static void getAccuracySineWaveState(STATE*,EQN_PARAMS*,double*,COMPONENT);
static void behind_state(int,double,double*,int,int,STATE*,STATE*);
static void assign_seeds(unsigned short int*);
static void Fracr(double,double,double,double,int,double,double*);
static void pertCoeff(double**,double**,double**,double**,int,double*,double*);
static void init_gun_and_bullet(Front*);
static double intfcPertHeight(FOURIER_POLY*,double*);
static double getStationaryVelocity(EQN_PARAMS*);

static void setRayleighTaylorParams(EQN_PARAMS*,char*);
static void setRichtmyerMeshkovParams(EQN_PARAMS*,char*);
static void setBubbleParams(EQN_PARAMS*,char*);
static void setImplosionParams(EQN_PARAMS*,char*);
static void setMTFusionParams(EQN_PARAMS*,char*);
static void setProjectileParams(EQN_PARAMS*,char*);
static void setRiemProbParams(EQN_PARAMS*,char*);
static void setRiemProbParams1d(EQN_PARAMS*,char*);
static void setRiemProbParams2d(EQN_PARAMS*,char*);
static void setOnedParams(EQN_PARAMS*,char*);

static void initSinePertIntfc(Front*,LEVEL_FUNC_PACK*);
static void initRandPertIntfc(Front*,LEVEL_FUNC_PACK*);
static void initCirclePlaneIntfc(Front*,LEVEL_FUNC_PACK*);
static void initImplosionIntfc(Front*,LEVEL_FUNC_PACK*);
static void initMTFusionIntfc(Front*,LEVEL_FUNC_PACK*);
static void initProjectileIntfc(Front*,LEVEL_FUNC_PACK*);
static void initProjectileIntfc2d(Front*,LEVEL_FUNC_PACK*);
static void initRectPlaneIntfc(Front*,LEVEL_FUNC_PACK*);
static void initTrianglePlaneIntfc(Front*,LEVEL_FUNC_PACK*);
static void initCylinderPlaneIntfc(Front*,LEVEL_FUNC_PACK*);
static void initRiemannProb(Front*,LEVEL_FUNC_PACK*);
static void initObliqueIntfc(Front*,LEVEL_FUNC_PACK*);

static void cFluid_setProbParams(char*,EQN_PARAMS*);



void read_cFluid_params(char *inname, EQN_PARAMS *eqn_params)
{
	char string[100];
	FILE *infile = fopen(inname,"r");
	int dim = eqn_params->dim;

	eqn_params->tracked = YES;  // Default
	eqn_params->prob_type = ERROR_TYPE;

	CursorAfterString(infile,"Enter problem type:");
	fscanf(infile,"%s",string);
	(void) printf("%s\n",string);
	if (string[0] == 'T' || string[0] == 't')
	{
	    if (string[10] == 'B' || string[10] == 'b')
	    	eqn_params->prob_type = TWO_FLUID_BUBBLE;
	    else if (string[10] == 'R' || string[10] == 'r')
	    {
		if (string[11] == 'T' || string[11] == 't')
	    	    eqn_params->prob_type = TWO_FLUID_RT;
		else if (string[11] == 'M' || string[11] == 'm')
	    	    eqn_params->prob_type = TWO_FLUID_RM;
	    }
	} 
	else if (string[0] == 'F' || string[0] == 'f')
	{
            if (string[12] == 'C' || string[12] =='c')
            {
                if (string[13] == 'I' || string[13] =='i')
                    eqn_params->prob_type = FLUID_SOLID_CIRCLE;
                else if (string[13] == 'Y' || string[13] =='y')
                    eqn_params->prob_type = FLUID_SOLID_CYLINDER;
            }
            else if (string[12] == 'R' || string[12] =='r')
                eqn_params->prob_type = FLUID_SOLID_RECT;
            else if (string[12] == 'T' || string[12] =='t')
                eqn_params->prob_type = FLUID_SOLID_TRIANGLE;
        }
	else if (string[0] == 'B' || string[0] == 'b')
	    eqn_params->prob_type = BUBBLE_SURFACE;
	else if (string[0] == 'I' || string[0] == 'i')
	    eqn_params->prob_type = IMPLOSION;
	else if (string[0] == 'P' || string[0] == 'p')
	    eqn_params->prob_type = PROJECTILE;
	else if (string[0] == 'R' || string[0] == 'r')
	    eqn_params->prob_type = RIEMANN_PROB;
	else if (string[0] == 'M' || string[0] == 'm')
	    eqn_params->prob_type = MT_FUSION;
	else if (string[0] == 'O' || string[0] == 'o')
	{
	    if (string[1] == 'B' || string[1] == 'b')
		eqn_params->prob_type = OBLIQUE_SHOCK_REFLECT;
	    else if (string[5] == 'S' || string[5] == 's')
	    	eqn_params->prob_type = ONED_SSINE;
	    else if (string[5] == 'B' || string[5] == 'b')
	    	eqn_params->prob_type = ONED_BLAST;
	    else if (string[5] == 'A' || string[5] == 'a')
	    	eqn_params->prob_type = ONED_ASINE;
	}
	CursorAfterString(infile,"Enter numerical scheme for interior solver:");
	fscanf(infile,"%s",string);
	(void) printf("%s\n",string);
	switch (string[0])
	{
	case 'T':
	case 't':
	    switch (string[4])
	    {
	    case '1':
		eqn_params->num_scheme = TVD_FIRST_ORDER;
		break;
	    case '2':
		eqn_params->num_scheme = TVD_SECOND_ORDER;
		break;
	    case '4':
		eqn_params->num_scheme = TVD_FOURTH_ORDER;
		break;
	    default:
		printf("Numerical scheme %s not implemented!\n",string);
		clean_up(ERROR);
	    }
	    break;
	case 'W':
	case 'w':
	    switch (string[5])
	    {
	    case '1':
		eqn_params->num_scheme = WENO_FIRST_ORDER;
		break;
	    case '2':
		eqn_params->num_scheme = WENO_SECOND_ORDER;
		break;
	    case '4':
		eqn_params->num_scheme = WENO_FOURTH_ORDER;
		break;
	    default:
		printf("Numerical scheme %s not implemented!\n",string);
		clean_up(ERROR);
	    }
	    eqn_params->articomp = NO;
	    if (CursorAfterStringOpt(infile,
		"Enter yes to use artificial compression:"))
	    {
	    	fscanf(infile,"%s",string);
            	(void) printf("%s\n",string);
	    	if (string[0] == 'y' || string[0] == 'Y')
	            eqn_params->articomp = YES;
	    }
	    break;
	default:
	    printf("Numerical scheme %s not implemented!\n",string);
	    clean_up(ERROR);
	}
	CursorAfterString(infile,"Enter order of point propagator:");
	fscanf(infile,"%s",string);
	(void) printf("%s\n",string);
	switch (string[0])
	{
	case '1':
	    eqn_params->point_prop_scheme = FIRST_ORDER;
	    break;
	case '2':
	    eqn_params->point_prop_scheme = SECOND_ORDER;
	    break;
	case '4':
	    eqn_params->point_prop_scheme = FOURTH_ORDER;
	    break;
	default:
	    printf("Point propagator order %s not implemented!\n",string);
	    clean_up(ERROR);
	}

	eqn_params->use_base_soln = NO;
	if (CursorAfterStringOpt(infile,
		"Enter yes for comparison with base data:"))
	{
            fscanf(infile,"%s",string);
            (void) printf("%s\n",string);
            if (string[0] == 'y' || string[0] == 'Y')
            	eqn_params->use_base_soln = YES;
	}

	if (eqn_params->use_base_soln == YES)
        {
	    CursorAfterString(infile,"Enter base directory name:");
            fscanf(infile,"%s",eqn_params->base_dir_name);
            (void) printf("%s\n",eqn_params->base_dir_name);
            CursorAfterString(infile,"Enter number of comparing steps:");
            fscanf(infile,"%d",&eqn_params->num_step);
            (void) printf("%d\n",eqn_params->num_step);
            FT_VectorMemoryAlloc((POINTER*)&eqn_params->steps,
                                eqn_params->num_step,sizeof(int));
            for (int i = 0; i < eqn_params->num_step; ++i)
            {
                sprintf(string,"Enter index of step %d:",i+1);
                CursorAfterString(infile,string);
                fscanf(infile,"%d",&eqn_params->steps[i]);
                (void) printf("%d\n",eqn_params->steps[i]);
            }
            FT_ScalarMemoryAlloc((POINTER*)&eqn_params->f_basic,
                                sizeof(F_BASIC_DATA));
            eqn_params->f_basic->dim = dim;
	}

	assert(eqn_params->prob_type != ERROR_TYPE);
	fclose(infile);

	if (eqn_params->use_base_soln == YES)
	    FT_ReadComparisonDomain(inname,eqn_params->f_basic);


    cFluid_setProbParams(inname, eqn_params);

}	/* end read_cFluid_params */


static void cFluid_setProbParams(char* inname,
        EQN_PARAMS* eqn_params)
{
	switch (eqn_params->prob_type)
	{
	case TWO_FLUID_RT:
	    setRayleighTaylorParams(eqn_params,inname);
	    break;
	case TWO_FLUID_RM:
	case TWO_FLUID_RM_RAND:
	    setRichtmyerMeshkovParams(eqn_params,inname);
	    break;
	case TWO_FLUID_BUBBLE:
	    setBubbleParams(eqn_params,inname);
	    break;
	case IMPLOSION:
	    setImplosionParams(eqn_params,inname);
	    break;
	case MT_FUSION:
	    setMTFusionParams(eqn_params,inname);
	    break;
	case PROJECTILE:
	case FLUID_SOLID_CIRCLE:
	case FLUID_SOLID_RECT:
	case FLUID_SOLID_TRIANGLE:
	case FLUID_SOLID_CYLINDER:
	    setProjectileParams(eqn_params,inname);
	    break;
	case RIEMANN_PROB:
	    setRiemProbParams(eqn_params,inname);
	    break;
	case ONED_BLAST:
	case ONED_SSINE:
	case ONED_ASINE:
	    setOnedParams(eqn_params,inname);
	    break;
	case OBLIQUE_SHOCK_REFLECT:
	    setRichtmyerMeshkovParams(eqn_params,inname);
	    break;
    default:
	    printf("In setProbParams(), unknown problem type!\n");
	    clean_up(ERROR);
	}
}

void cFluid_InitIntfc(Front* front,
        LEVEL_FUNC_PACK* level_func_pack)
{
    PROB_TYPE prob_type =
        (PROB_TYPE)((EQN_PARAMS*)front->extra1)->prob_type;

	switch( prob_type )
	{
	case TWO_FLUID_RT:
	case TWO_FLUID_RM:
	    initSinePertIntfc(front,level_func_pack);
	    break;
	case TWO_FLUID_RM_RAND:
	    initRandPertIntfc(front,level_func_pack);
	    break;
	case TWO_FLUID_BUBBLE:
	case FLUID_SOLID_CIRCLE:
	    initCirclePlaneIntfc(front,level_func_pack);
	    break;
	case IMPLOSION:
	    initImplosionIntfc(front,level_func_pack);
	    break;
	case MT_FUSION:
	    initMTFusionIntfc(front,level_func_pack);
	    break;
	case PROJECTILE:
	    initProjectileIntfc(front,level_func_pack);
	    break;
	case FLUID_SOLID_RECT:
	    initRectPlaneIntfc(front,level_func_pack);
	    break;
	case FLUID_SOLID_TRIANGLE:
	    initTrianglePlaneIntfc(front,level_func_pack);
	    break;
	case FLUID_SOLID_CYLINDER:
        initCylinderPlaneIntfc(front,level_func_pack);
        break;
	case RIEMANN_PROB:
	case ONED_BLAST:
	case ONED_SSINE:
	case ONED_ASINE:
	    initRiemannProb(front,level_func_pack);
	    break;
	case OBLIQUE_SHOCK_REFLECT:
	    initObliqueIntfc(front,level_func_pack);
	    break;
	default:
	    (void) printf("Problem type not implemented, code needed!\n");
	    clean_up(ERROR);
	}
}

void initSinePertIntfc(Front* front,
        LEVEL_FUNC_PACK* level_func_pack)
{
    char* inname = InName(front);
	FILE *infile = fopen(inname,"r");
	char mesg[100];

	level_func_pack->func = level_wave_func;
	level_func_pack->wave_type = FIRST_PHYSICS_WAVE_TYPE;
	level_func_pack->neg_component = GAS_COMP1;
	level_func_pack->pos_component = GAS_COMP2;

	static double L[3], U[3];
	ft_assign(L, front->rect_grid->L, 3*DOUBLE);
	ft_assign(U, front->rect_grid->U, 3*DOUBLE);

    static FOURIER_POLY *level_func_params;
	FT_ScalarMemoryAlloc((POINTER*)&level_func_params,sizeof(FOURIER_POLY));

    int dim = front->rect_grid->dim;    
    level_func_params->dim = dim;

	level_func_params->L = L;
	level_func_params->U = U;

	CursorAfterString(infile,"Enter mean position of fluid interface:");
	fscanf(infile,"%lf",&level_func_params->z0);
	(void) printf("%f\n",level_func_params->z0);

	int num_modes;
	CursorAfterString(infile,"Enter number of sine modes:");
	fscanf(infile,"%d",&num_modes);
	(void) printf("%d\n",num_modes);
	
    level_func_params->num_modes = num_modes;
	
    FT_MatrixMemoryAlloc((POINTER*)&level_func_params->nu,
            num_modes,	dim-1, sizeof(double));
	FT_VectorMemoryAlloc((POINTER*)&level_func_params->phase,
            num_modes, sizeof(double));
	FT_VectorMemoryAlloc((POINTER*)&level_func_params->A,
            num_modes, sizeof(double));
	
    for (int i = 0; i < num_modes; ++i)
	{
	    sprintf(mesg,"Enter frequency of mode %d:",i+1);
	    CursorAfterString(infile,mesg);
	    for (int j = 0; j < dim-1; ++j)
	    {
	    	fscanf(infile,"%lf",&level_func_params->nu[i][j]);
    		(void) printf("%f ",level_func_params->nu[i][j]);
	    }
	    (void) printf("\n");

	    sprintf(mesg,"Enter amplitude of mode %d:",i+1);
	    CursorAfterString(infile,mesg);
	    fscanf(infile,"%lf",&level_func_params->A[i]);
	    (void) printf("%f\n",level_func_params->A[i]);
	    
        sprintf(mesg,"Enter phase of mode %d:",i+1);
	    CursorAfterString(infile,mesg);
	    fscanf(infile,"%lf",&level_func_params->phase[i]);
	    (void) printf("%f\n",level_func_params->phase[i]);
	}

	level_func_pack->func_params = (POINTER)level_func_params;
	
    EQN_PARAMS* eqn_params = (EQN_PARAMS*)front->extra1;
	eqn_params->level_func_params = (POINTER)level_func_params;
	
    fclose(infile);

}	/* end initSinePertIntfc */




static double intfcPertHeight(
	FOURIER_POLY *wave_params,
	double *coords)
{
	double arg,z,k,phase;
	int i,j,num_modes,dim;
	double *L = wave_params->L;
	double *U = wave_params->U;

	dim = wave_params->dim;
	num_modes = wave_params->num_modes;
	z = wave_params->z0;

	for (i = 0; i < num_modes; ++i)
	{
	    arg = 0.0;
	    for (j = 0; j < dim-1; ++j)
	    {
		k = wave_params->nu[i][j]*2.0*PI/(U[j]-L[j]);
		arg += k*coords[j];
	    }
	    phase = wave_params->phase[i]*PI/180.0;
	    arg -= phase;
	    z += wave_params->A[i]*sin(arg);
	}
	return z;
}	/* end intfcPertHeight */


void setRayleighTaylorParams(
        EQN_PARAMS* eqn_params, char *inname)
{
	int i;
	FILE *infile = fopen(inname,"r");
	double pinf,einf,gamma;
	char s[100], str[256];

	sprintf(str, "Enter gamma, pinf, einf of the fluid with comp %d:", 
	             GAS_COMP1);
	CursorAfterString(infile,str);
	fscanf(infile,"%lf %lf %lf",&gamma,&pinf,&einf);
	(void) printf("%f %f %f\n",gamma,pinf,einf);

	(eqn_params->eos[GAS_COMP1]).gamma = gamma;
	(eqn_params->eos[GAS_COMP1]).pinf = pinf;
	(eqn_params->eos[GAS_COMP1]).einf = einf;
	
	sprintf(str, "Enter gamma, pinf, einf of the fluid with comp %d:", 
		     GAS_COMP2);
	CursorAfterString(infile,str);
	fscanf(infile,"%lf %lf %lf",&gamma,&pinf,&einf);
	(void) printf("%f %f %f\n",gamma,pinf,einf);
	
    (eqn_params->eos[GAS_COMP2]).gamma = gamma;
	(eqn_params->eos[GAS_COMP2]).pinf = pinf;
	(eqn_params->eos[GAS_COMP2]).einf = einf;

	CursorAfterString(infile,"Enter density of top fluid:");
	fscanf(infile,"%lf",&eqn_params->rho2);
	(void) printf("%f\n",eqn_params->rho2);
	
    CursorAfterString(infile,"Enter density of bottom fluid:");
	fscanf(infile,"%lf",&eqn_params->rho1);
	(void) printf("%f\n",eqn_params->rho1);
	
    int dim = eqn_params->dim;    
    CursorAfterString(infile,"Enter gravity:");
	for (i = 0; i < dim; ++i)
	{
	    fscanf(infile,"%lf",&eqn_params->gravity[i]);
	    (void) printf("%f ",eqn_params->gravity[i]);
	}
	(void) printf("\n");
	
    CursorAfterString(infile,"Enter pressure at interface:");
	fscanf(infile,"%lf",&eqn_params->p0);
	(void) printf("%f\n",eqn_params->p0);
	
    CursorAfterString(infile,"Type yes to track the interface:");
	fscanf(infile,"%s",s);
	(void) printf("%s\n",s);
	
    if (s[0] == 'y' || s[0] == 'Y')
	    eqn_params->tracked = YES;
	else
	    eqn_params->tracked = NO;
	
    fclose(infile);
}

void G_CARTESIAN::initRayleighTaylorStates()
{
	int i,j,k,l,index;
	//EQN_PARAMS *eqn_params = (EQN_PARAMS*)front->extra1;
	double coords[MAXD];
	COMPONENT comp;
	STATE *sl,*sr,state;
    POINT *p;
    HYPER_SURF *hs;
    HYPER_SURF_ELEMENT *hse;
	INTERFACE *intfc = front->interf;

	double *dens = field.dens;
	double *engy = field.engy;
	double *pres = field.pres;
	double **momn = field.momn;

    next_point(intfc,NULL,NULL,NULL);

    while (next_point(intfc,&p,&hse,&hs))
    {
	    FT_GetStatesAtPoint(p,hse,hs,(POINTER*)&sl,(POINTER*)&sr);
	    getRTState(sl,eqn_params,Coords(p),negative_component(hs));
	    getRTState(sr,eqn_params,Coords(p),positive_component(hs));
	}

	FT_MakeGridIntfc(front);
	setDomain();

	switch (dim)
	{
	case 2:
	    for (j = imin[1]; j <= imax[1]; ++j)
	    for (i = imin[0]; i <= imax[0]; ++i)
	    {
		index = d_index2d(i,j,top_gmax);
		comp = top_comp[index];
		getRectangleCenter(index,coords);
        getRTState(&state,eqn_params,coords,comp);
		dens[index] = state.dens;
		pres[index] = state.pres;
		engy[index] = state.engy;
		for (l = 0; l < dim; ++l)
		    momn[l][index] = state.momn[l];
	    }
	    break;
	case 3:
	    for (k = imin[2]; k <= imax[2]; ++k)
	    for (j = imin[1]; j <= imax[1]; ++j)
	    for (i = imin[0]; i <= imax[0]; ++i)
	    {
		index = d_index3d(i,j,k,top_gmax);
		comp = top_comp[index];
		getRectangleCenter(index,coords);
        getRTState(&state,eqn_params,coords,comp);
		dens[index] = state.dens;
		pres[index] = state.pres;
		engy[index] = state.engy;
		for (l = 0; l < dim; ++l)
		    momn[l][index] = state.momn[l];
	    }
	    break;
	}
	scatMeshStates();
}	/* end initRayleiTaylorStates */

void setRichtmyerMeshkovParams(
        EQN_PARAMS* eqn_params, char *inname)
{
	FILE *infile = fopen(inname,"r");
	char s[100], str[256];
	double	pinf, einf, gamma;
	int i, dim = FT_Dimension();

	sprintf(str,"Enter gamma, pinf, einf of the fluid with comp %d:", 
	             GAS_COMP1);
	CursorAfterString(infile,str);
	fscanf(infile,"%lf %lf %lf",&gamma,&pinf,&einf);
	(void) printf("%f %f %f\n",gamma,pinf,einf);
	(eqn_params->eos[GAS_COMP1]).gamma = gamma;
	(eqn_params->eos[GAS_COMP1]).pinf = pinf;
	(eqn_params->eos[GAS_COMP1]).einf = einf;
	
	sprintf(str,"Enter gamma, pinf, einf of the fluid with comp %d:", 
		     GAS_COMP2);
	CursorAfterString(infile,str);
	fscanf(infile,"%lf %lf %lf",&gamma,&pinf,&einf);
	(void) printf("%f %f %f\n",gamma,pinf,einf);
	(eqn_params->eos[GAS_COMP2]).gamma = gamma;
	(eqn_params->eos[GAS_COMP2]).pinf = pinf;
	(eqn_params->eos[GAS_COMP2]).einf = einf;

	CursorAfterString(infile,"Enter density of top fluid:");
	fscanf(infile,"%lf",&eqn_params->rho2);
	(void) printf("%f\n",eqn_params->rho2);
	CursorAfterString(infile,"Enter density of bottom fluid:");
	fscanf(infile,"%lf",&eqn_params->rho1);
	(void) printf("%f\n",eqn_params->rho1);
	for (int i = 0; i < dim; ++i)
	    eqn_params->v2[i] = 0.0; //default initial velocity
	if (CursorAfterStringOpt(infile,"Enter velocity of top fluid:"))
	{
	    for (i = 0; i < dim; ++i)
	    {
		fscanf(infile,"%lf",&eqn_params->v2[i]);
		(void) printf("%f ",eqn_params->v2[i]);
	    }
	    (void) printf("\n");
	}
	for (int i = 0; i < dim; ++i) //default initial velocity
	    eqn_params->v1[i] = 0.0;
	if(CursorAfterStringOpt(infile,"Enter velocity of bottom fluid:"))
	{
	    for (i = 0; i < dim; ++i)
	    {
		fscanf(infile,"%lf",&eqn_params->v1[i]);
		(void) printf("%f ",eqn_params->v1[i]);
	    }
	    (void) printf("\n");
	}
	CursorAfterString(infile,"Enter pressure at interface:");
	fscanf(infile,"%lf",&eqn_params->p0);
	(void) printf("%f\n",eqn_params->p0);
	CursorAfterString(infile,"Enter Mach number of shock:");
	fscanf(infile,"%lf",&eqn_params->Mach_number);
	(void) printf("%f\n",eqn_params->Mach_number);
	eqn_params->idir = dim - 1;	//default idir
	if (CursorAfterStringOpt(infile,"Enter idir of shock:"))
	{
	    fscanf(infile,"%d",&eqn_params->idir);
	    (void) printf("%d\n",eqn_params->idir);
	}
	CursorAfterString(infile,"Enter position of shock:");
	fscanf(infile,"%lf",&eqn_params->shock_position);
	(void) printf("%f\n",eqn_params->shock_position);
	CursorAfterString(infile,"Enter direction of shock:");
	fscanf(infile,"%d",&eqn_params->shock_side);
	(void) printf("%d\n",eqn_params->shock_side);
	CursorAfterString(infile,"Type yes to track the interface:");
	fscanf(infile,"%s",s);
	(void) printf("%s\n",s);
	if (s[0] == 'y' || s[0] == 'Y')
	    eqn_params->tracked = YES;
	else
	    eqn_params->tracked = NO;
	if (CursorAfterStringOpt(infile,"Type yes for stationary contact: "))
        {
	    fscanf(infile,"%s",s);
	    (void) printf("%s\n",s);
            if (s[0] == 'y' || s[0] == 'Y')
	    {
                eqn_params->contact_stationary = YES;
		eqn_params->contact_vel = getStationaryVelocity(eqn_params);
	    }
            else
                eqn_params->contact_stationary = NO;
        }

    fclose(infile);
}

void G_CARTESIAN::initRichtmyerMeshkovStates()
{
	int i,j,k,l,index;
	//EQN_PARAMS *eqn_params = (EQN_PARAMS*)front->extra1;
	double coords[MAXD];
	COMPONENT comp;
	STATE *sl,*sr,state;
        POINT *p;
        HYPER_SURF *hs;
        HYPER_SURF_ELEMENT *hse;
	INTERFACE *intfc = front->interf;
	double *dens = field.dens;
	double *engy = field.engy;
	double *pres = field.pres;
	double **momn = field.momn;

        next_point(intfc,NULL,NULL,NULL);
        while (next_point(intfc,&p,&hse,&hs))
        {
	    FT_GetStatesAtPoint(p,hse,hs,(POINTER*)&sl,(POINTER*)&sr);
	    getRMState(sl,eqn_params,Coords(p),negative_component(hs));
	    getRMState(sr,eqn_params,Coords(p),positive_component(hs));
	}

	FT_MakeGridIntfc(front);
	setDomain();

	switch (dim)
	{
	case 2:
	    for (j = imin[1]; j <= imax[1]; ++j)
	    for (i = imin[0]; i <= imax[0]; ++i)
	    {
		index = d_index2d(i,j,top_gmax);
		comp = top_comp[index];
		if (gas_comp(comp))
		{
		    getRectangleCenter(index,coords);
	    	    getRMState(&state,eqn_params,coords,comp);
		    dens[index] = state.dens;
		    pres[index] = state.pres;
		    engy[index] = state.engy;
		    for (l = 0; l < dim; ++l)
		    	momn[l][index] = state.momn[l];
		}
		else
		{
		    dens[index] = 0.0;
		    pres[index] = 0.0;
		    engy[index] = 0.0;
		    for (l = 0; l < dim; ++l)
		    	momn[l][index] = 0.0;
		}
	    }
	    break;
	case 3:
	    for (k = imin[2]; k <= imax[2]; ++k)
	    for (j = imin[1]; j <= imax[1]; ++j)
	    for (i = imin[0]; i <= imax[0]; ++i)
	    {
		index = d_index3d(i,j,k,top_gmax);
		comp = top_comp[index];
		getRectangleCenter(index,coords);
	    	getRMState(&state,eqn_params,coords,comp);
		dens[index] = state.dens;
		pres[index] = state.pres;
		engy[index] = state.engy;
		for (l = 0; l < dim; ++l)
		    momn[l][index] = state.momn[l];
	    }
	    break;
	}
	scatMeshStates();
}	/* end initRichtmyerMeshkovStates */

//EOS dep
static void getRTState(
	STATE *state,
	EQN_PARAMS *eqn_params,
	double *coords,
	COMPONENT comp)
{
	FOURIER_POLY	*wave_params;
	EOS_PARAMS	*eos;
	double		z_intfc;
	double 		rho1 = eqn_params->rho1;
	double 		rho2 = eqn_params->rho2;
	double 		p0 = eqn_params->p0;
	double 		*g = eqn_params->gravity;
	double 		dz, gz, c2, gamma;
	int    		i,dim;
	double		tmp;

	eos = &(eqn_params->eos[comp]);
	state->eos = eos;
	gamma = eos->gamma;

	wave_params = (FOURIER_POLY*)eqn_params->level_func_params;
	dim = wave_params->dim;
	z_intfc = intfcPertHeight(wave_params,coords);
	dz = coords[dim-1] - z_intfc;
	gz = g[dim-1];

	/* Constant density */
	for (i = 0; i < dim; ++i)
	    state->momn[i] = 0.0;

	switch (comp)
	{
	case GAS_COMP1:
	    c2 = gamma*(p0+eos->pinf)/rho1;
	    tmp = exp(gamma*gz/c2*dz);
	    state->dens = rho1*tmp;
	    state->pres = state->dens*c2/gamma - eos->pinf;
	    state->engy = EosInternalEnergy(state);
	    break;
	case GAS_COMP2:
	    c2 = gamma*(p0+eos->pinf)/rho2;
	    tmp = exp(gamma*gz/c2*dz);
	    state->dens = rho2*tmp;
	    state->pres = state->dens*c2/gamma - eos->pinf;
	    state->engy = EosInternalEnergy(state);
	    break;
	case EXT_COMP:
	    state->dens = 0.0;
	    state->pres = 0.0;
	    state->engy = 0.0;
	    break;
	default:
	    printf("ERROR: Unknown component %d in getRTState()!\n",comp);
	    clean_up(ERROR);
	}
}	/* end getRTState */

static void getRMState(
	STATE *state,
	EQN_PARAMS *eqn_params,
	double *coords,
	COMPONENT comp)
{
	FOURIER_POLY *wave_params;
	EOS_PARAMS	*eos;
	double rho1 = eqn_params->rho1;
	double rho2 = eqn_params->rho2;
	double p0 = eqn_params->p0;
	double shock_position = eqn_params->shock_position;
	double Mach_number = eqn_params->Mach_number;
	double shock_speed;
	double csp = eqn_params->contact_vel;
	int idir = eqn_params->idir;
	int shock_side = eqn_params->shock_side;
	int i,dim;
	double v1[MAXD], v2[MAXD];
 
	if (debugging("rm_state"))
	    printf("Entering getRMState(), coords = %f %f\n",
				coords[0],coords[1]);
	wave_params = (FOURIER_POLY*)eqn_params->level_func_params;
	dim = wave_params->dim;
	for (i = 0; i < dim; ++i)
	{
	    v1[i] = eqn_params->v1[i];
	    v2[i] = eqn_params->v2[i];
	}

	/* Constant density */
	for (i = 0; i < dim; ++i)
	    state->vel[i] = state->momn[i] = 0.0;
	state->dim = dim;
	eos = &(eqn_params->eos[comp]);
	state->eos = eos;
	
	switch (comp)
	{
	case GAS_COMP1:
	    state->dens = rho1;
	    state->pres = p0;
	    for (i = 0; i < dim; ++i)
	    {
		state->vel[i] = v1[i];
		state->momn[i] = rho1 * v1[i];
	    }
	    state->engy = EosInternalEnergy(state);
	    break;
	case GAS_COMP2:
	    state->dens = rho2;
	    state->pres = p0;
	    for (i = 0; i < dim; ++i)
	    {
		state->vel[i] = v2[i];
		state->momn[i] = rho2 * v2[i];
	    }
	    state->engy = EosInternalEnergy(state);
	    break;
	case EXT_COMP:
	    state->dens = 0.0;
	    state->pres = 0.0;
	    state->engy = 0.0;
	    break;
	default:
	    printf("ERROR: Unknown component %d in getRTState()!\n",comp);
	    clean_up(ERROR);
	}
	if (debugging("rm_state"))
	{
	    printf("Before calling behind_state()\n");
	    printf("state = %f %f %f\n",state->dens,state->pres,
					state->vel[0]);
	}
	if ((shock_side ==  1 && coords[idir] < shock_position) ||
	    (shock_side == -1 && coords[idir] > shock_position))
	{
	    behind_state(SHOCK_MACH_NUMBER,Mach_number,
			&shock_speed,idir,shock_side,state,state);
	    state->engy = EosEnergy(state);
	    if (debugging("rm_state"))
	    {
	    	printf("After calling behind_state()\n");
	    	printf("state = %f %f %f\n",state->dens,state->pres,
			state->vel[0]);
	    }
	}
	state->vel[idir] -= csp;
	state->momn[idir] = state->vel[idir]*state->dens;
	state->engy = EosEnergy(state);

}	/* end getRMState */

static void behind_state(
	int		which_parameter,
	double		parameter,
	double		*shock_speed,
	int		idir,
	int		shock_side,
	STATE		*ahead_state,
	STATE		*behind_state)
{
	double		r0, p0, u0;		/* ahead state */
	double		r1, p1, u1;		/* behind state */
	double		U;			/* shock speed */
	double		M0n;			/* shock mack number,
						   relative to ahead flow */
	double		M0nsq;			/* steady normal ahead Mach
						   number squared */
	int		dim;

	dim = ahead_state->dim;
	r0  = ahead_state->dens;
	p0  = ahead_state->pres;
	u0  = ahead_state->vel[idir]*shock_side;

	switch(which_parameter)
	{
	case SHOCK_MACH_NUMBER:
	    M0n = parameter;
	    *shock_speed = U = u0 + M0n*EosSoundSpeed(ahead_state);
	    M0nsq = sqr(M0n);
	    p1 = EosMaxBehindShockPres(M0nsq,ahead_state);
	    u1 =  u0 + (p0 - p1) / (r0*(u0 - U)); 
	    r1 = r0*((u0 - U)/(u1 - U));
	    if (debugging("rm_state"))
	    {
		printf("M0n = %f  shock_speed = %f\n",M0n,*shock_speed);
		printf("p1 = %f  u1 = %f  r1 = %f\n",p1,u1,r1);
	    }
	    break;
	default:
	    screen("ERROR in behind_state(), "
	           "unknown parameter %d\n",which_parameter);
	    clean_up(ERROR);
	}	
	behind_state->dens = r1;
	behind_state->pres = p1;
	behind_state->vel[idir] = u1*shock_side;
	behind_state->momn[idir] = r1*u1*shock_side;
}		/*end behind_state */

void initCirclePlaneIntfc(Front* front,
        LEVEL_FUNC_PACK *level_func_pack)
{
	FILE *infile = fopen(InName(front),"r");
	static CIRCLE_PARAMS *circle_params;
	int i,dim;

	EQN_PARAMS *eqn_params = (EQN_PARAMS*)front->extra1;
	PROB_TYPE prob_type = eqn_params->prob_type;

	FT_ScalarMemoryAlloc((POINTER*)&circle_params,sizeof(CIRCLE_PARAMS));
        circle_params->dim = dim = front->rect_grid->dim;
        circle_params->add_plan_surf = NO;
        CursorAfterString(infile,"Enter the center of the circle:");
        for (i = 0; i < dim; ++i)
	{
            fscanf(infile,"%lf",&circle_params->cen[i]);
            (void) printf("%f ",circle_params->cen[i]);
	}
	(void) printf("\n");
        CursorAfterString(infile,"Enter radius of the circle:");
        fscanf(infile,"%lf",&circle_params->R);
        (void) printf("%f\n",circle_params->R);

        if (prob_type == BUBBLE_SURFACE)
        {
            CursorAfterString(infile,"Enter height of the surface:");
            fscanf(infile,"%lf",&circle_params->H);
            (void) printf("%f\n",circle_params->H);
            circle_params->add_plan_surf = YES;
        }
	level_func_pack->func_params = (POINTER)circle_params;
	eqn_params->level_func_params = (POINTER)circle_params;

	switch (prob_type)
	{
	case TWO_FLUID_BUBBLE:
        case BUBBLE_SURFACE:
            level_func_pack->neg_component = GAS_COMP1;
            level_func_pack->pos_component = GAS_COMP2;
            level_func_pack->func = level_circle_func;
            level_func_pack->wave_type = FIRST_PHYSICS_WAVE_TYPE;
            break;
        case FLUID_SOLID_CIRCLE:
            level_func_pack->neg_component = SOLID_COMP;
            level_func_pack->pos_component = GAS_COMP1;
            level_func_pack->func = level_circle_func;
	    level_func_pack->wave_type = MOVABLE_BODY_BOUNDARY;
            break;
	default:
	    (void) printf("ERROR Wrong type in initCirclePlaneIntfc()\n");
	    clean_up(ERROR);
	}
	fclose(infile);	
}

void G_CARTESIAN::initBubbleStates()
{
	int i,j,k,l,index;
	//EQN_PARAMS *eqn_params = (EQN_PARAMS*)front->extra1;
	double coords[MAXD];
	COMPONENT comp;
	STATE *sl,*sr,state;
        POINT *p;
        HYPER_SURF *hs;
        HYPER_SURF_ELEMENT *hse;
	INTERFACE *intfc = front->interf;
	double *dens = field.dens;
	double *engy = field.engy;
	double *pres = field.pres;
	double **momn = field.momn;

        next_point(intfc,NULL,NULL,NULL);
        while (next_point(intfc,&p,&hse,&hs))
        {
	    FT_GetStatesAtPoint(p,hse,hs,(POINTER*)&sl,(POINTER*)&sr);
	    getBubbleState(sl,eqn_params,Coords(p),negative_component(hs));
	    getBubbleState(sr,eqn_params,Coords(p),positive_component(hs));
	}

	FT_MakeGridIntfc(front);
	setDomain();

	switch (dim)
	{
	case 2:
	    for (j = imin[1]; j <= imax[1]; ++j)
	    for (i = imin[0]; i <= imax[0]; ++i)
	    {
		index = d_index2d(i,j,top_gmax);
		comp = top_comp[index];
		getRectangleCenter(index,coords);
	    	getBubbleState(&state,eqn_params,coords,comp);
		dens[index] = state.dens;
		pres[index] = state.pres;
		engy[index] = state.engy;
		for (l = 0; l < dim; ++l)
		    momn[l][index] = state.momn[l];
	    }
	    break;
	case 3:
	    for (k = imin[2]; k <= imax[2]; ++k)
	    for (j = imin[1]; j <= imax[1]; ++j)
	    for (i = imin[0]; i <= imax[0]; ++i)
	    {
		index = d_index3d(i,j,k,top_gmax);
		comp = top_comp[index];
		getRectangleCenter(index,coords);
	    	getBubbleState(&state,eqn_params,coords,comp);
		dens[index] = state.dens;
		pres[index] = state.pres;
		engy[index] = state.engy;
		for (l = 0; l < dim; ++l)
		    momn[l][index] = state.momn[l];
	    }
	    break;
	}
	scatMeshStates();
}	/* end initRayleiTaylorStates */

static void getBubbleState(
	STATE *state,
	EQN_PARAMS *eqn_params,
	double *coords,
	COMPONENT comp)
{
	static CIRCLE_PARAMS	*circle_params;
	EOS_PARAMS	*eos;
	double		z0;
	double 		p1 = eqn_params->p1;
	double 		p2 = eqn_params->p2;
	double 		rho1 = eqn_params->rho1;
	double 		rho2 = eqn_params->rho2;
	double 		*g = eqn_params->gravity;
	double 		dz, gz, c2, gamma;
	int    		i,dim;
	double		tmp;

	eos = &(eqn_params->eos[comp]);
	state->eos = eos;
	gamma = eos->gamma;

	circle_params = (CIRCLE_PARAMS*)eqn_params->level_func_params;
	dim = circle_params->dim;
	z0 = circle_params->cen[dim-1];
	dz = coords[dim-1] - z0;
	gz = g[dim-1];

	/* Constant density */
	for (i = 0; i < dim; ++i)
	    state->momn[i] = 0.0;
	switch (comp)
	{
	case GAS_COMP1:
	    c2 = gamma*(p1+eos->pinf)/rho1;
	    tmp = exp(gamma*gz/c2*dz);
	    state->dens = rho1*tmp;
	    state->pres = state->dens*c2/gamma - eos->pinf;
	    state->engy = EosInternalEnergy(state);
	    break;
	case GAS_COMP2:
	    c2 = gamma*(p2+eos->pinf)/rho2;
	    tmp = exp(gamma*gz/c2*dz);
	    state->dens = rho2*tmp;
	    state->pres = state->dens*c2/gamma - eos->pinf;
	    state->engy = EosInternalEnergy(state);
	    break;
	case EXT_COMP:
	    state->dens = 0.0;
	    state->pres = 0.0;
	    state->engy = 0.0;
	    break;
	default:
	    printf("ERROR: Unknown component %d in getRTState()!\n",comp);
	    clean_up(ERROR);
	}
}	/* end getBubbleState */

void setBubbleParams(
        EQN_PARAMS* eqn_params,
        char *inname)
{
	int i;
	FILE *infile = fopen(inname,"r");
	double		pinf,einf,gamma;
	char s[100], str[256];

	sprintf(str, "Enter gamma, pinf, einf of the fluid with comp %d:", 
	             GAS_COMP1);
	CursorAfterString(infile,str);
	fscanf(infile,"%lf %lf %lf",&gamma,&pinf,&einf);
	(void) printf("%f %f %f\n",gamma,pinf,einf);
	(eqn_params->eos[GAS_COMP1]).gamma = gamma;
	(eqn_params->eos[GAS_COMP1]).pinf = pinf;
	(eqn_params->eos[GAS_COMP1]).einf = einf;
	
	sprintf(str, "Enter gamma, pinf, einf of the fluid with comp %d:", 
		     GAS_COMP2);
	CursorAfterString(infile,str);
	fscanf(infile,"%lf %lf %lf",&gamma,&pinf,&einf);
	(void) printf("%f %f %f\n",gamma,pinf,einf);
	(eqn_params->eos[GAS_COMP2]).gamma = gamma;
	(eqn_params->eos[GAS_COMP2]).pinf = pinf;
	(eqn_params->eos[GAS_COMP2]).einf = einf;

	CursorAfterString(infile,"Enter density and pressure inside bubble:");
	fscanf(infile,"%lf %lf",&eqn_params->rho2,&eqn_params->p2);
	(void) printf("%f %f\n",eqn_params->rho2,eqn_params->p2);
	CursorAfterString(infile,"Enter density and pressure outside bubble:");
	fscanf(infile,"%lf %lf",&eqn_params->rho1,&eqn_params->p1);
	(void) printf("%f %f\n",eqn_params->rho1,eqn_params->p1);
	CursorAfterString(infile,"Enter gravity:");
    int dim = eqn_params->dim;
	for (i = 0; i < dim; ++i)
	{
	    fscanf(infile,"%lf",&eqn_params->gravity[i]);
	    (void) printf("%f ",eqn_params->gravity[i]);
	}
	(void) printf("\n");
	CursorAfterString(infile,"Type yes to track the interface:");
	fscanf(infile,"%s",s);
	(void) printf("%s\n",s);
	if (s[0] == 'y' || s[0] == 'Y')
	    eqn_params->tracked = YES;
	else
	    eqn_params->tracked = NO;
	fclose(infile);

}

void setImplosionParams(
        EQN_PARAMS* eqn_params,
        char* inname)
{
	int i;
	FILE *infile = fopen(inname,"r");
	double		pinf,einf,gamma;
	char s[100], str[256];

	sprintf(str, "Enter gamma, pinf, einf of the fluid with comp %d:", 
	             GAS_COMP1);
	CursorAfterString(infile,str);
	fscanf(infile,"%lf %lf %lf",&gamma,&pinf,&einf);
	(void) printf("%f %f %f\n",gamma,pinf,einf);
	(eqn_params->eos[GAS_COMP1]).gamma = gamma;
	(eqn_params->eos[GAS_COMP1]).pinf = pinf;
	(eqn_params->eos[GAS_COMP1]).einf = einf;
	
	sprintf(str, "Enter gamma, pinf, einf of the fluid with comp %d:", 
		     GAS_COMP2);
	CursorAfterString(infile,str);
	fscanf(infile,"%lf %lf %lf",&gamma,&pinf,&einf);
	(void) printf("%f %f %f\n",gamma,pinf,einf);
	(eqn_params->eos[GAS_COMP2]).gamma = gamma;
	(eqn_params->eos[GAS_COMP2]).pinf = pinf;
	(eqn_params->eos[GAS_COMP2]).einf = einf;

	CursorAfterString(infile,"Enter density and pressure inside contact:");
	fscanf(infile,"%lf %lf",&eqn_params->rho2,&eqn_params->p2);
	(void) printf("%f %f\n",eqn_params->rho2,eqn_params->p2);
	CursorAfterString(infile,"Enter density and pressure outside contact:");
	fscanf(infile,"%lf %lf",&eqn_params->rho1,&eqn_params->p1);
	(void) printf("%f %f\n",eqn_params->rho1,eqn_params->p1);
	CursorAfterString(infile,"Enter density and pressure at the boundary:");
	fscanf(infile,"%lf %lf",&eqn_params->rho0,&eqn_params->p0);
	(void) printf("%f %f\n",eqn_params->rho0,eqn_params->p0);
	CursorAfterString(infile,"Enter gravity:");
    int dim = eqn_params->dim;
	for (i = 0; i < dim; ++i)
	{
	    fscanf(infile,"%lf",&eqn_params->gravity[i]);
	    (void) printf("%f ",eqn_params->gravity[i]);
	}
	(void) printf("\n");
	CursorAfterString(infile,"Type yes to track the interface:");
	fscanf(infile,"%s",s);
	(void) printf("%s\n",s);
	if (s[0] == 'y' || s[0] == 'Y')
	    eqn_params->tracked = YES;
	else
	    eqn_params->tracked = NO;
	fclose(infile);
}

void setMTFusionParams(
        EQN_PARAMS* eqn_params, char* inname)
{
	int i;
	FILE *infile = fopen(inname,"r");
	double		pinf,einf,gamma;
	char s[100], str[256];

	sprintf(str, "Enter gamma, pinf, einf of the fluid with comp %d:", 
	             GAS_COMP1);
	CursorAfterString(infile,str);
	fscanf(infile,"%lf %lf %lf",&gamma,&pinf,&einf);
	(void) printf("%f %f %f\n",gamma,pinf,einf);
	(eqn_params->eos[GAS_COMP1]).gamma = gamma;
	(eqn_params->eos[GAS_COMP1]).pinf = pinf;
	(eqn_params->eos[GAS_COMP1]).einf = einf;
	
	sprintf(str, "Enter gamma, pinf, einf of the fluid with comp %d:", 
		     GAS_COMP2);
	CursorAfterString(infile,str);
	fscanf(infile,"%lf %lf %lf",&gamma,&pinf,&einf);
	(void) printf("%f %f %f\n",gamma,pinf,einf);
	(eqn_params->eos[GAS_COMP2]).gamma = gamma;
	(eqn_params->eos[GAS_COMP2]).pinf = pinf;
	(eqn_params->eos[GAS_COMP2]).einf = einf;

	CursorAfterString(infile,"Enter density and pressure inside contact:");
	fscanf(infile,"%lf %lf",&eqn_params->rho2,&eqn_params->p2);
	(void) printf("%f %f\n",eqn_params->rho2,eqn_params->p2);
	CursorAfterString(infile,"Enter density and pressure outside contact:");
	fscanf(infile,"%lf %lf",&eqn_params->rho1,&eqn_params->p1);
	(void) printf("%f %f\n",eqn_params->rho1,eqn_params->p1);
	CursorAfterString(infile,"Enter gravity:");
    int dim = eqn_params->dim;
	for (i = 0; i < dim; ++i)
	{
	    fscanf(infile,"%lf",&eqn_params->gravity[i]);
	    (void) printf("%f ",eqn_params->gravity[i]);
	}
	(void) printf("\n");
	CursorAfterString(infile,"Type yes to track the interface:");
	fscanf(infile,"%s",s);
	(void) printf("%s\n",s);
	if (s[0] == 'y' || s[0] == 'Y')
	    eqn_params->tracked = YES;
	else
	    eqn_params->tracked = NO;
	fclose(infile);
}

void initImplosionIntfc(Front* front,
        LEVEL_FUNC_PACK *level_func_pack)
{
	FILE *infile = fopen(InName(front),"r");
	static CIRCLE_PARAMS *circle_params;
	COMPONENT neg_comp,pos_comp;
	int i,dim,num_segs;
	double  (*func)(POINTER,double*);
        POINTER func_params;
        CURVE **wall,**contact;
	char s[100];

	FT_ScalarMemoryAlloc((POINTER*)&circle_params,sizeof(CIRCLE_PARAMS));
        circle_params->dim = dim = front->rect_grid->dim;
        circle_params->add_plan_surf = NO;
        CursorAfterString(infile,"Enter the center of implosion:");
        for (i = 0; i < dim; ++i)
	{
            fscanf(infile,"%lf",&circle_params->cen[i]);
            (void) printf("%f ",circle_params->cen[i]);
	}
	(void) printf("\n");
        CursorAfterString(infile,"Enter radius of the wall:");
        fscanf(infile,"%lf",&circle_params->R);
        (void) printf("%f\n",circle_params->R);

	func_params = (POINTER)circle_params;
	func = level_circle_func;
	neg_comp = GAS_COMP1;
	pos_comp = SOLID_COMP;
	wall = (CURVE**)FT_CreateLevelHyperSurfs(front->rect_grid,front->interf,
                    	neg_comp,pos_comp,func,func_params,DIRICHLET_BOUNDARY,
			&num_segs);

	CursorAfterString(infile,
                        "Enter radius of the two fluid interface:");
        fscanf(infile,"%lf",&circle_params->R);
        (void) printf("%f\n",circle_params->R);
	CursorAfterString(infile,"Enter yes to add perturbation:");
        fscanf(infile,"%s",s);
	(void) printf("%s\n",s);
	if (s[0] == 'y' || s[0] == 'Y')
	{
	    FOURIER_POLY *fpoly;
            circle_params->add_perturbation = YES;
	    FT_ScalarMemoryAlloc((POINTER*)&circle_params->fpoly,
				sizeof(FOURIER_POLY));
	    fpoly = circle_params->fpoly;
	    CursorAfterString(infile,"Enter number of Fourier modes:");
            fscanf(infile,"%d",&fpoly->num_modes);
            (void) printf("%d\n",fpoly->num_modes);
	    FT_VectorMemoryAlloc((POINTER*)&fpoly->A,fpoly->num_modes,
				sizeof(double));
	    FT_VectorMemoryAlloc((POINTER*)&fpoly->phase,fpoly->num_modes,
				sizeof(double));
	    FT_MatrixMemoryAlloc((POINTER*)&fpoly->nu,1,fpoly->num_modes,
				sizeof(double));
	    for (i = 0; i < fpoly->num_modes; ++i)
	    {
		sprintf(s,"For %d-th mode",i);
	    	CursorAfterString(infile,s);
		sprintf(s,"enter frequency, amplitude and phase:");
	    	CursorAfterString(infile,s);
		fscanf(infile,"%lf %lf %lf",&fpoly->nu[0][i],&fpoly->A[i],
					&fpoly->phase[i]);
		(void) printf("%f %f %f\n",fpoly->nu[0][i],fpoly->A[i],
					fpoly->phase[i]);
	    }
	}
        neg_comp = GAS_COMP2;
        pos_comp = GAS_COMP1;
	contact = (CURVE**)FT_CreateLevelHyperSurfs(front->rect_grid,
			front->interf,neg_comp,pos_comp,func,func_params,
			FIRST_PHYSICS_WAVE_TYPE,&num_segs);

	level_func_pack->func = NULL;
	level_func_pack->point_array = NULL;

	fclose(infile);	
}

void initMTFusionIntfc(Front* front,
        LEVEL_FUNC_PACK *level_func_pack)
{
	FILE *infile = fopen(InName(front),"r");
	static CIRCLE_PARAMS *circle_params;
	COMPONENT neg_comp,pos_comp;
	int i,dim,num_segs;
	double  (*func)(POINTER,double*);
        POINTER func_params;
        CURVE **wall,**contact;
	char s[100];

	FT_ScalarMemoryAlloc((POINTER*)&circle_params,sizeof(CIRCLE_PARAMS));
        circle_params->dim = dim = front->rect_grid->dim;
        circle_params->add_plan_surf = NO;
	circle_params->add_perturbation = NO;
        CursorAfterString(infile,"Enter the center of the sphere:");
        for (i = 0; i < dim; ++i)
	{
            fscanf(infile,"%lf",&circle_params->cen[i]);
            (void) printf("%f ",circle_params->cen[i]);
	}
	(void) printf("\n");
        CursorAfterString(infile,"Enter radius of the sphere:");
        fscanf(infile,"%lf",&circle_params->R);
        (void) printf("%f\n",circle_params->R);

	func_params = (POINTER)circle_params;
	func = level_circle_func;
	neg_comp = GAS_COMP1;
	pos_comp = SOLID_COMP;
	wall = (CURVE**)FT_CreateLevelHyperSurfs(front->rect_grid,front->interf,
                    	neg_comp,pos_comp,func,func_params,DIRICHLET_BOUNDARY,
			&num_segs);

	CursorAfterString(infile,
                        "Enter radius of the two fluid interface:");
        fscanf(infile,"%lf",&circle_params->R);
        (void) printf("%f\n",circle_params->R);
	CursorAfterString(infile,"Enter yes to add perturbation:");
        fscanf(infile,"%s",s);
	(void) printf("%s\n",s);
	if (s[0] == 'y' || s[0] == 'Y')
	{
	    FOURIER_POLY *fpoly;
            circle_params->add_perturbation = YES;
	    FT_ScalarMemoryAlloc((POINTER*)&circle_params->fpoly,
				sizeof(FOURIER_POLY));
	    fpoly = circle_params->fpoly;
	    CursorAfterString(infile,"Enter number of Fourier modes:");
            fscanf(infile,"%d",&fpoly->num_modes);
            (void) printf("%d\n",fpoly->num_modes);
	    FT_VectorMemoryAlloc((POINTER*)&fpoly->A,fpoly->num_modes,
				sizeof(double));
	    FT_VectorMemoryAlloc((POINTER*)&fpoly->phase,fpoly->num_modes,
				sizeof(double));
	    FT_MatrixMemoryAlloc((POINTER*)&fpoly->nu,1,fpoly->num_modes,
				sizeof(double));
	    for (i = 0; i < fpoly->num_modes; ++i)
	    {
		sprintf(s,"For %d-th mode",i);
	    	CursorAfterString(infile,s);
		sprintf(s,"enter frequency, amplitude and phase:");
	    	CursorAfterString(infile,s);
		fscanf(infile,"%lf %lf %lf",&fpoly->nu[0][i],&fpoly->A[i],
					&fpoly->phase[i]);
		(void) printf("%f %f %f\n",fpoly->nu[0][i],fpoly->A[i],
					fpoly->phase[i]);
	    }
	}
        neg_comp = GAS_COMP2;
        pos_comp = GAS_COMP1;
	contact = (CURVE**)FT_CreateLevelHyperSurfs(front->rect_grid,
			front->interf,neg_comp,pos_comp,func,func_params,
			FIRST_PHYSICS_WAVE_TYPE,&num_segs);

	level_func_pack->func = NULL;
	level_func_pack->point_array = NULL;
	if (debugging("trace"))
	{
	    char dirname[100];
	    sprintf(dirname,"init_intfc-%d",pp_mynode());
	    if (dim == 2)
		xgraph_2d_intfc(dirname,front->interf);
	    else if (dim == 3)
		gview_plot_interface(dirname,front->interf);
	}

	fclose(infile);	
}

void initProjectileIntfc(Front* front,
        LEVEL_FUNC_PACK *level_func_pack)
{
	RECT_GRID *rgr = front->rect_grid;	
	int dim = rgr->dim;
	switch(dim)
	{
	case 2:
	    initProjectileIntfc2d(front,level_func_pack);
	    return;
	case 3:
	    level_func_pack->pos_component = GAS_COMP1;
	    level_func_pack->neg_component = SOLID_COMP;
	    return;
	}
}

void initProjectileIntfc2d(Front* front,
        LEVEL_FUNC_PACK *level_func_pack)
{
	FILE *infile = fopen(InName(front),"r");
	static PROJECTILE_PARAMS *proj_params;
	static RECTANGLE_PARAMS *rparams;
	COMPONENT neg_comp,pos_comp;
	int i,dim,num_segs;
	double  (*func)(POINTER,double*);
        POINTER func_params;
        CURVE **wall,**projectile;
	double gun_length,gun_thickness,gap;
	RECT_GRID *rgr = front->rect_grid;
	RG_PARAMS rgb_params;
	char string[100];

	FT_ScalarMemoryAlloc((POINTER*)&proj_params,sizeof(PROJECTILE_PARAMS));
	FT_ScalarMemoryAlloc((POINTER*)&rparams,sizeof(RECTANGLE_PARAMS));
        proj_params->dim = dim = rgr->dim;
	gap = 0.001*rgr->h[1];
        CursorAfterString(infile,"Enter the center of projectile:");
        for (i = 0; i < dim; ++i)
	{
            fscanf(infile,"%lf",&proj_params->cen[i]);
            (void) printf("%f ",proj_params->cen[i]);
	}
	(void) printf("\n");
        CursorAfterString(infile,"Enter the radius of the projectile:");
        fscanf(infile,"%lf",&proj_params->R);
        (void) printf("%f\n",proj_params->R);
        CursorAfterString(infile,"Enter the head height of the projectile:");
        fscanf(infile,"%lf",&proj_params->r);
        (void) printf("%f\n",proj_params->r);
        CursorAfterString(infile,"Enter the butt height of the projectile:");
        fscanf(infile,"%lf",&proj_params->h);
        (void) printf("%f\n",proj_params->h);
	proj_params->R -= 0.5*gap;

	func_params = (POINTER)proj_params;
	func = projectile_func;
	neg_comp = SOLID_COMP;
	pos_comp = GAS_COMP1;
	projectile = (CURVE**)FT_CreateLevelHyperSurfs(front->rect_grid,
			front->interf,neg_comp,pos_comp,func,func_params,
			MOVABLE_BODY_BOUNDARY,&num_segs);

	if (CursorAfterStringOpt(infile,"Enter yes if there is gun:"))
	{
	    fscanf(infile,"%s",string);
            (void) printf("%s\n",string);
	    if (string[0] == 'Y' || string[0] == 'y')
	    {
            	CursorAfterString(infile,"Enter the gun open position:");
            	fscanf(infile,"%lf",&gun_length);
            	(void) printf("%f\n",gun_length);
            	CursorAfterString(infile,"Enter the gun wall thickness:");
            	fscanf(infile,"%lf",&gun_thickness);
            	(void) printf("%f\n",gun_thickness);
	    	gun_thickness -= 0.5*gap;

	    	rparams->x0 = rgr->L[0] - rgr->h[0];
	    	rparams->y0 = proj_params->cen[1] + 
			(proj_params->R + gap);
	    	rparams->a = gun_length*2.0;
	    	rparams->b = gun_thickness;
		adjust_rectangle_params(rparams,front->rect_grid);

	    	func = rectangle_func;
	    	func_params = (POINTER)rparams;

	    	neg_comp = SOLID_COMP;
	    	pos_comp = GAS_COMP1;
	    	wall = (CURVE**)FT_CreateLevelHyperSurfs(front->rect_grid,
			front->interf,
                    	neg_comp,pos_comp,func,func_params,NEUMANN_BOUNDARY,
			&num_segs);

	    	rparams->y0 = proj_params->cen[1] - 
			(proj_params->R + gap) - gun_thickness;
		adjust_rectangle_params(rparams,front->rect_grid);
	    	wall = (CURVE**)FT_CreateLevelHyperSurfs(front->rect_grid,
			front->interf,
                    	neg_comp,pos_comp,func,func_params,NEUMANN_BOUNDARY,
			&num_segs);
	    }
	}

	level_func_pack->func = NULL;
	level_func_pack->point_array = NULL;

	fclose(infile);	
}

void G_CARTESIAN::initImplosionStates()
{
	int i,j,k,l,index;
	//EQN_PARAMS *eqn_params = (EQN_PARAMS*)front->extra1;
	double coords[MAXD];
	COMPONENT comp;
	STATE *sl,*sr,state;
        POINT *p;
        HYPER_SURF *hs;
        HYPER_SURF_ELEMENT *hse;
	INTERFACE *intfc = front->interf;
	double *dens = field.dens;
	double *engy = field.engy;
	double *pres = field.pres;
	double **momn = field.momn;

        next_point(intfc,NULL,NULL,NULL);
        while (next_point(intfc,&p,&hse,&hs))
        {
	    FT_GetStatesAtPoint(p,hse,hs,(POINTER*)&sl,(POINTER*)&sr);
	    getAmbientState(sl,eqn_params,Coords(p),negative_component(hs));
	    getAmbientState(sr,eqn_params,Coords(p),positive_component(hs));
	}

	FT_MakeGridIntfc(front);
	setDomain();

	switch (dim)
	{
	case 2:
	    for (j = imin[1]; j <= imax[1]; ++j)
	    for (i = imin[0]; i <= imax[0]; ++i)
	    {
		index = d_index2d(i,j,top_gmax);
		comp = top_comp[index];
		getRectangleCenter(index,coords);
	    	getAmbientState(&state,eqn_params,coords,comp);
		dens[index] = state.dens;
		pres[index] = state.pres;
		engy[index] = state.engy;
		for (l = 0; l < dim; ++l)
		    momn[l][index] = state.momn[l];
	    }
	    break;
	case 3:
	    for (k = imin[2]; k <= imax[2]; ++k)
	    for (j = imin[1]; j <= imax[1]; ++j)
	    for (i = imin[0]; i <= imax[0]; ++i)
	    {
		index = d_index3d(i,j,k,top_gmax);
		comp = top_comp[index];
		getRectangleCenter(index,coords);
	    	getAmbientState(&state,eqn_params,coords,comp);
		dens[index] = state.dens;
		pres[index] = state.pres;
		engy[index] = state.engy;
		for (l = 0; l < dim; ++l)
		    momn[l][index] = state.momn[l];
	    }
	    break;
	}
	scatMeshStates();
}	/* end initImplosionStates */

void G_CARTESIAN::initMTFusionStates()
{
	int i,j,k,l,index;
	//EQN_PARAMS *eqn_params = (EQN_PARAMS*)front->extra1;
	double coords[MAXD];
	COMPONENT comp;
	STATE *sl,*sr,state;
        POINT *p;
        HYPER_SURF *hs;
        HYPER_SURF_ELEMENT *hse;
	INTERFACE *intfc = front->interf;
	double *dens = field.dens;
	double *engy = field.engy;
	double *pres = field.pres;
	double **momn = field.momn;

        next_point(intfc,NULL,NULL,NULL);
        while (next_point(intfc,&p,&hse,&hs))
        {
	    FT_GetStatesAtPoint(p,hse,hs,(POINTER*)&sl,(POINTER*)&sr);
	    getAmbientState(sl,eqn_params,Coords(p),negative_component(hs));
	    getAmbientState(sr,eqn_params,Coords(p),positive_component(hs));
	}

	FT_MakeGridIntfc(front);
	setDomain();

	switch (dim)
	{
	case 2:
	    for (j = imin[1]; j <= imax[1]; ++j)
	    for (i = imin[0]; i <= imax[0]; ++i)
	    {
		index = d_index2d(i,j,top_gmax);
		comp = top_comp[index];
		getRectangleCenter(index,coords);
	    	getAmbientState(&state,eqn_params,coords,comp);
		dens[index] = state.dens;
		pres[index] = state.pres;
		engy[index] = state.engy;
		for (l = 0; l < dim; ++l)
		    momn[l][index] = state.momn[l];
	    }
	    break;
	case 3:
	    for (k = imin[2]; k <= imax[2]; ++k)
	    for (j = imin[1]; j <= imax[1]; ++j)
	    for (i = imin[0]; i <= imax[0]; ++i)
	    {
		index = d_index3d(i,j,k,top_gmax);
		comp = top_comp[index];
		getRectangleCenter(index,coords);
	    	getAmbientState(&state,eqn_params,coords,comp);
		dens[index] = state.dens;
		pres[index] = state.pres;
		engy[index] = state.engy;
		for (l = 0; l < dim; ++l)
		    momn[l][index] = state.momn[l];
	    }
	    break;
	}
	scatMeshStates();
}	/* end initMTFusionStates */

static void getAmbientState(
	STATE *state,
	EQN_PARAMS *eqn_params,
	double *coords,
	COMPONENT comp)
{
	EOS_PARAMS	*eos;
	double rho1 = eqn_params->rho1;
	double rho2 = eqn_params->rho2;
	double p1 = eqn_params->p1;
	double p2 = eqn_params->p2;
	double *v1 = eqn_params->v1;
	double *v2 = eqn_params->v2;
	int i,dim;
 
	if (debugging("ambient"))
	    printf("Entering getAmbientState(), coords = %f %f\n",
				coords[0],coords[1]);
	dim = eqn_params->dim;

	/* Constant density */
	for (i = 0; i < dim; ++i)
	    state->vel[i] = state->momn[i] = 0.0;
	state->dim = dim;
	eos = &(eqn_params->eos[comp]);
	state->eos = eos;
	
	switch (comp)
	{
	case GAS_COMP1:
	    state->dens = rho1;
	    state->pres = p1;
	    for (i = 0; i < dim; ++i)
	    {
	    	state->vel[i] = v1[i];
	    	state->momn[i] = rho1*v1[i];
	    }
	    state->engy = EosEnergy(state);
	    break;
	case GAS_COMP2:
	    state->dens = rho2;
	    state->pres = p2;
	    for (i = 0; i < dim; ++i)
	    {
	    	state->vel[i] = v2[i];
	    	state->momn[i] = rho2*v2[i];
	    }
	    state->engy = EosEnergy(state);
	    break;
	case SOLID_COMP:
	    state->dens = 0.0;
	    state->pres = 0.0;
	    state->engy = 0.0;
	    break;
	default:
	    printf("ERROR: Unknown component %d in getAmbientState()!\n",
				comp);
	    clean_up(ERROR);
	}
	if (debugging("ambient_state"))
	    (void) printf("Leaving getAmbientState(): state = %d %f %f %f\n",
			comp,state->dens,state->pres,state->engy);
}	/* end getAmbientState */

void setProjectileParams(
        EQN_PARAMS* eqn_params,
        char *inname)
{
	int i;
	FILE *infile = fopen(inname,"r");
	//EQN_PARAMS *eqn_params = (EQN_PARAMS*)front->extra1;
	double		pinf,einf,gamma;
	char str[256];

	sprintf(str, "Enter gamma, pinf, einf of ambient air:");
	CursorAfterString(infile,str);
	fscanf(infile,"%lf %lf %lf",&gamma,&pinf,&einf);
	(void) printf("%f %f %f\n",gamma,pinf,einf);
	(eqn_params->eos[GAS_COMP1]).gamma = gamma;
	(eqn_params->eos[GAS_COMP1]).pinf = pinf;
	(eqn_params->eos[GAS_COMP1]).einf = einf;
	
	CursorAfterString(infile,"Enter density and pressure of ambient air:");
	fscanf(infile,"%lf %lf",&eqn_params->rho1,&eqn_params->p1);
	(void) printf("%f %f\n",eqn_params->rho1,eqn_params->p1);
	CursorAfterString(infile,"Enter gravity:");
    int dim = eqn_params->dim;
	for (i = 0; i < dim; ++i)
	{
	    fscanf(infile,"%lf",&eqn_params->gravity[i]);
	    (void) printf("%f ",eqn_params->gravity[i]);
	}
	(void) printf("\n");
	fclose(infile);
}

void G_CARTESIAN::initProjectileStates()
{
	int i,j,k,l,index;
	//EQN_PARAMS *eqn_params = (EQN_PARAMS*)front->extra1;
	double coords[MAXD];
	COMPONENT comp;
	STATE *sl,*sr,state;
        POINT *p;
        HYPER_SURF *hs;
        HYPER_SURF_ELEMENT *hse;
	INTERFACE *intfc = front->interf;
	double *dens = field.dens;
	double *engy = field.engy;
	double *pres = field.pres;
	double **momn = field.momn;

        next_point(intfc,NULL,NULL,NULL);
        while (next_point(intfc,&p,&hse,&hs))
        {
	    FT_GetStatesAtPoint(p,hse,hs,(POINTER*)&sl,(POINTER*)&sr);
	    getAmbientState(sl,eqn_params,Coords(p),negative_component(hs));
	    getAmbientState(sr,eqn_params,Coords(p),positive_component(hs));
	}

	FT_MakeGridIntfc(front);
	setDomain();

	switch (dim)
	{
	case 2:
	    for (j = imin[1]; j <= imax[1]; ++j)
	    for (i = imin[0]; i <= imax[0]; ++i)
	    {
		index = d_index2d(i,j,top_gmax);
		comp = top_comp[index];
		getRectangleCenter(index,coords);
	    	getAmbientState(&state,eqn_params,coords,comp);
		dens[index] = state.dens;
		pres[index] = state.pres;
		engy[index] = state.engy;
		for (l = 0; l < dim; ++l)
		    momn[l][index] = state.momn[l];
	    }
	    break;
	case 3:
	    for (k = imin[2]; k <= imax[2]; ++k)
	    for (j = imin[1]; j <= imax[1]; ++j)
	    for (i = imin[0]; i <= imax[0]; ++i)
	    {
		index = d_index3d(i,j,k,top_gmax);
		comp = top_comp[index];
		getRectangleCenter(index,coords);
	    	getAmbientState(&state,eqn_params,coords,comp);
		dens[index] = state.dens;
		pres[index] = state.pres;
		engy[index] = state.engy;
		for (l = 0; l < dim; ++l)
		    momn[l][index] = state.momn[l];
	    }
	    break;
	}
	scatMeshStates();
}	/* end initProjectileStates */

void initRiemannProb(Front* front,
        LEVEL_FUNC_PACK *level_func_pack)
{
	FILE *infile = fopen(InName(front),"r");
	double x;
	char string[100];

	level_func_pack->neg_component = GAS_COMP1;
	level_func_pack->pos_component = GAS_COMP2;
	level_func_pack->func = NULL;
	level_func_pack->func_params = NULL;
	level_func_pack->point_array = NULL;

    int dim = front->rect_grid->dim;
	if (dim == 2)
	{
	    fclose(infile);
	    return;
	}

	EQN_PARAMS* eqn_params = (EQN_PARAMS*)front->extra1;
	eqn_params = (EQN_PARAMS*)front->extra1;

	switch (eqn_params->prob_type)
	{
	case RIEMANN_PROB:
	    level_func_pack->num_points = 1;
	    FT_MatrixMemoryAlloc((POINTER*)&level_func_pack->point_array,
			1,MAXD,sizeof(double));
	    CursorAfterString(infile,"Enter position of interface:");
	    fscanf(infile,"%lf",&x);
	    (void) printf("%f\n",x);
	    level_func_pack->point_array[0][0] = x;
	    level_func_pack->wave_type = FIRST_PHYSICS_WAVE_TYPE;
	    break;
	case ONED_BLAST:
	case ONED_SSINE:
	case ONED_ASINE:
	    break;
	default:
	    (void) printf("Unknow problem type\n");
	    clean_up(ERROR);
	}
	CursorAfterString(infile,"Type yes to track the interface:");
	fscanf(infile,"%s",string);
	(void) printf("%s\n",string);
	if (string[0] == 'y' || string[0] == 'Y')
            eqn_params->tracked = YES;
        else
            eqn_params->tracked = NO;
	fclose(infile);
}

void setRiemProbParams(
        EQN_PARAMS* eqn_params,
        char *inname)
{
	switch (eqn_params->dim)
	{
	case 1:
	    setRiemProbParams1d(eqn_params,inname);
	    return;
	case 2:
	    setRiemProbParams2d(eqn_params,inname);
	    return;
	}
}

void setRiemProbParams2d(
        EQN_PARAMS* eqn_params,char *inname)
{
	FILE *infile = fopen(inname,"r");
	//EQN_PARAMS *eqn_params = (EQN_PARAMS*)front->extra1;
	double		pinf,einf,gamma;
	char str[256];

	sprintf(str, "Enter gamma, pinf, einf of ambient air:");
	CursorAfterString(infile,str);
	fscanf(infile,"%lf %lf %lf",&gamma,&pinf,&einf);
	(void) printf("%f %f %f\n",gamma,pinf,einf);
	(eqn_params->eos[GAS_COMP1]).gamma = gamma;
	(eqn_params->eos[GAS_COMP1]).pinf = pinf;
	(eqn_params->eos[GAS_COMP1]).einf = einf;
	(eqn_params->eos[GAS_COMP2]).gamma = gamma;
	(eqn_params->eos[GAS_COMP2]).pinf = pinf;
	(eqn_params->eos[GAS_COMP2]).einf = einf;
	
	(void) printf("There four states in four quadrants\n");
	(void) printf("Enter states counter-clock-wise, "
			"starting from first quadrant\n");
	CursorAfterString(infile,"Enter density of four states:");
	fscanf(infile,"%lf %lf %lf %lf",
			&eqn_params->rho0,&eqn_params->rho1,
			&eqn_params->rho2,&eqn_params->rho3);
	(void) printf("%f %f %f %f\n",
			eqn_params->rho0,eqn_params->rho1,
			eqn_params->rho2,eqn_params->rho3);
	CursorAfterString(infile,"Enter pressure of four states:");
	fscanf(infile,"%lf %lf %lf %lf",
			&eqn_params->p0,&eqn_params->p1,
			&eqn_params->p2,&eqn_params->p3);
	(void) printf("%f %f %f %f\n",
			eqn_params->p0,eqn_params->p1,
			eqn_params->p2,eqn_params->p3);
	CursorAfterString(infile,"Enter x-velocity of four states:");
	fscanf(infile,"%lf %lf %lf %lf",
			&eqn_params->v0[0],&eqn_params->v1[0],
			&eqn_params->v2[0],&eqn_params->v3[0]);
	(void) printf("%f %f %f %f\n",
			eqn_params->v0[0],eqn_params->v1[0],
			eqn_params->v2[0],eqn_params->v3[0]);
	CursorAfterString(infile,"Enter y-velocity of four states:");
	fscanf(infile,"%lf %lf %lf %lf",
			&eqn_params->v0[1],&eqn_params->v1[1],
			&eqn_params->v2[1],&eqn_params->v3[1]);
	(void) printf("%f %f %f %f\n",
			eqn_params->v0[1],eqn_params->v1[1],
			eqn_params->v2[1],eqn_params->v3[1]);
	fclose(infile);
}

void setRiemProbParams1d(
        EQN_PARAMS* eqn_params,char *inname)
{
	FILE *infile = fopen(inname,"r");
	//EQN_PARAMS *eqn_params = (EQN_PARAMS*)front->extra1;
	double		pinf,einf,gamma;
	char str[256];

	sprintf(str, "Enter gamma, pinf, einf of ambient air:");
	CursorAfterString(infile,str);
	fscanf(infile,"%lf %lf %lf",&gamma,&pinf,&einf);
	(void) printf("%f %f %f\n",gamma,pinf,einf);
	(eqn_params->eos[GAS_COMP1]).gamma = gamma;
	(eqn_params->eos[GAS_COMP1]).pinf = pinf;
	(eqn_params->eos[GAS_COMP1]).einf = einf;
	(eqn_params->eos[GAS_COMP2]).gamma = gamma;
	(eqn_params->eos[GAS_COMP2]).pinf = pinf;
	(eqn_params->eos[GAS_COMP2]).einf = einf;
	
	CursorAfterString(infile,"Enter left and right density:");
	fscanf(infile,"%lf %lf",&eqn_params->rho1,&eqn_params->rho2);
	(void) printf("%f %f\n",eqn_params->rho1,eqn_params->rho2);
	CursorAfterString(infile,"Enter left and right pressure:");
	fscanf(infile,"%lf %lf",&eqn_params->p1,&eqn_params->p2);
	(void) printf("%f %f\n",eqn_params->p1,eqn_params->p2);
	CursorAfterString(infile,"Enter left and right velocity:");
	fscanf(infile,"%lf %lf",&eqn_params->v1[0],&eqn_params->v2[0]);
	(void) printf("%f %f\n",eqn_params->v1[0],eqn_params->v2[0]);
	fclose(infile);
}

void setOnedParams(
        EQN_PARAMS* eqn_params,char *inname)
{
	FILE *infile = fopen(inname,"r");
	//EQN_PARAMS *eqn_params = (EQN_PARAMS*)front->extra1;
	double		pinf,einf,gamma;
	char str[256];

	sprintf(str, "Enter gamma, pinf, einf of ambient air:");
        CursorAfterString(infile,str);
        fscanf(infile,"%lf %lf %lf",&gamma,&pinf,&einf);
        (void) printf("%f %f %f\n",gamma,pinf,einf);
        (eqn_params->eos[GAS_COMP1]).gamma = gamma;
        (eqn_params->eos[GAS_COMP1]).pinf = pinf;
        (eqn_params->eos[GAS_COMP1]).einf = einf;
        (eqn_params->eos[GAS_COMP2]).gamma = gamma;
        (eqn_params->eos[GAS_COMP2]).pinf = pinf;
        (eqn_params->eos[GAS_COMP2]).einf = einf;
}

void G_CARTESIAN::initRiemProbStates()
{
	int i,j,l,index;
	//EQN_PARAMS *eqn_params = (EQN_PARAMS*)front->extra1;
	double coords[MAXD];
	double center[MAXD];
	COMPONENT comp;
	STATE *sl,*sr,state;
        POINT *p;
        HYPER_SURF *hs;
        HYPER_SURF_ELEMENT *hse;
	INTERFACE *intfc = front->interf;
	double *dens = field.dens;
	double *engy = field.engy;
	double *pres = field.pres;
	double **momn = field.momn;

        next_point(intfc,NULL,NULL,NULL);
        while (next_point(intfc,&p,&hse,&hs))
        {
	    FT_GetStatesAtPoint(p,hse,hs,(POINTER*)&sl,(POINTER*)&sr);
	    getAmbientState(sl,eqn_params,Coords(p),negative_component(hs));
	    getAmbientState(sr,eqn_params,Coords(p),positive_component(hs));
	}

	FT_MakeGridIntfc(front);
	setDomain();

	switch (dim)
	{
	case 1:
	    for (i = imin[0]; i <= imax[0]; ++i)
	    {
		index = d_index1d(i,top_gmax);
		comp = top_comp[index];
		getRectangleCenter(index,coords);
	    	getAmbientState(&state,eqn_params,coords,comp);
		dens[index] = state.dens;
		pres[index] = state.pres;
		engy[index] = state.engy;
		for (l = 0; l < dim; ++l)
		    momn[l][index] = state.momn[l];
	    }
	    break;
	case 2:
	    for (l = 0; l < dim; ++l)
		center[l] = 0.5*(top_grid->GL[l] + top_grid->GU[l]);
	    state.dim = 2;
	    for (i = imin[0]; i <= imax[0]; ++i)
	    for (j = imin[1]; j <= imax[1]; ++j)
	    {
		index = d_index2d(i,j,top_gmax);
		comp = top_comp[index];
        	state.eos = &(eqn_params->eos[comp]);
		getRectangleCenter(index,coords);
		if (coords[0] > center[0] && coords[1] > center[1])
		{
		    /* First quadrant */
		    state.dens = dens[index] = eqn_params->rho0;
		    state.pres = pres[index] = eqn_params->p0;
		    for (l = 0; l < dim; ++l)
		    {
		    	momn[l][index] = eqn_params->rho0*eqn_params->v0[l];
			state.momn[l] = eqn_params->rho0*eqn_params->v0[l];
		    }
		}
		else if (coords[0] <= center[0] && coords[1] > center[1])
		{
		    /* Second quadrant */
		    state.dens = dens[index] = eqn_params->rho1;
		    state.pres = pres[index] = eqn_params->p1;
		    for (l = 0; l < dim; ++l)
		    {
		    	momn[l][index] = eqn_params->rho1*eqn_params->v1[l];
			state.momn[l] = eqn_params->rho1*eqn_params->v1[l];
		    }
		}
		else if (coords[0] <= center[0] && coords[1] <= center[1])
		{
		    /* Third quadrant */
		    state.dens = dens[index] = eqn_params->rho2;
		    state.pres = pres[index] = eqn_params->p2;
		    for (l = 0; l < dim; ++l)
		    {
		    	momn[l][index] = eqn_params->rho2*eqn_params->v2[l];
			state.momn[l] = eqn_params->rho2*eqn_params->v2[l];
		    }
		}
		else if (coords[0] > center[0] && coords[1] <= center[1])
		{
		    /* Fourth quadrant */
		    state.dens = dens[index] = eqn_params->rho3;
		    state.pres = pres[index] = eqn_params->p3;
		    for (l = 0; l < dim; ++l)
		    {
		    	momn[l][index] = eqn_params->rho3*eqn_params->v3[l];
			state.momn[l] = eqn_params->rho3*eqn_params->v3[l];
		    }
		}
		engy[index] = EosEnergy(&state);
	    }
	    break;
	case 3:
	default:
	    (void) printf("initRiemProbStates() not for case dim = %d\n",dim);
	    clean_up(ERROR);
	}
	scatMeshStates();
}	/* end initRiemProbStates */

void G_CARTESIAN::initBlastWaveStates()
{
	int i,l,index;
	//EQN_PARAMS *eqn_params = (EQN_PARAMS*)front->extra1;
	double coords[MAXD];
	COMPONENT comp;
	STATE *sl,*sr,state;
        POINT *p;
        HYPER_SURF *hs;
        HYPER_SURF_ELEMENT *hse;
	INTERFACE *intfc = front->interf;
	double *dens = field.dens;
	double *engy = field.engy;
	double *pres = field.pres;
	double **momn = field.momn;

        next_point(intfc,NULL,NULL,NULL);
        while (next_point(intfc,&p,&hse,&hs))
        {
	    FT_GetStatesAtPoint(p,hse,hs,(POINTER*)&sl,(POINTER*)&sr);
	    getBlastState(sl,eqn_params,Coords(p),negative_component(hs));
	    getBlastState(sr,eqn_params,Coords(p),positive_component(hs));
	}

	FT_MakeGridIntfc(front);
	setDomain();

	switch (dim)
	{
	case 1:
	    for (i = imin[0]; i <= imax[0]; ++i)
	    {
		index = d_index1d(i,top_gmax);
		comp = top_comp[index];
		getRectangleCenter(index,coords);
	    	getBlastState(&state,eqn_params,coords,comp);
		dens[index] = state.dens;
		pres[index] = state.pres;
		engy[index] = state.engy;
		for (l = 0; l < dim; ++l)
		    momn[l][index] = state.momn[l];
	    }
	    break;
	case 2:
	case 3:
	default:
	    (void) printf("initBlastWaveStates() not for case dim = %d\n",dim);
	    clean_up(ERROR);
	}
	scatMeshStates();
}	/* end initBlastWaveStates */

void G_CARTESIAN::initShockSineWaveStates()
{
	int i,l,index;
	//EQN_PARAMS *eqn_params = (EQN_PARAMS*)front->extra1;
	double coords[MAXD];
	COMPONENT comp;
	STATE *sl,*sr,state;
        POINT *p;
        HYPER_SURF *hs;
        HYPER_SURF_ELEMENT *hse;
	INTERFACE *intfc = front->interf;
	double *dens = field.dens;
	double *engy = field.engy;
	double *pres = field.pres;
	double **momn = field.momn;

        next_point(intfc,NULL,NULL,NULL);
        while (next_point(intfc,&p,&hse,&hs))
        {
	    FT_GetStatesAtPoint(p,hse,hs,(POINTER*)&sl,(POINTER*)&sr);
	    getShockSineWaveState(sl,eqn_params,Coords(p),
				negative_component(hs));
	    getShockSineWaveState(sr,eqn_params,Coords(p),
				positive_component(hs));
	}

	FT_MakeGridIntfc(front);
	setDomain();

	switch (dim)
	{
	case 1:
	    for (i = imin[0]; i <= imax[0]; ++i)
	    {
		index = d_index1d(i,top_gmax);
		comp = top_comp[index];
		getRectangleCenter(index,coords);
	    	getShockSineWaveState(&state,eqn_params,coords,comp);
		dens[index] = state.dens;
		pres[index] = state.pres;
		engy[index] = state.engy;
		for (l = 0; l < dim; ++l)
		    momn[l][index] = state.momn[l];
	    }
	    break;
	case 2:
	case 3:
	default:
	    (void) printf("initShockSineWaveStates() not for case dim = %d\n",
				dim);
	    clean_up(ERROR);
	}
	scatMeshStates();
}	/* end initShockSineWaveStates */

void G_CARTESIAN::initAccuracySineWaveStates()
{
	int i,l,index;
	//EQN_PARAMS *eqn_params = (EQN_PARAMS*)front->extra1;
	double coords[MAXD];
	COMPONENT comp;
	STATE *sl,*sr,state;
        POINT *p;
        HYPER_SURF *hs;
        HYPER_SURF_ELEMENT *hse;
	INTERFACE *intfc = front->interf;
	double *dens = field.dens;
	double *engy = field.engy;
	double *pres = field.pres;
	double **momn = field.momn;

        next_point(intfc,NULL,NULL,NULL);
        while (next_point(intfc,&p,&hse,&hs))
        {
	    FT_GetStatesAtPoint(p,hse,hs,(POINTER*)&sl,(POINTER*)&sr);
	    getAccuracySineWaveState(sl,eqn_params,Coords(p),
				negative_component(hs));
	    getAccuracySineWaveState(sr,eqn_params,Coords(p),
				positive_component(hs));
	}

	FT_MakeGridIntfc(front);
	setDomain();

	switch (dim)
	{
	case 1:
	    for (i = imin[0]; i <= imax[0]; ++i)
	    {
		index = d_index1d(i,top_gmax);
		comp = top_comp[index];
		getRectangleCenter(index,coords);
	    	getAccuracySineWaveState(&state,eqn_params,coords,comp);
		dens[index] = state.dens;
		pres[index] = state.pres;
		engy[index] = state.engy;
		for (l = 0; l < dim; ++l)
		    momn[l][index] = state.momn[l];
	    }
	    break;
	case 2:
	case 3:
	default:
	    (void) printf("initAcuracySineWaveStates() not for case dim = %d\n",
				dim);
	    clean_up(ERROR);
	}
	scatMeshStates();
}	/* end initAccuracySineWaveStates */

static void getBlastState(
	STATE *state,
	EQN_PARAMS *eqn_params,
	double *coords,
	COMPONENT comp)
{
	EOS_PARAMS	*eos;
	eos = &(eqn_params->eos[comp]);
	state->eos = eos;
	state->dim = 1;
	if (coords[0] < 0.1)
	{
	    state->dens = 1.0;
	    state->vel[0] = 0.0;
	    state->pres = 1000.0;
	}
	else if (coords[0] > 0.9)
	{
	    state->dens = 1.0;
	    state->vel[0] = 0.0;
	    state->pres = 100.0;
	}
	else
	{
	    state->dens = 1.0;
	    state->vel[0] = 0.0;
	    state->pres = 0.01;
	}
	state->momn[0] = state->vel[0] * state->dens;
	state->engy = EosInternalEnergy(state);
}	/* end getBlastState */

static void getShockSineWaveState(
	STATE *state,
	EQN_PARAMS *eqn_params,
	double *coords,
	COMPONENT comp)
{
	EOS_PARAMS	*eos;
	eos = &(eqn_params->eos[comp]);
	state->eos = eos;
	state->dim = 1;
	if (coords[0] < -4.0)
	{
	    state->dens = 3.857143;
	    state->momn[0] = 2.629329;
	    state->pres = 10.333333;
	}
	else
	{
	    state->dens = 1.0 + 0.2*sin(5.0*coords[0]);
	    state->momn[0] = 0.0;
	    state->pres = 1.0;
	}
	state->vel[0] = state->momn[0]/state->dens;
	state->engy = EosEnergy(state);
}	/* end getShockSineWaveState */

static void getAccuracySineWaveState(
	STATE *state,
	EQN_PARAMS *eqn_params,
	double *coords,
	COMPONENT comp)
{
	EOS_PARAMS	*eos;
	eos = &(eqn_params->eos[comp]);
	state->eos = eos;
	state->dim = 1;
	state->dens = 1.0 + 0.2 * sin(PI*coords[0]);
	state->vel[0] = 0.7;
	state->pres = 1.0;
	state->momn[0]= state->dens * state->vel[0];
	state->engy = EosEnergy(state);
}	/* end getAccuracySineWaveState */

void initRectPlaneIntfc(Front* front,
        LEVEL_FUNC_PACK *level_func_pack)
{
    FILE *infile = fopen(InName(front),"r");
    static RECT_BOX_PARAMS *rect_params;
    int i,dim;

	EQN_PARAMS* eqn_params = (EQN_PARAMS*)front->extra1;
    PROB_TYPE prob_type = eqn_params->prob_type;

    FT_ScalarMemoryAlloc((POINTER*)&rect_params,sizeof(RECT_BOX_PARAMS));
    rect_params->dim = dim = front->rect_grid->dim;
    CursorAfterString(infile,"Enter the center of the rectangle:");
    for (i = 0; i < dim; ++i)
    {
        fscanf(infile,"%lf",&rect_params->center[i]);
        (void) printf("%f ",rect_params->center[i]);
    }
    (void) printf("\n");
    CursorAfterString(infile,"Enter lengths of the rectangle:");
    for (i = 0; i < dim; ++i)
    {
        fscanf(infile,"%lf",&rect_params->length[i]);
        (void) printf("%f\n",rect_params->length[i]);
    }
    (void) printf("\n");

    level_func_pack->func_params = (POINTER)rect_params;

    switch (prob_type)
    {
    case FLUID_SOLID_RECT:
        level_func_pack->neg_component = SOLID_COMP;
        level_func_pack->pos_component = GAS_COMP1;
        level_func_pack->func = rect_box_func;
        level_func_pack->wave_type = MOVABLE_BODY_BOUNDARY;
        break;
    default:
        (void) printf("ERROR: entering wrong initialization function\n");
        clean_up(ERROR);
    }

    fclose(infile);
}

void initTrianglePlaneIntfc(Front* front,
        LEVEL_FUNC_PACK *level_func_pack)
{
    FILE *infile = fopen(InName(front),"r");
    static TRIANGLE_PARAMS *tri_params;
    int i,dim;
    char msg[100];
    
    EQN_PARAMS* eqn_params = (EQN_PARAMS*)front->extra1;
    PROB_TYPE prob_type = eqn_params->prob_type;

    FT_ScalarMemoryAlloc((POINTER*)&tri_params,sizeof(TRIANGLE_PARAMS));

    CursorAfterString(infile,"Triangle is specified by three vertices");
    (void) printf("\n");
    for (i = 0; i < 3; ++i)
    {
        sprintf(msg,"Enter coordinates of point %d:",i+1);
        CursorAfterString(infile,msg);
        fscanf(infile,"%lf %lf",&tri_params->x[i],&tri_params->y[i]);
        (void) printf("%f %f\n",tri_params->x[i],tri_params->y[i]);
    }

    level_func_pack->func_params = (POINTER)tri_params;

    switch (prob_type)
    {
    case FLUID_SOLID_TRIANGLE:
        level_func_pack->neg_component = SOLID_COMP;
        level_func_pack->pos_component = GAS_COMP1;
        level_func_pack->func = triangle_func;
        //level_func_pack->wave_type = NEUMANN_BOUNDARY;
        level_func_pack->wave_type = MOVABLE_BODY_BOUNDARY;
        break;
    default:
        (void) printf("ERROR: entering wrong initialization function\n");
        clean_up(ERROR);
    }
    fclose(infile);
}

void initCylinderPlaneIntfc(Front* front,
        LEVEL_FUNC_PACK *level_func_pack)
{
        FILE *infile = fopen(InName(front),"r");
        static CYLINDER_PARAMS *cylinder_params;
        int i;

        EQN_PARAMS* eqn_params = (EQN_PARAMS*)front->extra1;
        PROB_TYPE prob_type = eqn_params->prob_type;

        FT_ScalarMemoryAlloc((POINTER*)&cylinder_params,sizeof(CYLINDER_PARAMS));
        CursorAfterString(infile,"Enter the center of the cylinder:");
        for (i = 0; i < 3; ++i)
        {
            fscanf(infile,"%lf",&cylinder_params->center[i]);
            (void) printf("%f ",cylinder_params->center[i]);
        }
        (void) printf("\n");
        CursorAfterString(infile,"Enter radius of the cylinder:");
        fscanf(infile,"%lf",&cylinder_params->radius);
        (void) printf("%f\n",cylinder_params->radius);
        CursorAfterString(infile,"Enter height of the cylinder:");
        fscanf(infile,"%lf",&cylinder_params->height);
        (void) printf("%f\n",cylinder_params->height);

        level_func_pack->func_params = (POINTER)cylinder_params;

        switch (prob_type)
        {
        case FLUID_SOLID_CYLINDER:
            level_func_pack->neg_component = SOLID_COMP;
            level_func_pack->pos_component = GAS_COMP1;
            level_func_pack->func = cylinder_func;
            //level_func_pack->wave_type = NEUMANN_BOUNDARY;
            level_func_pack->wave_type = MOVABLE_BODY_BOUNDARY;
            break;
        default:
            (void) printf("ERROR: entering wrong initialization function\n");
            clean_up(ERROR);
        }
        fclose(infile);
}

extern  void prompt_for_rigid_body_params(
        int dim,
        char *inname,
        RG_PARAMS *rgb_params)
{
        int i;
        char msg[100],s[100],ss[100];
        FILE *infile = fopen(inname,"r");
        boolean is_preset_motion = NO;
        double mag_dir;

        if (debugging("rgbody"))
            (void) printf("Enter prompt_for_rigid_body_params()\n");

        rgb_params->dim = dim;
        CursorAfterString(infile,"Type yes if motion is preset: ");
        fscanf(infile,"%s",s);
        (void) printf("%s\n",s);
        if (s[0] == 'y' || s[0] == 'Y')
        {
            (void) printf("Available preset motion types are:\n");
            (void) printf("\tPRESET_TRANSLATION\n");
            (void) printf("\tPRESET_ROTATION\n");
            (void) printf("\tPRESET_MOTION (general)\n");
            CursorAfterString(infile,"Enter type of preset motion: ");
            fscanf(infile,"%s",s);
            (void) printf("%s\n",s);
            switch(s[7])
            {
            case 'M':
                rgb_params->motion_type = PRESET_MOTION;
                break;
            case 'C':
                rgb_params->motion_type = PRESET_COM_MOTION;
                break;
            case 'T':
                rgb_params->motion_type = PRESET_TRANSLATION;
                break;
            case 'R':
                rgb_params->motion_type = PRESET_ROTATION;
                break;
            default:
                (void) printf("Unknow type of preset motion!\n");
                clean_up(ERROR);
            }
        }
        else
        {
            (void) printf("Available dynamic motion types are:\n");
            (void) printf("\tFREE_MOTION:\n");
            (void) printf("\tCOM_MOTION (center of mass):\n");
            (void) printf("\tTRANSLATION:\n");
            (void) printf("\tROTATION:\n");
            CursorAfterString(infile,"Enter type of dynamic motion: ");
            fscanf(infile,"%s",s);
            (void) printf("%s\n",s);
            switch(s[0])
            {
            case 'F':
                rgb_params->motion_type = FREE_MOTION;
                break;
            case 'C':
                rgb_params->motion_type = COM_MOTION;
                break;
            case 'T':
                rgb_params->motion_type = TRANSLATION;
                break;
            case 'R':
                rgb_params->motion_type = ROTATION;
                break;
            default:
                (void) printf("Unknow type of motion!\n");
                clean_up(ERROR);
            }
        }
        if (rgb_params->motion_type == TRANSLATION ||
            rgb_params->motion_type == PRESET_TRANSLATION)
        {
            mag_dir = 0.0;
            CursorAfterString(infile,"Enter the direction of motion:");
            for (i = 0; i < dim; ++i)
            {
                fscanf(infile,"%lf",&rgb_params->translation_dir[i]);
                (void) printf("%f ",rgb_params->translation_dir[i]);
                mag_dir += sqr(rgb_params->translation_dir[i]);
            }
            (void) printf("\n");
            mag_dir = sqrt(mag_dir);
            for (i = 0; i < dim; ++i)
                rgb_params->translation_dir[i] /= mag_dir;
        }
        if (rgb_params->motion_type == FREE_MOTION ||
            rgb_params->motion_type == COM_MOTION ||
            rgb_params->motion_type == TRANSLATION)
        {
            sprintf(msg,"Enter the total mass for rigid body:");
            CursorAfterString(infile,msg);
            fscanf(infile,"%lf",&rgb_params->total_mass);
            (void) printf("%f\n",rgb_params->total_mass);
        }
        if (rgb_params->motion_type == FREE_MOTION ||
            rgb_params->motion_type == COM_MOTION ||
            rgb_params->motion_type == TRANSLATION ||
            rgb_params->motion_type == PRESET_MOTION ||
            rgb_params->motion_type == PRESET_TRANSLATION)
        {
            sprintf(msg,"Enter the initial center of mass for rigid body:");
            CursorAfterString(infile,msg);
            for (i = 0; i < dim; ++i)
            {
                fscanf(infile,"%lf",&rgb_params->center_of_mass[i]);
                (void) printf("%f ",rgb_params->center_of_mass[i]);
            }
            (void) printf("\n");
            sprintf(msg,"Enter the initial center of mass velocity:");
            CursorAfterString(infile,msg);
            for (i = 0; i < dim; ++i)
            {
                fscanf(infile,"%lf",&rgb_params->cen_of_mass_velo[i]);
                (void) printf("%f ",rgb_params->cen_of_mass_velo[i]);
            }
            (void) printf("\n");
        }
        if (rgb_params->motion_type == PRESET_ROTATION)
        {
            /* 2D preset rotation is always about the z-axis */
            /* 3D preset rotation axis along rotation_dir */
            if (dim == 3)
            {
                mag_dir = 0.0;
                CursorAfterString(infile,"Enter the direction of rotation:");
                for (i = 0; i < dim; ++i)
                {
                    fscanf(infile,"%lf",&rgb_params->rotation_dir[i]);
                    (void) printf("%f ",rgb_params->rotation_dir[i]);
                    mag_dir += sqr(rgb_params->rotation_dir[i]);
                }
                (void) printf("\n");
                mag_dir = sqrt(mag_dir);
                for (i = 0; i < dim; ++i)
                    rgb_params->rotation_dir[i] /= mag_dir;
                /* initialize the euler parameters */
                rgb_params->euler_params[0] = 1.0;
                for (i = 1; i < 4; ++i)
                    rgb_params->euler_params[i] = 0.0;
            }
            /* Center of axis is the coordinate of a point on the axis */
            CursorAfterString(infile,"Enter rotation center:");
            for (i = 0; i < dim; ++i)
            {
                fscanf(infile,"%lf",&rgb_params->rotation_cen[i]);
                (void) printf("%f ",rgb_params->rotation_cen[i]);
            }
            (void) printf("\n");
            CursorAfterString(infile,"Enter preset angular velocity:");
            fscanf(infile,"%lf",&rgb_params->angular_velo);
            (void) printf("%f\n",rgb_params->angular_velo);
            if (dim == 3)
            {
                /* used to update the maximum speed in 3D cases */
                for (i = 0; i < dim; ++i)
                    rgb_params->p_angular_velo[i] = rgb_params->angular_velo
                                        * rgb_params->rotation_dir[i];
            }
        }
        if (rgb_params->motion_type == ROTATION)
        {
            if (CursorAfterStringOpt(infile,
                "Type yes if rigid body will rotate about an point:"))
            {
                fscanf(infile,"%s",s);
                (void) printf("%s\n",s);
                if (s[0] == 'y' || s[0] == 'Y')
                {
                    sprintf(msg,"Enter rotation center:");
                    CursorAfterString(infile,msg);
                    for (i = 0; i < dim; ++i)
                    {
                        fscanf(infile,"%lf",&rgb_params->rotation_cen[i]);
                        (void) printf("%f ",rgb_params->rotation_cen[i]);
                    }
                    (void) printf("\n");
                }
            }
            if (CursorAfterStringOpt(infile,
                "Type yes if rigid body will rotate about an axis:"))
            {
                fscanf(infile,"%s",s);
                (void) printf("%s\n",s);
                if (s[0] == 'y' || s[0] == 'Y')
                {
                    /* For 2D, it is always about the z-axis */
                    if (dim == 3)
                    {
                        sprintf(msg,"Enter direction of the axis:");
                        CursorAfterString(infile,msg);
                        for (i = 0; i < dim; ++i)
                        {
                            fscanf(infile,"%lf",&rgb_params->rotation_dir[i]);
                            (void) printf("%f ",rgb_params->rotation_dir[i]);
                            mag_dir += sqr(rgb_params->rotation_dir[i]);
                        }
                        mag_dir = sqrt(mag_dir);
                        for (i = 0; i < dim; ++i)
                            rgb_params->rotation_dir[i] /= mag_dir;
                        (void) printf("\n");
                    }
                }
            }
        }
        if (rgb_params->motion_type == FREE_MOTION ||
            rgb_params->motion_type == ROTATION)
        {
            CursorAfterString(infile,"Enter the moment of inertia: ");
            if (dim == 2)
            {
                fscanf(infile,"%lf",&rgb_params->moment_of_inertia);
                (void) printf("%f\n",rgb_params->moment_of_inertia);
            }
            else if (dim == 3)
            {
                for (i = 0; i < dim; ++i)
                {
                    fscanf(infile,"%lf",&rgb_params->p_moment_of_inertia[i]);
                    (void) printf("%f ",rgb_params->p_moment_of_inertia[i]);
                }
                (void) printf("\n");
            }
            CursorAfterString(infile,"Enter initial angular velocity: ");
            if (dim == 2)
            {
                fscanf(infile,"%lf",&rgb_params->angular_velo);
                (void) printf("%f\n",rgb_params->angular_velo);
            }
            else if (dim == 3)
            {
                for (i = 0; i < dim; ++i)
                {
                    fscanf(infile,"%lf",&rgb_params->p_angular_velo[i]);
                    (void) printf("%f ",rgb_params->p_angular_velo[i]);
                }
                (void) printf("\n");
                /* initialize the euler parameters */
                rgb_params->euler_params[0] = 1.0;
                for (i = 1; i < 4; ++i)
                    rgb_params->euler_params[i] = 0.0;
            }
        }

        if (debugging("rgbody"))
            (void) printf("Leaving prompt_for_rigid_body_params()\n");
}       /* end prompt_for_rigid_body_params */

extern void set_rgbody_params(
        RG_PARAMS rg_params,
        HYPER_SURF *hs)
{
        int i,dim = rg_params.dim;
        total_mass(hs) = rg_params.total_mass;
        mom_inertia(hs) = rg_params.moment_of_inertia;
        angular_velo(hs) = rg_params.angular_velo;
        motion_type(hs) = rg_params.motion_type;
        surface_tension(hs) = 0.0;
        for (i = 0; i < dim; ++i)
        {
            center_of_mass(hs)[i] = rg_params.center_of_mass[i];
            center_of_mass_velo(hs)[i] =
                                rg_params.cen_of_mass_velo[i];
            rotation_center(hs)[i] =
                                rg_params.rotation_cen[i];
            translation_dir(hs)[i] = rg_params.translation_dir[i];
            if (dim == 3)
            {
                rotation_direction(hs)[i] = rg_params.rotation_dir[i];
                p_mom_inertia(hs)[i] = rg_params.p_moment_of_inertia[i];
                p_angular_velo(hs)[i] = rg_params.p_angular_velo[i];
            }
        }
        if (dim == 3)
        {
            for (i = 0; i < 4; i++)
                euler_params(hs)[i] = rg_params.euler_params[i];
        }
}       /* end set_rgbody_params */

static double getStationaryVelocity(
	EQN_PARAMS *eqn_params)
{
	RIEMANN_INPUT input;
	RIEMANN_SOLN riem_soln;
	int idir = eqn_params->idir;
	int shock_side = eqn_params->shock_side;
	STATE state;
	double shock_speed,contact_speed;
	double Mach_number = eqn_params->Mach_number;
	int dim = eqn_params->dim;

	input.left_state.d = eqn_params->rho1;
	input.right_state.d = eqn_params->rho2;
	input.left_state.gamma = eqn_params->eos[GAS_COMP1].gamma;
	input.right_state.gamma = eqn_params->eos[GAS_COMP2].gamma;

	input.left_state.p = input.right_state.p = eqn_params->p0;
	input.left_state.u = input.right_state.u = 0.0;
	state.dim = dim;
	if (shock_side == 1)
	{
	    state.dens = eqn_params->rho1;
	    state.pres = eqn_params->p1;
	    state.eos = &(eqn_params->eos[GAS_COMP1]);
	    state.engy = EosInternalEnergy(&state);
	}
	else if (shock_side == -1)
	{
	    state.dens = eqn_params->rho2;
	    state.pres = eqn_params->p2;
	    state.eos = &(eqn_params->eos[GAS_COMP2]);
	    state.engy = EosInternalEnergy(&state);
	}
	behind_state(SHOCK_MACH_NUMBER,Mach_number,&shock_speed,
			idir,shock_side,&state,&state);

	if (shock_side == 1)
        {
	    input.left_state.d = state.dens;
	    input.left_state.p = state.pres;
	    input.left_state.u = state.vel[idir];
	}
        else if (shock_side == -1)
        {
	    input.right_state.d = state.dens;
	    input.right_state.p = state.pres;
	    input.right_state.u = state.vel[idir];
	}
	RiemannSolution(input,&riem_soln);
	contact_speed = riem_soln.contact.speed_contact;
	return contact_speed;
}	/* end getStationaryVelocity */

void initRandPertIntfc(Front* front,
        LEVEL_FUNC_PACK *level_func_pack)
{
    static FOURIER_POLY_MD *level_func_params;
    FILE *infile = fopen(InName(front),"r");
    int i,j,num_modes;
    char mesg[100], s[10];
    static double   L[3], U[3];

    FT_ScalarMemoryAlloc((POINTER*)&level_func_params,
            sizeof(FOURIER_POLY_MD));

    int dim = level_func_params->dim = front->rect_grid->dim;
    ft_assign(L, front->rect_grid->L, 3*DOUBLE);
    ft_assign(U, front->rect_grid->U, 3*DOUBLE);

    level_func_params->L = L;
    level_func_params->U = U;

    level_func_pack->neg_component = GAS_COMP1;
    level_func_pack->pos_component = GAS_COMP2;
    CursorAfterString(infile,"Enter mean position of fluid interface:");
    fscanf(infile,"%lf",&level_func_params->z0);
    (void) printf("%f\n",level_func_params->z0);
    CursorAfterString(infile,"Enter number of Fourier modes:");
    fscanf(infile,"%d",&num_modes);
    (void) printf("%d\n",num_modes);
    CursorAfterString(infile, "Enter position of cross section:");
    fscanf(infile,"%lf",&level_func_params->croSectAtz);
    (void) printf("%f\n",level_func_params->croSectAtz);
    level_func_params->num_modes = num_modes;

    FT_MatrixMemoryAlloc((POINTER*)&level_func_params->A,num_modes+1,
                            num_modes+1,sizeof(double));
    FT_MatrixMemoryAlloc((POINTER*)&level_func_params->B,num_modes+1,                                       num_modes+1,sizeof(double));
    FT_MatrixMemoryAlloc((POINTER*)&level_func_params->C,num_modes+1,                                       num_modes+1,sizeof(double));
    FT_MatrixMemoryAlloc((POINTER*)&level_func_params->D,num_modes+1,                                       num_modes+1,sizeof(double));

    pertCoeff(level_func_params->A,level_func_params->B,
                    level_func_params->C,level_func_params->D,
                    num_modes,level_func_params->L,level_func_params->U);

    CursorAfterString(infile,"Type yes to print coefficient:");
    fscanf(infile,"%s",s);
    (void) printf("%s\n",s);
    if (s[0] == 'y' || s[0] == 'Y')
    {
        for (i = 0; i <= num_modes; i++)
        {
             for (j = 0; j <= num_modes; j++)
             {
                  printf("a%d%d = %f\n",i,j,level_func_params->A[i][j]);
                  printf("b%d%d = %f\n",i,j,level_func_params->B[i][j]);
                  printf("c%d%d = %f\n",i,j,level_func_params->C[i][j]);
                  printf("d%d%d = %f\n",i,j,level_func_params->D[i][j]);
                  printf("\n");
             }
        }
    }

    EQN_PARAMS* eqn_params = (EQN_PARAMS*)front->extra1;
    eqn_params->level_func_params = (POINTER)level_func_params;
    level_func_pack->func_params = (POINTER)level_func_params;
    level_func_pack->func = level_cwave_func;
    level_func_pack->wave_type = FIRST_PHYSICS_WAVE_TYPE;
    fclose(infile);
}

static void assign_seeds(

	unsigned short int *seeds)
{
	seeds[0] = 271;
	seeds[1] = 6253;
	seeds[2] = 176;
}

static void Fracr(
        double xk1,
        double xk2,
        double zk1,
        double zk2,
        int n,
        double xk,
        double *frc)
{

        /* frc is a fraction of box (zk1,zk2) x (xk1,xk2) 
           within circle of radius xk 
         */
        double r1,r2,r3,r4,rmax,rmin;
        double f,ff,dxx,dzz,xx,zz,rr;
        int i,j;

        r1 = sqrt(xk1*xk1+zk1*zk1);
        r2 = sqrt(xk2*xk2+zk1*zk1);
        r3 = sqrt(xk1*xk1+zk2*zk2);
        r4 = sqrt(xk2*xk2+zk2*zk2);
        rmax = std::max(std::max(r1,r2),std::max(r3,r4));
        rmin = std::min(std::min(r1,r2),std::max(r3,r4));

        if (xk >= rmax)
        {
            f = 1.0;
        }
        else
        {
            if (xk <= rmin)
            {
                f = 0.0;
            }
            else
            {
                f = 0.0;
                ff = 1.0/(n*n);
                dxx = (xk2-xk1)/n;
                dzz = (zk2-zk1)/n;
                for (i =1; i <= n; i++)
                {
                     for (j = 1; j <= n; j++)
                     {
                          xx = xk1+(i-0.5)*dxx;
                          zz = zk1+(j-0.5)*dzz;
                          rr = sqrt(xx*xx+zz*zz);
                          if (rr <= xk) f += ff;
                     }
                 }
             }
        }

        *frc = f;
}       /* end fracr */

static void pertCoeff(
        double **A,
        double **B,
        double **C,
        double **D,
        int num_modes,
        double *L,
        double *U)
{
        GAUSS_PARAMS *gauss_params; 
	double rescaled_sd,width,lambdaMin,lambdaMax;
        double sdModCoeff; /* Scaling parameter for stand derivation*/
        double r,ra,rb,rc,rd,calculated_sd;
        double ylimit,xmn,zmn,desired_sd;
        double zk1,zk2,xk1,xk2,kMin,kMax,frc1,frc2,frc;
        static unsigned short int seeds[3];
        int i,j;

        width = U[0]-L[0];
        lambdaMin = width/8;
        lambdaMax = width/4;
        rescaled_sd = 0.1*lambdaMin;

	assign_seeds(seeds);
	
	FT_ScalarMemoryAlloc((POINTER*)&gauss_params,sizeof(GAUSS_PARAMS));
	gauss_params->mu = 0.0;
	gauss_params->sigma = 1.0;

        for (j = 0; j <= num_modes; j++)
        {
             for (i = 0; i <= num_modes; i++)
             {
                  zmn = (double)j;
                  xmn = (double)i;
                  zk2 = (zmn+0.5)/width;
                  xk2 = (xmn+0.5)/width;
                  zk1 = (std::max(zmn-0.5,0.0))/width;
                  xk1 = (std::max(xmn-0.5,0.0))/width;
                  kMin = 1.0/lambdaMax;
                  kMax = 1.0/lambdaMin;

                  Fracr(xk1, xk2, zk1, zk2, 10, kMin, &frc1);
                  Fracr(xk1, xk2, zk1, zk2, 10, kMax, &frc2);
                  frc = frc2-frc1;

                  ylimit = 1.0/sqrt(pow(xmn/width,2)+pow(zmn/width,2)+1.0E-50);
                  desired_sd = pow((lambdaMin/ylimit),(0.5*(-1.0)))*sqrt(frc);
                  
		  if (frc >= 1.0E-10)
                  {

                      /* choose ra,rb,rc,rd from unit Gaussian distribution 
                         limit to range -4.0 to +4.0 */
		      r = gauss_box_muller((POINTER)gauss_params,seeds);
                      if (r > 4.0) r = 4.0;
                      if (r < -4.0) r = -4.0;
                      ra = r;

                      r = gauss_box_muller((POINTER)gauss_params,seeds);
		      if (r > 4.0) r = 4.0;
                      if (r < -4.0) r = -4.0;
                      rb = r;

                      r = gauss_box_muller((POINTER)gauss_params,seeds);
		      if (r > 4.0) r = 4.0;
                      if (r < -4.0) r = -4.0;
                      rc = r;

                      r = gauss_box_muller((POINTER)gauss_params,seeds);
		      if (r > 4.0) r = 4.0;
                      if (r < -4.0) r = -4.0;
                      rd = r;

                      A[i][j] = ra*desired_sd;
                      B[i][j] = rb*desired_sd;
                      C[i][j] = rc*desired_sd;
                      D[i][j] = rd*desired_sd;
                  }
                  else
                  {
                      A[i][j] = 0;
                      B[i][j] = 0;
                      C[i][j] = 0;
                      D[i][j] = 0;
                  }
             }
        }

        A[0][0] = 0;
        B[0][0] = 0;
        C[0][0] = 0;
        D[0][0] = 0;

        for(i = 1; i <= num_modes; i++)
        {
            A[i][0] = 1.0/sqrt(2.0)*A[i][0];
            B[i][0] = 0.0;
            C[i][0] = 1.0/sqrt(2.0)*C[i][0];
            D[i][0] = 0.0;
        }

        for (j = 1; j <= num_modes; j++)
        {
             A[0][j] = 1.0/sqrt(2.0)*A[0][j];
             B[0][j] = 1.0/sqrt(2.0)*B[0][j];
             C[0][j] = 0.0;
             D[0][j] = 0.0;
        }

        /* Calculate current SD */
        calculated_sd = 0.0;
        for (i = 0; i <= num_modes; i++)
        {
             for (j = 0; j <= num_modes; j++)
             {
                  if ((i > 0) && (j > 0))
                  {
                      calculated_sd += (1./4.)*(pow(A[i][j],2)+pow(B[i][j],2)+
                                pow(C[i][j],2)+pow(D[i][j],2));
                  }
                  else if ((i > 0) && (j == 0))
                  {
                      calculated_sd += (1./2.)*(pow(A[i][j],2)+pow(C[i][j],2));
                  }
                  else if ((j > 0) && (i == 0))
                  {
                      calculated_sd += (1./2.)*(pow(A[i][j],2)+pow(B[i][j],2));
                  }
                  else
                  {
                      calculated_sd += pow(A[i][j],2)+pow(B[i][j],2)+pow(C[i][j],2)+pow(D[i][j],2);
                  }
             }
        }

        calculated_sd = sqrt(calculated_sd);
        sdModCoeff = rescaled_sd/calculated_sd;

        /* Scale all coefficients to give the desired Standard derivation */
        for (i = 0; i <= num_modes; i++)
        {
             for (j = 0; j <= num_modes; j++)
             {
                  A[i][j] *= sdModCoeff;
                  B[i][j] *= sdModCoeff;
                  C[i][j] *= sdModCoeff;
                  D[i][j] *= sdModCoeff;
             }
        }
} /* end pertCoeff */


void initObliqueIntfc(Front* front,
        LEVEL_FUNC_PACK *level_func_pack)
{
	static FOURIER_POLY *level_func_params;
    FILE *infile = fopen(InName(front),"r");
    int i, num_nodes;
	COMPONENT neg_comp,pos_comp;
	double **node_coords;
	CURVE *oblique;
    EQN_PARAMS* eqn_params = (EQN_PARAMS*)front->extra1;
    PROB_TYPE prob_type = eqn_params->prob_type;
	char test_name[100];

	FT_ScalarMemoryAlloc((POINTER*)&level_func_params,sizeof(FOURIER_POLY));
    level_func_params->dim = front->rect_grid->dim;
	eqn_params->level_func_params = (POINTER)level_func_params;

	CursorAfterString(infile,"Enter number of node points of oblique: ");
	fscanf(infile,"%d",&num_nodes);
	(void) printf("%d\n",num_nodes);
	FT_MatrixMemoryAlloc((POINTER*)&node_coords,num_nodes,MAXD,
				sizeof(double));
	CursorAfterString(infile,"Enter coordinates of node points: ");
	(void) printf("\n");
	for (i = 0; i < num_nodes; ++i)
	{
	    fscanf(infile,"%lf %lf",&node_coords[i][0],&node_coords[i][1]);
	    (void) printf("%f %f\n",node_coords[i][0],node_coords[i][1]);
	}
	neg_comp = SOLID_COMP;
	pos_comp = GAS_COMP2;
	oblique = FT_MakeNodeArrayCurve(front,num_nodes,node_coords,neg_comp,
				pos_comp,NO,0.75,NEUMANN_BOUNDARY);
        fclose(infile);
}

extern void insert_objects(
	Front *front)
{
	int dim = front->rect_grid->dim;
    EQN_PARAMS *eqn_params = (EQN_PARAMS*)front->extra1;
	switch (eqn_params->prob_type)
    {
    	case PROJECTILE:
            init_gun_and_bullet(front);
            break;
	}

}	/* end insert_objects */

static void init_gun_and_bullet(
	Front *front)
{
	char *inname = InName(front);
	FILE *infile = fopen(inname,"r");
	COMPONENT neg_comp,pos_comp;
	int i,dim,num_segs;
	RECT_GRID *rgr = front->rect_grid;
	char string[100];
	int idir;
	double center[MAXD],r_in,r_out,h,gap;
	SURFACE *gun,**projectile;
	static PROJECTILE_PARAMS *proj_params;
	double  (*func)(POINTER,double*);
	POINTER func_params;

	dim = rgr->dim;
	if (dim != 3) return;

	neg_comp = SOLID_COMP;
	pos_comp = GAS_COMP1;
	gap = 0.01*rgr->h[1];
	if (CursorAfterStringOpt(infile,"Enter yes if there is a gun:"))
	{
	    fscanf(infile,"%s",string);
            (void) printf("%s\n",string);
	    if (string[0] == 'Y' || string[0] == 'y')
	    {
		CursorAfterString(infile,"Enter direction of the axis: ");
		fscanf(infile,"%d",&idir);
		(void) printf("%d\n",idir);
		CursorAfterString(infile,"Enter center coordinates: ");
		fscanf(infile,"%lf %lf %lf",center,center+1,center+2);
		(void) printf("%f %f %f\n",center[0],center[1],center[2]);
		CursorAfterString(infile,"Enter inner and outer radii: ");
		fscanf(infile,"%lf %lf",&r_in,&r_out);
		(void) printf("%f %f\n",r_in,r_out);
		CursorAfterString(infile,"Enter height of the cylinder: ");
		fscanf(infile,"%lf",&h);
		(void) printf("%f\n",h);
		r_in += 0.5*gap;
		FT_MakeAnnularCylinderSurf(front,center,r_out,r_in,h,
				neg_comp,pos_comp,idir,NEUMANN_BOUNDARY,
				2,YES,&gun);
	    }
	}

	if (CursorAfterStringOpt(infile,"Enter yes to remove the projectile:"))	
	{
	    fscanf(infile,"%s",string);
            (void) printf("%s\n",string);
            if (string[0] == 'Y' || string[0] == 'y')
		return;
	}

	FT_ScalarMemoryAlloc((POINTER*)&proj_params,sizeof(PROJECTILE_PARAMS));
	proj_params->dim = dim = rgr->dim;
	CursorAfterString(infile,"Enter the center of projectile:");
        for (i = 0; i < dim; ++i)
        {
            fscanf(infile,"%lf",&proj_params->cen[i]);
            (void) printf("%f ",proj_params->cen[i]);
        }
        (void) printf("\n");
        CursorAfterString(infile,"Enter the radius of the projectile:");
        fscanf(infile,"%lf",&proj_params->R);
        (void) printf("%f\n",proj_params->R);
        CursorAfterString(infile,"Enter the head height of the projectile:");
        fscanf(infile,"%lf",&proj_params->r);
        (void) printf("%f\n",proj_params->r);
        CursorAfterString(infile,"Enter the butt height of the projectile:");
        fscanf(infile,"%lf",&proj_params->h);
        (void) printf("%f\n",proj_params->h);
        proj_params->R -= 0.5*gap;
	func_params = (POINTER)proj_params;
        func = projectile_func;
        projectile = (SURFACE**)FT_CreateLevelHyperSurfs(front->rect_grid,
                        front->interf,neg_comp,pos_comp,func,func_params,
                        MOVABLE_BODY_BOUNDARY,&num_segs);

}

