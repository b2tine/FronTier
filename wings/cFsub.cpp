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


#include "cFluid.h"

	/*  Function Declarations */
static void neumann_point_propagate(Front*,POINTER,POINT*,POINT*,
                        HYPER_SURF_ELEMENT*,HYPER_SURF*,double,double*);
static void dirichlet_point_propagate(Front*,POINTER,POINT*,POINT*,
                        HYPER_SURF_ELEMENT*,HYPER_SURF*,double,double*);
static void contact_point_propagate(Front*,POINTER,POINT*,POINT*,
                        HYPER_SURF_ELEMENT*,HYPER_SURF*,double,double*);
static void promptForDirichletBdryState(FILE*,Front*,HYPER_SURF**,int,int);

static double intrp_between(double,double,double,double,double);

static double (*getStateVel[MAXD])(Locstate) =
               {getStateXvel,getStateYvel,getStateZvel};
static double (*getStateMom[MAXD])(Locstate) =
               {getStateXmom,getStateYmom,getStateZmom};
static void set_state_max_speed(Front*,STATE*,double*);
static void get_variable_bdry_params(int,FILE*,POINTER*);
static void cF_variableBoundaryState2d(double*,HYPER_SURF*,Front*,
					POINTER,POINTER);
static void cF_variableBoundaryState3d(double*,HYPER_SURF*,Front*,
					POINTER,POINTER);

double getStateMach(POINTER state)
{
	STATE *fstate = (STATE*)state;
	double c = EosSoundSpeed(fstate);
	double speed;
	int i,dim = fstate->dim;
	speed = 0.0;
	for (i = 0; i < dim; ++i)
	    speed += sqr(fstate->momn[i]/fstate->dens);
	return sqrt(speed)/c;
}	/* end getStateDens */

double getStateDens(POINTER state)
{
	STATE *fstate = (STATE*)state;
	return fstate->dens;
}	/* end getStateDens */

double getStateEngy(POINTER state)
{
	STATE *fstate = (STATE*)state;
	return fstate->engy;
}	/* end getStateEngy */

double getStatePres(POINTER state)
{
	STATE *fstate = (STATE*)state;
	return fstate->pres;
}	/* end getStatePres */

double getStateVort(POINTER state)
{
	STATE *fstate = (STATE*)state;
	return fstate->vort;
}	/* end getStateVort */

double getStateXmom(POINTER state)
{
	STATE *fstate = (STATE*)state;
	return fstate->momn[0];
}	/* end getStateXmom */

double getStateYmom(POINTER state)
{
	STATE *fstate = (STATE*)state;
	return fstate->momn[1];
}	/* end getStateYmom */

double getStateZmom(POINTER state)
{
	STATE *fstate = (STATE*)state;
	return fstate->momn[2];
}	/* end getStateZmom */

double getStateXvel(POINTER state)
{
	STATE *fstate = (STATE*)state;
	return fstate->vel[0];
}	/* end getStateXvel */

double getStateYvel(POINTER state)
{
	STATE *fstate = (STATE*)state;
	return fstate->vel[1];
}	/* end getStateYvel */

double getStateZvel(POINTER state)
{
	STATE *fstate = (STATE*)state;
	return fstate->vel[2];
}	/* end getStateZvel */

void read_dirichlet_bdry_data(
	char *inname,
	Front *front)
{
	char msg[100];
	int i,j,k,nhs,dim = front->rect_grid->dim;
	FILE *infile = fopen(inname,"r");
	HYPER_SURF *hs,**hss;
	INTERFACE *intfc = front->interf;
	int i_hs = 0;

	for (i = 0; i < dim; ++i)
	for (j = 0; j < 2; ++j)
	{
	    if (rect_boundary_type(intfc,i,j) == DIRICHLET_BOUNDARY)
	    {
		hs = FT_RectBoundaryHypSurf(intfc,DIRICHLET_BOUNDARY,i,j);
		if (hs == NULL)
		{
		    printf("ERROR: cannot find Dirichlet boundary"
			   " in dimension %d direction %d\n",i,j);
		    clean_up(ERROR);
		}
		if (j == 0)
		    sprintf(msg,"For lower boundary in %d-th dimension",i);
		else
		    sprintf(msg,"For upper boundary in %d-th dimension",i);
		CursorAfterString(infile,msg);
		(void) printf("\n");
		promptForDirichletBdryState(infile,front,&hs,1,i_hs);
		i_hs++;
	    }
	    else if (rect_boundary_type(intfc,i,j) == MIXED_TYPE_BOUNDARY)
	    {
		hss = FT_MixedBoundaryHypSurfs(intfc,i,j,DIRICHLET_BOUNDARY,
					&nhs);
		printf("Number of Dirichlet boundaries on dir %d side %d: %d\n",
					i,j,nhs);
		if (dim == 2)
		{
		    for (k = 0; k < nhs; ++k)
		    {
			CURVE *c = Curve_of_hs(hss[k]);
		    	(void) printf("Curve %d start and end at: ",k+1);
		    	(void) printf("(%f %f)->(%f %f)\n",
				  Coords(c->start->posn)[0],
				  Coords(c->start->posn)[1],
				  Coords(c->end->posn)[0],
				  Coords(c->end->posn)[1]);
			promptForDirichletBdryState(infile,front,hss+k,1,i_hs);
			i_hs++;
		    }
		}
	    }
	}
	hss = FT_InteriorHypSurfs(intfc,DIRICHLET_BOUNDARY,&nhs);
	if (nhs == 0 || hss == NULL) return;

	sprintf(msg,"For interior Dirichlet boundary:");
	CursorAfterString(infile,msg);
	(void) printf("\n");
	promptForDirichletBdryState(infile,front,hss,nhs,i_hs);
	i_hs++;
}	/* end read_dirichlet_bdry_data */

void cF_variableBoundaryState(
        double          *p0,
        HYPER_SURF      *hs,
        Front           *front,
        POINTER         params,
        POINTER         state)
{
	switch (front->rect_grid->dim)
	{
	case 2:
	    cF_variableBoundaryState2d(p0,hs,front,params,state);
	    return;
	case 3:
	    cF_variableBoundaryState3d(p0,hs,front,params,state);
	    return;
	}
}	/* end cF_variableBoundaryState */

void cF_variableBoundaryState2d(
        double          *p0,
        HYPER_SURF      *hs,
        Front           *front,
        POINTER         params,
        POINTER         state)
{
	VAR_BDRY_PARAMS *bdry_params;
	STATE *newst = (STATE*) state;
	int i=0, dim, nbr_pist;
	double *angles,half_angular_width,*center,jet_duration_time;
	double radius,theta,vec[MAXD];
	boolean within_piston = NO;
	FLOW_THROUGH_PARAMS *ft_params = (FLOW_THROUGH_PARAMS*)params;

	bdry_params = (VAR_BDRY_PARAMS*)boundary_state_function_params(hs);
	jet_duration_time = bdry_params->jet_duration_time;
	if (front->time > jet_duration_time)
	{
	    cF_flowThroughBoundaryState(p0,hs,front,params,state);
	    return;
	}

	dim = bdry_params->dim;
	center = bdry_params->center;
	nbr_pist = bdry_params->number_pistons;
	half_angular_width = bdry_params->half_angular_width;
	angles = bdry_params->angles_pistons;

	radius = 0.0;
	for (i = 0; i < dim; ++i) 
	{
	    vec[i] = p0[i] - center[i];
	    radius += sqr(vec[i]);
	}
	radius = sqrt(radius);
	for (i = 0; i < dim; ++i) 
	    vec[i] /= -radius;
	theta = asin(fabs(p0[1] - center[1])/radius);
	if (p0[0]-center[0] < 0 && p0[1]-center[1] > 0)
            theta = PI - theta;
	else if (p0[0]-center[0] < 0 && p0[1]-center[1] < 0)
            theta = PI + theta;
	else if (p0[0]-center[0] > 0 && p0[1]-center[1] < 0)
            theta = 2*PI - theta;
	for (i = 0; i < nbr_pist; ++i)
	{
	    if (theta > angles[i] - half_angular_width &&
		theta < angles[i] + half_angular_width)
	    {
		within_piston = YES;
	    }
	}
	if (within_piston)
	{
	    POINT *oldp = ft_params->oldp;
	    HYPER_SURF *oldhs = oldp->hs;
	    HYPER_SURF_ELEMENT *oldhse = oldp->hse;
	    STATE *sl,*sr;
	    COMPONENT comp;
	    slsr(oldp,oldhse,oldhs,(POINTER*)&sl,(POINTER*)&sr);
	    if (gas_comp(negative_component(oldhs)))
	    {
	        newst = (STATE*)state;
	        comp = negative_component(oldhs);
	        newst->eos = sl->eos;
	    }
	    else if (gas_comp(positive_component(oldhs)))
	    {
	        newst = (STATE*)state;
	        comp = positive_component(oldhs);
	        newst->eos = sr->eos;
	    }

	    newst->dens = bdry_params->bdry_dens;
	    newst->pres = bdry_params->bdry_pres;
	    for (i = 0; i < dim; ++i)
	    {
		newst->vel[i] = bdry_params->bdry_vel*vec[i];
		newst->momn[i] = (newst->dens)*(newst->vel[i]);
	    }
	    newst->engy = EosEnergy(newst);
	    set_state_max_speed(front,newst,p0);
	}
	else
	{
	    cF_flowThroughBoundaryState(p0,hs,front,params,state);
	}
}	/* end cF_variableBoundaryState2d */

void cF_variableBoundaryState3d(
        double          *p0,
        HYPER_SURF      *hs,
        Front           *front,
        POINTER         params,
        POINTER         state)
{
}	/* end cF_variableBoundaryState3d */

void cF_flowThroughBoundaryState(
        double          *p0,
        HYPER_SURF      *hs,
        Front           *front,
        POINTER         params,
        POINTER         state)
{
	Tan_stencil **tsten;
	Nor_stencil *nsten;
	FLOW_THROUGH_PARAMS *ft_params = (FLOW_THROUGH_PARAMS*)params;
	POINT *oldp = ft_params->oldp;
	COMPONENT comp = ft_params->comp;
	EQN_PARAMS *eqn_params = ft_params->eqn_params;
	static SWEEP *st_stencil;
	static FSWEEP *st_flux;
	double dir[MAXD];
	double u[3];		/* velocity in the sweeping direction */
	double v[3][MAXD];	/* velocity in the orthogonal direction */
	double vort[3];		/* vorticity stencil */
	double pres[3];		/* pressure stencil */
	double dens[3];		/* pressure stencil */
	double f_u;		/* u flux in the sweeping direction */
	double f_v[MAXD];	/* v flux in the orthogonal direction */
	double f_vort;		/* vort flux */
	double f_pres;		/* pressure flux */
	double f_dens;		/* density flux */
	double dn,dt = front->dt;
	STATE *newst = (STATE*)state;
	STATE  *s0,*sl,*sr,**sts;
	static STATE *s1;
	int i,j,dim = front->rect_grid->dim;
	int nrad = 3;
	int size = 2*nrad + 1;
	
	if (debugging("flow_through"))
	    printf("Entering cF_flowThroughBoundaryState()\n");
	if (s1 == NULL)
	{
	    FT_ScalarMemoryAlloc((POINTER*)&s1,sizeof(STATE));
	    FT_ScalarMemoryAlloc((POINTER*)&st_stencil,sizeof(SWEEP));
	    FT_ScalarMemoryAlloc((POINTER*)&st_flux,sizeof(FSWEEP));
	    FT_VectorMemoryAlloc((POINTER*)&st_stencil->dens,size,
					sizeof(double));
	    FT_VectorMemoryAlloc((POINTER*)&st_stencil->engy,size,
					sizeof(double));
	    FT_VectorMemoryAlloc((POINTER*)&st_stencil->pres,size,
					sizeof(double));
	    FT_MatrixMemoryAlloc((POINTER*)&st_stencil->momn,MAXD,size,
					sizeof(double));
	    FT_VectorMemoryAlloc((POINTER*)&st_flux->dens_flux,size,
					sizeof(double));
	    FT_VectorMemoryAlloc((POINTER*)&st_flux->engy_flux,size,
					sizeof(double));
	    FT_MatrixMemoryAlloc((POINTER*)&st_flux->momn_flux,MAXD,size,
					sizeof(double));
	}

	tsten = FrontGetTanStencils(front,oldp,nrad);
	if (tsten != NULL)
	{
	    for (i = 0; i < dim; ++i)
	    	dir[i] = tsten[0]->dir[i];
	    dn = FT_GridSizeInDir(dir,front);

	    if (comp == negative_component(hs))  
	    {
	        sts = (STATE**)tsten[0]->leftst;
		s0 = sts[0];
	    }
	    else 
	    {
	        sts = (STATE**)tsten[0]->rightst;
		s0 = sts[0];
	    }

	    if (debugging("flow_through"))
	    {
	    	(void) printf("Ambient component: %d\n",comp);
	    	(void) printf("hs = %p  oldp->hs = %p\n",
					(POINTER)hs,(POINTER)oldp->hs);
	    	(void) printf("Time step = %f  Tangential grid size = %f\n",
					dt,dn);
	    	(void) printf("Tangential direction: ");
	    	for (j = 0; j < dim; ++j)
		    (void) printf("%f ",tsten[0]->dir[j]);
	    	(void) printf("\n");
	    	(void) printf("Tan_stencil at point p(%f %f)\n",Coords(oldp)[0],
				Coords(oldp)[1]);
	    	(void) printf("Left points:\n");
	    	for (i = 0; i < nrad; ++i)
	    	{
		    for (j = 0; j < dim; ++j)
	    	    	(void) printf("%f ",Coords(tsten[0]->p[-i])[j]);
		    (void) printf("\n");
	    	}
	    	(void) printf("Right points:\n");
	    	for (i = 0; i < nrad; ++i)
	    	{
		    for (j = 0; j < dim; ++j)
	    	    	(void) printf("%f ",Coords(tsten[0]->p[i])[j]);
		    (void) printf("\n");
	    	}
	    }

	    for (j = 0; j < 3; ++j)
	    	u[j] = 0.0;
	    for (j = 0; j < 3; ++j)
	    {
	    	vort[j] = sts[j-1]->vort;
	    	pres[j] = sts[j-1]->pres;
	    	dens[j] = sts[j-1]->dens;
	    	for (i = 0; i < dim; ++i)
	    	{
		    u[j] += sts[j-1]->vel[i]*dir[i];
		    v[j][i] = sts[j-1]->vel[i]*(1.0 - dir[i]);
	    	}
	    }

	    f_u = burger_flux(u[0],u[1],u[2]);
	    for (i = 0; i < dim; ++i)
	    	f_v[i] = linear_flux(u[1],v[0][i],v[1][i],v[2][i]);
	    f_vort = linear_flux(u[1],vort[0],vort[1],vort[2]);
	    f_pres = linear_flux(u[1],pres[0],pres[1],pres[2]);
	    f_dens = linear_flux(u[1],dens[0],dens[1],dens[2]);

	    for (i = 0; i < dim; ++i)
	    	newst->vel[i] = sts[0]->vel[i] - dt/dn*(f_u*dir[i] + f_v[i]) ;
	    newst->vort = sts[0]->vort - dt/dn*f_vort;
	    newst->pres = sts[0]->pres - dt/dn*f_pres;
	    newst->dens = sts[0]->dens - dt/dn*f_dens;
	}
	else
	{
	    slsr(oldp,oldp->hse,oldp->hs,(POINTER*)&sl,(POINTER*)&sr);
	    if (comp == negative_component(hs))  
		s0 = sl;
	    else
		s0 = sr;
	}
	
	nsten = FT_CreateNormalStencil(front,oldp,comp,nrad);
	for (i = 0; i < dim; ++i)
	    dir[i] = nsten->nor[i];
	dn = FT_GridSizeInDir(dir,front);

	if (debugging("flow_through"))
	{
	    printf("Time step = %f  Normal grid size = %f\n",dt,dn);
	    printf("Normal direction: ");
	    for (j = 0; j < dim; ++j)
		printf("%f ",nsten->nor[j]);
	    printf("\n");
	    printf("Nor_stencil at point p(%f %f)\n",Coords(oldp)[0],
				Coords(oldp)[1]);
	    printf("Nor_stencil:\n");
	    for (i = 0; i < nrad; ++i)
	    {
		for (j = 0; j < dim; ++j)
	    	    printf("%f ",nsten->pts[i][j]);
		printf("\n");
	    }
	}

	for (j = 0; j < 3; ++j)
	    u[j] = 0.0;
	for (j = 0; j < 2; ++j)
	{
	    vort[j] = s0->vort;
	    pres[j] = s0->pres;
	    dens[j] = s0->dens;
	    for (i = 0; i < dim; ++i)
	    {
		u[j] += s0->vel[i]*dir[i];
		v[j][i] = s0->vel[i]*(1.0 - dir[i]);
	    }
	}
	for (i = 0; i < dim; ++i)
	{
	    double vtmp;
	    FT_IntrpStateVarAtCoords(front,comp,nsten->pts[1],
			eqn_params->vel[i],getStateVel[i],&vtmp,&s0->vel[i]);
	    s1->vel[i] = vtmp;
	}
	if (dim == 2)
	{
	    FT_IntrpStateVarAtCoords(front,comp,nsten->pts[1],
			eqn_params->vort,getStateVort,&vort[2],&s0->vort);
	    s1->vort = vort[2];
	}
	FT_IntrpStateVarAtCoords(front,comp,nsten->pts[1],eqn_params->pres,
                            getStatePres,&pres[2],&s0->pres);
	FT_IntrpStateVarAtCoords(front,comp,nsten->pts[1],eqn_params->dens,
                            getStateDens,&dens[2],&s0->dens);
	s1->pres = pres[2];
	s1->dens = dens[2];
	for (i = 0; i < dim; ++i)
	{
	    u[2] += s1->vel[i]*dir[i];
	    v[2][i] = s1->vel[i] - s1->vel[i]*dir[i];
	}

	f_u = burger_flux(u[0],u[1],u[2]);
	for (i = 0; i < dim; ++i)
	    f_v[i] = linear_flux(u[1],v[0][i],v[1][i],v[2][i]);
	f_vort = linear_flux(u[1],vort[0],vort[1],vort[2]);
	f_pres = linear_flux(u[1],pres[0],pres[1],pres[2]);
	f_dens = linear_flux(u[1],dens[0],dens[1],dens[2]);

	for (i = 0; i < dim; ++i)
	    newst->vel[i] += - dt/dn*(f_u*dir[i] + f_v[i]) ;
	newst->vort += - dt/dn*f_vort;
	newst->pres += - dt/dn*f_pres;
	newst->dens += - dt/dn*f_dens;
	set_state_max_speed(front,newst,p0);
	if (debugging("flow_through"))
	{
	    printf("flow through boundary state:\n");
	    print_general_vector("Velocity: ",newst->vel,dim,"\n");
	    printf("Pressure: %f\n",newst->pres);
	    printf("Vorticity: %f\n",newst->vort);
	}
}       /* end cF_flowThroughBoundaryState */

void cFluid_point_propagate(
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
	case NEUMANN_BOUNDARY:
	case MOVABLE_BODY_BOUNDARY:
	    return neumann_point_propagate(front,wave,oldp,newp,oldhse,
					oldhs,dt,V);
	case DIRICHLET_BOUNDARY:
	    return dirichlet_point_propagate(front,wave,oldp,newp,oldhse,
					oldhs,dt,V);
	case SUBDOMAIN_BOUNDARY:
	case PASSIVE_BOUNDARY:
            return;
	default:
	    return contact_point_propagate(front,wave,oldp,newp,oldhse,
					oldhs,dt,V);
	}
}       /* cFluid_point_propagate */

static  void neumann_point_propagate(
        Front *front,
        POINTER wave,
        POINT *oldp,
        POINT *newp,
        HYPER_SURF_ELEMENT *oldhse,
        HYPER_SURF         *oldhs,
        double              dt,
        double              *V)
{
	EQN_PARAMS *eqn_params = (EQN_PARAMS*)front->extra1;
        double vel[MAXD],s;
        int i, dim = front->rect_grid->dim;
	double dn,*h = front->rect_grid->h;
	double *m_pres = eqn_params->pres;
	double *m_dens = eqn_params->dens;
	double *m_engy = eqn_params->engy;
	double nor[MAXD],tan[MAXD],p1[MAXD];
	double *p0 = Coords(oldp);
	STATE *oldst,*newst;
	POINTER sl,sr;
	COMPONENT comp;

	FT_GetStatesAtPoint(oldp,oldhse,oldhs,&sl,&sr);
	if (gas_comp(negative_component(oldhs)))
	{
	    comp = negative_component(oldhs);
	    oldst = (STATE*)sl;
	    newst = (STATE*)left_state(newp);
	}
	else if (gas_comp(positive_component(oldhs)))
	{
	    comp = positive_component(oldhs);
	    oldst = (STATE*)sr;
	    newst = (STATE*)right_state(newp);
	}
	FT_NormalAtPoint(oldp,front,nor,comp);

	dn = grid_size_in_direction(nor,h,dim);
	for (i = 0; i < dim; ++i)
	    p1[i] = p0[i] + nor[i]*dn;
	tan[0] = -nor[1]; 	tan[1] = nor[0];

	if (wave_type(oldhs) == MOVABLE_BODY_BOUNDARY)
	{
            double omega_dt,crds_com[MAXD];
            omega_dt = angular_velo(oldhs)*dt;
            for (i = 0; i < dim; ++i)
            {
                vel[i] = center_of_mass_velo(oldhs)[i];
                crds_com[i] = Coords(oldp)[i] +dt*vel[i] - 
			center_of_mass(oldhs)[i];
            }
            vel[0] += -angular_velo(oldhs)*crds_com[1]*cos(omega_dt) -
                     angular_velo(oldhs)*crds_com[0]*sin(omega_dt);
            vel[1] +=  angular_velo(oldhs)*crds_com[0]*cos(omega_dt) -
                     angular_velo(oldhs)*crds_com[1]*sin(omega_dt);
	}
	else
	{
            for (i = 0; i < dim; ++i)
	    	vel[i] = 0.0;
	}
	for (i = 0; i < dim; ++i)
	{
            Coords(newp)[i] = Coords(oldp)[i] + dt*vel[i];
	    newst->vel[i] = vel[i];
            FT_RecordMaxFrontSpeed(i,fabs(vel[i]),NULL,Coords(newp),front);
	}
	FT_IntrpStateVarAtCoords(front,comp,p1,m_dens,
			getStateDens,&newst->dens,&oldst->dens);
	FT_IntrpStateVarAtCoords(front,comp,p1,m_engy,
			getStateEngy,&newst->engy,&oldst->engy);
        newst->eos = oldst->eos;
	newst->pres = EosPressure(newst);
	for (i = 0; i < dim; ++i)
	{
	    newst->vel[i] = vel[i];
	    newst->momn[i] = newst->dens*vel[i];
	}
	s = mag_vector(vel,dim);
	FT_RecordMaxFrontSpeed(dim,s,NULL,Coords(newp),front);
	set_state_max_speed(front,newst,Coords(newp));
}	/* end neumann_point_propagate */

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
	EQN_PARAMS *eqn_params = (EQN_PARAMS*)front->extra1;
        int i, dim = front->rect_grid->dim;
	STATE *sl,*sr,*newst = NULL;
	STATE *bstate;
	FLOW_THROUGH_PARAMS ft_params;
	COMPONENT comp;

	if (debugging("dirichlet_bdry"))
	{
	    printf("Entering dirichlet_point_propagate()\n");
	    print_general_vector("oldp:  ",Coords(oldp),dim,"\n");
	}

	slsr(oldp,oldhse,oldhs,(POINTER*)&sl,(POINTER*)&sr);
	if (gas_comp(negative_component(oldhs)))
	{
	    newst = (STATE*)left_state(newp);
	    comp = negative_component(oldhs);
	    newst->eos = sl->eos;
	}
	else if (gas_comp(positive_component(oldhs)))
	{
	    newst = (STATE*)right_state(newp);
	    comp = positive_component(oldhs);
	    newst->eos = sr->eos;
	}
	if (newst == NULL) return;	// node point

	if (boundary_state(oldhs) != NULL)
	{
	    bstate = (STATE*)boundary_state(oldhs);
	    newst->dens = bstate->dens;
	    newst->engy = bstate->engy;
            for (i = 0; i < dim; ++i)
	    {
	    	newst->vel[i] = bstate->vel[i];
	    	newst->momn[i] = bstate->vel[i]*bstate->dens;
	    }
            newst->pres = bstate->pres;
            newst->vort = 0.0;
	    set_state_max_speed(front,newst,Coords(oldp));

	    if (debugging("dirichlet_bdry"))
	    {
		printf("Preset boundary state:\n");
		print_general_vector("Velocity: ",newst->vel,dim,"\n");
		printf("Density: %f\n",newst->dens);
		printf("Energy: %f\n",newst->engy);
		printf("Pressure: %f\n",newst->pres);
		printf("Vorticity: %f\n",newst->vort);
	    }
	}
	else if (boundary_state_function(oldhs))
	{
	    oldp->hse = oldhse;
	    oldp->hs = oldhs;
	    ft_params.oldp = oldp;
	    ft_params.eqn_params = eqn_params;
	    ft_params.comp = comp;
	    (*boundary_state_function(oldhs))(Coords(oldp),oldhs,front,
			(POINTER)&ft_params,(POINTER)newst);	
	}
	if (debugging("dirichlet_bdry"))
	    printf("Leaving dirichlet_point_propagate()\n");
        return;
}	/* end dirichlet_point_propagate */

static  void contact_point_propagate(
        Front *front,
        POINTER wave,
        POINT *oldp,
        POINT *newp,
        HYPER_SURF_ELEMENT *oldhse,
        HYPER_SURF         *oldhs,
        double              dt,
        double              *V)
{
	EQN_PARAMS *eqn_params = (EQN_PARAMS*)front->extra1;
        int i, dim = front->rect_grid->dim;
	double **m_mom = eqn_params->mom;
	double *m_dens = eqn_params->dens;
	double *m_engy = eqn_params->engy;
	double *p0;
	STATE *oldst,*newst;
	POINTER sl,sr;
	EOS_PARAMS *eos = eqn_params->eos;
	double default_var;

	switch (eqn_params->point_prop_scheme)
	{
	case FIRST_ORDER:
	    first_order_point_propagate(front,wave,oldp,newp,
				oldhse,oldhs,dt,V);
	    break;
	case SECOND_ORDER:
	    second_order_point_propagate(front,wave,oldp,newp,
				oldhse,oldhs,dt,V);
	    break;
	case FOURTH_ORDER:
	    fourth_order_point_propagate(front,wave,oldp,newp,
				oldhse,oldhs,dt,V);
	    break;
	}
	if (debugging("point_propagate"))
	{
	    double hmin,dist;
	    hmin = front->rect_grid->h[0];
	    dist = 0.0;
	    for (i = 0; i < dim; ++i)
		dist += sqr(Coords(newp)[i] - Coords(oldp)[i]);
	    dist = sqrt(dist);
	    if (dist > 0.5*hmin || isnan(dist))
	    {
		printf("WARNING: propagating over half grid size!\n");
		printf("dist/hmin = %f\n",dist/hmin);
		if (dist > hmin || isnan(dist))
		    clean_up(ERROR);
	    }
	}

	FT_GetStatesAtPoint(oldp,oldhse,oldhs,&sl,&sr);
	p0 = Coords(newp);
	newst = (STATE*)left_state(newp);
	oldst = (STATE*)sl;
	*newst = *oldst;
	newst->eos = &eos[negative_component(oldhs)];
	FT_NearestRectGridVarInRange(front,negative_component(oldhs),p0,
			m_dens,2,&default_var);
	FT_IntrpStateVarAtCoords(front,negative_component(oldhs),p0,
		m_dens,getStateDens,&newst->dens,&default_var);
	FT_NearestRectGridVarInRange(front,negative_component(oldhs),p0,
			m_engy,2,&default_var);
	FT_IntrpStateVarAtCoords(front,negative_component(oldhs),p0,
		m_engy,getStateEngy,&newst->engy,&default_var);
	for (i = 0; i < dim; ++i)
	{
	    FT_NearestRectGridVarInRange(front,negative_component(oldhs),p0,
			m_mom[i],2,&default_var);
	    FT_IntrpStateVarAtCoords(front,negative_component(oldhs),p0,
		    m_mom[i],getStateMom[i],&newst->momn[i],&default_var);
	    newst->vel[i] = newst->momn[i]/newst->dens;
	}
	newst->pres = EosPressure(newst);
	set_state_max_speed(front,newst,Coords(oldp));

	newst = (STATE*)right_state(newp);
	oldst = (STATE*)sr;
	*newst = *oldst;
	newst->eos = &eos[positive_component(oldhs)];
	FT_NearestRectGridVarInRange(front,negative_component(oldhs),p0,
			m_dens,2,&default_var);
	FT_IntrpStateVarAtCoords(front,positive_component(oldhs),p0,
		m_dens,getStateDens,&newst->dens,&default_var);
	FT_NearestRectGridVarInRange(front,negative_component(oldhs),p0,
			m_engy,2,&default_var);
	FT_IntrpStateVarAtCoords(front,positive_component(oldhs),p0,
		m_engy,getStateEngy,&newst->engy,&default_var);
	for (i = 0; i < dim; ++i)
	{
	    FT_NearestRectGridVarInRange(front,negative_component(oldhs),p0,
			m_mom[i],2,&default_var);
	    FT_IntrpStateVarAtCoords(front,positive_component(oldhs),p0,
		    m_mom[i],getStateMom[i],&newst->momn[i],&default_var);
	    newst->vel[i] = newst->momn[i]/newst->dens;
	}
	newst->pres = EosPressure(newst);
	set_state_max_speed(front,newst,Coords(oldp));
}	/* end contact_point_propagate */

static void set_state_max_speed(
	Front *front,
	STATE *state,
	double *coords)
{
	int i,dim = front->rect_grid->dim;
	double c,s;
	s = 0.0;
	for (i = 0; i < dim; ++i)
	    s += sqr(state->momn[i]/state->dens);
	s = sqrt(s);
	c = EosSoundSpeed(state);
	set_max_front_speed(dim,s+c,NULL,coords,front);
}	/* end set_state_max_speed */

// Flux of Riemann solution of Burgers equation u_t + uu_x = 0

double burger_flux(	
	double ul,
	double um,
	double ur)
{
	double u_Rl,u_Rr;
	if (ul < um)
	{
	    if (ul > 0.0) u_Rl = ul;
	    else if (um < 0.0) u_Rl = um;
	    else u_Rl = 0.0;
	}
	else
	{
	    if (ul + um > 0.0) u_Rl = ul;
	    else u_Rl = um;
	}

	if (um < ur)
	{
	    if (um > 0.0) u_Rr = um;
	    else if (ur < 0.0) u_Rr = ur;
	    else u_Rr = 0.0;
	}
	else
	{
	    if (um + ur > 0.0) u_Rr = um;
	    else u_Rr = ur;
	}
	return 0.5*(u_Rr*u_Rr - u_Rl*u_Rl);
}	/* end flux */

// Flux of Riemann solution of linear equation u_t + au_x = 0

double linear_flux(	
	double a,
	double ul,
	double um,
	double ur)
{
	if (a > 0.0)
	    return a*(um - ul);
	else
	    return a*(ur - um);
}	/* end net_uwind_flux */

void readFrontStates(
	Front		*front,
	char		*restart_name)
{
	FILE 		*infile;
	EQN_PARAMS 	*eqn_params = (EQN_PARAMS*)front->extra1;
	INTERFACE 	*intfc = front->interf;
        STATE 		*sl,*sr;
        POINT 		*p;
        HYPER_SURF 	*hs;
        HYPER_SURF_ELEMENT *hse;
	STATE 		*lstate,*rstate;
	char 		fname[100];
	int 		i,dim = front->rect_grid->dim;
	int		comp;
	EOS_PARAMS	*eos = eqn_params->eos;

	sprintf(fname,"%s-gas",restart_name);
	infile = fopen(fname,"r");
	
	/* Initialize states at the interface */
        next_output_line_containing_string(infile,"Interface gas states:");
        next_point(intfc,NULL,NULL,NULL);
        while (next_point(intfc,&p,&hse,&hs))
        {
            FT_GetStatesAtPoint(p,hse,hs,(POINTER*)&sl,(POINTER*)&sr);
	    lstate = (STATE*)sl;	rstate = (STATE*)sr;
            fscanf(infile,"%lf %lf",&lstate->dens,&rstate->dens);
            fscanf(infile,"%lf %lf",&lstate->engy,&rstate->engy);
	    for (i = 0; i < dim; ++i)
            	fscanf(infile,"%lf %lf",&lstate->momn[i],&rstate->momn[i]);
	    
	    comp = negative_component(hs);
	    lstate->eos = &eos[comp];
	    lstate->dim = dim;
	    if(gas_comp(comp))
	    	lstate->pres = EosPressure(lstate);
		
	    comp = positive_component(hs);

	    rstate->eos = &eos[comp];
	    rstate->dim = dim;
	    if(gas_comp(comp))
	    	rstate->pres = EosPressure(rstate);
	    lstate->dim = rstate->dim = dim;
        }
	FT_MakeGridIntfc(front);
	fclose(infile);
}

extern void reflectVectorThroughPlane(
	double *vec,
	double *nor,
	double *vec_ref,
	int dim)
{
	int i;
	double vec_nor[MAXD];
	for (i = 0; i < dim; ++i)
	{
	    vec_nor[i] = vec[i]*fabs(nor[i]);
	    vec_ref[i] = vec[i] - 2.0*vec_nor[i];
	}	
}	/* end reflectVectorThroughPlane */

extern boolean reflectNeumannState(
	Front *front,
	HYPER_SURF *hs,
	double *coords,
	COMPONENT comp,
	SWEEP *m_vst,
	STATE *state)
{
	int i,dim = front->rect_grid->dim;
	double coordsref[MAXD],nor[MAXD];
	double momn[MAXD];

	if (!FrontReflectPointViaNeumannBdry(coords,coordsref,nor,comp,
				hs,front))
	{
	    printf("ERROR: in appendGhostBuffer(), cannot reflect point!\n");
	    return NO;
	}
	FT_IntrpStateVarAtCoords(front,comp,coordsref,m_vst->dens,getStateDens,
					&state->dens,NULL);
	FT_IntrpStateVarAtCoords(front,comp,coordsref,m_vst->engy,getStateEngy,
					&state->engy,NULL);
	FT_IntrpStateVarAtCoords(front,comp,coordsref,m_vst->pres,getStatePres,
					&state->pres,NULL);
	for (i = 0; i < dim; ++i)
	{
	    FT_IntrpStateVarAtCoords(front,comp,coordsref,m_vst->momn[i],
				getStateMom[i],&momn[i],NULL);
	}
        reflectVectorThroughPlane(momn,nor,state->momn,dim);
	return YES;
}	/* end reflectNeumannState */	

extern void findGhostState(
	STATE intfc_st,
	STATE inter_st,
	STATE *ghost_st)
{
	double vel[MAXD];
	vel[0] = inter_st.momn[0]/inter_st.dens;
	vel[1] = inter_st.momn[1]/inter_st.dens;
	vel[2] = inter_st.momn[2]/inter_st.dens;


	ghost_st->dens = intfc_st.dens;
	ghost_st->pres = intfc_st.pres;
	ghost_st->momn[0] = intfc_st.dens*vel[0];
	ghost_st->momn[1] = intfc_st.dens*vel[1];
	ghost_st->momn[2] = intfc_st.dens*vel[2];
	ghost_st->engy = EosEnergy(ghost_st);
}	/* end findGhostState */

static void promptForDirichletBdryState(
	FILE *infile,
	Front *front,
	HYPER_SURF **hs,
	int nhs,
	int i_hs)
{
	static STATE *state;
	EQN_PARAMS *eqn_params = (EQN_PARAMS*)front->extra1;
	char s[100];
	COMPONENT comp;
	int i,k,dim = front->rect_grid->dim;
	POINTER func_params;

	FT_ScalarMemoryAlloc((POINTER*)&state,sizeof(STATE));
	state->dim = dim;

	CursorAfterString(infile,"Enter type of Dirichlet boundary:");
	fscanf(infile,"%s",s);
	(void) printf("%s\n",s);
	switch (s[0])
	{
	case 'c':			// Constant state
	case 'C':
	    comp = gas_comp(positive_component(hs[0])) ? 
				positive_component(hs[0]) :
				negative_component(hs[0]);
	    state->eos = &(eqn_params->eos[comp]);
	    CursorAfterString(infile,"Enter velocity:");
	    for (k = 0; k < dim; ++k)
	    {
		fscanf(infile,"%lf",&state->vel[k]);
		(void) printf("%f ",state->vel[k]);
	    }
	    (void) printf("\n");
	    CursorAfterString(infile,"Enter pressure:");
	    fscanf(infile,"%lf",&state->pres);
	    (void) printf("%f\n",state->pres);
	    CursorAfterString(infile,"Enter density:");
	    fscanf(infile,"%lf",&state->dens);
	    (void) printf("%f\n",state->dens);
	    for (k = 0; k < dim; ++k)
                       state->momn[k] = state->dens*state->vel[k];
	    state->engy = EosEnergy(state);
	    FT_InsertDirichletBoundary(front,NULL,NULL,NULL,
			(POINTER)state,hs[0],i_hs);
	    for (i = 1; i < nhs; ++i)
		bstate_index(hs[i]) = bstate_index(hs[0]);
	    break;
	case 'f':			// Flow through state
	case 'F':
	    FT_InsertDirichletBoundary(front,cF_flowThroughBoundaryState,
			"cF_flowThroughBoundaryState",NULL,NULL,hs[0],i_hs);
	    for (i = 1; i < nhs; ++i)
		bstate_index(hs[i]) = bstate_index(hs[0]);
	    break;
	case 'v':			// Flow through state
	case 'V':
	    get_variable_bdry_params(dim,infile,&func_params);
	    FT_InsertDirichletBoundary(front,cF_variableBoundaryState,
			"cF_variableBoundaryState",func_params,NULL,hs[0],i_hs);
	    for (i = 1; i < nhs; ++i)
		bstate_index(hs[i]) = bstate_index(hs[0]);
	    break;
	}
} 	/* end  promptForDirichletBdryState */

static double intrp_between(
	double x1,
	double x2,
	double x,
	double y1,
	double y2)
{
	double y;
	if (x1 == x2) return y1;
	y = y1 + (y2 - y1)/(x2 - x1)*(x - x1);
	return y;
}

static void get_variable_bdry_params(
	int dim,
	FILE *infile,
	POINTER *func_params)
{
	static VAR_BDRY_PARAMS params;
	int nbr_pistons;
	double half_angular_width;
	int i;

	params.dim=dim;
	CursorAfterString(infile,"Enter the center of the circle:");
        for (i = 0; i < dim; ++i)
	{
            fscanf(infile,"%lf",&params.center[i]);
            (void) printf("%f ",params.center[i]);
	}
	(void) printf("\n");

	CursorAfterString(infile,"Enter number of pistons:");
	fscanf(infile,"%d",&nbr_pistons);
	(void) printf("%d\n",nbr_pistons);

	params.number_pistons = nbr_pistons;

	FT_VectorMemoryAlloc((POINTER*)&params.angles_pistons,nbr_pistons+1,
				sizeof(double));
	for (i = 0; i < nbr_pistons+1; ++i)
		params.angles_pistons[i] = 2*PI*i/nbr_pistons;
	params.number_pistons = nbr_pistons + 1;

	CursorAfterString(infile,"Enter angular width of pistons:");
	fscanf(infile,"%lf",&half_angular_width);
	(void) printf("%f\n",half_angular_width);
	params.half_angular_width = half_angular_width;

	CursorAfterString(infile,
		"Enter radial velocity, density and pressure at piston:");
	fscanf(infile,"%lf %lf %lf",&params.bdry_vel,&params.bdry_dens,
					&params.bdry_pres);
	(void) printf("%f %f %f\n",params.bdry_vel,params.bdry_dens,
					params.bdry_pres);
	CursorAfterString(infile,
		"Enter time duration of the piston:");
	fscanf(infile,"%lf",&params.jet_duration_time);
	(void) printf("%f\n",params.jet_duration_time);

	*func_params = (POINTER)&params;
}	/* end get_variable_bdry_params */

extern void restart_set_dirichlet_bdry_function(Front *front)
{
        INTERFACE *intfc = front->interf;
        int i;
        BOUNDARY_STATE  *bstate;
        const char *s;
        for (i = 0; i < num_bstates(intfc); ++i)
        {
            bstate = bstate_list(intfc)[i];
            if (bstate == NULL) continue;
            s = bstate->_boundary_state_function_name;
            if (s == NULL) continue;
            if (strcmp(s,"cF_flowThroughBoundaryState") == 0)
                bstate->_boundary_state_function = cF_flowThroughBoundaryState;
	    else if (strcmp(s,"cF_variableBoundaryState") == 0)
                bstate->_boundary_state_function = cF_variableBoundaryState;
        }
}       /* end restart_set_dirichlet_bdry_function */

extern	void preset_wing_motion(
	Front *fr,
	HYPER_SURF *hs,
	double dt,
	double *force,
	double *torque)
{
	CURVE *c = Curve_of_hs(hs);
	MOTION_PARAMS *motion_params = (MOTION_PARAMS*)c->extra;
	double A = motion_params->amplitude;
	double T = motion_params->period;
	double Omega = motion_params->angular_velo;

	switch (motion_params->wing_motion_type)
	{
	case OSCILLATORY_MOTION:
	    angular_velo(hs) = 0.01*cos(fr->time*2.0*PI/T);
	    break;
	case ROTATIONAL_MOTION:
	    angular_velo(hs) = Omega;
	    break;
	default:
	    (void) printf("Unknown type of wing motion\n");
	    clean_up(ERROR);
	}
}	/* end preset_motion */

/**********************************************************************
 * 	This is a special purpose function to print out the physical  *
 * 	variable (density, pressure, etc) on top and at bottom of a   *
 * 	wing. The assumption is that the wing is a counter-clock-wise *
 * 	curve with starting point from the left tip. Would need       *
 * 	adjustment if configuration is different.                     *
 * ********************************************************************/

extern void recordWingVar(
	Front *front)
{
	INTERFACE *intfc = front->interf;
	POINT *p;
	BOND *b;
	CURVE **c;
	int num_pts,max_npt,max_gindex;
	double *dens,*pres,*x,*tmp_dens,*tmp_pres,*tmp_x;
	int *gindex;
	int wing_tag = 888;
	int i,indx,n;
	char *out_name = OutName(front);
	int step = front->step;
	STATE *sl,*sr,*state;
	static char **sample_color;
	static int count = 0;

        if (sample_color == NULL && pp_mynode() == 0)
        {
            FT_MatrixMemoryAlloc((POINTER*)&sample_color,10,20,sizeof(char));
            sprintf(sample_color[0],"red");
            sprintf(sample_color[1],"blue");
            sprintf(sample_color[2],"green");
            sprintf(sample_color[3],"violet");
            sprintf(sample_color[4],"orange");
            sprintf(sample_color[5],"yellow");
            sprintf(sample_color[6],"pink");
            sprintf(sample_color[7],"cyan");
            sprintf(sample_color[8],"light-gray");
            sprintf(sample_color[9],"dark-gray");
        }

	if (debugging("trace"))
	    (void) printf("Entering recordWingVar()\n");
	num_pts = 0;
	max_gindex = 0;
	intfc_curve_loop(intfc,c)
	{
	    if (wave_type(*c) == MOVABLE_BODY_BOUNDARY)
	    {
		if (Gindex((*c)->first->start) != ERROR_INDEX)
		{
		    max_gindex = Gindex((*c)->first->start);
		    num_pts++;
		}
		curve_bond_loop(*c,b)
		{
		    if (Gindex(b->end) == ERROR_INDEX) continue;
		    num_pts++;
		    if (Gindex(b->end) > max_gindex)
			max_gindex = Gindex(b->end);
		}
	    }
	}
#if defined(__MPI__)
	max_npt = num_pts;
	pp_global_imax(&max_npt,1);
	pp_global_imax(&max_gindex,1);
#endif /* defined(__MPI__) */
	if (pp_mynode() == 0)
	{
	    FT_VectorMemoryAlloc((POINTER*)&x,max_gindex,sizeof(double));
	    FT_VectorMemoryAlloc((POINTER*)&dens,max_gindex,sizeof(double));
	    FT_VectorMemoryAlloc((POINTER*)&pres,max_gindex,sizeof(double));
	    intfc_curve_loop(intfc,c)
	    {
	    	if (wave_type(*c) == MOVABLE_BODY_BOUNDARY)
	    	{
		    b = (*c)->first;
		    p = b->start;
		    if (Gindex(p) != ERROR_INDEX)
		    {
		    	indx = Gindex(p);
		    	FT_GetStatesAtPoint(p,Hyper_surf_element(b),
				Hyper_surf(*c),(POINTER*)&sl,(POINTER*)&sr);
		    	state = gas_comp(negative_component(*c)) ? sl : sr;
		    	dens[indx] = state->dens;
		    	pres[indx] = state->pres;
			x[indx] = Coords(p)[0];
		    }
		    curve_bond_loop(*c,b)
		    {
		    	p = b->end;
		    	indx = Gindex(p);
		    	if (Gindex(p) == ERROR_INDEX) continue;
		    	FT_GetStatesAtPoint(p,Hyper_surf_element(b),
				Hyper_surf(*c),(POINTER*)&sl,(POINTER*)&sr);
		    	state = gas_comp(negative_component(*c)) ? sl : sr;
		    	dens[indx] = state->dens;
		    	pres[indx] = state->pres;
			x[indx] = Coords(p)[0];
		    }
		}
	    }
#if defined(__MPI__)
	    FT_VectorMemoryAlloc((POINTER*)&gindex,max_npt,sizeof(int));
	    FT_VectorMemoryAlloc((POINTER*)&tmp_x,max_npt,sizeof(double));
	    FT_VectorMemoryAlloc((POINTER*)&tmp_dens,max_npt,sizeof(double));
	    FT_VectorMemoryAlloc((POINTER*)&tmp_pres,max_npt,sizeof(double));
	    for (i = 1; i < pp_numnodes(); ++i)
	    {
		pp_recv(wing_tag,i,(POINTER)&num_pts,sizeof(int));
		if (num_pts != 0)
		{
		    pp_recv(wing_tag,i,(POINTER)gindex,num_pts*sizeof(int));
		    pp_recv(wing_tag,i,(POINTER)tmp_x,num_pts*sizeof(double));
		    pp_recv(wing_tag,i,(POINTER)tmp_dens,
					num_pts*sizeof(double));
		    pp_recv(wing_tag,i,(POINTER)tmp_pres,
					num_pts*sizeof(double));
		    for (n = 0; n < num_pts; ++n)
		    {
			indx = gindex[n];
			dens[indx] = tmp_dens[n];
			pres[indx] = tmp_pres[n];
			x[indx] = tmp_x[n];
		    }
		}
	    }
#endif /* defined(__MPI__) */
	}
#if defined(__MPI__)
	else
	{
	    FT_VectorMemoryAlloc((POINTER*)&gindex,num_pts,sizeof(int));
	    FT_VectorMemoryAlloc((POINTER*)&x,num_pts,sizeof(double));
	    FT_VectorMemoryAlloc((POINTER*)&dens,num_pts,sizeof(double));
	    FT_VectorMemoryAlloc((POINTER*)&pres,num_pts,sizeof(double));
	    n = 0;
	    intfc_curve_loop(intfc,c)
	    {
	    	if (wave_type(*c) == MOVABLE_BODY_BOUNDARY)
	    	{
		    b = (*c)->first;
		    p = b->start;
		    if (Gindex(p) != ERROR_INDEX)
		    {
		    	gindex[n] = Gindex(p);
		    	FT_GetStatesAtPoint(p,Hyper_surf_element(b),
				Hyper_surf(*c),(POINTER*)&sl,(POINTER*)&sr);
		    	state = gas_comp(negative_component(*c)) ? sl : sr;
			x[n] = Coords(p)[0];
		    	dens[n] = state->dens;
		    	pres[n] = state->pres;
		    	n++;
		    }
		    curve_bond_loop(*c,b)
		    {
		    	p = b->end;
		    	if (Gindex(p) == ERROR_INDEX) continue;
		    	gindex[n] = Gindex(p);
		    	FT_GetStatesAtPoint(p,Hyper_surf_element(b),
				Hyper_surf(*c),(POINTER*)&sl,(POINTER*)&sr);
		    	state = gas_comp(negative_component(*c)) ? sl : sr;
			x[n] = Coords(p)[0];
		    	dens[n] = state->dens;
		    	pres[n] = state->pres;
		    	n++;
		    }
		}
	    }
	    pp_send(wing_tag,(POINTER)&num_pts,sizeof(int),0);
	    if (num_pts != 0)
	    {
	    	pp_send(wing_tag,(POINTER)gindex,num_pts*sizeof(int),0);
	    	pp_send(wing_tag,(POINTER)x,num_pts*sizeof(double),0);
	    	pp_send(wing_tag,(POINTER)dens,num_pts*sizeof(double),0);
	    	pp_send(wing_tag,(POINTER)pres,num_pts*sizeof(double),0);
	    }
	}
#endif /* defined(__MPI__) */
	if (pp_mynode() == 0)
	{
	    char dirname[200],fname[200];
	    int itop_end = 0;
	    FILE *ofile;

	    for (i = 1; i < max_gindex; ++i)
	    {
		if (x[i-1] < x[i])
		    itop_end = i+1;
	    }

	    sprintf(dirname, "%s/wing-top/wing-top-%s",out_name,
			right_flush(front->step,7));
	    if (!create_directory(dirname,NO))
            {
            	screen("Cannot create directory %s\n",dirname);
            	clean_up(ERROR);
            }
	    /* Print top density */
	    sprintf(fname,"%s/density.xg",dirname);
	    ofile = fopen(fname,"w");
	    fprintf(ofile,"Next\n");
            fprintf(ofile,"color=%s\n",sample_color[count%10]);
            fprintf(ofile,"thickness=1.5\n");
	    for (i = 0; i < itop_end; ++i)
	    {
		fprintf(ofile,"%f  %f\n",x[i],dens[i]);
	    }
	    fclose(ofile);
	    /* Print top pressure */
	    sprintf(fname,"%s/pressure.xg",dirname);
	    ofile = fopen(fname,"w");
	    fprintf(ofile,"Next\n");
            fprintf(ofile,"color=%s\n",sample_color[count%10]);
            fprintf(ofile,"thickness=1.5\n");
	    for (i = 0; i < itop_end; ++i)
	    {
		fprintf(ofile,"%f  %f\n",x[i],pres[i]);
	    }
	    fclose(ofile);

	    sprintf(dirname, "%s/wing-bot/wing-bot-%s",out_name,
			right_flush(front->step,7));
	    if (!create_directory(dirname,NO))
            {
            	screen("Cannot create directory %s\n",dirname);
            	clean_up(ERROR);
            }
	    /* Print bottom density */
	    sprintf(fname,"%s/density.xg",dirname);
	    ofile = fopen(fname,"w");
	    fprintf(ofile,"Next\n");
            fprintf(ofile,"color=%s\n",sample_color[count%10]);
            fprintf(ofile,"thickness=1.5\n");
	    for (i = itop_end-1; i < max_gindex; ++i)
	    {
		fprintf(ofile,"%f  %f\n",x[i],dens[i]);
	    }
	    fclose(ofile);
	    /* Print top bottom */
	    sprintf(fname,"%s/pressure.xg",dirname);
	    ofile = fopen(fname,"w");
	    fprintf(ofile,"Next\n");
            fprintf(ofile,"color=%s\n",sample_color[count%10]);
            fprintf(ofile,"thickness=1.5\n");
	    for (i = itop_end-1; i < max_gindex; ++i)
	    {
		fprintf(ofile,"%f  %f\n",x[i],pres[i]);
	    }
	    fclose(ofile);
	    count++;
	}
	if (debugging("trace"))
	    (void) printf("Leaving recordWingVar()\n");
}	/* end recordWingVar */
