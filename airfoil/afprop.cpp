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

#include <iFluid.h>
#include <airfoil.h>

static double (*getStateVel[3])(POINTER) = {getStateXvel,getStateYvel,
                                        getStateZvel};
static SURFACE *canopy_of_string_node(NODE*);
static void string_curve_propagation(Front*,POINTER,CURVE*,CURVE*,double);
static void mono_curve_propagation(Front*,POINTER,CURVE*,CURVE*,double);
static void gore_curve_propagation(Front*,POINTER,CURVE*,CURVE*,double);
static void passive_curve_propagation(Front*,POINTER,CURVE*,CURVE*,double);
static void gore_point_propagate(Front*,POINTER,POINT*,POINT*,BOND*,double);
static void load_node_propagate(Front*,NODE*,NODE*,double);
static void rg_string_node_propagate(Front*,NODE*,NODE*,double);
static void coating_mono_hyper_surf2d(Front*);
static void coating_mono_hyper_surf3d(Front*);
static	int arrayOfMonoHsbdry(INTERFACE*,CURVE**);
static	int arrayOfGoreHsbdry(INTERFACE*,CURVE**);
static 	int getGoreNodes(INTERFACE*,NODE**);

extern void elastic_point_propagate(
        Front *front,
        POINTER wave,
        POINT *oldp,
        POINT *newp,
        HYPER_SURF_ELEMENT *oldhse,
        HYPER_SURF         *oldhs,
        double              dt,
        double              *V)
{
	STATE *newsl,*newsr;
	STATE *sl,*sr;
	IF_PARAMS *iFparams = (IF_PARAMS*)front->extra1;
	AF_PARAMS *af_params = (AF_PARAMS*)front->extra2;
        IF_FIELD *field = iFparams->field;
	int i, dim = front->rect_grid->dim;
	double *vort = field->vort;
	double **vel = field->vel;
	double *pres = field->pres;
	COMPONENT base_comp = positive_component(oldhs);
	double pp[MAXD],pm[MAXD],nor[MAXD],h;
	double area_dens = af_params->area_dens;
	double left_nor_speed,right_nor_speed;
	double dv[MAXD];

	if (af_params->no_fluid)
	{
	    fourth_order_point_propagate(front,wave,oldp,newp,oldhse,
				oldhs,dt,V);
	    ft_assign(left_state(newp),left_state(oldp),front->sizest);
	    ft_assign(right_state(newp),right_state(oldp),front->sizest);
	    return;
	}

	//FT_GetStatesAtPoint(oldp,oldhse,oldhs,(POINTER*)&sl,(POINTER*)&sr);
	sl = (STATE*)left_state(oldp);
	sr = (STATE*)right_state(oldp);
	newsl = (STATE*)left_state(newp);
	newsr = (STATE*)right_state(newp);

	FT_NormalAtPoint(oldp,front,nor,NO_COMP);
	h = FT_GridSizeInDir(nor,front);
	for (i = 0; i < dim; ++i)
	{
	    pm[i] = Coords(oldp)[i] - h*nor[i];
	    pp[i] = Coords(oldp)[i] + h*nor[i];
	}

	if (dim == 2 && wave_type(oldhs) == ELASTIC_STRING)
	{
            FT_IntrpStateVarAtCoords(front,base_comp,Coords(oldp),pres,
			getStatePres,&newsl->pres,&sl->pres);
            FT_IntrpStateVarAtCoords(front,base_comp,Coords(oldp),pres,
			getStatePres,&newsr->pres,&sr->pres);
	}
	else
	{
            FT_IntrpStateVarAtCoords(front,base_comp-1,pm,pres,
			getStatePres,&newsl->pres,&sl->pres);
            FT_IntrpStateVarAtCoords(front,base_comp+1,pp,pres,
			getStatePres,&newsr->pres,&sr->pres);
	}
	/* Impulse is incremented by the fluid pressure force */
	for (i = 0; i < dim; ++i)
	{
	    dv[i] = 0.0;

	    if (debugging("rigid_canopy"))
		dv[i] = 0.0;
	    else if (front->step > 5)
		dv[i] = (sl->pres - sr->pres)*nor[i]/area_dens;
	    newsr->fluid_accel[i] = newsl->fluid_accel[i] = dv[i];
	    newsr->other_accel[i] = newsl->other_accel[i] = 0.0;
	    newsr->impulse[i] = newsl->impulse[i] = sl->impulse[i];
	    newsr->vel[i] = newsl->vel[i] = sl->vel[i];
	}

	/* Interpolating vorticity for the hyper surface point */
	if (dim == 2)
	{
	    if (wave_type(oldhs) == ELASTIC_STRING)
	    {
	        FT_IntrpStateVarAtCoords(front,base_comp,Coords(oldp),pres,
                        getStateVort,&newsl->vort,&sl->vort);
                FT_IntrpStateVarAtCoords(front,base_comp,Coords(oldp),pres,
                        getStateVort,&newsr->vort,&sr->vort);
	        for (i = 0; i < dim; ++i)
	        {
	            newsr->impulse[i] = newsl->impulse[i] = sl->impulse[i];
		    FT_IntrpStateVarAtCoords(front,base_comp,Coords(oldp),
			vel[i],getStateVel[i],&newsl->vel[i],&sl->vel[i]);
		    FT_IntrpStateVarAtCoords(front,base_comp,Coords(oldp),
			vel[i],getStateVel[i],&newsr->vel[i],&sr->vel[i]);
	        }
	    }
	    else
            {
                FT_IntrpStateVarAtCoords(front,base_comp-1,pm,vort,
				getStateVort,&newsl->vort,&sl->vort);
                FT_IntrpStateVarAtCoords(front,base_comp+1,pp,vort,
				getStateVort,&newsr->vort,&sr->vort);
            }
	}
}       /* elastic_point_propagate */

/*	Given string node, the function finds the corresponding
*	canopy surface.
*/

static SURFACE *canopy_of_string_node(NODE *n)
{
	SURFACE *canopy,**s;
	CURVE *c,**curves;
	int i,nc;
	boolean canopy_found = NO;

	canopy = NULL;
	nc = I_NumOfNodeCurves(n);
	FT_VectorMemoryAlloc((POINTER*)&curves,nc,sizeof(CURVE*));
	I_ArrayOfNodeCurves(n,curves);

	for (i = 0; i < nc; ++i)
	{
	    c = curves[i];
	    for (s = c->pos_surfaces; s && *s; ++s)
	    {
		if (wave_type(*s) == ELASTIC_BOUNDARY)
		{
		    canopy_found = YES;
		    canopy = *s;
		    break;
		}
	    }
	   if (canopy_found) break;
	    for (s = c->neg_surfaces; s && *s; ++s)
	    {
		if (wave_type(*s) == ELASTIC_BOUNDARY)
		{
		    canopy_found = YES;
		    canopy = *s;
		    break;
		}
	    }
	}
	FT_FreeThese(1,curves);
	return (canopy_found == YES) ? canopy : NULL;
}	/* end canopy_of_string_node */

extern void airfoil_point_propagate(
        Front *front,
        POINTER wave,
        POINT *oldp,
        POINT *newp,
        HYPER_SURF_ELEMENT *oldhse,
        HYPER_SURF         *oldhs,
        double              dt,
        double              *V)
{
        if (wave_type(oldhs) == ELASTIC_BOUNDARY ||
	    wave_type(oldhs) == ELASTIC_STRING)
            return elastic_point_propagate(front,wave,oldp,newp,oldhse,oldhs,
                                        dt,V);
        else
            return ifluid_point_propagate(front,wave,oldp,newp,oldhse,oldhs,
                                        dt,V);
}       /* airfoil_point_propagate */

extern void airfoil_curve_propagate(
        Front *front,
        POINTER wave,
	CURVE *oldc,
	CURVE *newc,
        double dt)
{
	int dim = front->rect_grid->dim;

	if (dim == 3)
	{
	    switch (hsbdry_type(oldc))
	    {
		case STRING_HSBDRY:
	    	    return string_curve_propagation(front,wave,oldc,newc,dt);
		case MONO_COMP_HSBDRY:
	    	    return mono_curve_propagation(front,wave,oldc,newc,dt);
		case GORE_HSBDRY:
	    	    return gore_curve_propagation(front,wave,oldc,newc,dt);
		case PASSIVE_HSBDRY:
		    return passive_curve_propagation(front,wave,oldc,newc,dt);
		default:
	    	    return;
	    }
	}
	else if (dim == 2)
	{
	    if (wave_type(oldc) == ELASTIC_BOUNDARY)
		string_curve_propagation(front,wave,oldc,newc,dt);
	}
}	/* end airfoil_curve_propagate */

static void string_curve_propagation(
        Front *front,
        POINTER wave,
	CURVE *oldc,
	CURVE *newc,
        double dt)
{
	BOND *oldb,*newb;
	POINT *oldp,*newp;

	if (!is_load_node(oldc->start))
	{
	    oldp = oldc->start->posn;
	    newp = newc->start->posn;
	    ft_assign(left_state(newp),left_state(oldp),front->sizest);
	    ft_assign(right_state(newp),right_state(oldp),front->sizest);
	}

	if (!is_load_node(oldc->end))
	{
	    oldp = oldc->end->posn;
	    newp = newc->end->posn;
	    ft_assign(left_state(newp),left_state(oldp),front->sizest);
	    ft_assign(right_state(newp),right_state(oldp),front->sizest);
	}

	for (oldb = oldc->first, newb = newc->first; oldb != oldc->last;
		oldb = oldb->next, newb = newb->next)
	{
	    oldp = oldb->end;
	    newp = newb->end;
	    ft_assign(left_state(newp),left_state(oldp),front->sizest);
	    ft_assign(right_state(newp),right_state(oldp),front->sizest);
	}
}	/* end string_curve_propagation */

static void gore_curve_propagation(
        Front *front,
        POINTER wave,
	CURVE *oldc,
	CURVE *newc,
        double dt)
{
	BOND *oldb,*newb;
	POINT *oldp,*newp;

	if (debugging("interact_curve"))
	{
	    (void) printf("Entering gore_curve_propagation()\n");
	}
	oldp = oldc->start->posn;
	newp = newc->start->posn;
	ft_assign(left_state(newp),left_state(oldp),front->sizest);
	ft_assign(right_state(newp),right_state(oldp),front->sizest);

	oldp = oldc->end->posn;
	newp = newc->end->posn;
	ft_assign(left_state(newp),left_state(oldp),front->sizest);
	ft_assign(right_state(newp),right_state(oldp),front->sizest);

	for (oldb = oldc->first, newb = newc->first; oldb != oldc->last;
		oldb = oldb->next, newb = newb->next)
	{
	    oldp = oldb->end;
	    newp = newb->end;
	    gore_point_propagate(front,wave,oldp,newp,oldb,dt);
	}
	if (debugging("interact_curve"))
	{
	    (void) printf("Leaving gore_curve_propagation()\n");
	}
}	/* end gore_curve_propagation */

static void gore_point_propagate(
	Front *front,
        POINTER wave,
	POINT *oldp,
	POINT *newp,
	BOND *oldb,
	double dt)
{
	BOND_TRI **btris;
	HYPER_SURF_ELEMENT *oldhse;
	HYPER_SURF         *oldhs;
	STATE *sl,*sr,*newsl,*newsr;
	double mag_nor,branch_nor[MAXD],nor[MAXD];
	double pm[MAXD],pp[MAXD],h;
	IF_PARAMS *iFparams = (IF_PARAMS*)front->extra1;
	AF_PARAMS *af_params = (AF_PARAMS*)front->extra2;
	IF_FIELD *field = iFparams->field;
	double **vel = field->vel;
        double *pres = field->pres;
	double area_dens = af_params->area_dens;
	double left_nor_speed,right_nor_speed,dv;
	COMPONENT base_comp;
	double V[MAXD];
	int i;

	if (af_params->no_fluid)
	{
	    for (btris = Btris(oldb); btris && *btris; ++btris)
	    {
	    	oldhse = Hyper_surf_element((*btris)->tri);
	    	oldhs = Hyper_surf((*btris)->surface);
	    }
	    fourth_order_point_propagate(front,wave,oldp,newp,oldhse,
				oldhs,dt,V);
	    ft_assign(left_state(newp),left_state(oldp),front->sizest);
	    ft_assign(right_state(newp),right_state(oldp),front->sizest);
	    return;
	}
	sl = (STATE*)left_state(oldp);		
	sr = (STATE*)right_state(oldp);
	newsl = (STATE*)left_state(newp);	
	newsr = (STATE*)right_state(newp);

	for (i = 0; i < 3; ++i) nor[i] = 0.0;
	for (btris = Btris(oldb); btris && *btris; ++btris)
	{
	    oldp->hse = oldhse = Hyper_surf_element((*btris)->tri);
	    oldp->hs = oldhs = Hyper_surf((*btris)->surface);
	    FT_NormalAtPoint(oldp,front,branch_nor,NO_COMP);
	    base_comp = positive_component(oldhs);
	    for (i = 0; i < 3; ++i) nor[i] += branch_nor[i];
	}
	mag_nor = Mag3d(nor);
	for (i = 0; i < 3; ++i) nor[i] /= mag_nor;
	h = FT_GridSizeInDir(nor,front);
	for (i = 0; i < 3; ++i)
	{
	    pm[i] = Coords(oldp)[i] - h*nor[i];
            pp[i] = Coords(oldp)[i] + h*nor[i];
	}
	FT_IntrpStateVarAtCoords(front,base_comp-1,pm,pres,
                        getStatePres,&newsl->pres,&sl->pres);
        FT_IntrpStateVarAtCoords(front,base_comp+1,pp,pres,
                        getStatePres,&newsr->pres,&sr->pres);
	for (i = 0; i < 3; ++i)
        {
            FT_IntrpStateVarAtCoords(front,base_comp-1,pm,vel[i],
                        getStateVel[i],&newsl->vel[i],&sl->vel[i]);
            FT_IntrpStateVarAtCoords(front,base_comp+1,pp,vel[i],
                        getStateVel[i],&newsr->vel[i],&sr->vel[i]);
        }
        left_nor_speed = Dot3d(newsl->vel,nor);
        right_nor_speed = Dot3d(newsr->vel,nor);
	for (i = 0; i < 3; ++i)
        {
            newsl->vel[i] -= left_nor_speed*nor[i];
            newsr->vel[i] -= right_nor_speed*nor[i];
        }
	/* Impulse is incremented by the fluid pressure force */
        for (i = 0; i < 3; ++i)
        {
	    dv = 0.0;

	    if (front->step > 5)
		dv = (sl->pres - sr->pres)*nor[i]/area_dens;
	    if (debugging("rigid_canopy"))
	    	dv = 0.0;
	    newsr->fluid_accel[i] = newsl->fluid_accel[i] = dv;
	    newsr->other_accel[i] = newsl->other_accel[i] = 0.0;
	    newsr->impulse[i] = newsl->impulse[i] = sl->impulse[i];
	}
}	/* end gore_point_propagate */

extern int numOfMonoHsbdry(
	INTERFACE *intfc)
{
	CURVE **c;
	int nc = 0;
	intfc_curve_loop(intfc,c)
	{
	    if (hsbdry_type(*c) == MONO_COMP_HSBDRY) nc++;
	} 
	return nc;
}	/* end numOfMonoBdry */

extern int numOfGoreHsbdry(
	INTERFACE *intfc)
{
	CURVE **c;
	int nc = 0;
	intfc_curve_loop(intfc,c)
	{
	    if (hsbdry_type(*c) == GORE_HSBDRY) nc++;
	} 
	return nc;
}	/* end numOfMonoBdry */

static int arrayOfMonoHsbdry(
	INTERFACE *intfc,
	CURVE **mono_curves)
{
	CURVE **c;
	int nc = 0;
	intfc_curve_loop(intfc,c)
	{
	    if (hsbdry_type(*c) == MONO_COMP_HSBDRY) 
	    {
		mono_curves[nc] = *c;
		nc++;
	    }
	} 
	return nc;
}	/* end arrayOfMonoBdry */

static int arrayOfGoreHsbdry(
	INTERFACE *intfc,
	CURVE **gore_curves)
{
	CURVE **c;
	int nc = 0;
	intfc_curve_loop(intfc,c)
	{
	    if (hsbdry_type(*c) == GORE_HSBDRY) 
	    {
		gore_curves[nc] = *c;
		nc++;
	    }
	} 
	return nc;
}	/* end arrayOfGoreBdry */

extern int numOfGoreNodes(
	INTERFACE *intfc)
{
	NODE **n;
	CURVE **c;
	int num_gore_nodes = 0;
	AF_NODE_EXTRA *extra;
	boolean is_string_node;

	intfc_node_loop(intfc,n)
	{
	    if ((*n)->extra == NULL)
		continue;
	    is_string_node = NO;
	    for (c = (*n)->in_curves; c && *c; ++c)
		if (hsbdry_type(*c) == STRING_HSBDRY)
		    is_string_node = YES;
	    for (c = (*n)->out_curves; c && *c; ++c)
		if (hsbdry_type(*c) == STRING_HSBDRY)
		    is_string_node = YES;
	    if (is_string_node) continue;
	    extra = (AF_NODE_EXTRA*)(*n)->extra;
	    if (extra->af_node_type == GORE_NODE)
		num_gore_nodes++;
	}
	return num_gore_nodes;
}	/* numOfGoreNodes */

extern boolean is_bdry_node(
	NODE *node)
{
	CURVE **c;
	for (c = node->in_curves; c && *c; ++c)
	{
	    if (hsbdry_type(*c) == NEUMANN_HSBDRY ||
		hsbdry_type(*c) == DIRICHLET_HSBDRY ||
		hsbdry_type(*c) == SUBDOMAIN_HSBDRY) 
	    {
		return YES;
	    }
	} 
	for (c = node->out_curves; c && *c; ++c)
	{
	    if (hsbdry_type(*c) == NEUMANN_HSBDRY ||
		hsbdry_type(*c) == DIRICHLET_HSBDRY ||
		hsbdry_type(*c) == SUBDOMAIN_HSBDRY) 
	    {
		return YES;
	    }
	} 
	return NO;
}	/* is_bdry_node */

extern boolean is_gore_node(
	NODE *node)
{
	CURVE **c;
	AF_NODE_EXTRA *extra;

	if (node->extra == NULL)
	    return NO;
	for (c = node->in_curves; c && *c; ++c)
	    if (hsbdry_type(*c) == STRING_HSBDRY)
		return NO;
	for (c = node->out_curves; c && *c; ++c)
	    if (hsbdry_type(*c) == STRING_HSBDRY)
		return NO;
	extra = (AF_NODE_EXTRA*)(node)->extra;
	if (extra->af_node_type == GORE_NODE)
	    return YES;
	else 
	    return NO;
}	/* end is_gore_node */

extern boolean is_load_node(NODE *n)
{
        AF_NODE_EXTRA *af_node_extra;
        if (n->extra == NULL) return NO;
        af_node_extra = (AF_NODE_EXTRA*)n->extra;
        if (af_node_extra->af_node_type == LOAD_NODE) return YES;
        return NO;
}       /* end is_load_node */

extern boolean is_rg_string_node(NODE *n)
{
        AF_NODE_EXTRA *af_node_extra;
        if (n->extra == NULL) return NO;
        af_node_extra = (AF_NODE_EXTRA*)n->extra;
        if (af_node_extra->af_node_type == RG_STRING_NODE) return YES;
        return NO;
}       /* end is_rg_string_node */

static int getGoreNodes(
	INTERFACE *intfc,
	NODE **gore_nodes)
{
	NODE **n;
	int num_nodes = 0;

	intfc_node_loop(intfc,n)
	{
	    if (is_gore_node(*n))
		gore_nodes[num_nodes++] = *n;
	}
	return num_nodes;
}	/* getGoreNodes */

extern boolean goreInIntfc(
	INTERFACE *intfc)
{
	NODE **n;

	for (n = intfc->nodes; n && *n; ++n)
	{
	    if (is_gore_node(*n))
		return YES;
	}
	return NO;
}	/* end goreInIntfc */

static void mono_curve_propagation(
        Front *front,
        POINTER wave,
	CURVE *oldc,
	CURVE *newc,
        double dt)
{
	BOND *oldb,*newb;
	POINT *oldp,*newp;
	double V[MAXD];
	BOND_TRI **btris;
	HYPER_SURF_ELEMENT *oldhse;
	HYPER_SURF         *oldhs;

	if (debugging("interact_curve"))
	{
	    (void) printf("Entering mono_curve_propagation()\n");
	}

	/*
	oldb = oldc->first;
	newb = newc->first;
	oldp = oldb->start;
	newp = newb->start;
	for (btris = Btris(oldb); btris && *btris; ++btris)
	{
	    oldp->hse = oldhse = Hyper_surf_element((*btris)->tri);
	    oldp->hs = oldhs = Hyper_surf((*btris)->surface);
	    elastic_point_propagate(front,wave,oldp,newp,oldhse,oldhs,dt,V);
	}
	*/
	for (oldb = oldc->first, newb = newc->first; oldb != NULL;
		oldb = oldb->next, newb = newb->next)
	{
	    oldp = oldb->end;
	    newp = newb->end;
	    for (btris = Btris(oldb); btris && *btris; ++btris)
	    {
	    	oldp->hse = oldhse = Hyper_surf_element((*btris)->tri);
	    	oldp->hs = oldhs = Hyper_surf((*btris)->surface);
		elastic_point_propagate(front,wave,oldp,newp,oldhse,oldhs,dt,V);
	    }
	}
	if (debugging("interact_curve"))
	{
	    (void) printf("Leaving mono_curve_propagation()\n");
	}
}	/* end mono_curve_propagation */

extern double springCharTimeStep(
	Front *fr)
{
	AF_PARAMS *af_params = (AF_PARAMS*)fr->extra2;
	double dt_tol;
	dt_tol = sqrt((af_params->m_s)/(af_params->ks));
        if (af_params->m_l != 0.0 &&
            dt_tol > sqrt((af_params->m_l)/(af_params->kl)))
            dt_tol = sqrt((af_params->m_l)/(af_params->kl));
        if (af_params->m_g != 0.0 &&
            dt_tol > sqrt((af_params->m_g)/(af_params->kg)))
            dt_tol = sqrt((af_params->m_g)/(af_params->kg));
	return dt_tol;
}	/* end springCharTimeStep */

extern void coating_mono_hyper_surf(
	Front *front)
{
	int dim = front->rect_grid->dim;
	switch (dim)
	{
	case 2:
	    coating_mono_hyper_surf2d(front);
	    return;
	case 3:
	    coating_mono_hyper_surf3d(front);
	    return;
	}
}	/* end coating_mono_hyper_surf */

static void coating_mono_hyper_surf2d(
	Front *front)
{
	INTERFACE *grid_intfc = front->grid_intfc;
	RECT_GRID *top_grid = &topological_grid(grid_intfc);
	struct Table *T = table_of_interface(grid_intfc);

    COMPONENT *top_comp = T->components;
    COMPONENT          comp;
    INTERFACE          *intfc = front->interf;

    double 		   *L = top_grid->L;
	double 		   *h = top_grid->h;
	double             coords[MAXD];
    double             t[MAXD],p[MAXD],vec[MAXD];
	const double 	   *nor;

    CURVE **c,*immersed_curve;
    HYPER_SURF_ELEMENT *hse;
    HYPER_SURF         *hs;
	BOND *b;
	COMPONENT base_comp;

    int i,index,nb,index_nb,*top_gmax = top_grid->gmax;
	int dim = top_grid->dim;
	int icoords[MAXD],icn[MAXD],smin[MAXD],smax[MAXD];
	GRID_DIRECTION dir[4] = {WEST,EAST,SOUTH,NORTH};

	if (debugging("trace"))
	    (void) printf("Entering coating_mono_hyper_surf2d()\n");

	immersed_curve = NULL;
	for (c = grid_intfc->curves; c && *c; ++c)
	{
	    if (wave_type(*c) == ELASTIC_BOUNDARY)
	    {
            immersed_curve = *c;
            comp = base_comp = negative_component(*c);
	    	hs = Hyper_surf(immersed_curve);
		    break;
	    }
	}
	if (immersed_curve == NULL)
	    return;

	for (icoords[0] = 1; icoords[0] < top_gmax[0]; ++icoords[0])
	for (icoords[1] = 1; icoords[1] < top_gmax[1]; ++icoords[1])
	{
	    index = d_index(icoords,top_gmax,dim);
	    for (i = 0; i < dim; ++i)
            coords[i] = L[i] + icoords[i]*h[i];

	    if (nearest_interface_point_within_range(coords,comp,grid_intfc,
			NO_BOUNDARIES,hs,p,t,&hse,&hs,3))
	    {
		    if (wave_type(hs) != ELASTIC_BOUNDARY) 
                continue;

            b = Bond_of_hse(hse);
	    	t[1] = Coords(b->start)[0] - Coords(b->end)[0]; // t is normal
	    	t[0] = Coords(b->end)[1] - Coords(b->start)[1];
	    	
            for (i = 0; i < dim; ++i)
                vec[i] = coords[i] - p[i];
	
            if (scalar_product(vec,t,dim) > 0.0)
	    	    top_comp[index] = base_comp + 1;
	    	else
		        top_comp[index] = base_comp - 1;
	    }
	}
	negative_component(immersed_curve) = base_comp - 1;
	positive_component(immersed_curve) = base_comp + 1;
	if (debugging("coat_comp"))
	{
	    for (icoords[1] = 1; icoords[1] < top_gmax[1]; ++icoords[1])
	    {
	    	for (icoords[0] = 1; icoords[0] < top_gmax[0]; ++icoords[0])
	    	{
		    index = d_index(icoords,top_gmax,dim);
		    (void) printf("%d",top_comp[index]);
	    	}
	    	(void) printf("\n");
	    }
	}
	if (debugging("trace"))
	    (void) printf("Leaving coating_mono_hyper_surf2d()\n");
}	/* end coating_mono_hyper_surf2d */

static void coating_mono_hyper_surf3d(
	Front *front)
{
	INTERFACE *grid_intfc = front->grid_intfc;
	RECT_GRID *top_grid = &topological_grid(grid_intfc);
	struct Table *T = table_of_interface(grid_intfc);
	COMPONENT *top_comp = T->components;
        COMPONENT          comp;
        INTERFACE          *intfc = front->interf;
	double 		   *L = top_grid->L;
	double 		   *h = top_grid->h;
	double             coords[MAXD];
        double             t[MAXD],p[MAXD],vec[MAXD];
	const double 	   *nor;
	SURFACE **s,*immersed_surf;
        HYPER_SURF_ELEMENT *hse;
        HYPER_SURF         *hs;
	COMPONENT base_comp;
	int i,index,nb,index_nb,*top_gmax = top_grid->gmax;
	int dim = top_grid->dim;
	int icoords[MAXD],icn[MAXD],smin[MAXD],smax[MAXD];
	GRID_DIRECTION dir[6] = {WEST,EAST,SOUTH,NORTH,LOWER,UPPER};

	if (debugging("trace"))
	    (void) printf("Entering coating_mono_hyper_surf3d()\n");
	immersed_surf = NULL;
	for (s = grid_intfc->surfaces; s && *s; ++s)
	{
	    if (wave_type(*s) == ELASTIC_BOUNDARY)
	    {
		immersed_surf = *s;
		comp = base_comp = negative_component(*s);
		break;
	    }
	}
	if (immersed_surf == NULL)
	    return;

	for (icoords[0] = 1; icoords[0] < top_gmax[0]; ++icoords[0])
	for (icoords[1] = 1; icoords[1] < top_gmax[1]; ++icoords[1])
	for (icoords[2] = 1; icoords[2] < top_gmax[2]; ++icoords[2])
	{
	    index = d_index(icoords,top_gmax,dim);
	    for (i = 0; i < dim; ++i)
		coords[i] = L[i] + icoords[i]*h[i];
	    if (nearest_interface_point_within_range(coords,comp,grid_intfc,
			NO_BOUNDARIES,NULL,p,t,&hse,&hs,3))
	    {
		if (wave_type(hs) != ELASTIC_BOUNDARY) continue;
	    	nor = Tri_normal(Tri_of_hse(hse));
	    	for (i = 0; i < dim; ++i)
		    vec[i] = coords[i] - p[i];
	    	if (scalar_product(vec,nor,dim) > 0.0)
		    top_comp[index] = base_comp + 1;
	    	else
		    top_comp[index] = base_comp - 1;
	    }
	}
	for (s = grid_intfc->surfaces; s && *s; ++s)
	{
	    if (wave_type(*s) == ELASTIC_BOUNDARY)
	    {
		if (base_comp == negative_component(*s))
		{
		    negative_component(*s) = base_comp - 1;
		    positive_component(*s) = base_comp + 1;
		}
	    }
	}
	if (debugging("coat_comp"))
	{
	    icoords[0] = top_gmax[0]/2;
	    for (icoords[2] = 0; icoords[2] <= top_gmax[2]; ++icoords[2])
	    {
	    	for (icoords[1] = 0; icoords[1] <= top_gmax[1]; ++icoords[1])
	    	{
		    index = d_index(icoords,top_gmax,dim);
		    printf("%d",top_comp[index]);
	    	}
	    	printf("\n");
	    }
	}
	if (debugging("immersed_surf") && front->step%1 == 0)
	{
	    int icrd_nb[MAXD],index_nb,n;
	    POINTER l_state,u_state;
	    double crx_coords[MAXD],crx_nb[MAXD];
	    static double *pl,*pu,*vz,*x;
	    FILE *pfile;
	    char pname[200];

	    n = 0;
	    if (pu == NULL)
	    {
	    	FT_VectorMemoryAlloc((POINTER*)&pu,top_gmax[1],sizeof(double));
	    	FT_VectorMemoryAlloc((POINTER*)&pl,top_gmax[1],sizeof(double));
	    	FT_VectorMemoryAlloc((POINTER*)&vz,top_gmax[1],sizeof(double));
	    	FT_VectorMemoryAlloc((POINTER*)&x,top_gmax[1],sizeof(double));
	    }
	    icrd_nb[0] = icoords[0] = top_gmax[0]/2;
	    for (icoords[1] = 0; icoords[1] <= top_gmax[1]; ++icoords[1])
	    {
	    	icrd_nb[1] = icoords[1];
	        for (icoords[2] = 2; icoords[2] < top_gmax[2]-1; ++icoords[2])
		{
		    index = d_index(icoords,top_gmax,dim);
	    	    icrd_nb[2] = icoords[2] + 1;
		    index_nb = d_index(icrd_nb,top_gmax,dim);
		    if (top_comp[index] != top_comp[index_nb] &&
			FT_StateStructAtGridCrossing(front,grid_intfc,icoords,
                                UPPER,top_comp[index],&l_state,&hs,crx_coords)
                        &&
                        FT_StateStructAtGridCrossing(front,grid_intfc,icrd_nb,
                                LOWER,top_comp[index_nb],&u_state,&hs,
                                crx_coords))
		    {
			pl[n] = getStatePres(l_state);
			pu[n] = getStatePres(u_state);
			vz[n] = getStateZvel(l_state);
			x[n] = crx_coords[1];
			n++;
		    }
		}
	    }
	    sprintf(pname,"cpres-%d.xg",front->step);
	    pfile = fopen(pname,"w");
	    fprintf(pfile,"\"Lower pressure\"\n");
	    for (i = 0; i < n; ++i)
		fprintf(pfile,"%f %f\n",x[i],pl[i]);
	    fprintf(pfile,"\n\n\"Upper pressure\"\n");
	    for (i = 0; i < n; ++i)
		fprintf(pfile,"%f %f\n",x[i],pu[i]);
	    fprintf(pfile,"\n\n\"Pressure difference\"\n");
	    for (i = 0; i < n; ++i)
		fprintf(pfile,"%f %f\n",x[i],pl[i]-pu[i]);
	    fclose(pfile);
	    sprintf(pname,"cvelz-%d.xg",front->step);
	    pfile = fopen(pname,"w");
	    fprintf(pfile,"\"Z-velocity\"\n");
	    for (i = 0; i < n; ++i)
		fprintf(pfile,"%f %f\n",x[i],vz[i]);
	    fclose(pfile);
	}
	if (debugging("trace"))
	    (void) printf("Leaving coating_mono_hyper_surf3d()\n");
}	/* end coating_mono_hyper_surf3d */

static void load_node_propagate(
	Front *front,
	NODE *oldn,
	NODE *newn,
	double dt)
{
	IF_PARAMS *iFparams = (IF_PARAMS*)front->extra1;
	AF_PARAMS *af_params = (AF_PARAMS*)front->extra2;
	POINT *oldp,*newp;
	double *g = iFparams->gravity;
	double f[MAXD],accel[MAXD];
	double kl = af_params->kl;
	double mass = af_params->payload;
	CURVE **c;
	STATE *sl,*sr,*newsl,*newsr;
	double vec[MAXD],vec_mag;
	BOND *b;
	int i,dim = FT_Dimension();

	if (!is_load_node(oldn)) return;
	oldp = oldn->posn;
	newp = newn->posn;
	sl = (STATE*)left_state(oldp);
	sr = (STATE*)right_state(oldp);
	newsl = (STATE*)left_state(newp);
	newsr = (STATE*)right_state(newp);

	if (debugging("trace"))
	{
	    (void)printf("\nEntering load_node_propagate()\n");
	}
	for (i = 0; i < dim; ++i)
	    f[i] = 0.0;
	node_out_curve_loop(oldn,c)
	{
	    b = (*c)->first;
	    for (i = 0; i < dim; ++i)
	    {
		vec[i] = Coords(b->end)[i] - Coords(b->start)[i];
		vec[i] /= bond_length(b);
		f[i] += kl*(bond_length(b) - bond_length0(b))*vec[i];
	    }
	}
	node_in_curve_loop(oldn,c)
	{
	    b = (*c)->last;
	    for (i = 0; i < dim; ++i)
	    {
		vec[i] = Coords(b->start)[i] - Coords(b->end)[i];
		vec[i] /= bond_length(b);
		f[i] += kl*(bond_length(b) - bond_length0(b))*vec[i];
	    }
	}
	for (i = 0; i < dim; ++i)
	{
	    accel[i] = f[i]/mass;
	    newsl->fluid_accel[i] = newsr->fluid_accel[i] = 0.0;
	    newsr->other_accel[i] = newsl->other_accel[i] = accel[i];
	    newsl->impulse[i] = newsr->impulse[i] = sl->impulse[i];
	    newsl->vel[i] = newsr->vel[i] = sl->vel[i] + 
				(accel[i] + g[i]) * dt;
	}
	node_out_curve_loop(newn,c)
	{
	    b = (*c)->first;
	    set_bond_length(b,dim);
	}
	node_in_curve_loop(newn,c)
	{
	    b = (*c)->last;
	    set_bond_length(b,dim);
	}
	
	if (debugging("trace"))
	{
	    (void)printf("accel = %f %f %f\n",accel[0],accel[1],
				accel[2]);
	    (void)printf("Leaving load_node_propagate()\n\n");
	}
}	/* end load_node_propagate */

static void rg_string_node_propagate(
        Front *front,
        NODE *oldn,
        NODE *newn,
        double dt)
{
	IF_PARAMS *iFparams = (IF_PARAMS*)front->extra1;
        AF_PARAMS *af_params = (AF_PARAMS*)front->extra2;
	RECT_GRID *gr = computational_grid(front->interf);
        POINT *oldp,*newp;
        double *g = iFparams->gravity;
        double f[MAXD],accel[MAXD];
        double kl = af_params->kl;
        double mass = af_params->payload;
        CURVE **c;
        STATE *sl,*sr,*newsl,*newsr;
        double vec[MAXD],vec_mag;
        BOND *b;
        int i,dim = FT_Dimension();
	HYPER_SURF *hs;
	HYPER_SURF_ELEMENT *hse;
	POINTER wave;
	double V[MAXD];

        if (!is_rg_string_node(oldn)) return;
	for (i = 0; i < dim; ++i)
	{
	    if (Coords(oldn->posn)[i] <= gr->L[i] || 
		Coords(oldn->posn)[i] > gr->U[i])
		break;
	}
	if (i != dim || oldn->extra == NULL) return;

        if (debugging("trace"))
	{
            (void)printf("\nEntering rg_string_node_propagate()\n");
	}

        oldp = oldn->posn;
        newp = newn->posn;
        sl = (STATE*)left_state(oldp);
        sr = (STATE*)right_state(oldp);
        newsl = (STATE*)left_state(newp);
        newsr = (STATE*)right_state(newp);

        for (i = 0; i < dim; ++i)
            f[i] = 0.0;
	/* calculate the force from the string chords */
	node_out_curve_loop(oldn,c)
        {
            if (hsbdry_type(*c) == PASSIVE_HSBDRY) continue;
            b = (*c)->first;
            for (i = 0; i < dim; ++i)
            {
                vec[i] = Coords(b->end)[i] - Coords(b->start)[i];
                vec[i] /= bond_length(b);
                f[i] += kl*(bond_length(b) - bond_length0(b))*vec[i];
            }
        }
        node_in_curve_loop(oldn,c)
        {
            if (hsbdry_type(*c) == PASSIVE_HSBDRY) continue;
            b = (*c)->last;
            for (i = 0; i < dim; ++i)
            {
                vec[i] = Coords(b->start)[i] - Coords(b->end)[i];
                vec[i] /= bond_length(b);
                f[i] += kl*(bond_length(b) - bond_length0(b))*vec[i];
            }
        }
	/* propagate the nodes along with the rigid body */
	node_out_curve_loop(oldn,c)
        {
            if (hsbdry_type(*c) == PASSIVE_HSBDRY)
	    {
		b = (*c)->first;
		hs = Hyper_surf(b->_btris[0]->surface);
		hse = Hyper_surf_element(b->_btris[0]->tri);
		break;
	    }
        }
        node_in_curve_loop(oldn,c)
        {
            if (hsbdry_type(*c) == PASSIVE_HSBDRY)
	    {
		b = (*c)->last;
		hs = Hyper_surf(b->_btris[0]->surface);
		hse = Hyper_surf_element(b->_btris[0]->tri);
		break;
	    }
        }
	if (hs == NULL || hse == NULL)
	{
	    printf("ERROR in rg_string_node_propagate \n");
	    printf("No related hs or hse found");
	    clean_up(ERROR);
	}
	ifluid_point_propagate(front,wave,oldp,newp,hse,hs,dt,V);
	if (dt > 0.0)
	{
	    for (i = 0; i < dim; ++i)
		accel[i] = (Coords(newp)[i] - Coords(oldp)[i] - 
				oldp->vel[i] * dt) * 2.0 / dt / dt;
	}
	else
	{
	    for (i = 0; i < dim; ++i)
		accel[i] = 0.0;
	}
	for (i = 0; i < dim; ++i)
	    accel[i] -= g[i];

        if (debugging("rigid_body"))
        {
	    (void)printf("accel = %f %f %f\n", accel[0], accel[1], accel[2]);
	    (void)printf("old coords = %f %f %f\n", Coords(oldp)[0], 
				Coords(oldp)[1], Coords(oldp)[2]);
	    (void)printf("oldsl velo = %f %f %f\n", sl->vel[0], 
				sl->vel[1], sl->vel[2]);
	    (void)printf("oldsr velo = %f %f %f\n", sr->vel[0], 
				sr->vel[1], sr->vel[2]);
	    (void)printf("new coords = %f %f %f\n", Coords(newp)[0], 
				Coords(newp)[1], Coords(newp)[2]);
	    (void)printf("newsl velo = %f %f %f\n", newsl->vel[0], 
				newsl->vel[1], newsl->vel[2]);
	    (void)printf("newsr velo = %f %f %f\n", newsr->vel[0], 
				newsr->vel[1], newsr->vel[2]);
	}

	/* Do not change coords, but record the force */
	for (i = 0; i < dim; ++i)
	{
	    Coords(newp)[i] = Coords(oldp)[i];
	    newp->force[i] = f[i];
	    newsl->fluid_accel[i] = newsr->fluid_accel[i] = 
					accel[i] - f[i]/mass;
	    newsr->other_accel[i] = newsl->other_accel[i] = f[i]/mass;
	    newsl->impulse[i] = newsr->impulse[i] = sl->impulse[i];
	}

        if (debugging("trace"))
        {
            (void)printf("Leaving rg_string_node_propagate()\n\n");
        }
}	/* end rg_string_node_propagate */

extern int airfoil_node_propagate(
	Front *front,
	POINTER wave,
	NODE *oldn,
	NODE *newn,
	RPROBLEM        **rp,
        double          dt,
        double          *dt_frac,
        NODE_FLAG       flag,
        POINTER         user)
{
	if (is_load_node(oldn))
	    load_node_propagate(front,oldn,newn,dt);
	else if (is_rg_string_node(oldn))
	    rg_string_node_propagate(front,oldn,newn,dt);
	else
	    return GOOD_NODE;
	return GOOD_NODE;
}	/* end airfoil_node_propagate */

static void passive_curve_propagation(
        Front *front,
        POINTER wave,
        CURVE *oldc,
        CURVE *newc,
        double dt)
{
        BOND *oldb,*newb;
        POINT *oldp,*newp;
	BOND_TRI **btris;
	HYPER_SURF_ELEMENT *hse;
        HYPER_SURF         *hs;
	STATE *sl,*sr;
	STATE *oldst,*btrist;
	double V[MAXD];
	int i,dim = FT_Dimension();

        if (debugging("trace"))
        {
            (void) printf("Entering passive_curve_propagation()\n");
        }
        for (oldb = oldc->first, newb = newc->first; oldb != NULL;
                oldb = oldb->next, newb = newb->next)
        {
            oldp = oldb->end;
            newp = newb->end;
	    for (btris = Btris(oldb); btris && *btris; ++btris)
	    {
		hse = Hyper_surf_element((*btris)->tri);
		hs = Hyper_surf((*btris)->surface);
		FT_GetStatesAtPoint(oldp,hse,hs,(POINTER*)&sl,(POINTER*)&sr);
		if (ifluid_comp(negative_component(hs)))
		{
		    oldst = (STATE*)left_state(oldp);
		    btrist = (STATE*)sl;
		}
		else if (ifluid_comp(positive_component(hs)))
		{
		    oldst = (STATE*)right_state(oldp);
		    btrist = (STATE*)sr;
		}
		for (i = 0; i < dim; ++i)
		    btrist->vel[i] = oldst->vel[i];
	    }
            ifluid_point_propagate(front,wave,oldp,newp,
                        Hyper_surf_element(oldb->_btris[0]->tri),
                        Hyper_surf(oldb->_btris[0]->surface),dt,V);
        }
        if (debugging("trace"))
        {
            (void) printf("Leaving passive_curve_propagation()\n\n");
        }
}       /* end passive_curve_propagation */
