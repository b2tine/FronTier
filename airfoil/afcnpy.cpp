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
#include "solver.h"
#ifdef __COLLISION__
	#include "collid.h"
#endif

static void spring_force_at_point1(double*,POINT*,TRI*,SURFACE*,double);
static void spring_force_at_point2(double*,POINT*,TRI*,SURFACE*,double);
static boolean is_pore(Front*,HYPER_SURF_ELEMENT*,double*);
static void compute_total_canopy_force2d(Front*,double*,double*);
static void compute_total_canopy_force3d(Front*,double*,double*);
static void compute_center_of_mass_velo(ELASTIC_SET*);
static void reduce_high_freq_vel(Front*,SURFACE*);
static void smooth_vel(double*,POINT*,TRI*,SURFACE*);
static boolean curve_in_pointer_list(CURVE*,CURVE**);
static void setNodeVelocity(ELASTIC_SET*,NODE*,double**,GLOBAL_POINT**);
static void setCurveVelocity(ELASTIC_SET*,CURVE*,double**,GLOBAL_POINT**);
static void setSurfVelocity(ELASTIC_SET*,SURFACE*,double**,GLOBAL_POINT**);
static void new_setNodeVelocity(ELASTIC_SET*,NODE*,GLOBAL_POINT**);
static void new_setNodeVelocity2d(ELASTIC_SET*,NODE*,GLOBAL_POINT**);
static void new_setNodeVelocity3d(ELASTIC_SET*,NODE*,GLOBAL_POINT**);
static void new_setCurveVelocity(ELASTIC_SET*,CURVE*,double**,GLOBAL_POINT**);
static void new_setSurfVelocity(ELASTIC_SET*,SURFACE*,double**,GLOBAL_POINT**);
static void setCollisionFreePoints3d(INTERFACE*);
static void break_string_curve(CURVE*,double);

#define 	MAX_NUM_RING1		30

static void spring_force_at_point1(
	double *f,
	POINT *p,
	TRI *tri,
	SURFACE *surf,
	double ks)
{
	TRI *tris[MAX_NUM_RING1];
	int i,j,k,nt;
	POINT *p_nb;
	double length0,length,dir[3];
	
	if (is_registered_point(surf,p))
	{
	    for (i = 0; i < 3; ++i)
		f[i] = 0.0;
	    return;
	}
	PointAndFirstRingTris(p,Hyper_surf_element(tri),Hyper_surf(surf),
				&nt,tris);
	for (k = 0; k < 3; ++k) f[k] = 0.0;
	for (i = 0; i < nt; ++i)
	{
	    for (j = 0; j < 3; ++j)
	    {
		if (Point_of_tri(tris[i])[j] == p)
		{
		    length0 = tris[i]->side_length0[j];
		    p_nb = Point_of_tri(tris[i])[(j+1)%3];
		    length = separation(p,p_nb,3);
	    	    for (k = 0; k < 3; ++k)
		    {
			dir[k] = (Coords(p_nb)[k] - Coords(p)[k])/length;
			f[k] += ks*(length - length0)*dir[k];
		    }
		    if (is_side_bdry(tris[i],(j+2)%3))
		    {
			(void) printf("Detect boundary "
				"in spring_force_at_point1()\n");
			clean_up(ERROR);
		    }
		}
	    }
	}
}	/* end spring_force_at_point1 */

extern void compute_total_canopy_force(
	Front *front,
	double *pos_force,
	double *neg_force)
{
	int dim = front->rect_grid->dim;
	switch (dim)
	{
	case 2:
	    compute_total_canopy_force2d(front,pos_force,neg_force);
	    return;
	case 3:
	    compute_total_canopy_force3d(front,pos_force,neg_force);
	    return;
	}
}	/* end compute_total_canopy_force */

static void compute_total_canopy_force2d(
	Front *front,
        double *pos_force,
        double *neg_force)
{
		
}	/* end compute_total_canopy_force2d */

static void compute_total_canopy_force3d(
	Front *front,
        double *pos_force,
        double *neg_force)
{
	TRI *tri;
	SURFACE *surf;
	INTERFACE *intfc = front->interf;
	POINT *p;
	STATE *sl,*sr;
	double pres_p,pres_m;
	double area[MAXD];
	int i;
	static FILE *pfile;

	if (debugging("trace"))
	    (void) printf("Entering compute_total_canopy_force3d()\n");
	if (pfile == NULL)
	{
	    pfile = fopen("payload","w");
	    fprintf(pfile,"\"Net lift vs time\"\n");
	}
	for (i = 0; i < 3; ++i)
	    pos_force[i] = neg_force[i] = 0.0;
	next_tri(intfc,NULL,NULL);
	while (next_tri(intfc,&tri,&surf))
	{
	    if (wave_type(surf) != ELASTIC_BOUNDARY)
		continue; 
	    pres_p = pres_m = 0.0;
	    for (i = 0; i < 3; ++i)
	    {
		p = Point_of_tri(tri)[i];
		FT_GetStatesAtPoint(p,Hyper_surf_element(tri),Hyper_surf(surf),
				(POINTER*)&sl,(POINTER*)&sr);
		pres_m += sl->pres;
		pres_p += sr->pres;
		area[i] = Tri_normal(tri)[i];
	    }
	    for (i = 0; i < 3; ++i)
	    {
		pos_force[i] -= pres_p*area[i]/3.0;
		neg_force[i] += pres_m*area[i]/3.0;
	    }
	}
	if (debugging("trace"))
	    (void) printf("Leaving compute_total_canopy_force3d()\n");
}	/* end compute_total_canopy_force3d */

extern int airfoil_velo(
        POINTER params,
        Front *front,
        POINT *p,
        HYPER_SURF_ELEMENT *hse,
        HYPER_SURF *hs,
        double *vel)
{
	/*TMP to be written*/
        vel[0] = vel[1] = vel[2] = 0.0;
	return YES;
}       /* end airfoil_velo */

extern int af_find_state_at_crossing(
        Front *front,
        int *icoords,
        GRID_DIRECTION dir,
        int comp,
        POINTER *state,
        HYPER_SURF **hs,
        double *crx_coords)
{
        boolean status;
	HYPER_SURF_ELEMENT *hse;
	AF_PARAMS *af_params = (AF_PARAMS*)front->extra2;

        status = FT_StateStructAtGridCrossing2(front,icoords,dir,comp,state,hs,
                                        &hse,crx_coords);
        if (status == NO) 
	    return NO_PDE_BOUNDARY;
        if (wave_type(*hs) == FIRST_PHYSICS_WAVE_TYPE) 
	    return NO_PDE_BOUNDARY;
	if (wave_type(*hs) == ELASTIC_STRING)//2d string 
	    return NO_PDE_BOUNDARY;
	if (wave_type(*hs) == ELASTIC_BOUNDARY && 
		af_params->with_porosity)
	    return NO_PDE_BOUNDARY;
	if (wave_type(*hs) == ELASTIC_BOUNDARY) //TMP
	    return CONST_V_PDE_BOUNDARY;
	if (wave_type(*hs) == DIRICHLET_BOUNDARY)
	{
	    if (boundary_state_function(*hs) &&
             strcmp(boundary_state_function_name(*hs),
                "flowThroughBoundaryState") == 0)
            {
                return CONST_P_PDE_BOUNDARY;
            }
	    else
	    {
		return CONST_V_PDE_BOUNDARY;
	    }
	}
        return NEUMANN_PDE_BOUNDARY;
}       /* af_find_state_at_crossing */

static boolean is_pore(
	Front *front,
	HYPER_SURF_ELEMENT *hse,
	double *crx_coords)
{
	INTERFACE *intfc = front->grid_intfc;
	SURFACE **s;
	TRI *tri;
	AF_PARAMS *af_params = (AF_PARAMS*)front->extra2;
	double gamma = af_params->gamma;
	static int current_step = -1;

	if (front->rect_grid->dim != 3) return NO;
	if (gamma == 0.0) return NO;
	if (front->step != current_step)
	{
	    double R;
	    current_step = front->step;
	    for (s = intfc->surfaces; s && *s; ++s)
	    {
		if (wave_type(*s) != ELASTIC_BOUNDARY) continue;
		for (tri = first_tri(*s); !at_end_of_tri_list(tri,*s);
				tri = tri->next)
		{
		    R = drand48();
		    if (R < gamma) Tri_index(tri) = 0;
		    else Tri_index(tri) = 1;
		}
	    }
	}
	tri = Tri_of_hse(hse);
	return (Tri_index(tri) == 0) ? YES : NO;
}	/* end is_pore */

extern void assign_node_field(
	NODE *node,
	double **x,
	double **v,
	int *n)
{
	int i,dim = Dimension(node->interface);
	CURVE **c;

	for (i = 0; i < dim; ++i)
	{
	    Coords(node->posn)[i] = x[*n][i];
	    node->posn->vel[i] = v[*n][i];
	}
	for (c = node->out_curves; c && *c; ++c)
	    set_bond_length((*c)->first,dim);
	for (c = node->in_curves; c && *c; ++c)
	    set_bond_length((*c)->last,dim);
	(*n)++;
}	/* end assign_node_field */
	
extern void assign_curve_field(
	CURVE *curve,
	double **x,
	double **v,
	int *n)
{
	int i,j,dim = Dimension(curve->interface);
	BOND *b;

	i = *n;
	for (b = curve->first; b != curve->last; b = b->next)
	{
	    for (j = 0; j < dim; ++j)
	    {
	    	Coords(b->end)[j] = x[i][j];
	    	b->end->vel[j] = v[i][j];
	    }
	    set_bond_length(b,dim);
	    i++;
	}
	set_bond_length(curve->first,dim);
	set_bond_length(curve->last,dim);
	*n = i;
}	/* end assign_curve_field */
	
extern void assign_surf_field(
	SURFACE *surf,
	double **x,
	double **v,
	int *n)
{
	int i,j,k;
	TRI *tri;
	POINT *p;
	STATE *sl,*sr;

	unsort_surf_point(surf);
	i = *n;
	for (tri = first_tri(surf); !at_end_of_tri_list(tri,surf); 
			tri = tri->next)
	{
	    for (j = 0; j < 3; ++j)
	    {
		p = Point_of_tri(tri)[j];
		if (sorted(p) || Boundary_point(p)) continue;
		sl = (STATE*)left_state(p);
		sr = (STATE*)right_state(p);
		for (k = 0; k < 3; ++k)
		{
		    Coords(p)[k] = x[i][k];
		    p->vel[k] = v[i][k];
		}
		sorted(p) = YES;
	    	++i;
	    }
	}
	*n = i;
}	/* end assign_surf_field */
	
extern void compute_surf_accel1(
	ELASTIC_SET *geom_set,
	SURFACE *surf,
	double **f,
	double **x,
	double **v,
	int *n)
{
	int j,k;
	TRI *tri;
	POINT *p;
	int dim = 3;
	double ks = geom_set->ks;
	double m_s = geom_set->m_s;
	double lambda_s = geom_set->lambda_s;

	unsort_surf_point(surf);
	for (tri = first_tri(surf); !at_end_of_tri_list(tri,surf); 
			tri = tri->next)
	{
	    for (j = 0; j < 3; ++j)
	    {
		p = Point_of_tri(tri)[j];
		if (sorted(p) || Boundary_point(p)) continue;
		for (k = 0; k < dim; ++k)
		{
		    x[*n][k] = Coords(p)[k];
		    v[*n][k] = p->vel[k];
		}
	    	spring_force_at_point1(f[*n],p,tri,surf,ks);
		for (k = 0; k < dim; ++k)
		{
		    f[*n][k] -= lambda_s*(v[*n][k]);
		    f[*n][k] /= m_s;
		}
		sorted(p) = YES;
	    	++(*n);
	    }
	}
}	/* end compute_surf_accel1 */

extern void compute_curve_accel1(
	ELASTIC_SET *geom_set,
	CURVE *curve,
	double **f,
	double **x,
	double **v,
	int *n)
{
	int i,j;
	double x_diff;
	BOND *b;
	double dir[MAXD],len,len0,vect[MAXD];
	int dim = Dimension(curve->interface);
	double kl,m_l,lambda_l;

	if (dim == 3)
	{
	    if (hsbdry_type(curve) == STRING_HSBDRY)
	    {
	    	kl = geom_set->kl;
	    	m_l = geom_set->m_l;
	    	lambda_l = geom_set->lambda_l;
	    }
	    else if (hsbdry_type(curve) == GORE_HSBDRY)
	    {
	    	kl = geom_set->kg;
	    	m_l = geom_set->m_g;
	    	lambda_l = geom_set->lambda_g;
	    }
	    else
	    {
	    	kl = geom_set->ks;
	    	m_l = geom_set->m_s;
	    	lambda_l = geom_set->lambda_s;
	    }
	}
	else
	{
	    kl = geom_set->kl;
	    m_l = geom_set->m_l;
	    lambda_l = geom_set->lambda_l;
	}
	i = *n;
	for (b = curve->first; b != curve->last; b = b->next)
	{
	    x[i] = Coords(b->end);
	    v[i] = b->end->vel;
	    for (j = 0; j < dim; ++j)
	    {
		f[i][j] = -lambda_l*v[i][j]/m_l;
	    }
	    i++;
	}

	i = *n;
	for (b = curve->first; b != NULL; b = b->next)
	{
	    len = separation(b->start,b->end,dim);
	    len0 = bond_length0(b);
	    x_diff = len - len0; 
	    for (j = 0; j < dim; ++j)
	    {
		dir[j] = (Coords(b->end)[j] - Coords(b->start)[j])
				/len;
		vect[j] = x_diff*dir[j];
		if (b != curve->first)
		{
	    	    f[i-1][j]   += kl*vect[j]/m_l;
		}
		if (b != curve->last)
		{
	    	    f[i][j] -= kl*vect[j]/m_l;
		}
	    }
	    if (b != curve->last) i++;
	}

	if (dim == 3)
	{
	    POINT *p,*p_nb;
	    BOND_TRI **btris;
	    TRI **tris;
	    int j,k,side,nt;
	    double length0,length;
	    double ks = geom_set->ks;
	    i = *n;
	    for (b = curve->first; b != curve->last; b = b->next)
	    {
		p = b->end;
		for (btris = Btris(b); btris && *btris; ++btris)
		{
		    nt = I_FirstRingTrisAroundPoint(p,(*btris)->tri,&tris);
		    for (j = 0; j < nt; ++j)
		    {
			for (side = 0; side < 3; ++side)
			{
			    if (p == Point_of_tri(tris[j])[side])
			    {
				if (is_side_bdry(tris[j],side))
				    continue;
				p_nb = Point_of_tri(tris[j])[(side+1)%3];
				length0 = tris[j]->side_length0[side];
				length = separation(p,p_nb,3);
				for (k = 0; k < 3; ++k)
                        	{
                            	    dir[k] = (Coords(p_nb)[k] - 
						Coords(p)[k])/length;
                            	    f[i][k] += ks*(length - length0)*
					dir[k]/m_l;
                        	}
			    }
			}
		    }
		}
		i++;
	    }
	}
	*n = i;
}	/* end compute_curve_accel1 */

extern void compute_node_accel1(
	ELASTIC_SET *geom_set,
	NODE *node,
	double **f,
	double **x,
	double **v,
	int *n)
{
	CURVE **c;
	BOND *b;
	double x_diff,len0,len,dir[MAXD],vect[MAXD];
	POINT *p,*p_nb;
	INTERFACE *intfc = geom_set->front->interf;
	int i,j,dim = Dimension(intfc);
	double ks = geom_set->ks;
	double kl = geom_set->kl;
	double kg = geom_set->kg;
	double mass;
	double lambda_s = geom_set->lambda_s;
	double lambda_l = geom_set->lambda_l;
	double lambda_g = geom_set->lambda_g;

	if (dim == 3)
	{
	    AF_NODE_EXTRA *extra = (AF_NODE_EXTRA*)node->extra;
	    if (extra != NULL)
	    {
		if (extra->af_node_type == LOAD_NODE)
		{
	    	    Front *front = geom_set->front;
	    	    AF_PARAMS *af_params = (AF_PARAMS*)front->extra2;
		    mass = af_params->payload;
		}
		else if (extra->af_node_type == GORE_NODE)
                    mass = geom_set->m_g;
		else if (extra->af_node_type == STRING_NODE)
                    mass = geom_set->m_s;
	    }
	    else
                mass = geom_set->m_s;
	}
	else if (dim == 2)
	{
	    AF_NODE_EXTRA *extra = (AF_NODE_EXTRA*)node->extra;
            mass = geom_set->m_l;
	    if (extra->af_node_type == LOAD_NODE)
            {
                Front *front = geom_set->front;
                AF_PARAMS *af_params = (AF_PARAMS*)front->extra2;
                mass = af_params->payload;
            }
	}

	x[*n] = Coords(node->posn);
	v[*n] = node->posn->vel;
	for (i = 0; i < dim; ++i)
	{
	    f[*n][i] = 0.0;
	}
	for (c = node->out_curves; c && *c; ++c)
	{
	    b = (*c)->first;
	    len = separation(b->start,b->end,dim);
	    len0 = bond_length0(b);
	    x_diff = len - len0; 
	    for (j = 0; j < dim; ++j)
	    {
		dir[j] = (Coords(b->end)[j] - Coords(b->start)[j])
				/len;
		vect[j] = x_diff*dir[j];
		if (dim == 3)
		{
		    if (is_load_node(node) == YES)
		    	f[*n][j]   += kl*vect[j]/mass;
		    else if (hsbdry_type(*c) == STRING_HSBDRY)
	    	    	f[*n][j]   += kl*vect[j]/mass;
		    else if (hsbdry_type(*c) == MONO_COMP_HSBDRY)
	    	    	f[*n][j]   += ks*vect[j]/mass;
		    else if (hsbdry_type(*c) == GORE_HSBDRY)
	    	    	f[*n][j]   += kg*vect[j]/mass;
		}
		else
		    f[*n][j]   += kl*vect[j]/mass;
	    }
	}
	for (c = node->in_curves; c && *c; ++c)
	{
	    if (curve_in_pointer_list(*c,node->out_curves) && 
		!is_closed_curve(*c)) 
		continue;
	    b = (*c)->last;
	    len = separation(b->start,b->end,dim);
	    len0 = bond_length0(b);
	    x_diff = len - len0; 
	    for (j = 0; j < dim; ++j)
	    {
		dir[j] = (Coords(b->end)[j] - Coords(b->start)[j])
				/len;
		vect[j] = x_diff*dir[j];
		if (dim == 3)
		{
		    if (is_load_node(node) == YES)
		    	f[*n][j]   -= kl*vect[j]/mass;
		    else if (hsbdry_type(*c) == STRING_HSBDRY)
	    	    	f[*n][j]   -= kl*vect[j]/mass;
		    else if (hsbdry_type(*c) == MONO_COMP_HSBDRY)
	    	    	f[*n][j]   -= ks*vect[j]/mass;
		    else if (hsbdry_type(*c) == GORE_HSBDRY)
	    	    	f[*n][j]   -= kg*vect[j]/mass;
		}
		else
		    f[*n][j]   -= kl*vect[j]/mass;
	    }
	}
	if (dim == 3)
	{
	    BOND_TRI **btris;
	    TRI **tris,*tri_list[500];
	    int k,side,nt,num_tris;
	    TRI *tri;

	    num_tris = 0;
	    p = node->posn;
	    for (c = node->out_curves; c && *c; ++c)
	    {
		b = (*c)->first;
		for (btris = Btris(b); btris && *btris; ++btris)
		{
		    nt = I_FirstRingTrisAroundPoint(p,(*btris)->tri,&tris);
		    for (j = 0; j < nt; ++j)
		    {
			if (!pointer_in_list((POINTER)tris[j],num_tris,
					(POINTER*)tri_list))
			    tri_list[num_tris++] = tris[j];
		    }
		}
	    }
	    for (c = node->in_curves; c && *c; ++c)
	    {
		b = (*c)->last;
		for (btris = Btris(b); btris && *btris; ++btris)
		{
		    nt = I_FirstRingTrisAroundPoint(p,(*btris)->tri,&tris);
		    for (j = 0; j < nt; ++j)
		    {
			if (!pointer_in_list((POINTER)tris[j],num_tris,
					(POINTER*)tri_list))
			    tri_list[num_tris++] = tris[j];
		    }
		}
	    }
	    for (i = 0; i < num_tris; ++i)
	    {
		tri = tri_list[i];
		for (side = 0; side < 3; ++side)
		{
		    if (p == Point_of_tri(tri)[side])
		    {
			if (is_side_bdry(tri,side))
			    continue;
			p_nb = Point_of_tri(tri)[(side+1)%3];
			len0 = tri->side_length0[side];
			len = separation(p,p_nb,3);
    			x_diff = len - len0; 
			for (k = 0; k < 3; ++k)
                       	{
                       	    dir[k] = (Coords(p_nb)[k] - 
					Coords(p)[k])/len;
                       	    f[*n][k] += ks*x_diff*dir[k]/mass;
                       	}
		    }
		}
	    }
	    if (!is_load_node(node))
	    {
	    	for (i = 0; i < 3; ++i)
	    	    f[*n][i] -= lambda_s*v[*n][i]/mass;
	    }
	}
	else
	{
	    for (i = 0; i < 3; ++i)
	    	f[*n][i] -= lambda_l*v[*n][i]/mass;
	}
	(*n)++;
}	/* end compute_node_accel1 */

static void compute_center_of_mass_velo(
	ELASTIC_SET *geom_set)
{
	int i,j,n;
	TRI *tri;
	POINT *p;
	STATE *state;
	Front *front = geom_set->front;
	SURFACE **surfs = geom_set->surfs;
	SURFACE *canopy;
	NODE *node;
	AF_PARAMS *af_params = (AF_PARAMS*)front->extra2;
	double area_dens = af_params->area_dens;
	double xt[MAXD],vt[MAXD],xcan[MAXD],vcan[MAXD],xload[MAXD],vload[MAXD];
	double area,mass_canopy,payload;
	double *xcom,*vcom;

	if (debugging("canopy"))
	    (void) printf("Entering compute_center_of_mass_velo()\n");

	for (n = 0; n < geom_set->num_surfs; ++n)
	{
	    canopy = geom_set->surfs[n];
	    for (j = 0; j < 3; ++j)
	    	vcan[j] = 0.0;
	    area = mass_canopy = 0.0;
	    surf_tri_loop(canopy,tri)
	    {
	    	for (j = 0; j < 3; ++j)
	    	{
		    vt[j] = 0.0;
		    xt[j] = 0.0;
	    	}
	    	for (i = 0; i < 3; ++i)
	    	{
		    p = Point_of_tri(tri)[i];
		    state = (STATE*)left_state(p);
		    for (j = 0; j < 3; ++j)
		    {
		    	vt[j] += state->vel[j]/3.0;
		    	xt[j] += Coords(p)[j]/3.0;
		    }
	    	}
	    	for (j = 0; j < 3; ++j)
	    	{
		    vcan[j] += vt[j]*tri_area(tri);
		    xcan[j] += xt[j]*tri_area(tri);
	    	}
	    	area += tri_area(tri);
	    }
	    mass_canopy += area_dens*area;
	    for (j = 0; j < 3; ++j)
	    {
	    	vcan[j] /= area;
	    	xcan[j] /= area;
	    }

	    if (NULL != geom_set->load_node)
	    {
	  	node = geom_set->load_node;
	  	state = (STATE*)left_state(node->posn);
	  	for (j = 0; j < 3; ++j)
	  	{
	    	    vload[j] = state->vel[j];
	    	    xload[j] = Coords(node->posn)[j];
	  	}
	  	payload = af_params->payload;
	  	xcom = center_of_mass(Hyper_surf(canopy));
	  	vcom = center_of_mass_velo(Hyper_surf(canopy));
	  	for (j = 0; j < 3; ++j)
	  	{
	    	    vcom[j] = (vcan[j]*mass_canopy + vload[j]*payload)/
				(mass_canopy + payload);
	    	    xcom[j] = (xcan[j]*mass_canopy + xload[j]*payload)/
				(mass_canopy + payload);
	  	}
	    }
	    else
	    {
	  	xcom = center_of_mass(Hyper_surf(canopy));
	  	vcom = center_of_mass_velo(Hyper_surf(canopy));
	    }	
	}
	if (debugging("canopy"))
	    (void) printf("Leaving compute_center_of_mass_velo()\n");
}	/* end compute_center_of_mass_velo */

static void smooth_vel(
	double *vel,
	POINT *p,
	TRI *tri,
	SURFACE *surf)
{
	TRI *tris[20];
	HYPER_SURF_ELEMENT *hse = Hyper_surf_element(tri);
        HYPER_SURF         *hs = Hyper_surf(surf);
	int i,j,k,nt,np;
	POINT *pt_list[20],*pt;
	STATE *sl,*sr;
	static double max_speed = 0.0;
	
	PointAndFirstRingTris(p,hse,hs,&nt,tris);
	np = 0;
	for (k = 0; k < 3; ++k)
	    vel[k] = 0.0;
	for (i = 0; i < nt; ++i)
	{
	    for (j = 0; j < 3; ++j)
	    {
		pt = Point_of_tri(tris[i])[j];
		if (pointer_in_list((POINTER)pt,np,(POINTER*)pt_list))
		    continue;
		pt_list[np++] = pt;	
		hse = Hyper_surf_element(tris[i]);
		FT_GetStatesAtPoint(pt,hse,hs,(POINTER*)&sl,(POINTER*)&sr);
		for (k = 0; k < 3; ++k)
	    	    vel[k] += sl->vel[k];
		
	    }
	}
	for (k = 0; k < 3; ++k)
	    vel[k] /= (double)np;
}	/* end smooth_vel */

extern void compute_node_accel2(
	ELASTIC_SET *geom_set,
	NODE *node,
	double **f,
	double **x,
	double **v,
	int *n)
{
	CURVE **c;
	BOND *b;
	double x_diff,len0,len,dir[MAXD],vect[MAXD];
	POINT *p,*p_nb;
	int i,j,dim = Dimension(node->interface);
	double ks = geom_set->ks;
	double kl = geom_set->kl;
	double kg = geom_set->kg;
	double mass;
	double lambda_s = geom_set->lambda_s;
	double lambda_l = geom_set->lambda_l;
	double lambda_g = geom_set->lambda_g;

	if (dim == 3)
	{
	    AF_NODE_EXTRA *extra = (AF_NODE_EXTRA*)node->extra;
	    if (extra != NULL)
	    {
		if (extra->af_node_type == LOAD_NODE)
		{
	    	    Front *front = geom_set->front;
	    	    AF_PARAMS *af_params = (AF_PARAMS*)front->extra2;
	    	    mass = af_params->payload;
		}
		else if (extra->af_node_type == GORE_NODE)
		    mass = geom_set->m_g;
		else if (extra->af_node_type == STRING_NODE)
		    mass = geom_set->m_s;
	    }
	    else
		mass = geom_set->m_s;
	}
	else
	    mass = geom_set->m_l;

	for (i = 0; i < dim; ++i)
	{
	    x[*n][i] = Coords(node->posn)[i];
	    v[*n][i] = node->posn->vel[i];
	    f[*n][i] = 0.0;
	}
	for (c = node->out_curves; c && *c; ++c)
	{
	    b = (*c)->first;
	    len0 = bond_length0(b);
	    for (j = 0; j < dim; ++j)
	    {
		vect[j] = Coords(b->end)[j] - Coords(b->start)[j] -
				len0*b->dir0[j];
		if (dim == 3)
		{
		    if (is_load_node(node) == YES)
		    	f[*n][j]   += kl*vect[j]/mass;
		    else if (hsbdry_type(*c) == STRING_HSBDRY)
	    	    	f[*n][j]   += kl*vect[j]/mass;
		    else if (hsbdry_type(*c) == MONO_COMP_HSBDRY)
	    	    	f[*n][j]   += ks*vect[j]/mass;
		    else if (hsbdry_type(*c) == GORE_HSBDRY)
	    	    	f[*n][j]   += kg*vect[j]/mass;
		}
		else
		    f[*n][j]   += kl*vect[j]/mass;
	    }
	}
	for (c = node->in_curves; c && *c; ++c)
	{
	    if (curve_in_pointer_list(*c,node->out_curves) && 
		!is_closed_curve(*c)) 
		continue;
	    b = (*c)->last;
	    len0 = bond_length0(b);
	    for (j = 0; j < dim; ++j)
	    {
		vect[j] = Coords(b->end)[j] - Coords(b->start)[j] -
				len0*b->dir0[j];
		if (dim == 3)
		{
		    if (is_load_node(node) == YES)
		    	f[*n][j]   -= kl*vect[j]/mass;
		    else if (hsbdry_type(*c) == STRING_HSBDRY)
	    	    	f[*n][j]   -= kl*vect[j]/mass;
		    else if (hsbdry_type(*c) == MONO_COMP_HSBDRY)
	    	    	f[*n][j]   -= ks*vect[j]/mass;
		    else if (hsbdry_type(*c) == GORE_HSBDRY)
	    	    	f[*n][j]   -= kg*vect[j]/mass;
		}
		else
		    f[*n][j]   -= kl*vect[j]/mass;
	    }
	}
	if (dim == 3)
	{
	    BOND_TRI **btris;
	    TRI **tris,*tri_list[500];
	    int k,side,nt,num_tris;
	    TRI *tri;

	    num_tris = 0;
	    p = node->posn;
	    for (c = node->out_curves; c && *c; ++c)
	    {
		b = (*c)->first;
		for (btris = Btris(b); btris && *btris; ++btris)
		{
		    nt = I_FirstRingTrisAroundPoint(p,(*btris)->tri,&tris);
		    for (j = 0; j < nt; ++j)
		    {
			if (!pointer_in_list((POINTER)tris[j],num_tris,
					(POINTER*)tri_list))
			    tri_list[num_tris++] = tris[j];
		    }
		}
	    }
	    for (c = node->in_curves; c && *c; ++c)
	    {
		b = (*c)->last;
		for (btris = Btris(b); btris && *btris; ++btris)
		{
		    nt = I_FirstRingTrisAroundPoint(p,(*btris)->tri,&tris);
		    for (j = 0; j < nt; ++j)
		    {
			if (!pointer_in_list((POINTER)tris[j],num_tris,
					(POINTER*)tri_list))
			    tri_list[num_tris++] = tris[j];
		    }
		}
	    }
	    for (i = 0; i < num_tris; ++i)
	    {
		tri = tri_list[i];
		for (side = 0; side < 3; ++side)
		{
		    if (p == Point_of_tri(tri)[side])
		    {
			if (is_side_bdry(tri,side))
			    continue;
			p_nb = Point_of_tri(tri)[(side+1)%3];
			len0 = tri->side_length0[side];
			len = separation(p,p_nb,3);
			x_diff = len - len0;
			for (k = 0; k < 3; ++k)
                       	{
                       	    dir[k] = tri->side_dir0[side][k]; 
                       	    vect[k] = (Coords(p_nb)[k] - Coords(p)[k])
					- len0*dir[k];
                       	    f[*n][k] += ks*vect[k]/mass;
                       	}
		    }
		}
	    }
	    if (!is_load_node(node))
	    {
	    	for (i = 0; i < 3; ++i)
	    	    f[*n][i] -= lambda_s*v[*n][i]/mass;
	    }
	}
	else
	{
	    for (i = 0; i < 3; ++i)
	    	f[*n][i] -= lambda_l*v[*n][i]/mass;
	}
	(*n)++;
}	/* end compute_node_accel2 */

extern void compute_curve_accel2(
	ELASTIC_SET *geom_set,
	CURVE *curve,
	double **f,
	double **x,
	double **v,
	int *n)
{
	int i,j;
	double x_diff;
	BOND *b;
	double dir[MAXD],len0,vect[MAXD];
	int dim = Dimension(curve->interface);
	double kl,m_l,lambda_l,ks,lambda_s;

	if (dim == 3)
	{
	    if (hsbdry_type(curve) == STRING_HSBDRY)
	    {
	    	kl = geom_set->kl;
	    	m_l = geom_set->m_l;
	    	lambda_l = geom_set->lambda_l;
	    }
	    else if (hsbdry_type(curve) == GORE_HSBDRY)
	    {
	    	kl = geom_set->kg;
	    	m_l = geom_set->m_g;
	    	lambda_l = geom_set->lambda_g;
	    }
	    else
	    {
	    	kl = geom_set->ks;
	    	m_l = geom_set->m_s;
	    	lambda_l = geom_set->lambda_s;
	    }
	    ks = geom_set->ks;
	    lambda_s = geom_set->lambda_s;
	}
	else
	{
	    kl = geom_set->kl;
	    m_l = geom_set->m_l;
	    lambda_l = geom_set->lambda_l;
	}
	i = *n;
	for (b = curve->first; b != curve->last; b = b->next)
	{
	    for (j = 0; j < dim; ++j)
	    {
	    	x[i][j] = Coords(b->end)[j];
	    	v[i][j] = b->end->vel[j];
		f[i][j] = -lambda_l*v[i][j]/m_l;
	    }
	    i++;
	}

	i = *n;
	for (b = curve->first; b != NULL; b = b->next)
	{
	    len0 = bond_length0(b);
	    for (j = 0; j < dim; ++j)
	    {
		vect[j] = Coords(b->end)[j] - Coords(b->start)[j] -
				len0*b->dir0[j];
		if (b != curve->first)
		{
	    	    f[i-1][j]   += kl*vect[j]/m_l;
		}
		if (b != curve->last)
		{
	    	    f[i][j] -= kl*vect[j]/m_l;
		}
	    }
	    if (b != curve->last) i++;
	}

	if (dim == 3)
	{
	    POINT *p,*p_nb;
	    BOND_TRI **btris;
	    TRI **tris;
	    int j,k,side,nt;
	    double length0,length;
	    i = *n;
	    for (b = curve->first; b != curve->last; b = b->next)
	    {
		p = b->end;
		for (btris = Btris(b); btris && *btris; ++btris)
		{
		    nt = I_FirstRingTrisAroundPoint(p,(*btris)->tri,&tris);
		    for (j = 0; j < nt; ++j)
		    {
			for (side = 0; side < 3; ++side)
			{
			    if (p == Point_of_tri(tris[j])[side])
			    {
				if (is_side_bdry(tris[j],side))
				    continue;
				p_nb = Point_of_tri(tris[j])[(side+1)%3];
				length0 = tris[j]->side_length0[side];
				for (k = 0; k < 3; ++k)
                        	{
                            	    dir[k] = tris[j]->side_dir0[side][k]; 
				    vect[k] = Coords(p_nb)[k] - Coords(p)[k]
					- length0*dir[k];
                            	    f[i][k] += ks*vect[k]/m_l;
                        	}
			    }
			}
		    }
		}
		i++;
	    }
	}
	*n = i;
}	/* end compute_curve_accel2 */

extern void compute_surf_accel2(
	ELASTIC_SET *geom_set,
	SURFACE *surf,
	double **f,
	double **x,
	double **v,
	int *n)
{
	int j,k;
	TRI *tri;
	POINT *p;
	int dim = 3;
	double ks = geom_set->ks;
	double m_s = geom_set->m_s;
	double lambda_s = geom_set->lambda_s;

	unsort_surf_point(surf);
	for (tri = first_tri(surf); !at_end_of_tri_list(tri,surf); 
			tri = tri->next)
	{
	    for (j = 0; j < 3; ++j)
	    {
		p = Point_of_tri(tri)[j];
		if (sorted(p) || Boundary_point(p)) continue;
		for (k = 0; k < dim; ++k)
		{
		    x[*n][k] = Coords(p)[k];
		    v[*n][k] = p->vel[k];
		}
	    	spring_force_at_point2(f[*n],p,tri,surf,ks);
		for (k = 0; k < dim; ++k)
		{
		    f[*n][k] -= lambda_s*(v[*n][k]);
		    f[*n][k] /= m_s;
		}
		sorted(p) = YES;
	    	++(*n);
	    }
	}
}	/* end compute_surf_accel2 */

static void spring_force_at_point2(
	double *f,
	POINT *p,
	TRI *tri,
	SURFACE *surf,
	double ks)
{
	TRI *tris[MAX_NUM_RING1];
	int i,j,k,nt;
	POINT *p_nb;
	double length0,length,dir[3],vect[3];
	
	PointAndFirstRingTris(p,Hyper_surf_element(tri),Hyper_surf(surf),
				&nt,tris);
	for (k = 0; k < 3; ++k) f[k] = 0.0;
	for (i = 0; i < nt; ++i)
	{
	    for (j = 0; j < 3; ++j)
	    {
		if (Point_of_tri(tris[i])[j] == p)
		{
		    length0 = tris[i]->side_length0[j];
		    p_nb = Point_of_tri(tris[i])[(j+1)%3];
		    length = separation(p,p_nb,3);
	    	    for (k = 0; k < 3; ++k)
		    {
			dir[k] = tris[i]->side_dir0[j][k];
			vect[k] = (Coords(p_nb)[k] - Coords(p)[k]) -
				//length0*tris[i]->side_dir0[j][k];
				length0*dir[k];
			f[k] += ks*vect[k];
		    }
		    if (is_side_bdry(tris[i],(j+2)%3))
		    {
			(void) printf("Detect boundary "
				"in spring_force_at_point2()\n");
			clean_up(ERROR);
		    }
		}
	    }
	}
}	/* end spring_force_at_point2 */

extern void compute_node_accel3(
	ELASTIC_SET *geom_set,
	NODE *node,
	double **f,
	double **x,
	double **v,
	int *n)
{
	CURVE **c;
	BOND *b;
	double x_diff,len0,len,dir[MAXD],vect[MAXD];
	POINT *p,*p_nb;
	int i,j,dim = Dimension(node->interface);
	double ks = geom_set->ks;
	double kl = geom_set->kl;
	double m_s = geom_set->m_s;
	double m_l = geom_set->m_l;
	double lambda_s = geom_set->lambda_s;
	double lambda_l = geom_set->lambda_l;
	double payload;

	if (dim == 3)
	{
	    if (is_load_node(node))
	    {
	    	Front *front = geom_set->front;
	    	AF_PARAMS *af_params = (AF_PARAMS*)front->extra2;

	    	payload = af_params->payload;
	    }
	}

	for (i = 0; i < dim; ++i)
	{
	    x[*n][i] = Coords(node->posn)[i];
	    v[*n][i] = node->posn->vel[i];
	    f[*n][i] = 0.0;
	}
	for (c = node->out_curves; c && *c; ++c)
	{
	    b = (*c)->first;
	    len0 = bond_length0(b);
	    for (j = 0; j < dim; ++j)
	    {
		vect[j] = len0*((Coords(b->end)[j] - Coords(b->start)[j])
				/bond_length(b) - b->dir0[j]);
		if (dim == 3)
		{
		    if (is_load_node(node) == YES)
		    	f[*n][j]   += kl*vect[j]/payload;
		    else if (hsbdry_type(*c) == STRING_HSBDRY)
	    	    	f[*n][j]   += kl*vect[j]/m_s;
		    else if (hsbdry_type(*c) == MONO_COMP_HSBDRY)
	    	    	f[*n][j]   += ks*vect[j]/m_s;
		}
		else
		    f[*n][j]   += kl*vect[j]/m_l;
	    }
	}
	for (c = node->in_curves; c && *c; ++c)
	{
	    b = (*c)->last;
	    len0 = bond_length0(b);
	    for (j = 0; j < dim; ++j)
	    {
		vect[j] = len0*((Coords(b->end)[j] - Coords(b->start)[j])
				/bond_length(b) - b->dir0[j]);
		if (dim == 3)
		{
		    if (is_load_node(node) == YES)
		    	f[*n][j]   -= kl*vect[j]/payload;
		    else if (hsbdry_type(*c) == STRING_HSBDRY)
	    	    	f[*n][j]   -= kl*vect[j]/m_s;
		    else if (hsbdry_type(*c) == MONO_COMP_HSBDRY)
	    	    	f[*n][j]   -= ks*vect[j]/m_s;
		}
		else
		    f[*n][j]   -= kl*vect[j]/m_l;
	    }
	}
	if (dim == 3)
	{
	    BOND_TRI **btris;
	    TRI **tris;
	    int k,side,nt,ns;
	    SURFACE *out_surfs[10];

	    ns = 0;
	    for (c = node->out_curves; c && *c; ++c)
	    {
		b = (*c)->first;
		p = b->start;
		for (btris = Btris(b); btris && *btris; ++btris)
		{
		    out_surfs[ns++] = (*btris)->surface;
		    nt = I_FirstRingTrisAroundPoint(p,(*btris)->tri,&tris);
		    for (j = 0; j < nt; ++j)
		    {
			for (side = 0; side < 3; ++side)
			{
			    if (p == Point_of_tri(tris[j])[side])
			    {
				if (is_side_bdry(tris[j],side))
				    continue;
				p_nb = Point_of_tri(tris[j])[(side+1)%3];
				len0 = tris[j]->side_length0[side];
				len = separation(p,p_nb,3);
	    			x_diff = len - len0; 
				for (k = 0; k < 3; ++k)
                        	{
                            	    dir[k] = (Coords(p_nb)[k] - 
						Coords(p)[k])/len;
                            	    vect[k] = len0*(dir[k]
					- tris[j]->side_dir0[side][k]);
                            	    f[*n][k] += ks*vect[k]/m_s;
                        	}
			    }
			}
		    }
		}
	    }
	    for (c = node->in_curves; c && *c; ++c)
	    {
		if (is_closed_curve(*c)) continue;
		b = (*c)->last;
		p = b->end;
		for (btris = Btris(b); btris && *btris; ++btris)
		{
		    boolean duplicate_surf = NO;
		    for (j = 0; j < ns; ++j)
			if ((*btris)->surface == out_surfs[j])
			    duplicate_surf = YES;
		    if (duplicate_surf == YES) continue;
		    nt = I_FirstRingTrisAroundPoint(p,(*btris)->tri,&tris);
		    for (j = 0; j < nt; ++j)
		    {
			for (side = 0; side < 3; ++side)
			{
			    if (p == Point_of_tri(tris[j])[side])
			    {
				if (is_side_bdry(tris[j],side))
				    continue;
				p_nb = Point_of_tri(tris[j])[(side+1)%3];
				len0 = tris[j]->side_length0[side];
				len = separation(p,p_nb,3);
				x_diff = len - len0;
				for (k = 0; k < 3; ++k)
                        	{
                            	    dir[k] = (Coords(p_nb)[k] - 
						Coords(p)[k])/len;
                            	    vect[k] = len0*(dir[k]
					- tris[j]->side_dir0[side][k]);
                            	    f[*n][k] += ks*vect[k]/m_s;
                        	}
			    }
			}
		    }
		}
	    }
	    if (!is_load_node(node))
	    {
	    	for (i = 0; i < 3; ++i)
	    	    f[*n][i] -= lambda_s*v[*n][i]/m_s;
	    }
	}
	else
	{
	    for (i = 0; i < 3; ++i)
	    	f[*n][i] -= lambda_l*v[*n][i]/m_l;
	}
	(*n)++;
}	/* end compute_node_accel3 */

extern void compute_curve_accel3(
	ELASTIC_SET *geom_set,
	CURVE *curve,
	double **f,
	double **x,
	double **v,
	int *n)
{
	int i,j;
	double x_diff;
	BOND *b;
	double dir[MAXD],len0,vect[MAXD];
	int dim = Dimension(curve->interface);
	double kl,m_l,lambda_l;

	if (dim == 3)
	{
	    if (hsbdry_type(curve) == STRING_HSBDRY)
	    {
	    	kl = geom_set->kl;
	    	m_l = geom_set->m_l;
	    	lambda_l = geom_set->lambda_l;
	    }
	    else
	    {
	    	kl = geom_set->ks;
	    	m_l = geom_set->m_s;
	    	lambda_l = geom_set->lambda_s;
	    }
	}
	else
	{
	    kl = geom_set->kl;
	    m_l = geom_set->m_l;
	    lambda_l = geom_set->lambda_l;
	}
	i = *n;
	for (b = curve->first; b != curve->last; b = b->next)
	{
	    for (j = 0; j < dim; ++j)
	    {
	    	x[i][j] = Coords(b->end)[j];
	    	v[i][j] = b->end->vel[j];
		f[i][j] = -lambda_l*v[i][j]/m_l;
	    }
	    i++;
	}

	i = *n;
	for (b = curve->first; b != NULL; b = b->next)
	{
	    len0 = bond_length0(b);
	    for (j = 0; j < dim; ++j)
	    {
		vect[j] = len0*((Coords(b->end)[j] - Coords(b->start)[j])
				/bond_length(b) - b->dir0[j]);
		if (b != curve->first)
		{
	    	    f[i-1][j]   += kl*vect[j]/m_l;
		}
		if (b != curve->last)
		{
	    	    f[i][j] -= kl*vect[j]/m_l;
		}
	    }
	    if (b != curve->last) i++;
	}
	printf("f[20] = %f %f\n",f[20][0],f[20][1]);

	if (dim == 3)
	{
	    POINT *p,*p_nb;
	    BOND_TRI **btris;
	    TRI **tris;
	    int j,k,side,nt;
	    double length0,length;
	    i = *n;
	    for (b = curve->first; b != curve->last; b = b->next)
	    {
		p = b->end;
		for (btris = Btris(b); btris && *btris; ++btris)
		{
		    nt = I_FirstRingTrisAroundPoint(p,(*btris)->tri,&tris);
		    for (j = 0; j < nt; ++j)
		    {
			for (side = 0; side < 3; ++side)
			{
			    if (p == Point_of_tri(tris[j])[side])
			    {
				if (is_side_bdry(tris[j],side))
				    continue;
				p_nb = Point_of_tri(tris[j])[(side+1)%3];
				length0 = tris[j]->side_length0[side];
				length = separation(p,p_nb,3);
				for (k = 0; k < 3; ++k)
                        	{
                            	    dir[k] = (Coords(p_nb)[k] - 
						Coords(p)[k])/length;
				    vect[k] = length0*(dir[k]
					- tris[j]->side_dir0[side][k]);
                            	    f[i][k] += kl*vect[k]/m_l;
                        	}
			    }
			}
		    }
		}
		i++;
	    }
	}
	*n = i;
}	/* end compute_curve_accel3 */

static boolean curve_in_pointer_list(
	CURVE *c,
	CURVE **c_list)
{
	CURVE **pc;
	if (c_list == NULL) return NO;
	for (pc = c_list; pc && *pc; ++pc)
	{
	    if (c == *pc) return YES;
	}
	return NO;
}	/* end curve_in_pointer_list */

extern boolean is_registered_point(
	SURFACE *surf,
	POINT *p)
{
	REGISTERED_PTS *rgp = (REGISTERED_PTS*)surf->extra;
	int i,num_pts;
	int *global_ids;
	
	if (rgp == NULL) return NO;

	num_pts = rgp->num_pts;
	global_ids = rgp->global_ids;
	for (i = 0; i < num_pts; ++i)
	{
	    if (Gindex(p) == global_ids[i])
		return YES;
	}
	return NO;
}	/* end is_registered_point */

extern void propagate_surface(
        ELASTIC_SET *geom_set,
        SURFACE *surf,
        double **x,
        int *n)
{
        int i,j;
        TRI *tri;
        POINT *p;
        STATE *sl,*sr;
        HYPER_SURF_ELEMENT *hse;
        HYPER_SURF         *hs;
	double dt = geom_set->dt;
	Front *front = geom_set->front;
	IF_PARAMS *iFparams = (IF_PARAMS*)front->extra1;
	double *g = iFparams->gravity;

	hs = Hyper_surf(surf);
	unsort_surf_point(surf);
        for (tri = first_tri(surf); !at_end_of_tri_list(tri,surf);
                        tri = tri->next)
        {
            hse = Hyper_surf_element(tri);
            for (i = 0; i < 3; ++i)
            {
                p = Point_of_tri(tri)[i];
                if (sorted(p) || Boundary_point(p)) continue;
                FT_GetStatesAtPoint(p,hse,hs,(POINTER*)&sl,(POINTER*)&sr);
		if (is_registered_point(surf,p))
		{
                    for (j = 0; j < 3; ++j)
                    {
                        x[*n][j] += sl->impulse[j]*dt;
                        sr->impulse[j] = sl->impulse[j] = sl->impulse[j];
                    }
		}
		else
                {
                    for (j = 0; j < 3; ++j)
                    {
                        x[*n][j] += (sl->impulse[j] + 0.5*g[j]*dt)*dt;
                        sr->impulse[j] = sl->impulse[j] = 
					sl->impulse[j] + g[j]*dt;
                    }
                }
                sorted(p) = YES;
                ++(*n);
            }
        }
}       /* propagate_surface */

extern void propagate_node(
        ELASTIC_SET *geom_set,
	NODE *node,
        double **x,
        int *n)
{
        int i,j;
        POINT *p;
        STATE *sl,*sr;
	double dt = geom_set->dt;
	Front *front = geom_set->front;
	IF_PARAMS *iFparams = (IF_PARAMS*)front->extra1;
	double *g = iFparams->gravity;
	int dim = front->rect_grid->dim;

        sl = (STATE*)left_state(node->posn);
        sr = (STATE*)right_state(node->posn);
        for (j = 0; j < dim; ++j)
        {
            x[*n][j] += (sl->impulse[j] + 0.5*g[j]*dt)*dt;
            sr->impulse[j] = sl->impulse[j] = sl->impulse[j] + g[j]*dt;
        }
        ++(*n);
}	/* end propagate_node */

extern void propagate_curve(
        ELASTIC_SET *geom_set,
	CURVE *curve,
        double **x,
        int *n)
{
        int i,j;
        POINT *p;
	BOND *b;
        STATE *sl,*sr;
	double dt = geom_set->dt;
	Front *front = geom_set->front;
	IF_PARAMS *iFparams = (IF_PARAMS*)front->extra1;
	double *g = iFparams->gravity;
	int dim = front->rect_grid->dim;

	for (b = curve->first; b != curve->last; b = b->next)
        {
            p = b->end;
            sl = (STATE*)left_state(p);
            sr = (STATE*)right_state(p);
            for (j = 0; j < dim; ++j)
            {
                x[*n][j] += (sl->impulse[j] + 0.5*g[j]*dt)*dt;
                sr->impulse[j] = sl->impulse[j] = sl->impulse[j] + g[j]*dt;
            }
            ++(*n);
        }
}	/* end propagate_curve */

static void reduce_high_freq_vel(
	Front *front,
	SURFACE *canopy)
{
	AF_PARAMS *af_params = (AF_PARAMS*)front->extra2;
	double **vv;
	int ncan;
	POINT *p;
	TRI *tri;
	HYPER_SURF *hs;
	HYPER_SURF_ELEMENT *hse;
	STATE *sl,*sr;
	double max_speed;
	double crds_max[MAXD];
	int i,j,gindex_max;
	int l,num_layers = af_params->num_smooth_layers;

	hs = Hyper_surf(canopy);

	ncan = 0;
	unsort_surf_point(canopy);
	for (tri = first_tri(canopy); !at_end_of_tri_list(tri,canopy); 
			tri = tri->next)
	{
	    hse = Hyper_surf_element(tri);
	    for (i = 0; i < 3; ++i)
	    {
		p = Point_of_tri(tri)[i];
		if (sorted(p) || Boundary_point(p)) continue;
		sorted(p) = YES;
		ncan++;
	    }
	}

	FT_MatrixMemoryAlloc((POINTER*)&vv,ncan,3,sizeof(double));

	l = 0;
	while (l < num_layers)
	{
	    ncan = 0;
	    unsort_surf_point(canopy);
	    for (tri = first_tri(canopy); !at_end_of_tri_list(tri,canopy); 
			tri = tri->next)
	    {
	    	for (i = 0; i < 3; ++i)
	    	{
		    p = Point_of_tri(tri)[i];
		    if (sorted(p) || Boundary_point(p)) continue;
		    smooth_vel(vv[ncan],p,tri,canopy);
		    ncan++;
		    sorted(p) = YES;
	    	}
	    }

	    ncan = 0;
	    max_speed = 0.0;
	    unsort_surf_point(canopy);
	    for (tri = first_tri(canopy); !at_end_of_tri_list(tri,canopy); 
			tri = tri->next)
	    {
	    	hse = Hyper_surf_element(tri);
	    	for (i = 0; i < 3; ++i)
	    	{
		    p = Point_of_tri(tri)[i];
		    if (sorted(p) || Boundary_point(p)) continue;
		    FT_GetStatesAtPoint(p,hse,hs,(POINTER*)&sl,(POINTER*)&sr);
		    for (j = 0; j < 3; ++j)
		    {
		    	sl->vel[j] = vv[ncan][j];
		    	sr->vel[j] = vv[ncan][j];
		    }
		    if (max_speed < Mag3d(sl->vel)) 
		    	max_speed = Mag3d(sl->vel);
		    ncan++;
		    sorted(p) = YES;
	    	}
	    }
	    if (debugging("smooth_canopy_vel"))
		(void) printf("Max speed after smoothing round %d: %f\n",
					l,max_speed);
	    l++;
	}

	unsort_surf_point(canopy);
	max_speed = 0.0;
	for (tri = first_tri(canopy); !at_end_of_tri_list(tri,canopy); 
			tri = tri->next)
	{
	    hse = Hyper_surf_element(tri);
	    for (i = 0; i < 3; ++i)
	    {
		p = Point_of_tri(tri)[i];
		if (sorted(p) || Boundary_point(p)) continue;
		FT_GetStatesAtPoint(p,hse,hs,(POINTER*)&sl,(POINTER*)&sr);
		for (j = 0; j < 3; ++j)
	    	    FT_RecordMaxFrontSpeed(j,sl->vel[j],NULL,Coords(p),front);
	    	FT_RecordMaxFrontSpeed(3,Mag3d(sl->vel),NULL,Coords(p),front);
		if (max_speed < Mag3d(sl->vel)) 
		{
		    	max_speed = Mag3d(sl->vel);
		    	gindex_max = Gindex(p);
		    	for (j = 0; j < 3; ++j)
			    crds_max[j] = Coords(p)[j];
		}
		sorted(p) = YES;
	    }
	}
	FT_FreeThese(1,vv);
}	/* end reduce_high_freq_vel */

static void print_elastic_params(
	ELASTIC_SET geom_set)
{
	int i;
	double *spfr;
	Front *fr = geom_set.front;

	spfr = Spfr(fr);
        for (i = 0; i <= 3; ++i)
            printf("Max front speed(%d) = %f\n",i,spfr[i]);
        (void) printf("Input surface parameters:\n");
        (void) printf("ks = %f  m_s = %f  lambda_s = %f\n",
                    	geom_set.ks,
                        geom_set.m_s,
                        geom_set.lambda_s);
        (void) printf("Input string parameters:\n");
        (void) printf("kl = %f  m_l = %f  lambda_l = %f\n",
                        geom_set.kl,
                        geom_set.m_l,
                        geom_set.lambda_l);
        (void) printf("Input gore parameters:\n");
        (void) printf("kg = %f  m_g = %f  lambda_g = %f\n",
                        geom_set.kg,
                        geom_set.m_g,
                        geom_set.lambda_g);
	(void) printf("\ndt_tol = %20.14f  dt = %20.14f\n",
                        geom_set.dt_tol,geom_set.dt);
}	/* end print_elastic_params */

extern void fourth_order_elastic_set_propagate(
	Front           *fr,
        double           fr_dt)
{
	static ELASTIC_SET geom_set;
	static int size = 0,owner_size,client_size;
	static int *client_size_old, *client_size_new;
        AF_PARAMS *af_params = (AF_PARAMS*)fr->extra2;
        int i,j,k,n_sub;
        double dt;
        static SPRING_VERTEX *sv;
        static boolean first = YES;
        static GLOBAL_POINT **point_set;
        static GLOBAL_POINT *point_set_store;
	static GLOBAL_POINT **client_point_set_store;
        int dim = FT_Dimension();
        long max_point_gindex = fr->interf->max_point_gindex;
	int owner[MAXD];
	int owner_id = af_params->node_id[0];
        int myid = pp_mynode();
	int gindex;
        INTERFACE *elastic_intfc = NULL;
	double *L = fr->rect_grid->L;
	double *U = fr->rect_grid->U;
	double client_L[MAXD],client_U[MAXD];
	static boolean first_break_strings = YES;
	static double break_strings_time = af_params->break_strings_time;
	static int break_strings_num = af_params->break_strings_num;
#ifdef __COLLISION__
	static CollisionSolver* collision_solver = new CollisionSolver3d();
#endif

	if (debugging("trace"))
	    (void) printf("Entering fourth_order_elastic_set_propagate()\n");
	geom_set.front = fr;

	if (first_break_strings && break_strings_num > 0 &&
	    break_strings_time >= 0.0 && 
	    fr->time + fr->dt >= break_strings_time)
	{
	    printf("Some strings break! Count and set spring vertex again.\n");
	    first_break_strings = NO;
	    first = YES;
	}
	if (first)
        {
            set_elastic_params(&geom_set,fr_dt);
            if (debugging("step_size"))
                print_elastic_params(geom_set);
        }
        if (fr_dt > geom_set.dt_tol)
        {
            n_sub = (int)(fr_dt/geom_set.dt_tol);
            dt = fr_dt/n_sub;
        }
	else
        {
            n_sub = af_params->n_sub;
            dt = fr_dt/n_sub;
        }
	if (first)
	{
            owner[0] = 0;
            owner[1] = 0;
            owner[2] = 0;
	    if (point_set != NULL)
		FT_FreeThese(1, point_set);
	    FT_VectorMemoryAlloc((POINTER*)&point_set,max_point_gindex,
					sizeof(GLOBAL_POINT*));
	    for (i = 0; i < max_point_gindex; ++i)
		point_set[i] = NULL;

	    if (pp_numnodes() > 1)
	    {
            	elastic_intfc = FT_CollectHypersurfFromSubdomains(fr,owner,
				ELASTIC_BOUNDARY);
		collectNodeExtra(fr,elastic_intfc,owner_id);
	    }
	    else
		elastic_intfc = fr->interf;
	    start_clock("set_data");
	    if (myid == owner_id)
            {
		if (client_size_old != NULL)
		    FT_FreeThese(3, client_size_old, client_size_new, 
					client_point_set_store);
		FT_VectorMemoryAlloc((POINTER*)&client_size_old,pp_numnodes(),
                                        sizeof(int));
                FT_VectorMemoryAlloc((POINTER*)&client_size_new,pp_numnodes(),
                                        sizeof(int));
                FT_VectorMemoryAlloc((POINTER*)&client_point_set_store,
                                        pp_numnodes(),sizeof(GLOBAL_POINT*));
                for (i = 0; i < pp_numnodes(); i++)
                    client_size_old[i] = client_size_new[i] = 0;

		assembleParachuteSet(elastic_intfc,&geom_set);
		owner_size = geom_set.num_verts;
		if (point_set_store != NULL) 
		    FT_FreeThese(2,point_set_store, sv);
		FT_VectorMemoryAlloc((POINTER*)&point_set_store,owner_size,
                                        sizeof(GLOBAL_POINT));
                FT_VectorMemoryAlloc((POINTER*)&sv,owner_size,
                                        sizeof(SPRING_VERTEX));
		link_point_set(&geom_set,point_set,point_set_store);
	    	count_vertex_neighbors(&geom_set,sv);
	    	set_spring_vertex_memory(sv,owner_size);
	    	set_vertex_neighbors(&geom_set,sv,point_set);
		if (elastic_intfc != fr->interf)
		    delete_interface(elastic_intfc);
	    }
	    stop_clock("set_data");
	    first = NO;
	}

	elastic_intfc = fr->interf;
	assembleParachuteSet(elastic_intfc,&geom_set);
	if (myid != owner_id)
	{
	    client_size = geom_set.num_verts;
	    if (size < client_size)
	    {
	    	size = client_size;
	    	if (point_set_store != NULL)
		{
		    FT_FreeThese(2,point_set_store,sv);
		}
	    	FT_VectorMemoryAlloc((POINTER*)&point_set_store,size,
                                        sizeof(GLOBAL_POINT));
                FT_VectorMemoryAlloc((POINTER*)&sv,size,sizeof(SPRING_VERTEX));
	    }
	    for (i = 0; i < max_point_gindex; ++i)
                point_set[i] = NULL;
	    link_point_set(&geom_set,point_set,point_set_store);
	    count_vertex_neighbors(&geom_set,sv);
	    set_spring_vertex_memory(sv,client_size);
	    set_vertex_neighbors(&geom_set,sv,point_set);
	    get_point_set_from(&geom_set,point_set);
	    pp_send(5,L,MAXD*sizeof(double),owner_id);
	    pp_send(6,U,MAXD*sizeof(double),owner_id);
	    pp_send(1,&(client_size),sizeof(int),owner_id);
            pp_send(2,point_set_store,client_size*sizeof(GLOBAL_POINT),
					owner_id);
	}
	else
	    size = owner_size;

	if (myid == owner_id)
	{
#ifdef __COLLISION__
	    if (FT_Dimension() == 3)
            {
                // collision setup
                setCollisionFreePoints3d(fr->interf);
                collision_solver->assembleFromInterface(fr->interf,fr->dt);
                collision_solver->recordOriginPosition();
                collision_solver->setSpringConstant(af_params->ks);
                collision_solver->setFrictionConstant(0.0);
                collision_solver->setPointMass(af_params->m_s);
                collision_solver->setFabricThickness(1.0e-4);
		collision_solver->setRestitutionCoef(0.0);
            }
#endif

	    get_point_set_from(&geom_set,point_set);
	    for (i = 0; i < pp_numnodes(); i++)
	    {
		if (i == myid) continue;
		pp_recv(5,i,client_L,MAXD*sizeof(double));
		pp_recv(6,i,client_U,MAXD*sizeof(double));
		pp_recv(1,i,client_size_new+i,sizeof(int));
		if (client_size_new[i] > client_size_old[i])
		{
		    client_size_old[i] = client_size_new[i];
		    if (client_point_set_store[i] != NULL)
		    	FT_FreeThese(1,client_point_set_store[i]);
	    	    FT_VectorMemoryAlloc((POINTER*)&client_point_set_store[i],
				client_size_new[i], sizeof(GLOBAL_POINT));
		}
		pp_recv(2,i,client_point_set_store[i],
		    client_size_new[i]*sizeof(GLOBAL_POINT));
		copy_from_client_point_set(point_set,client_point_set_store[i],
				client_size_new[i],client_L,client_U);
	    } 

	    start_clock("spring_model");
#if defined(__GPU__)
            if (af_params->use_gpu)
            {
            	if (debugging("trace"))
                    (void) printf("Enter gpu_spring_solver()\n");
            	gpu_spring_solver(sv,dim,size,n_sub,dt);
            	if (debugging("trace"))
                    (void) printf("Left gpu_spring_solver()\n");
            }
            else
#endif
            	generic_spring_solver(sv,dim,size,n_sub,dt);
	    stop_clock("spring_model");

	    for (i = 0; i < pp_numnodes(); i++)
	    {
		if (i == myid) continue;
		copy_to_client_point_set(point_set,client_point_set_store[i],
				client_size_new[i]);
		pp_send(3,client_point_set_store[i],
                        client_size_new[i]*sizeof(GLOBAL_POINT),i);
	    }
	}
	if (myid != owner_id)
        {
            pp_recv(3,owner_id,point_set_store,
				client_size*sizeof(GLOBAL_POINT));
        }
	/* Owner send and patch point_set_store from other processors */	
	put_point_set_to(&geom_set,point_set);
	/* Calculate the real force on load_node and rg_string_node */
	setSpecialNodeForce(fr, geom_set.kl);

	set_vertex_impulse(&geom_set,point_set);
	set_geomset_velocity(&geom_set,point_set);
	compute_center_of_mass_velo(&geom_set);

#ifdef __COLLISION__
	if (myid == owner_id)
        {
            if (FT_Dimension() == 3)
                // resolve collision
                collision_solver->resolveCollision();
        }
#endif

	if (debugging("trace"))
	    (void) printf("Leaving fourth_order_elastic_set_propagate()\n");
}	/* end fourth_order_elastic_set_propagate() */

static void setSurfVelocity(
	ELASTIC_SET *geom_set,
	SURFACE *surf,
	GLOBAL_POINT **point_set)
{
	int i,j;
	TRI *tri;
	POINT *p;
	STATE *sl,*sr;
	HYPER_SURF_ELEMENT *hse;
        HYPER_SURF         *hs;
	Front *front = geom_set->front;
	double nor[MAXD],nor_speed;
	double *vel;
	int gindex_max;
	long gindex;

	unsort_surf_point(surf);
	hs = Hyper_surf(surf);
	for (tri = first_tri(surf); !at_end_of_tri_list(tri,surf); 
			tri = tri->next)
	{
	    hse = Hyper_surf_element(tri);
	    for (i = 0; i < 3; ++i)
	    {
		p = Point_of_tri(tri)[i];
		if (sorted(p) || Boundary_point(p)) continue;
		gindex = Gindex(p);
		FT_NormalAtPoint(p,front,nor,NO_COMP);
		FT_GetStatesAtPoint(p,hse,hs,(POINTER*)&sl,(POINTER*)&sr);
		vel = point_set[gindex]->v;
		nor_speed = scalar_product(vel,nor,3);
		for (j = 0; j < 3; ++j)
		{
		    sl->vel[j] = nor_speed*nor[j];
		    sr->vel[j] = nor_speed*nor[j];
		}
		sorted(p) = YES;
	    }
	}
	//reduce_high_freq_vel(front,surf);
}	/* end setSurfVelocity */

static void setCurveVelocity(
	ELASTIC_SET *geom_set,
	CURVE *curve,
	GLOBAL_POINT **point_set)
{
	int i,j;
	BOND *b;
	POINT *p;
	BOND_TRI **btris;
	STATE *sl,*sr;
	HYPER_SURF_ELEMENT *hse;
        HYPER_SURF         *hs;
	Front *front = geom_set->front;
	double nor[MAXD],nor_speed;
	double *vel;
	double crds_max[MAXD];
	int gindex_max;
	long gindex;
	int dim = FT_Dimension();

	for (b = curve->first; b != curve->last; b = b->next)
        {
            p = b->end;
	    for (btris = Btris(b); btris && *btris; ++btris)
            {
                p->hse = hse = Hyper_surf_element((*btris)->tri);
                p->hs = hs = Hyper_surf((*btris)->surface);
		gindex = Gindex(p);
                FT_GetStatesAtPoint(p,hse,hs,(POINTER*)&sl,(POINTER*)&sr);
		FT_NormalAtPoint(p,front,nor,NO_COMP);
		vel = point_set[gindex]->v;
		nor_speed = scalar_product(vel,nor,3);
                for (j = 0; j < 3; ++j)
		    sl->vel[j] = sr->vel[j] = nor_speed*nor[j];
            }
        }
	for (b = curve->first; b != NULL; b = b->next)
	    set_bond_length(b,dim);
}	/* end setCurveVelocity */

static void setNodeVelocity(
	ELASTIC_SET *geom_set,
	NODE *node,
	GLOBAL_POINT **point_set)
{
	int dim = FT_Dimension();
	switch(dim)
	{
	case 2:
	    new_setNodeVelocity2d(geom_set,node,point_set);
	    return;
	case 3:
	    new_setNodeVelocity3d(geom_set,node,point_set);
	    return;
	}
}	/* end setNodeVelocity */

static void new_setNodeVelocity2d(
	ELASTIC_SET *geom_set,
	NODE *node,
	GLOBAL_POINT **point_set)
{
	int i,j;
	BOND *b;
	POINT *p;
	STATE *sl,*sr;
	double *vel;
	long gindex;

	if (is_load_node(node))
	{
	    sl = (STATE*)left_state(node->posn);
            sr = (STATE*)right_state(node->posn);
	    gindex = Gindex(node->posn);
	    vel = point_set[gindex]->v;
            for (j = 0; j < 3; ++j)
            {
            	sl->vel[j] = vel[j];
            	sr->vel[j] = vel[j];
            }
	}
}	/* end setNodeVelocity2d */

static void new_setNodeVelocity3d(
	ELASTIC_SET *geom_set,
	NODE *node,
	GLOBAL_POINT **point_set)
{
	int i,j;
	BOND *b;
	POINT *p;
	BOND_TRI **btris;
	STATE *sl,*sr;
	HYPER_SURF_ELEMENT *hse;
        HYPER_SURF         *hs;
	Front *front = geom_set->front;
	CURVE **c;
	double nor[MAXD],nor_speed,max_speed;
	double *vel;
	double crds_max[MAXD];
	int gindex_max;
	long gindex;

	for (c = node->out_curves; c && *c; ++c)
        {
		if (hsbdry_type(*c) != MONO_COMP_HSBDRY &&
		    hsbdry_type(*c) != GORE_HSBDRY && 
		    hsbdry_type(*c) != PASSIVE_HSBDRY) 
		    continue;
                b = (*c)->first;
                p = b->start;
		for (btris = Btris(b); btris && *btris; ++btris)
		{
                    p->hse = hse = Hyper_surf_element((*btris)->tri);
                    p->hs = hs = Hyper_surf((*btris)->surface);
                    FT_GetStatesAtPoint(p,hse,hs,(POINTER*)&sl,(POINTER*)&sr);
		    FT_NormalAtPoint(p,front,nor,NO_COMP);
		    gindex = Gindex(p);
		    vel = point_set[gindex]->v;
		    if (hsbdry_type(*c) == PASSIVE_HSBDRY)
		    {
			for (j = 0; j < 3; ++j)
			    sl->vel[j] = sr->vel[j] = vel[j];
			continue;
		    }
		    nor_speed = scalar_product(vel,nor,3);
		    if (max_speed < fabs(nor_speed)) 
		    {
		    	max_speed = fabs(nor_speed);
		    	gindex_max = Gindex(p);
		    	for (j = 0; j < 3; ++j)
			    crds_max[j] = Coords(p)[j];
		    }
                    for (j = 0; j < 3; ++j)
		    	sl->vel[j] = sr->vel[j] =  nor_speed*nor[j];
		}
        }
        for (c = node->in_curves; c && *c; ++c)
        {
		if (hsbdry_type(*c) != MONO_COMP_HSBDRY &&
		    hsbdry_type(*c) != GORE_HSBDRY && 
		    hsbdry_type(*c) != PASSIVE_HSBDRY) 
		    continue;
                b = (*c)->last;
                p = b->end;
		for (btris = Btris(b); btris && *btris; ++btris)
		{
                    p->hse = hse = Hyper_surf_element((*btris)->tri);
                    p->hs = hs = Hyper_surf((*btris)->surface);
                    FT_GetStatesAtPoint(p,hse,hs,(POINTER*)&sl,(POINTER*)&sr);
		    FT_NormalAtPoint(p,front,nor,NO_COMP);
		    gindex = Gindex(p);
		    vel = point_set[gindex]->v;
		    if (hsbdry_type(*c) == PASSIVE_HSBDRY)
		    {
			for (j = 0; j < 3; ++j)
			    sl->vel[j] = sr->vel[j] = vel[j];
			continue;
		    }
		    nor_speed = scalar_product(vel,nor,3);
		    if (max_speed < fabs(nor_speed)) 
		    {
		    	max_speed = fabs(nor_speed);
		    	gindex_max = Gindex(p);
		    	for (j = 0; j < 3; ++j)
			    crds_max[j] = Coords(p)[j];
		    }
                    for (j = 0; j < 3; ++j)
		    	sl->vel[j] = sr->vel[j] = nor_speed*nor[j];
		}
        }
}	/* end setNodeVelocity3d */

extern void set_geomset_velocity(
	ELASTIC_SET *geom_set,
	GLOBAL_POINT **point_set)
{
	int i,ns,nc,nn;

	ns = geom_set->num_surfs;
	nc = geom_set->num_curves;
	nn = geom_set->num_nodes;
	for (i = 0; i < ns; ++i)
	    setSurfVelocity(geom_set,geom_set->surfs[i],point_set);
	for (i = 0; i < nc; ++i)
	    setCurveVelocity(geom_set,geom_set->curves[i],point_set);
	for (i = 0; i < nn; ++i)
	{
	    if (is_load_node(geom_set->nodes[i])) continue;
	    setNodeVelocity(geom_set,geom_set->nodes[i],point_set);
	}

}	/* end set_geomset_velocity */

extern void collectNodeExtra(
	Front *front,
	INTERFACE *host_intfc,
	int owner_id)
{
	NODE **n;
	NODE **node_with_extra;
	int i,j,k,num_nodes;
	INTERFACE *intfc = front->interf;
	long global_index;
	AF_NODE_EXTRA extra_recv, *extra;
	RECT_GRID *gr = front->rect_grid;
        double *L = gr->L;
        double *U = gr->U;
        int dim = gr->dim;

	if (pp_mynode() != owner_id)
	{
	    num_nodes = 0;
	    intfc_node_loop(intfc,n)
	    {
		for (k = 0; k < dim; ++k)
		{
		    if (Coords((*n)->posn)[k] <= L[k] ||
		        Coords((*n)->posn)[k] > U[k])
		        break;
		}
		if (k != dim || (*n)->extra == NULL) continue;
		num_nodes++;
	    }
	    pp_send(10,&num_nodes,sizeof(int),owner_id);
	    intfc_node_loop(intfc,n)
	    {
		for (k = 0; k < dim; ++k)
		{
		    if (Coords((*n)->posn)[k] <= L[k] ||
		        Coords((*n)->posn)[k] > U[k])
		        break;
		}
		if (k != dim || (*n)->extra == NULL) continue;
	    	pp_send(11,&(Gindex((*n)->posn)),sizeof(long),owner_id);
	    	pp_send(12,(*n)->extra,sizeof(AF_NODE_EXTRA),owner_id);
	    }
	}
	else
	{
	    for (i = 0; i < pp_numnodes(); ++i)
	    {
		if (i == owner_id) continue;
		pp_recv(10,i,&num_nodes,sizeof(int));
		for (j = 0; j < num_nodes; ++j)
		{
		    pp_recv(11,i,&global_index,sizeof(long));
		    pp_recv(12,i,&extra_recv,sizeof(AF_NODE_EXTRA));
		    intfc_node_loop(host_intfc,n)
		    {
			if (Gindex((*n)->posn) != global_index)
			    continue;
		        FT_ScalarMemoryAlloc((POINTER*)&extra,
					sizeof(AF_NODE_EXTRA));
			*extra = extra_recv;
			(*n)->extra = (POINTER)extra;
			(*n)->size_of_extra = sizeof(AF_NODE_EXTRA);
		    }
		}
	    }
	}
}	/* end collectNodeExtra */

static void setCollisionFreePoints3d(INTERFACE* intfc)
{
        POINT *p;
        HYPER_SURF *hs;
        HYPER_SURF_ELEMENT *hse;
        SURFACE* surf;
        if (intfc->dim == 2){
            printf("ERROR dim = %d\n",intfc->dim);
            clean_up(ERROR);
        }
        next_point(intfc,NULL,NULL,NULL);
        while(next_point(intfc,&p,&hse,&hs)){
            STATE* sl = (STATE*)left_state(p);
            sl->is_fixed = false;
	    sl->is_movableRG = false;
            if ((surf = Surface_of_hs(hs)) &&
                (is_registered_point(surf,p) ||
                 wave_type(hs) == NEUMANN_BOUNDARY))
            {
                sl->is_fixed = true;
            }
	    if ((surf = Surface_of_hs(hs)) &&
		(wave_type(hs) == MOVABLE_BODY_BOUNDARY))
	    {
		sl->is_movableRG = true;
	    }
        }

        CURVE **c;
        BOND* b;
        intfc_curve_loop(intfc,c){
            if (hsbdry_type(*c) != FIXED_HSBDRY)
                continue;
            for (b = (*c)->first; b != (*c)->last; b = b->next)
            {
                STATE* sl = (STATE*)left_state(b->end);
                sl->is_fixed = true;
            }
        }

        NODE** n;
        intfc_node_loop(intfc,n){
            STATE* sl = (STATE*)left_state((*n)->posn);
            sl->is_fixed = false;
            AF_NODE_EXTRA* extra;
            if ((extra = (AF_NODE_EXTRA*)(*n)->extra) &&
                (extra->af_node_type == PRESET_NODE ||
                 is_load_node(*n) ||
                 is_rg_string_node(*n)))
            {
                sl->is_fixed = true;
            }
            else if ((*n)->hsb && is_fixed_node(*n))
            {
                sl->is_fixed = true;
            }
        }
}       /* setCollisionFreePoints3d() */

extern void scatterAirfoilExtra(
	Front *front)
{
	NODE **n;
	NODE **node_with_extra;
	int i,j,k,num_nodes;
	INTERFACE *intfc = front->interf;
	long global_index;
	AF_NODE_EXTRA extra_recv,*extra;
	RECT_GRID *gr = front->rect_grid;
	double *L = gr->L;
	double *U = gr->U;
	int dim = gr->dim;

	num_nodes = 0;
        intfc_node_loop(intfc,n)
        {
            for (k = 0; k < dim; ++k)
            {
                if (Coords((*n)->posn)[k] <= L[k] ||
                    Coords((*n)->posn)[k] > U[k])
                    break;
            }
            if (k != dim || (*n)->extra == NULL) continue;
            num_nodes++;
        }
        for (i = 0; i < pp_numnodes(); ++i)
        {
            if (i == pp_mynode()) continue;
            pp_send(30,&num_nodes,sizeof(int),i);
            intfc_node_loop(intfc,n)
            {
                for (k = 0; k < dim; ++k)
                {
                    if (Coords((*n)->posn)[k] <= L[k] ||
                        Coords((*n)->posn)[k] > U[k])
                        break;
                }
                if (k != dim || (*n)->extra == NULL) continue;
                pp_send(31,&(Gindex((*n)->posn)),sizeof(long),i);
                pp_send(32,(*n)->extra,sizeof(AF_NODE_EXTRA),i);
            }
        }
	pp_gsync();
        for (i = 0; i < pp_numnodes(); ++i)
        {
            if (i == pp_mynode()) continue;
            pp_recv(30,i,&num_nodes,sizeof(int));
            for (j = 0; j < num_nodes; ++j)
            {
                pp_recv(31,i,&global_index,sizeof(long));
                pp_recv(32,i,&extra_recv,sizeof(AF_NODE_EXTRA));
                intfc_node_loop(intfc,n)
                {
                    if (Gindex((*n)->posn) != global_index)
                        continue;
                    FT_ScalarMemoryAlloc((POINTER*)&extra,
                                sizeof(AF_NODE_EXTRA));
                    *extra = extra_recv;
                    (*n)->extra = (POINTER)extra;
                    (*n)->size_of_extra = sizeof(AF_NODE_EXTRA);
                }
            }
	}
}	/* end scatterAirfoilExtra */

extern void setSpecialNodeForce(
	Front *front,
	double kl)
{
	INTERFACE *intfc = front->interf;
	int i, k;
	double f[MAXD], vec[MAXD];
	NODE **n;
	CURVE **c;
	BOND *b;
	RECT_GRID *gr = front->rect_grid;
	double *L = gr->L;
	double *U = gr->U;
	int dim = gr->dim;
	
	if (debugging("trace"))
	    printf("Entering setSpecialNodeForce() \n");

	intfc_node_loop(intfc, n)
	{
	    if ( (!is_load_node(*n)) && (!is_rg_string_node(*n)) )
		continue;
            for (k = 0; k < dim; ++k)
            {
                if (Coords((*n)->posn)[k] <= L[k] ||
                    Coords((*n)->posn)[k] > U[k])
                    break;
            }
            if (k != dim) continue;

	    for (i = 0; i < dim; ++i)
		f[i] = 0.0;
	    node_out_curve_loop(*n,c)
	    {
		b = (*c)->first;
		set_bond_length(b,dim);
	    }
	    node_in_curve_loop(*n,c)
	    {
		b = (*c)->last;
		set_bond_length(b,dim);
	    }
	    node_out_curve_loop(*n,c)
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
	    node_in_curve_loop(*n,c)
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
	    for (i = 0; i < dim; ++i)
		(*n)->posn->force[i] = f[i];

	    if (debugging("rigid_body"))
	    {
		printf("special node coords = %f %f %f \n", 
			Coords((*n)->posn)[0], Coords((*n)->posn)[1], 
			Coords((*n)->posn)[2]);
		printf("force on the node = %f %f %f \n", f[0], f[1], f[2]);
		printf("velo of the node = %f %f %f \n", (*n)->posn->vel[0], 
				(*n)->posn->vel[1], (*n)->posn->vel[2]);
	    }
	}

	if (debugging("trace"))
	    printf("Leaving setSpecialNodeForce() \n");
}	/* end setSpecialNodeForce */

extern void break_strings(Front *front)
{
        INTERFACE *intfc = front->interf;
	AF_PARAMS *af_params = (AF_PARAMS*)front->extra2;
        static boolean first = YES;
	static double break_strings_time = af_params->break_strings_time;
	static int break_strings_num = af_params->break_strings_num;
	static int *bs_gindex = af_params->break_strings_gindex;
        CURVE **c, **curves;
        int i = 0;

	if (FT_Dimension() != 3) return;
	if (break_strings_time < 0.0 || break_strings_num <= 0) return;
        if ((!first) || (front->time + front->dt < break_strings_time)) return;

	if (debugging("trace"))
	    printf("Entering break_strings() \n");

        first = NO;
	double z_max = -HUGE, z_min = HUGE;
	intfc_curve_loop(intfc, c)
	{
	    if (hsbdry_type(*c) != STRING_HSBDRY) continue;
	    if (z_max < Coords((*c)->start->posn)[2])
		z_max = Coords((*c)->start->posn)[2];
	    if (z_max < Coords((*c)->end->posn)[2])
		z_max = Coords((*c)->end->posn)[2];
	    if (z_min > Coords((*c)->start->posn)[2])
		z_min = Coords((*c)->start->posn)[2];
	    if (z_min > Coords((*c)->end->posn)[2])
		z_min = Coords((*c)->end->posn)[2];
	}
	pp_gsync();
	pp_global_max(&z_max, 1);
	pp_global_min(&z_min, 1);
	double break_height = 0.5 * (z_max + z_min);
	for (int i = 0; i < break_strings_num; ++i)
	{
	    intfc_curve_loop(intfc, c)
	    {
		if (Gindex(*c) == bs_gindex[i]) 
		    break;
	    }
	    if (c && *c)
		break_string_curve(*c, break_height);
	    pp_gsync();
	    pp_global_lmax(&(current_interface()->max_point_gindex),1);
	    pp_global_lmax(&(current_interface()->max_curve_gindex),1);
	    exchange_curve_gindex(front);
	    scatter_front(front);
	}

	if (debugging("trace"))
	    printf("Leaving break_strings() \n");
}	/* end break_strings */

static void break_string_curve(CURVE* c, double height)
{
	static CURVE *curves[2] = {NULL, NULL};

	if (debugging("trace"))
	    printf("Entering split_string_curve()\n");
	INTERFACE *cur_intfc = current_interface();
	RECT_GRID *gr = computational_grid(cur_intfc);
	BOND *b = NULL;

	set_current_interface(c->interface);
	curve_bond_loop(c, b)
	{
	    double tmp = (Coords(b->start)[2] - height) * 
				(Coords(b->end)[2] - height);
	    if (tmp <= 0.0) break;
	}
	if (b == c->first) b = b->next;
	if (b == NULL || c->num_points <= 2)
	{
	    set_current_interface(cur_intfc);
	    return;
	}

	POINT *p[2];
	p[0] = b->start;
	p[1] = Point(Coords(p[0]));
	/* for consistency and uniqueness of global index */
	Gindex(p[1]) = cur_intfc->max_point_gindex + 1;
	cur_intfc->max_point_gindex += 1;

	NODE *n[2];
	n[0] = make_node(p[0]);
	n[1] = make_node(p[1]);
	if (n[0] == NULL || n[1] == NULL)
	{
	    printf("ERROR: split_string_curve returning NULL, "
	    	"cannot split at point p\n");
	    clean_up(ERROR);
	}
	curves[0] = make_curve(0,0,c->start,n[0]);
	curves[1] = make_curve(0,0,n[1],c->end);
	if (curves[0] == NULL || curves[1] == NULL)
	{
	    printf("ERROR: split_string_curve returning NULL, "
	    	"cannot make curve\n");
	    clean_up(ERROR);
	}
	curves[0]->first = c->first;
	curves[0]->last = b->prev;
	curves[0]->first->prev = NULL;
	curves[0]->last->next = NULL;
	curves[1]->first = b;
	curves[1]->last = c->last;
	curves[1]->first->prev = NULL;
	curves[1]->last->next = NULL;

	/* for consistency and uniqueness of global index */
	Gindex(curves[1]) = Gindex(c);
	Gindex(curves[0]) = cur_intfc->max_curve_gindex + 1;
	cur_intfc->max_curve_gindex += 1;

	hsbdry_type(curves[0]) = hsbdry_type(c);
	hsbdry_type(curves[1]) = hsbdry_type(c);
	curve_tangent_function(curves[0]) = curve_tangent_function(c);
	curve_tangent_function(curves[1]) = curve_tangent_function(c);

	delete_curve(c);

	curves[0]->num_points = num_points_on_curve(curves[0]);
	curves[1]->num_points = num_points_on_curve(curves[1]);

	AF_NODE_EXTRA *extra[2];
	FT_ScalarMemoryAlloc((POINTER*)&extra[0],sizeof(AF_NODE_EXTRA));
	extra[0]->af_node_type = STRING_NODE;
	n[0]->extra = (POINTER)extra[0];
	n[0]->size_of_extra = sizeof(AF_NODE_EXTRA);
	FT_ScalarMemoryAlloc((POINTER*)&extra[1],sizeof(AF_NODE_EXTRA));
	extra[1]->af_node_type = STRING_NODE;
	n[1]->extra = (POINTER)extra[1];
	n[1]->size_of_extra = sizeof(AF_NODE_EXTRA);

	set_current_interface(cur_intfc);
	if (debugging("trace"))
	    printf("Leaving split_string_curve()\n");
}	/* end split_string_curve */

extern void record_break_strings_gindex(Front *front)
{
	INTERFACE *intfc = front->interf;
	AF_PARAMS *af_params = (AF_PARAMS*)front->extra2;
	int *bs_gindex = af_params->break_strings_gindex;
	CURVE **c;

	if(FT_Dimension() != 3) return;

	/* build up the map from string gindex to string id */
	/* to make the input and output more intuitive to users */
	std::vector<CURVE*> string_curves = af_params->string_curves;
	int *array, n = string_curves.size();
	FT_ScalarMemoryAlloc((POINTER*)&array, sizeof(int) * n);
	for (int i = 0; i < n; ++i)
	    array[i] = Gindex(string_curves[i]);
	pp_gsync();
	pp_global_imax(array, n);
	for (int i = 0; i < n; ++i)
	    af_params->string_hash[array[i]] = i;

	if (af_params->break_strings_num <= 0) return;

	if (debugging("trace"))
	    printf("Entering record_break_strings_gindex()\n");

	for (int i = 0; i < af_params->break_strings_num; ++i)
	    bs_gindex[i] = array[bs_gindex[i]];

	FT_FreeThese(1,array);
	if (debugging("trace"))
	    printf("Leaving record_break_strings_gindex()\n");
}	/* end record_break_strings_gindex */

extern void set_unequal_strings(Front *front)
{
	INTERFACE *intfc = front->interf;
	AF_PARAMS *af_params = (AF_PARAMS*)front->extra2;
	int *us_gindex = af_params->unequal_strings_gindex;
	double tol = MACH_EPS;
	CURVE **c;

	if(FT_Dimension() != 3) return;
	if (af_params->unequal_strings_num <= 0) return;
	if (fabs(af_params->unequal_coeff - 1.0) < tol) return;

	if (debugging("trace"))
	    printf("Entering set_unequal_strings\n");

	std::vector<CURVE*> string_curves = af_params->string_curves;
	for (int i = 0; i < af_params->unequal_strings_num; ++i)
	{
	    CURVE *curve = string_curves[us_gindex[i]];
	    BOND *bond;
	    curve_bond_loop(curve, bond)
		bond->length0 *= af_params->unequal_coeff;
	}

	/* become invalid after first step and no longer used */
	af_params->string_curves.clear();

	if (debugging("trace"))
	    printf("Leaving record_break_strings_gindex()\n");
}	/* end set_unequal_strings */
