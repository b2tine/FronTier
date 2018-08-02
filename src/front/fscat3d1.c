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
*			fscat3d1.c
*
*	Copyright 1999 by The University at Stony Brook, All rights reserved.
*
*/


#define DEBUG_STRING    "fscatter"
#include <front/fdecs.h>

struct _POINT_LIST {
	POINT              *p;
	HYPER_SURF         *hs;
	HYPER_SURF_ELEMENT *hse;
	struct _POINT_LIST *prev, *next;
};
typedef struct _POINT_LIST POINT_LIST;

struct _POINT_LIST_STORE {
	POINT_LIST *pl;
	int        len;
};
typedef struct _POINT_LIST_STORE POINT_LIST_STORE;

LOCAL   boolean    open_fg = NO;

void    set_open_bdry_flag(boolean);
void    set_open_bdry_flag(boolean fg)
{
        open_fg = fg;
}

	/* LOCAL Function Declarations */
LOCAL	POINT_LIST	*set_point_list(TRI**,int,HYPER_SURF*,
					POINT_LIST_STORE*);
LOCAL	TRI*	find_following_null_tri(TRI*,POINT**,int*,ANGLE_DIRECTION);
LOCAL	boolean	add_matching_pt_to_hash_table(TRI**,TRI**,int,int,SURFACE*,
					      SURFACE*,P_LINK*,int);
LOCAL	boolean	buffer_extension3d1(INTERFACE*,INTERFACE*,int,int,boolean);
LOCAL	boolean    f_intfc_communication3d1(Front*);
LOCAL	boolean	match_tris_at_subdomain_bdry(SURFACE*,SURFACE*,TRI**,TRI**,
	                                     int,int);
LOCAL	boolean	match_two_tris(TRI*,TRI*);
LOCAL	boolean	null_side_on_surface(SURFACE*,TRI**,int*);
LOCAL	boolean	tri_out_domain(TRI*,double*,double*,int,int);
LOCAL	boolean	bond_out_domain1(BOND*,double*,double*,int,int);
LOCAL	int	append_adj_intfc_to_buffer1(INTERFACE*,INTERFACE*,
					   RECT_GRID*,int,int);
LOCAL	int	append_buffer_surface1(SURFACE*,SURFACE*,RECT_GRID*,int,int,
				      P_LINK*,int);
LOCAL	boolean	append_other_curves1(INTERFACE*,INTERFACE*,RECT_GRID*,int,int,
				      P_LINK*,int);
LOCAL	void	clip_intfc_at_grid_bdry1(INTERFACE*);
LOCAL	boolean	tri_set_out_domain(TRI*,double*,double*,int,int);
LOCAL	boolean	tri_cross_line(TRI*,double,int);
LOCAL	void	copy_tri_state_to_btri(BOND_TRI*,BOND*,ORIENTATION,INTERFACE *);
LOCAL	void	merge_point_pointers_at_subdomain_bdry(TRI**,TRI**,
						       int,P_LINK*,int);
LOCAL	void	strip_curve_from_surf(CURVE*,SURFACE*,ORIENTATION);
LOCAL	void	synchronize_tris_at_subdomain_bdry(TRI**,TRI**,int,P_LINK*,int);
	void 	sep_common_pt_for_open_bdry(INTERFACE*);
LOCAL	void 	incremental_alloc_tris(TRI***,int*);
LOCAL   boolean    find_null_side_loop(INTERFACE*,TRI*,int,TRI**,int*);
LOCAL	COMPONENT buffer_component(INTERFACE*,int,int);
LOCAL  void    test_curve_link(CURVE *);
LOCAL   double   *constr_position(double*,double*,
                        boolean (*constr_func)(POINTER,double*),POINTER);
LOCAL	void merge_gindex_of_tris(TRI**,int,RECT_GRID*,int,int,P_LINK*,int);
LOCAL	void merge_gindex_of_curve(CURVE*,RECT_GRID*,int,int,P_LINK*,int);
LOCAL	void merge_point_pair_gindex(POINT*,POINT*,RECT_GRID*,int,int);
LOCAL	boolean bond_match1(BOND*,BOND*);

LOCAL	double	tol1[MAXD]; /*TOLERANCE*/
#define MAX_NULL_SIDE      8000

#define MAX_SUBDOMAIN_TRIS      3000
LOCAL   boolean find_ending_null_side(TRI*,int,TRI**,int*);

LOCAL  void    test_curve_link(CURVE *c)
{
BOND    *b, *bn;

	if(c->first->prev != NULL || c->last->next != NULL)
	    printf("#test_curve_link, first or last != NULL.\n");
	
	for(b=c->first; b; b=b->next)
	{
	    bn = b->next;
	    if(bn != NULL)
	    {
	        if(bn->start != b->end)
	    	    printf("#test_curve_link, point error %p %p\n", 
		        (void*)bn->start, (void*)b->end);
		if(bn->prev != b)
	    	    printf("#test_curve_link, link error %p %p %p\n", 
		        (void*)b, (void*)bn, (void*)bn->prev);
	    }
	}
}

LOCAL   boolean    deact_default_comp = NO;
EXPORT  void    set_default_comp(boolean flag)
{
        deact_default_comp = flag;
}

/*ARGSUSED*/
EXPORT boolean f_intfc_communication3d(
	Front		*fr)
{
#if defined __MPI__
	if (pp_numnodes() > 1 && deact_default_comp == NO)
	{
	    communicate_default_comp(fr);
	    set_default_comp(YES);
	}
#endif /* defined __MPI__ */
	if (static_mesh(fr->interf))
	{
	    return f_intfc_communication3d3(fr);
	}
	else if (interface_reconstructed(fr->interf))
	{
	    return f_intfc_communication3d2(fr);
	}
	else
	{
	    return f_intfc_communication3d1(fr);
	}

}	/* end f_intfc_communication3d */

/*void sep_common_pt_for_open_bdry(INTERFACE*);*/
void     cut_buffer_tris_from_intfc(INTERFACE*, int, int);
int     cut_buffer_tris(SURFACE*, RECT_GRID*, int, int);
void    sep_common_pt_for_ref_bdry(INTERFACE *intfc, int, int);
void    make_ref_strip(double*, double*, INTERFACE*);
void    extend_long_ref_side(double*,double*,int, INTERFACE*,Front*);
void    cut_intfc_along_bdry(int, int, INTERFACE*);
void    flag_buffer_tris(int, int, INTERFACE*);
void   clip_intfc_in_dir(INTERFACE*, int, int);
void    reflect_curves_on_intfc(INTERFACE*,double*,double*);

void    set_open_bdry_flag(boolean);

LOCAL   boolean ref_append = NO;
void    set_ref_append(boolean);

void    set_ref_append(boolean fg)
{
        ref_append = fg;
}

void    construct_reflect_bdry(
        Front           *fr)
{
        INTERFACE    *intfc = fr->interf, *buf_intfc;
        INTERFACE    *cur_intfc = current_interface();
        RECT_GRID    *gr = fr->rect_grid;
        SURFACE      **s;
        double        *U = gr->U, *L = gr->L;
        double        *nor, p[3];
        int          me[MAXD], him[MAXD];
        int          i, j, k;
        int          dim = intfc->dim;
        boolean         status;
        static double nors[] = { 1.0,  0.0,  0.0,
                                 0.0,  1.0,  0.0,
                                 0.0,  0.0,  1.0,
                                -1.0,  0.0,  0.0,
                                 0.0, -1.0,  0.0,
                                 0.0,  0.0, -1.0};

        DEBUG_ENTER(construct_reflect_bdry)

        set_current_interface(intfc);
        strip_subdomain_bdry_curves(intfc);

        for (i = 0; i < dim; ++i)
        {
            for (j = 0; j < 2; ++j)
            {
                for (k = 0; k < dim; ++k)
                    him[k] = me[k] = pp_mynode();

                if (rect_boundary_type(intfc,i,j) != REFLECTION_BOUNDARY)
                    continue;

                set_open_bdry_flag(NO);

                clip_intfc_in_dir(intfc,i,j);
                buf_intfc = cut_buf_interface1(intfc,i,j,me,him);

                nor = nors + 3*i + 9*j;
                p[i] = (j > 0) ? U[i] : L[i];
                for(k = 1; k < dim; ++k)
                    p[(k+i)%dim] = 0.5*(U[(k+i)%dim] + L[(k+i)%dim]);

                if(debugging("cut_ref"))
                {
                    cut_intfc_along_bdry(i,j,intfc);
                    cut_intfc_along_bdry(i,j,buf_intfc);

                    flag_buffer_tris(i,j,intfc);
                    printf("#cut_ref  buf_intfc\n");
                    flag_buffer_tris(i,j,buf_intfc);
                }
                else
                {
                    cut_buffer_tris_from_intfc(intfc,i,j);
                    cut_buffer_tris_from_intfc(buf_intfc,i,j);
                }

                sep_common_pt_for_ref_bdry(intfc,i,j);
                sep_common_pt_for_ref_bdry(buf_intfc,i,j);

                reflect_interface(buf_intfc,p,nor);
                reflect_curves_on_intfc(buf_intfc,p,nor);

                make_ref_strip(p,nor,intfc);
                make_ref_strip(p,nor,buf_intfc);

                set_current_interface(intfc);

                status = FUNCTION_SUCCEEDED;
                set_ref_append(YES);
                status = buffer_extension3d1(intfc,buf_intfc,i,j,status);
                set_ref_append(NO);

                if (!status)
                {
                    (void) printf("WARNING in construct_reflect_bdry "
                              "buffer_extension3d1 failed for "
                              "i = %d, j = %d \n",i,j);
                    clean_up(ERROR);
                }

                (void) delete_interface(buf_intfc);
                set_current_interface(intfc);

                for (s = intfc->surfaces; s && *s; ++s)
                {
                    if (no_tris_on_surface(*s))
                    {
                        (void) delete_surface(*s);
                        --s;
                    }
                }
                if (intfc->surfaces != NULL && *(intfc->surfaces) == NULL)
                    intfc->surfaces = NULL;
                if (intfc->surfaces == NULL)
                {
                    if(debugging("cut_ref"))
                    	(void) printf("no tri, exit construct_reflect_bdry.\n");
                    goto exit_ref;
                }

            }   /*for j=0, 1, two direction communication.*/
            reset_intfc_num_points(intfc);
        }
exit_ref:
        set_open_bdry_flag(NO);
        set_current_interface(cur_intfc);
        DEBUG_LEAVE(construct_reflect_bdry)
}

int     cut_buffer_tris(SURFACE*, RECT_GRID*, int, int);
void    cut_buffer_tris_from_intfc(INTERFACE*, int, int);
void    make_ref_strip_for_surf(double*, double*, SURFACE*, INTERFACE*);
void    make_ref_strip(double*, double*, INTERFACE*);

#define  cut_case(a)  ((a[0].fg == 0) + (a[1].fg == 0) + (a[2].fg == 0) == 2 ? YES : NO)
/* if line cut two sides of the tri, it returns YES.*/
#define  out_side(a)  (a[0].fg != 0 ? 0 : (a[1].fg !=0 ? 1 : 2))
/* line without intersection with the side.*/
#define  cut_side(a)  (cut_case(a) ? Next_m3(out_side(a)) : (a[0].fg == 0 ? 0 : (a[1].fg == 0 ? 1 : 2)))
/* line cut the side*/

typedef  struct {
        int     fg;
        /*edge flag  0: side across the line;
                    -1: side inside the line;
                    -2: side outside the line;
        */

        POINT   *pt;
        /* corssing point if fg == 0 */
}       EDGE_CUT;

#define MAX_EDGE_CUT    3000

int     neighbor_tri_side(TRI*, TRI*);

/* nbtri is on which side of tri */
int     neighbor_tri_side(
        TRI     *tri,
        TRI     *nbtri)
{
        int     i;

        for(i=0; i<3; i++)
            if(!is_side_bdry(tri,i) && Tri_on_side(tri,i) == nbtri)
                return i;
        return -1;
}

int     cut_split_tris(
        TRI             **new_tris,
        TRI             *tri,
        SURFACE         *s,
        EDGE_CUT        **edge_cut)
{
        EDGE_CUT        *A;
        int             i, j, k, side, nside, pside;
        int             nnb, nt, nnt, nnt1;
        POINT           **p, *pa, *pb, *pte[2];
        TRI             *new_tri1, *new_tri2, *nbtri;
        TRI             *nbtris[2], *new_nbtris[3], **ptris;

        DEBUG_ENTER(cut_split_tris)

        /* tri is inside or outside, return. */
        if(Tri_order(tri) < 0)
        {
            DEBUG_LEAVE(cut_split_tris)
            return 0;
        }

        A = edge_cut[Tri_order(tri)];

        nnb = 0;
        if(cut_case(A))
        {
            p = Point_of_tri(tri);

            side = out_side(A);
            pside = Prev_m3(side);
            nside = Next_m3(side);

            pa = A[nside].pt;
            pb = A[pside].pt;

            /*fact: pside and nside has cut point on them. if a neighbor 
	      tri is on the side Tri_order must have a value
            */
            nbtris[0] = Tri_on_side(tri,pside);
            pte[0] = pb;
            if(nbtris[0] != NULL)
                nnb++;

            nbtris[nnb] = Tri_on_side(tri,nside);
            pte[nnb] = pa;
            if(nbtris[nnb] != NULL)
                nnb++;

            /* disconnect the neighbor tris */
            for(i=0; i<nnb; i++)
                Tri_on_side(nbtris[i], neighbor_tri_side(nbtris[i],tri)) = NULL;
            Tri_on_side(tri,pside) = NULL;
            Tri_on_side(tri,nside) = NULL;

            new_tri1 = make_tri(p[nside],pa,pb, NULL, NULL, NULL, NO);
            new_tri2 = make_tri(pa,p[pside],pb, NULL, NULL, NULL, NO);
            insert_tri_at_tail_of_list(new_tri1,s);
            insert_tri_at_tail_of_list(new_tri2,s);
            p[pside] = pb;

            link_neighbor_null_side_tris(new_tri1, new_tri2);
            link_neighbor_null_side_tris(new_tri1, tri);

            /* out side is inside the line */
            if(A[side].fg == -1)
            {
                Tri_order(tri) = -1;
                Tri_order(new_tri1) = -1;
                Tri_order(new_tri2) = -2;
            }
            else
            {
                Tri_order(tri) = -2;
                Tri_order(new_tri1) = -2;
                Tri_order(new_tri2) = -1;
            }

            new_tris[0] = tri;
            new_tris[1] = new_tri1;
            new_tris[2] = new_tri2;
            nt = 3;
        }
        else
        {
            p = Point_of_tri(tri);

            side = cut_side(A);
            pside = Prev_m3(side);
            nside = Next_m3(side);

            pa = A[side].pt;

            nbtris[0] = Tri_on_side(tri,side);
            pte[0] = pa;
            if(nbtris[0] != NULL)
                nnb++;

            /* disconnect the neighbor tris */
            for(i=0; i<nnb; i++)
                Tri_on_side(nbtris[i], neighbor_tri_side(nbtris[i],tri)) = NULL;

            new_tri1 = make_tri(p[pside], pa, p[nside], NULL,NULL,NULL,NO);
            insert_tri_at_tail_of_list(new_tri1,s);

            /* reconnect nside neighbor */
            Neighbor_on_side(new_tri1,2) = Neighbor_on_side(tri,nside);
            if(is_side_bdry(tri,nside))
            {
                set_20_bdry(Boundary_tri(new_tri1), YES);
                set_side_bdry(Boundary_tri(tri), nside, NO);
            }
            if(Neighbor_on_side(tri,nside) != NULL)
                if(is_side_bdry(tri,nside))
                    Bond_tri_on_side(tri,nside)->tri = new_tri1;
                else
                {
                    nbtri = Tri_on_side(tri,nside);
                    i = neighbor_tri_side(nbtri,tri);
                    Tri_on_side(nbtri,i) = new_tri1;
                }

            p[nside] = pa;
            Tri_on_side(tri,side) = NULL;
            Tri_on_side(tri,nside) = NULL;

            link_neighbor_null_side_tris(new_tri1, tri);

            if(A[pside].fg == -1)
            {
                Tri_order(tri) = -1;
                Tri_order(new_tri1) = -2;
            }
            else
            {
                Tri_order(tri) = -2;
                Tri_order(new_tri1) = -1;
            }

            new_tris[0] = tri;
            new_tris[1] = new_tri1;
            nt = 2;
        }

        /* split the neighbor tris and linking them. */
        for(i=0; i<nnb; i++)
        {
            nnt = cut_split_tris(new_nbtris, nbtris[i], s, edge_cut);
            if(nnt == 0)
            {
                printf("#cut_split_tris, enter nnt == 0 case.\n");

                /* check if pte[i] is a vertex of nbtri */
                for(j=0; j<3; j++)
                    if(Point_of_tri(nbtris[i])[j] == pte[i])
                    {
                        nbtri = nbtris[i];
                        break;
                    }

                /* find nbtri which has pte[i] as its vertex. */
                if(j == 3)
                {
                    for(j=0; j<3; j++)
                    {
                        if(is_side_bdry(nbtris[i], j))
                            continue;
                        nbtri = Tri_on_side(nbtris[i],j);
                        if(nbtri == NULL)
                            continue;
                        for(k=0; k<3; k++)
                            if(Point_of_tri(nbtri)[k] == pte[i])
                                break;
                        if(k < 3)
                            break;
                    }
                    if(j == 3)
                    {
                        printf("ERROR cut_split_tris, unable to split nbtris[%d].\n", i);
                        clean_up(ERROR);
                    }
                }

                /* nbtri has vertex pte[i] */
                nnt1 = set_tri_list_around_point(pte[i],
                    nbtri, &ptris, s->interface);
                for(j=0; j<nnt1; j++)
                    new_nbtris[j] = ptris[j];
                nnt = nnt1;

                for(j=0; j<nnt; j++)
                    print_tri(new_nbtris[j], s->interface);

                for(j=0; j<nt; j++)
                    print_tri(new_tris[j], s->interface);
            }

            for(j=0; j<nnt; j++)
                for(k=0; k<nt; k++)
                    link_neighbor_null_side_tris(new_tris[k], new_nbtris[j]);
        }

        DEBUG_LEAVE(cut_split_tris)

        return nt;
}

void	flag_point_pos(
	double		plane,
	int		dir,
	int		nb,
	SURFACE		*s,
	INTERFACE	*intfc)
{
	RECT_GRID	*gr = computational_grid(intfc);
	POINT		**p;
	TRI		*tri;
	double		crds, tol;
	int		i;

	tol = 0.01*min3(gr->h[0], gr->h[1], gr->h[2]);

	/* 1 not calculated */
	for(tri=first_tri(s); !at_end_of_tri_list(tri,s); tri=tri->next)
	    for(i=0; i<3; i++)
		Point_order(Point_of_tri(tri)[i]) = 1;

	/* check points 0 on line -1 insdie line  -2 outside line */
	for(tri=first_tri(s); !at_end_of_tri_list(tri,s); tri=tri->next)
	    for(i=0; i<3; i++)
	    {
		p = Point_of_tri(tri);
		if(Point_order(p[i]) != 1)
		    continue;
		
		crds = Coords(p[i])[dir];
		if(crds > plane + tol)
		    Point_order(p[i]) = nb == 0 ? -1 : -2;
		else if(crds < plane - tol)
		    Point_order(p[i]) = nb == 0 ? -2 : -1 ;
		else
		    Point_order(p[i]) = 0;
	    }
}

void	cut_curve_along_line(
	int		dir,
	int		nb,
	double		plane,
	INTERFACE	*intfc)
{
	RECT_GRID	*gr = computational_grid(intfc);
	CURVE	**c;
	BOND	*bd;
	POINT	*p, *p1, *p2, *newp;
	double	crds, a, b, tol;

	tol = 0.01*min3(gr->h[0], gr->h[1], gr->h[2]);
	
	for (c = intfc->curves; c && *c; ++c)
	{
	    for (bd = (*c)->first; bd != NULL; bd = bd->next)
	    {
		p = bd->start;
		crds = Coords(p)[dir];
		if(crds > plane + tol)
		    Point_order(p) = nb == 0 ? -1 : -2;
		else if(crds < plane - tol)
		    Point_order(p) = nb == 0 ? -2 : -1 ;
		else
		    Point_order(p) = 0;
	
		p = bd->end;
		crds = Coords(p)[dir];
		if(crds > plane + tol)
		    Point_order(p) = nb == 0 ? -1 : -2;
		else if(crds < plane - tol)
		    Point_order(p) = nb == 0 ? -2 : -1 ;
		else
		    Point_order(p) = 0;
	    }
cut_curve:
	    for (bd = (*c)->first; bd != NULL; bd = bd->next)
	    {
		p1 = bd->start;
		p2 = bd->end;
		if(Point_order(p1) == Point_order(p2))
		    continue;
		if(Point_order(p1) == 0 || Point_order(p2) == 0)
		    continue;
		
		a = fabs(Coords(p1)[dir]-plane);
		b = fabs(Coords(p2)[dir]-plane);
		newp = copy_point(p1);

		Coords(newp)[dir] = plane;
		Coords(newp)[Next_m3(dir)] = 
		    b*Coords(p1)[Next_m3(dir)] + a*Coords(p2)[Next_m3(dir)];
		Coords(newp)[Next_m3(dir)] /= (a+b);
		    
		Coords(newp)[Prev_m3(dir)] = 
		    b*Coords(p1)[Prev_m3(dir)] + a*Coords(p2)[Prev_m3(dir)];
		Coords(newp)[Prev_m3(dir)] /= (a+b);
		Point_order(newp) = 0;
	
		insert_point_in_bond(newp, bd, *c);
		goto cut_curve;
	    }
	}
}

void	flag_null_on_line_tris(
	double		plane,
	int		dir,
	int		nb,
	SURFACE		*s,
	INTERFACE	*intfc)
{
	int	i;
	TRI	*tri;
	POINT	**p;

	flag_point_pos(plane, dir, nb, s, intfc);
	
	for(tri=first_tri(s); !at_end_of_tri_list(tri,s); tri=tri->next)
	    Tri_order(tri) = 0;
	
	for(tri=first_tri(s); !at_end_of_tri_list(tri,s); tri=tri->next)
	    for(i=0; i<3; i++)
	    {
		p = Point_of_tri(tri);
		if(Point_order(p[i]) == 0 && Point_order(p[Next_m3(i)]) == 0)
		{
		    if(Tri_on_side(tri, i) != NULL)
		    {
			continue;
			printf("ERROR flag_null_on_line_tris, bdry tri side is not null.\n");
			printf("side %d\n", i);
			print_tri(tri, intfc);
			clean_up(ERROR);
		    }
		    Tri_order(tri) |= Bin_side(i);
		    break;
		}
	    }
}


void	set_surf_edge_cut(
	EDGE_CUT	**edge_cut,
	double		plane,
	int		dir,
	int		nb,
	SURFACE		*s,
	INTERFACE	*intfc)
{
	int	i, j, k, nt, pf[3];
	double	a, b;
	POINT	**p, *newp;
	TRI	*tri, *nbtri;

	flag_point_pos(plane, dir, nb, s, intfc);

	for(tri=first_tri(s); !at_end_of_tri_list(tri,s); tri=tri->next)
	    Tri_order(tri) = -3;

	nt = 0;
	for(tri=first_tri(s); !at_end_of_tri_list(tri,s); tri=tri->next)
	{
	    p = Point_of_tri(tri);
	    for(i=0; i<3; i++)
		pf[i] = Point_order(p[i]);

	    /* no inside point, remove tri */
	    if(pf[0] != -1 && pf[1] != -1 && pf[2] != -1)
	    {
		Tri_order(tri) = -2;
		continue;
	    }
	    
	    /* no outside point, keep tri */
	    if(pf[0] != -2 && pf[1] != -2 && pf[2] != -2)
	    {
		Tri_order(tri) = -1;
		continue;
	    }

	    /* there are edge cuts and at least 1 inside 1 outside */
	    Tri_order(tri) = nt;
	    for(i=0; i<3; i++)
	    {
		j = Next_m3(i);
		if(pf[i] == pf[j])  /* edge inside or outside */
		    edge_cut[nt][i].fg = pf[i];
		else  if(pf[i] == 0)  /*1 point on */
		    edge_cut[nt][i].fg = pf[j];
		else  if(pf[j] == 0)  /*1 point on */
		    edge_cut[nt][i].fg = pf[i];
		else /*1 in and 1 out */
		{
		    edge_cut[nt][i].fg = 0;
		   
		    /* check if edge crx is calculated by the neighbor tri */
		    nbtri = Tri_on_side(tri, i);
		    if(nbtri != NULL && Tri_order(nbtri) >= 0)
		    {
			k = Tri_order(nbtri);
			edge_cut[nt][i].pt = 
			    edge_cut[k][neighbor_tri_side(nbtri,tri)].pt;
		    }
		    else
		    {
			a = fabs(Coords(p[i])[dir]-plane);
			b = fabs(Coords(p[j])[dir]-plane);
			newp = copy_point(p[i]);

			Coords(newp)[dir] = plane;
			Coords(newp)[Next_m3(dir)] = 
			    b*Coords(p[i])[Next_m3(dir)] + a*Coords(p[j])[Next_m3(dir)];
			Coords(newp)[Next_m3(dir)] /= (a+b);
		    
			Coords(newp)[Prev_m3(dir)] = 
			    b*Coords(p[i])[Prev_m3(dir)] + a*Coords(p[j])[Prev_m3(dir)];
			Coords(newp)[Prev_m3(dir)] /= (a+b);
		   
			interpolate_crx_pt_states_on_edge(intfc,newp,tri,s,i);
			
			edge_cut[nt][i].pt = newp; 
		    }
		}
	    }

	    /* crx tri is found */
	    nt++;
	    if(nt >= MAX_EDGE_CUT)
	    {
		printf("ERROR in set_surf_edge_cut, too many cut edges.\n");
		clean_up(ERROR);
	    }
	}
}

void	cut_surf_along_line(
	SURFACE		*s,
	EDGE_CUT	**edge_cut)
{
	TRI	*tri, *tris[MAX_EDGE_CUT], *new_tris[3];
	int	i, nt;

	nt = 0;
	for(tri=first_tri(s); !at_end_of_tri_list(tri,s); tri=tri->next)
	    if(Tri_order(tri) >= 0)
		tris[nt++] = tri;
	
	for(i=0; i<nt; i++)
	    cut_split_tris(new_tris, tris[i], s, edge_cut);
	
	for(tri=first_tri(s); !at_end_of_tri_list(tri,s); tri=tri->next)
	    if(Tri_order(tri) == -2)
		remove_tri_from_surface(tri,s,NO);
	
}
		
void	cut_intfc_along_bdry(int, int, INTERFACE*);

void	cut_intfc_along_bdry(
	int		dir,
	int		nb,
	INTERFACE	*intfc)
{
	static EDGE_CUT	**edge_cut = NULL;
	double		plane, tol;
	RECT_GRID	*gr = computational_grid(intfc);
	SURFACE		**s;
	INTERFACE	*sav_intfc = current_interface();

	DEBUG_ENTER(cut_intfc_along_bdry)
	
	set_current_interface(intfc);

	if(edge_cut == NULL)
	    bi_array(&edge_cut, MAX_EDGE_CUT,3,sizeof(EDGE_CUT));

	if(nb == 0)
	    plane = gr->L[dir] + 0.5*gr->h[dir];
	else
	    plane = gr->U[dir] - 0.5*gr->h[dir];
       
	cut_curve_along_line(dir, nb, plane, intfc);
	
	for(s = intfc->surfaces; s && *s; ++s)
	{
	    set_surf_edge_cut(edge_cut,plane,dir,nb,*s,intfc);
	    cut_surf_along_line(*s, edge_cut);
	}
	
	{
	CURVE	**c;
	BOND	*bs;
	BOND_TRI	**btris;

	for (c = intfc->curves; c && *c; ++c)
	{
	    for (bs = (*c)->first; bs != NULL; bs = bs->next)
	    {
		btris = Btris(bs);
		if(btris != NULL && *btris == NULL)
		{
		    printf(" #  ");
		    Btris(bs) = NULL;
		}
	    }
	}
	}

	cut_out_curves_in_buffer(intfc);
	
	set_current_interface(sav_intfc);

	DEBUG_LEAVE(cut_intfc_along_bdry)
}

void	flag_buffer_tris(int, int, INTERFACE*);
void	flag_buffer_tris(
	int		dir,
	int		nb,
	INTERFACE	*intfc)
{
	double		plane;
	RECT_GRID	*gr = computational_grid(intfc);
	SURFACE		**s;

	if(nb == 0)
	    plane = gr->L[dir] + 0.5*gr->h[dir];
	else
	    plane = gr->U[dir] - 0.5*gr->h[dir];
        
	for(s = intfc->surfaces; s && *s; ++s)
	    flag_null_on_line_tris(plane, dir, nb, *s, intfc);
	
	DEBUG_LEAVE(cut_intfc_along_bdry)
}

void	make_ref_strip(
	double		*pt,
	double		*nor,
	INTERFACE	*intfc)
{
	SURFACE		**s;
 	int		i;
	POINT		*p;
	TRI		*tri;
	INTERFACE	*sav_intfc;

	sav_intfc = current_interface();
	set_current_interface(intfc);

	for (s = intfc->surfaces; s && *s; ++s)
	    for(tri=first_tri(*s); !at_end_of_tri_list(tri,*s); tri=tri->next)
	        for(i=0; i<3; i++)
	        {
		    p = Point_of_tri(tri)[i];
		    Cor_point(p) = NULL;
	        }
      
	for (s = intfc->surfaces; s && *s; ++s)
	    make_ref_strip_for_surf(pt, nor, *s, intfc);

	/*do not know the order of btris constructed in the above function
	  is consistent with the previous ones, must order interface. 
	*/
	order_interface(intfc);
	set_current_interface(sav_intfc);
}

#define	 swap_pointers(a, b)   {st = a; a = b; b = st;}
/* assume the curve c is already reflected. */

void	reflect_btris_on_curve(
	CURVE		*c,
	INTERFACE	*intfc,
	double		*p,
	double		*nor)
{
	BOND_TRI	**btris;
	BOND		*b;
	Locstate	st;

	/* reflect btri states  */
	for(b=c->first; b != NULL; b=b->next)
	{
	    for(btris=Btris(b); btris && *btris; btris++)
	    {
		/* after reverse curve, b->start is the end of btri. */
		reflect_state(left_end_btri_state(*btris), intfc,
			Coords(b->start), p, nor);
		reflect_state(right_end_btri_state(*btris), intfc,
			Coords(b->start), p, nor);

		/* the last point on curve */
		if(!is_closed_curve(c) && b->next == NULL)
		{
		    reflect_state(left_start_btri_state(*btris), intfc,
			    Coords(b->end), p, nor);
		    reflect_state(right_start_btri_state(*btris), intfc,
			    Coords(b->end), p, nor);
		}
		
		swap_pointers(left_start_btri_state(*btris),
					left_end_btri_state(*btris));
		swap_pointers(right_start_btri_state(*btris),
					right_end_btri_state(*btris));
	    }
	}

	/*relink btri considering the orientation of c.
	  surface, c and b are reflected
	*/
	for(b=c->first; b != NULL; b=b->next)
	{
	    for(btris=Btris(b); btris && *btris; btris++)
	    {
		link_tri_to_bond(*btris,(*btris)->tri,(*btris)->surface,b,c);
	    }
	}

}

void	reflect_curves_on_intfc(INTERFACE*,double*,double*);

void	reflect_curves_on_intfc(
	INTERFACE	*intfc,
	double		*p,
	double		*nor)
{
	CURVE		**c;
	
	for(c=intfc->curves; c && *c; c++)
	{
	    reverse_curve(*c);
	    reflect_btris_on_curve(*c, intfc, p, nor);
	}
}

BOND	*reflect_curve_bond(
	CURVE		**ref_cur,
	POINT		*p,
	POINT		*pf,
	INTERFACE	*intfc)
{
	BOND		*bnew, *b;
	CURVE		**c, *curve;
	int		dir;

	b = NULL;
	for(c=intfc->curves; c && *c; c++)	
	{
	    if((*c)->first->start == p)
	    {
		b = (*c)->first;
		dir = NEGATIVE_ORIENTATION;
		curve = *c;
		break;
	    }
	    if((*c)->last->end == p)
	    {
		b = (*c)->last;
		dir = POSITIVE_ORIENTATION;
		curve = *c;
		break;
	    }
	}
	if(b == NULL)
	{
	    printf("ERROR reflect_curve_bond, no curve has point.\n");
	    clean_up(ERROR);
	}

	Boundary_point(pf) = YES;
	if(dir == POSITIVE_ORIENTATION)
	{
	    bnew = Bond(p, pf);
	    bnew->prev = b;
	    curve->last = bnew;
	    bnew->next = NULL;
	    b->next = bnew;
	    ++curve->num_points;
	    curve->end->posn = pf;
	}
	else
	{
	    bnew = Bond(pf, p);
	    bnew->next = b;
	    curve->first = bnew;
	    bnew->prev = NULL;
	    b->prev = bnew;
	    ++curve->num_points;
	    curve->start->posn = pf;
	}

	*ref_cur = curve;
	return bnew;
}

void	contruct_ref_btri(
	TRI		*tri,
	SURFACE		*s,
	POINT		*p1,
	POINT		*p2,
	double		*pt,
	double		*nor,
	INTERFACE	*intfc)
{
	BOND_TRI	*btri;
	int		sizest;
	CURVE		**c, *curve;
	BOND		*b;

	for(c=intfc->curves; c && *c; c++)
	{
	    /* find the corresponding bond. */
	    for(b=(*c)->first; b != NULL; b=b->next)
	    {
		if((b->start == p1 && b->end == p2) ||
		   (b->start == p2 && b->end == p1))
		    break;
	    }
	    if(b != NULL)
	    {
		curve = *c;
		break;
	    }
	}
	if(b == NULL)
	{
	    printf("ERROR contruct_ref_btri, the bond is not found.\n");
	    clean_up(ERROR);
	}

	sizest = size_of_state(intfc);
	btri = link_tri_to_bond(NULL, tri, s, b, curve);

	if(b == curve->first)
	{
	    ft_assign(left_start_btri_state(btri),left_end_btri_state(btri), 
					sizest);
	    ft_assign(right_start_btri_state(btri),right_end_btri_state(btri), 
					sizest);
	    reflect_state(left_start_btri_state(btri),intfc,Coords(b->start), 
					pt,nor);
	    reflect_state(right_start_btri_state(btri),intfc,Coords(b->start), 
					pt,nor);
	}
	else
	{
	    ft_assign(left_end_btri_state(btri),left_start_btri_state(btri), 
					sizest);
	    ft_assign(right_end_btri_state(btri),right_start_btri_state(btri), 
					sizest);
	    reflect_state(left_end_btri_state(btri),intfc,Coords(b->end),pt, 
					nor);
	    reflect_state(right_end_btri_state(btri),intfc,Coords(b->end),pt, 
					nor);
	}
}

void	make_ref_strip_for_surf(
	double		*pt,
	double		*nor,
	SURFACE		*s,
	INTERFACE	*intfc)
{
	TRI	*tri, **ptris, *new_tri1, *new_tri2;
	int	i,j, nt;
	POINT	*p, *p1, *p2, *p1f, *p2f;
 	BOND	*b1, *b2;
	CURVE	*c1, *c2;

	/* set curve ref point */
	for(tri=first_tri(s); !at_end_of_tri_list(tri,s); tri=tri->next)
	    for(i=0; i<3; i++)
	    {
		if((Tri_order(tri) & Bin_side(i)) == 0)
		    continue;
	
		p1 = Point_of_tri(tri)[i];
		if(Cor_point(p1) == NULL)
		{
		    p1f = copy_point(p1);
		    reflect_point(p1f, pt, nor, intfc);

		    /* if a point on curve is reflected, a bond will be made. */
		    if(Boundary_point(p1))
			b1 = reflect_curve_bond(&c1, p1, p1f, intfc);
		    Cor_point(p1) = p1f;
		}
		else
		    p1f = (POINT*)Cor_point(p1);

		p2 = Point_of_tri(tri)[Next_m3(i)];
		if(Cor_point(p2) == NULL)
		{
		    p2f = copy_point(p2);
		    reflect_point(p2f, pt, nor, intfc);
		    if(Boundary_point(p2))
			b2 = reflect_curve_bond(&c2, p2, p2f, intfc);
		    Cor_point(p2) = p2f;
		}
		else
		    p2f = (POINT*)Cor_point(p2);

		new_tri1 = make_tri(p2,p1,p2f,NULL,NULL,NULL,NO);

		insert_tri_at_tail_of_list(new_tri1,s);
		link_neighbor_null_side_tris(new_tri1,tri);

		new_tri2 = make_tri(p2f,p1,p1f,NULL,NULL,NULL,NO);
		insert_tri_at_tail_of_list(new_tri2,s);
		link_neighbor_null_side_tris(new_tri1,new_tri2);

		nt = set_tri_list_around_point(p1, tri, &ptris, intfc);
		link_neighbor_null_side_tris(ptris[nt-1], new_tri2);
		
		if(Boundary_point(p2))
		    contruct_ref_btri(new_tri1, s, p2, p2f, pt, nor, intfc);
		
		nt = set_tri_list_around_point(p2, tri, &ptris, intfc);
		link_neighbor_null_side_tris(ptris[0], new_tri1);
		
		if(Boundary_point(p1))
		    contruct_ref_btri(new_tri2, s, p1, p1f, pt, nor, intfc);
	    }

}

void   clip_intfc_in_dir(INTERFACE*, int, int);

void   clip_intfc_in_dir(
	INTERFACE	*intfc,
	int		dir,
	int		nb)
{
	INTERFACE	*cur_intfc = current_interface();
	RECT_GRID	*gr = computational_grid(intfc);
	int		dim = intfc->dim;
	double		L[MAXD],U[MAXD];

	DEBUG_ENTER(clip_intfc_in_dir)
	
	set_current_interface(intfc);
	
	L[dir] = gr->L[dir];
	U[dir] = gr->U[dir];
	if (gr->lbuf[dir] == 0) L[dir] -= 0.5*gr->h[dir];
	if (gr->ubuf[dir] == 0) U[dir] += 0.5*gr->h[dir];
	    
	if(rect_boundary_type(intfc,dir,0) == OPEN_BOUNDARY)
	{
	    L[dir] = gr->VL[dir];
            if(debugging("ref_scat"))
		L[dir] -= 4.0*gr->h[dir];
	}
	if(rect_boundary_type(intfc,dir,1) == OPEN_BOUNDARY)
	{
	    U[dir] = gr->VU[dir];
	    if(debugging("ref_scat"))
		U[dir] += 4.0*gr->h[dir];
	}
	
	open_null_sides1(intfc,L,U,dir,nb);
	
	reset_intfc_num_points(intfc);
	set_current_interface(cur_intfc);
	
	DEBUG_LEAVE(clip_intfc_in_dir)
}



typedef  struct {
	TRI		*tri;
	int		side;
	double		len;
}    TRI_REF;

int compare_tri_ref(const void *a, const void *b)
{
	TRI_REF  *c1=(TRI_REF*)a, *c2=(TRI_REF*)b;

	if(c1->len > c2->len)
	    return -1;
	else
	    return 1;
}
#define		MAX_TRI_REF	1000

void	extend_long_ref_side_for_surf(double*,double*,int, SURFACE*,INTERFACE*,Front*);
void	extend_long_ref_side(double*,double*,int, INTERFACE*,Front*);

void	extend_long_ref_side(
	double		*pt,
	double		*nor,
	int		dir,
	INTERFACE	*intfc, 
	Front		*fr)
{
	SURFACE		**s;
        
	for (s = intfc->surfaces; s && *s; ++s)
	    extend_long_ref_side_for_surf(pt, nor, dir, *s, intfc, fr);
}

void	extend_long_ref_side_for_surf(
	double		*pt,
	double		*nor,
	int		dir,
	SURFACE		*s,
	INTERFACE	*intfc, 
	Front		*fr)
{
	TRI		*tri, **ptris, *new_tri1, *new_tri2;
	int		i,j, nt, nt1, side;
	POINT		*p, *p1, *pf;
	RECT_GRID	*gr = computational_grid(intfc);
	double		tol, len, v[3], hmin;
	TRI_REF		tri_ref[MAX_TRI_REF];
	static Locstate	ans = NULL;

	if(ans == NULL)
	    scalar(&ans, fr->sizest);
	
	hmin = min3(gr->h[0], gr->h[1], gr->h[2]);
	tol = 0.7;

	/* set curve ref point */
	nt = 0;
	for(tri=first_tri(s); !at_end_of_tri_list(tri,s); tri=tri->next)
	    for(i=0; i<3; i++)
	    {
		if((Tri_order(tri) & Bin_side(i)) == 0)
		    continue;
		
		p1 = Point_of_tri(tri)[i];
		difference(Coords(p1), pt, v, 3);
		len = fabs(Dot3d(v, nor));
		if(len < tol*hmin)
		    continue;

		tri_ref[nt].tri = tri;
		tri_ref[nt].side = i;
		tri_ref[nt].len = len;
		nt++;
		if(nt >= MAX_TRI_REF)
		{
		    printf("ERROR in extend_long_ref_side, nt is too large.\n");
		    clean_up(ERROR);
		}
	    }

	if(nt == 0)
	    return;

	qsort((POINTER)tri_ref, nt, sizeof(TRI_REF), compare_tri_ref);

	printf("#ref_extend %d\n", nt);
	for(i=0; i<nt; i++)
	{
	    tri = tri_ref[i].tri;
	    side = tri_ref[i].side;

	    if((Tri_order(tri) & Bin_side(side)) == 0)
		continue;
	    
	    p = Point_of_tri(tri)[side];
	    pf = copy_point(p);
	    reflect_point(pf, pt, nor, intfc);
	    
	    Coords(pf)[dir] = Coords(pf)[dir]/3.0 + Coords(p)[dir]*2.0/3.0;
	    
	    bi_interpolate_intfc_states(intfc, 1.0/3.0, 2.0/3.0, 
		Coords(pf), left_state(pf), Coords(p), left_state(p), ans);
	    ft_assign(left_state(pf), ans, fr->sizest);
	    
	    bi_interpolate_intfc_states(intfc, 1.0/3.0, 2.0/3.0, 
		Coords(pf), right_state(pf), Coords(p), right_state(p), ans);
	    ft_assign(right_state(pf), ans, fr->sizest);

	    p1 = Point_of_tri(tri)[Next_m3(side)];
	    new_tri1 = make_tri(p1, p ,pf, NULL,NULL,NULL,NO);
	    insert_tri_at_tail_of_list(new_tri1,s);

	    link_neighbor_null_side_tris(new_tri1,tri);
            
	    j = Vertex_of_point(new_tri1,pf);
	    Tri_order(new_tri1) = Bin_side(j);
	    
	    nt1 = set_tri_list_around_point(p, tri, &ptris, intfc);
	    j = Vertex_of_point(ptris[nt1-1], p);
	    
	    p1 = Point_of_tri(ptris[nt1-1])[Prev_m3(j)];
	    new_tri2 = make_tri(p, p1, pf, NULL,NULL,NULL,NO);
	    insert_tri_at_tail_of_list(new_tri2,s);
	    
	    link_neighbor_null_side_tris(new_tri2,ptris[nt1-1]);
	    link_neighbor_null_side_tris(new_tri1,new_tri2);
	    
	    j = Vertex_of_point(new_tri2,p1);
	    Tri_order(new_tri2) = Bin_side(j);

	    set_side_bdry(Tri_order(tri), side, 0);
	    set_side_bdry(Tri_order(ptris[nt1-1]), Vertex_of_point(ptris[nt1-1],p1), 0);
	}
}

void sep_common_pt_for_ref_bdry(
	INTERFACE	*intfc,
	int		dir,
	int		nb)
{
	SURFACE			**s;
	TRI			*tri;
	POINT			**p, **p1;
	int			i, j, nt, num;
	double			tmp_fac, cut_plane, crds;
	double			tol = tol1[dir];
	TRI	 **tris = NULL;
	int	 num_tri_allocated = 0;
	
	cut_plane = nb==0 ? -HUGE_VAL : HUGE_VAL;

	nt = 0;
	incremental_alloc_tris(&tris,&num_tri_allocated);
	for(s=intfc->surfaces; s && *s; s++)
	{
	    for(tri=first_tri(*s); !at_end_of_tri_list(tri,*s); tri=tri->next)
	    {
		if(Tri_order(tri) != 0)
		{
		    tris[nt++] = tri;
		    if(nt >= num_tri_allocated)
		    {
			incremental_alloc_tris(&tris,&num_tri_allocated);
		    }

		    /*tri is the triangle needed to connect with the strip
		      cut_plane is the highest/lowest coord
		    */
		    p = Point_of_tri(tri);
		    for(i=0; i<3; i++)
		    {
			crds = Coords(p[i])[dir];
			if((nb == 0 && cut_plane < crds) ||
			   (nb == 1 && cut_plane > crds) )
			    cut_plane = crds;
		    }
		}
	    }
	}

	num = nt;
	for(s=intfc->surfaces; s && *s; s++)
	{
	    for(tri=first_tri(*s); !at_end_of_tri_list(tri,*s); tri=tri->next)
	    {
		if(Tri_order(tri) != 0)
		    continue;
		p = Point_of_tri(tri);
		for(i=0; i<3; i++)
		{
		    crds = Coords(p[i])[dir];
		    if((nb == 0 && cut_plane + tol >= crds) ||
		       (nb == 1 && cut_plane - tol <= crds) )
			break;
		}

		/* tri is in the range of ref tris */
		if(i < 3)
		{
		    for(i=0; i<3; i++)
		    {
			p = Point_of_tri(tri);
			for(j=0; j<num; j++)
			{
			    p1 = Point_of_tri(tris[j]);
			    if(p1[0] == p[i] || p1[1] == p[i] || p1[2] == p[i])
				break;
			}
		
			/* p is a point of tri on ref bdry */
			if(j<num)
			{
			    tris[nt++] = tri;
			    if(nt >= num_tri_allocated)
			    {
				incremental_alloc_tris(&tris,
						&num_tri_allocated);
			    }
			    break;
			}
		    }
		}
	    }
	}
	
	tmp_fac = get_comm_pt_fac();
	set_comm_pt_fac(0.999);
	if(nt > 0)
	    sep_common_point_from_loop(tris, nt, NULL, NULL, intfc);
	set_comm_pt_fac(tmp_fac);
	free_these(1,tris);
}

void cut_buffer_tris_from_intfc(
	INTERFACE	*intfc,
	int		dir,
	int		nb)
{
	SURFACE		**s;
	RECT_GRID	*gr = computational_grid(intfc);
	INTERFACE	*sav_intfc;
	
	sav_intfc = current_interface();

	set_current_interface(intfc);

	for (s = intfc->surfaces; s && *s; ++s)
	    cut_buffer_tris(*s, gr, dir, nb);
	
	{
	CURVE	**c;
	BOND	*bs;
	BOND_TRI	**btris;

	for (c = intfc->curves; c && *c; ++c)
	{
	    for (bs = (*c)->first; bs != NULL; bs = bs->next)
	    {
		btris = Btris(bs);
		if(btris != NULL && *btris == NULL)
		{
		    printf(" #  ");
		    Btris(bs) = NULL;
		}
	    }
	}
	}

	cut_out_curves_in_buffer(intfc);

	set_current_interface(sav_intfc);
}

int cut_buffer_tris(
	SURFACE		*surf,
	RECT_GRID	*grid,
	int		dir,
	int		nb)
{
	TRI	   *tri, *nbtri;
	int	   ns, i, j;
	double	   crx_coord;
	static TRI **tris_s = NULL;
	static int len_tris_s = 0;
	POINT	   **p;

	DEBUG_ENTER(cut_bubber_tris)
	
	if (len_tris_s < surf->num_tri)
	{
	    len_tris_s = surf->num_tri;
	    if (tris_s != NULL)
		free(tris_s);
	    uni_array(&tris_s,len_tris_s,sizeof(TRI *));
	}

	crx_coord = (nb == 0) ? grid->L[dir] : grid->U[dir];

	for (tri=first_tri(surf); !at_end_of_tri_list(tri,surf); tri=tri->next)
	    Tri_order(tri) = 0;

	flag_point_pos(crx_coord, dir, nb, surf, surf->interface);

	ns = 0;
	for (tri=first_tri(surf); !at_end_of_tri_list(tri,surf); tri=tri->next)
	{
	    /*on reflect bdry and 3comp curve is perpendicular to the ref plane,
	      this condition is enough 
	    */
	    if (tri_cross_line(tri,crx_coord,dir) == YES)
	    {
		tris_s[ns++] = tri;
		for(i=0; i<3; i++)
		{
		    if(!is_side_bdry(tri,i))
		    {
		        nbtri = Tri_on_side(tri,i);
			if(nbtri != NULL)
			{
			    for(j=0; j<3; j++)
				if(!is_side_bdry(nbtri,j) && 
				    Tri_on_side(nbtri,j) == tri)
				{
				    Tri_order(nbtri) |= Bin_side(j);
		        	}
			}
		    }
		}
	    }

	    p = Point_of_tri(tri);
	    if(Point_order(p[0]) == -2 && 
	       Point_order(p[1]) == -2 && Point_order(p[2]) == -2)
		tris_s[ns++] = tri;
	}
	
	for(i=0; i<ns; i++)
	    remove_tri_from_surface(tris_s[i], surf, NO);

	DEBUG_LEAVE(cut_buffer_tris)
	return YES;
}		/*end append_buffer_surface1*/

void tecplot_surface_curve(FILE*, SURFACE*);
	char *  get_directory();
	char *  get_directory()
{       
	static char def_dir[200];
        
        sprintf(def_dir,"/gpfs/home2/hklim/FronTier_Contactor/iFluid/tecout/");
        return  def_dir;
}       

LOCAL	boolean f_intfc_communication3d1(
	Front		*fr)
{
	INTERFACE    *intfc = fr->interf;
	INTERFACE    *adj_intfc[2], *sav_intfc, *buf_intfc;
	PP_GRID	     *pp_grid = fr->pp_grid;
	RECT_GRID    *gr = fr->rect_grid;
	boolean	     sav_copy;
	boolean      status = FUNCTION_SUCCEEDED;
	double        *U = gr->U, *L = gr->L;
	double        *nor, p[3];
	double        T;
	int	     me[MAXD], him[MAXD];
	int	     myid, dst_id;
	int	     *G;
	int	     i, j, k;
	int	     dim = intfc->dim;
	static double nors[] = { 1.0,  0.0,  0.0,
			         0.0,  1.0,  0.0,
			         0.0,  0.0,  1.0,
			        -1.0,  0.0,  0.0,
			         0.0, -1.0,  0.0,
			         0.0,  0.0, -1.0};

	if (debugging("trace"))
	    (void) printf("Entering f_intfc_communication3d1()\n");

	set_floating_point_tolerance1(fr->rect_grid->h);
	sav_copy = copy_intfc_states();
	sav_intfc = current_interface();
	set_copy_intfc_states(YES);

	myid = pp_mynode();
	G = pp_grid->gmax;
	find_Cartesian_coordinates(myid,pp_grid,me);

	if (DEBUG)
	{
	    (void) printf("myid = %d, ",myid);
	    print_int_vector("me = ",me,dim,"\n");
	    print_PP_GRID_structure(pp_grid);
	    (void) printf("Input interface:\n");
	    print_interface(intfc);
	}

	/* Extend interface in three directions */

        construct_reflect_bdry(fr);
	clip_intfc_at_grid_bdry1(intfc);

	for (i = 0; i < dim; ++i)
	{
	    adj_intfc[0] = adj_intfc[1] = NULL;
	    for (j = 0; j < 2; ++j)
	    {
	    	pp_gsync();

		for (k = 0; k < dim; ++k)
		    him[k] = me[k];

		if (rect_boundary_type(intfc,i,j) == SUBDOMAIN_BOUNDARY)
		{
		    him[i] = (me[i] + 2*j - 1 + G[i])%G[i];
		    dst_id = domain_id(him,G,dim);
		    buf_intfc = cut_buf_interface1(intfc,i,j,me,him);
		    if ((j == 0) && (me[i] == 0))
		    {
			T = gr->GU[i] - gr->GL[i];
	                shift_interface(buf_intfc,T,i);
		    }
		    else if ((j == 1) && (me[i] == (G[i]-1)))
		    {
			T = gr->GL[i] - gr->GU[i];
	                shift_interface(buf_intfc,T,i);
		    }
		    if (dst_id != myid)
		    {
		    	send_interface(buf_intfc,dst_id);
		        (void) delete_interface(buf_intfc);
		    }
		    else
		        adj_intfc[(j+1)%2] = buf_intfc;
		}
		else if (rect_boundary_type(intfc,i,j) == REFLECTION_BOUNDARY)
		{
		    nor = nors + 3*i + 9*j;
		    buf_intfc = cut_buf_interface1(intfc,i,j,me,him);
		    p[i] = (j > 0) ? U[i] : L[i];
		    for (k = 1; k < dim; ++k)
			p[(k+i)%dim] = 0.5*(U[(k+i)%dim] + L[(k+i)%dim]);
		    reflect_interface(buf_intfc,p,nor);
		    adj_intfc[j] = buf_intfc;
		}
		if (rect_boundary_type(intfc,i,(j+1)%2) == SUBDOMAIN_BOUNDARY)
		{
		    him[i] = (me[i] - 2*j + 1 + G[i])%G[i];
		    dst_id = domain_id(him,G,dim);
		    if (dst_id != myid)
		    {
			adj_intfc[(j+1)%2] = receive_interface(dst_id);
		    }
		}

	    }
	    for (j = 0; j < 2; ++j)
	    {
		
		status = FUNCTION_SUCCEEDED;
		if (adj_intfc[j] != NULL)
		{
		    status = buffer_extension3d1(intfc,adj_intfc[j],i,j,status);
                   
		    if (!status)
                    {
		        (void) printf("WARNING in f_intfc_communication3d1 "
				      "buffer_extension3d1 failed for "
                                        "i = %d, j = %d \n", i, j);
                    }

		    (void) delete_interface(adj_intfc[j]);
		    set_current_interface(intfc);
		}    /*if (adj_intfc[j] != NULL) */
		
		/*surface in this direction can not mathch, return */
		status = pp_min_status(status);
		if(!status)
		    goto comm_exit;

	    }	
	    reset_intfc_num_points(intfc);
	}

	if (status == FUNCTION_SUCCEEDED)
	{
	    sep_common_pt_for_open_bdry(intfc);

	    install_subdomain_bdry_curves(intfc);

	    reset_intfc_num_points(intfc);
	    if (DEBUG)
	    {
	    	(void) printf("Final intfc:\n");
	    	print_interface(intfc);
	    }
	}

comm_exit:
	set_copy_intfc_states(sav_copy);
	set_current_interface(sav_intfc);
	if (debugging("trace"))
	    (void) printf("Leaving f_intfc_communication3d1()\n");
	return status;
}	/*end f_intfc_communication3d1 */

#define	MAX_SUBDOMAIN_TRIS	3000
EXPORT void install_subdomain_bdry_curves(
	INTERFACE	*intfc)
{
	BOND	 *b;
	BOND_TRI *btri;
	CURVE	 *curve, **c;
	NODE	 *ns, *ne, **n, **n1;
	POINT	 *p, *ps, *pe;
	TRI	 *tri_start, *start_tri;
	SURFACE	 **s;
	Locstate sl, sr;
	boolean  	 sav_intrp, dup_nodes;
	int	 side_start, start_side, ntris, i;
	int	 dim = intfc->dim;
	TRI	 **null_tris = NULL;
	int	 num_tri_allocated = 0;

	DEBUG_ENTER(install_subdomain_bdry_curves)
	
	sav_intrp = interpolate_intfc_states(intfc);
	interpolate_intfc_states(intfc) = NO;
	
	incremental_alloc_tris(&null_tris,&num_tri_allocated);
	for (s = intfc->surfaces; s && *s; ++s)
	{
	    while (null_side_on_surface(*s,&start_tri,&start_side))
	    {
		find_ending_null_side(start_tri,start_side,
				&tri_start,&side_start);
		
		ps = Point_of_tri(tri_start)[side_start];
	    	ns = make_node(ps);
		pe = Point_of_tri(tri_start)[Next_m3(side_start)];
		ne = make_node(pe);
		curve = make_curve(NO_COMP,NO_COMP,ns,ne);

		hsbdry_type(curve) = SUBDOMAIN_HSBDRY;
	    	install_curve_in_surface_bdry(*s,curve,POSITIVE_ORIENTATION);
	    	
		null_tris[0] = tri_start;
		ntris = 1;

	    	b = curve->first;
		while((tri_start = find_following_null_tri(tri_start,&p,
					&side_start,CLOCKWISE)) != NULL)
		{
		    if (insert_point_in_bond(ns->posn,b,curve) !=
			FUNCTION_SUCCEEDED)
		    {
			screen("ERROR in install_subdomain_bdry_curves(), "
		               "insert_point_in_bond() failed\n");
			clean_up(ERROR);
		    }
		    
		    ns->posn = b->start = p;
		    set_bond_length(b,dim);
	
		    null_tris[ntris] = tri_start;
		    ntris++;
		    if(ntris >= num_tri_allocated)
			incremental_alloc_tris(&null_tris,&num_tri_allocated);
	    
		    if (ns->posn == ne->posn) /*Closed loop formed */
		    {
			change_node_of_curve(curve,POSITIVE_ORIENTATION,ne);
			if(!delete_node(ns))
			{
			    printf("ERROR install_subdomain_bdry_curves "
			    	   "can not delete ns node for a loop.\n");
			    print_node(ns);
			    clean_up(ERROR);
			}
			break;
		    }
		}
		
		for(i=0, b=curve->first; b!=NULL; b=b->next, i++)
		{
		    btri = link_tri_to_bond(NULL,null_tris[ntris-1-i],*s,b,
						curve);
		    copy_tri_state_to_btri(btri,b,NEGATIVE_ORIENTATION,intfc);
		    copy_tri_state_to_btri(btri,b,POSITIVE_ORIENTATION,intfc);
		}
	    }   /*null_side_on_surface.*/
	}

	dup_nodes = YES;
	while(dup_nodes)
	{
	    dup_nodes = NO;
	    for(n=intfc->nodes; n && *n; n++)
	    {
		for(n1=n+1; n1 && *n1; n1++)
		{
		    if((*n)->posn != (*n1)->posn)
			continue;
		    dup_nodes = YES;

		    /*merge nodes *n1 with node *n */
		    for(c=(*n1)->in_curves; c && *c; c++)
			change_node_of_curve(*c,NEGATIVE_ORIENTATION,*n);
		    for(c=(*n1)->out_curves; c && *c; c++)
			change_node_of_curve(*c,POSITIVE_ORIENTATION,*n);
		    
		    if(!delete_node(*n1))
		    {
			printf("ERROR install_subdomain_bdry_curves "
		    	       "can not delete node.\n");
			print_node(*n1);
			clean_up(ERROR);
		    }
		    break;
		}
		if(dup_nodes)
		    break;
	    }
	}
	free_these(1,null_tris);
	reset_intfc_num_points(intfc);

	if (debugging("consistency"))
	{
	    if (!consistent_interface(intfc))
	    {
		screen("ERROR in install_subdomain_bdry_curves(), "
		       "intfc is inconsistent\n");
		clean_up(ERROR);
	    }
	}
	interpolate_intfc_states(intfc) = sav_intrp;
	DEBUG_LEAVE(install_subdomain_bdry_curves)
}		/*end install_subdomain_bdry_curves*/

void sep_common_pt_for_open_bdry(
	INTERFACE	*intfc)
{
	SURFACE	**s;
	POINT	*p;
	TRI	*tri;
	int	i, k, nt;
	double	nor[4], tmp_fac;
	TRI	**tris = NULL;
	int	num_tri_allocated = 0;

	if(!debugging("sep_for_open"))
	    return;

	nt = 0;
	incremental_alloc_tris(&tris,&num_tri_allocated);
	for(s=intfc->surfaces; s && *s; s++)
	{
	    for(tri=first_tri(*s); !at_end_of_tri_list(tri,*s); tri=tri->next)
	    {
		for(i=0; i<3; i++)
		{
		    p = Point_of_tri(tri)[i];
		    if(point_outside_open_bdry(&k, nor, p, intfc))
		    {
		        tris[nt++] = tri;
			if(nt >= num_tri_allocated)
			{
			     incremental_alloc_tris(&tris,&num_tri_allocated);
			}
		    }
		}
	    }
	}
	
	tmp_fac = get_comm_pt_fac();
	set_comm_pt_fac(0.999);
	if(nt > 0)
	    sep_common_point_from_loop(tris, nt, NULL, NULL, intfc);
	set_comm_pt_fac(tmp_fac);
	free_these(1,tris);
}

LOCAL boolean null_side_on_surface(
	SURFACE *s,
	TRI     **tri_start,
	int     *side_start)
{
	TRI *tri;
	int i;

	for (tri = first_tri(s); !at_end_of_tri_list(tri,s); tri = tri->next)
	{
	    for (i = 0; i < 3; ++i)
	    {
		if (Tri_on_side(tri,i) == NULL)
		{
		    *side_start = i;
		    *tri_start = tri;
		    return YES;
		}
	    }
	}
	return NO;
}	/* end null_side_on_surface */

LOCAL boolean buffer_extension3d1(
	INTERFACE	*intfc,
	INTERFACE	*adj_intfc,
	int		dir,
	int		nb,
	boolean		status)
{
	RECT_GRID	*gr = computational_grid(intfc);

	if (debugging("trace"))
	    (void) printf("Entering buffer_extension3d1()\n");

	set_current_interface(intfc);

		/* Patch tris from adj_intfc to intfc */

	if (!append_adj_intfc_to_buffer1(intfc,adj_intfc,gr,dir,nb))
	{
	    status = FUNCTION_FAILED;
	    (void) printf("WARNING in buffer_extension3d1(), "
	                  "append_adj_intfc_to_buffer1() failed\n");
	}

	if (debugging("trace"))
	    (void) printf("Leaving buffer_extension3d1()\n");
	return status;
}	/*end buffer_extension3d1*/

EXPORT	void shift_interface(
	INTERFACE *intfc,
	double     T,
	int       dir)
{
	POINT   *p;
	SURFACE **s;
	NODE	**n;
	CURVE	**c;
	TRI     *t;
	BOND	*b;
	int     i;

	if (debugging("trace"))
	    (void) printf("Entering shift_interface()\n");

	/* Reset sort status of all points on interface */
	intfc_node_loop(intfc,n)
	    sorted((*n)->posn) = NO;

	intfc_curve_loop(intfc,c)
	    curve_bond_loop(*c,b)
	    	sorted(b->end) = NO;

	intfc_surface_loop(intfc,s)
	    surf_tri_loop(*s,t)
		for (i = 0; i < 3; ++i)
		    sorted(Point_of_tri(t)[i]) = NO;

	/* Shift points on interface */
	intfc_node_loop(intfc,n)
	{
	    p = (*n)->posn;
	    Coords(p)[dir] += T;
	    p->pshift[dir] += T;
	    sorted(p) = YES;	    
	}
	intfc_curve_loop(intfc,c)
	{
	    curve_bond_loop(*c,b)
	    {
		p = b->end;
		if (sorted(p)) continue;
	    	Coords(p)[dir] += T;
	    	p->pshift[dir] += T;
	    	sorted(p) = YES;	    
	    }
	}
	intfc_surface_loop(intfc,s)
	{
	    for (t = first_tri(*s); !at_end_of_tri_list(t,*s); t = t->next)
	    {
		for (i = 0; i < 3; ++i)
		{
		    p = Point_of_tri(t)[i];
		    if (sorted(p)) continue;
	            Coords(p)[dir] += T;
		    p->pshift[dir] += T;
		    sorted(p) = YES;
		}
	    }
	}

	if (debugging("trace"))
	    (void) printf("Leaving shift_interface()\n");
}		/*end shift_interface*/

LOCAL 	int append_adj_intfc_to_buffer1(
	INTERFACE	*intfc,		/* my interface 	       */
	INTERFACE	*adj_intfc,	/* received interface	       */
	RECT_GRID	*grid,		/* Rectangular grid for region */
	int		dir,
	int		nb)
{
	INTERFACE	*cur_intfc;
	SURFACE		**s, **as;
	CURVE		**ac;
	int		p_size;		/*Size of space allocated for p_table*/
	static P_LINK	*p_table = NULL;/* Table of matching points on intfc
					 * and adj_intfc*/
	static int      len_p_table = 0;
	boolean      	corr_surf_found;

	if (debugging("append_curves"))
	{
	    (void) printf("Entering append_adj_intfc_to_buffer1()\n");
            (void) printf("Number of curves = %d\n",I_NumOfIntfcCurves(intfc));
	}

	cur_intfc = current_interface();
	set_current_interface(intfc);

	p_size = 4*(adj_intfc->num_points) + 1;
	if (len_p_table < p_size)
	{
	    len_p_table = p_size;
	    if (p_table != NULL)
		free(p_table);
	    uni_array(&p_table,len_p_table,sizeof(P_LINK));
	}
	
	reset_hash_table(p_table,p_size);
	
	/* Begin patching adj_intfc to current interface */
	for (as = adj_intfc->surfaces; as && *as; ++as)
	{
	    corr_surf_found = NO;
	    for (s = intfc->surfaces; s && *s; ++s)
	    {
		/*
		*  COMMENT -
		*  The Hyper_surf_index() function is not
		*  fully supported.  This will fail in the
		*  presence of interface changes in topology
		*  TODO: FULLY SUPPORT THIS OBJECT
		*/
                
		/*  assume one neg and pos comp gives ONLY ONE surf */
		if (surfaces_matched(*as,*s))
		{
		    corr_surf_found = YES;
		    if (!append_buffer_surface1(*s,*as,grid,dir,nb,p_table,
						   p_size))
		    {
			set_current_interface(cur_intfc);
			
			(void) printf("WARNING in "
			              "append_adj_intfc_to_buffer1(), "
			              "append surface failed\n");
			
			DEBUG_LEAVE(append_adj_intfc_to_buffer1)
			return NO;
		    }
		}
	    }
	    if (!corr_surf_found && as && *as)
	    {
	    	SURFACE *surf;
		if (as==NULL)
		    continue;
		surf = copy_buffer_surface(*as,p_table,p_size);
		Hyper_surf_index(surf) = Hyper_surf_index((*as));
	    }
	}
	append_other_curves1(intfc,adj_intfc,grid,dir,nb,p_table,p_size);
	/* append curves not on surfaces */
        for (ac = adj_intfc->curves; ac && *ac; ++ac)
	{
            matching_curve(*ac,p_table,p_size);
	    merge_gindex_of_curve(*ac,grid,dir,nb,p_table,p_size);
	}
	
	merge_curves(intfc,adj_intfc,NO);
	
	set_current_interface(cur_intfc);
	
	if (debugging("append_curves"))
	{
	    (void) printf("Leaving append_adj_intfc_to_buffer1()\n");
            (void) printf("Number of curves = %d\n",I_NumOfIntfcCurves(intfc));
	}
	DEBUG_LEAVE(append_adj_intfc_to_buffer1)
	return YES;
}		/*end append_adj_intfc_to_buffer1*/

EXPORT void set_floating_point_tolerance1(
	double		*h)
{
	double eps;
	static const double mfac = 10.0;/*TOLERANCE*/
	static const double gfac = 0.00004;/*TOLERANCE*/
	int   i;

	eps = mfac*MACH_EPS;
	for (i = 0; i < 3; ++i)
	{
	    tol1[i] = gfac*h[i];/*TOLERANCE*/
	    tol1[i] = max(tol1[i],eps);
	}
}		/*end set_floating_point_tolerance1*/

LOCAL int append_buffer_surface1(
	SURFACE		*surf,
	SURFACE		*adj_surf,
	RECT_GRID	*grid,
	int		dir,
	int		nb,
	P_LINK		*p_table,
	int		p_size)
{
	SURFACE    *new_adj_surf;
	TRI	   *tri;
	CURVE	   **c;
	int	   ns, na, new_na;
	double	   crx_coord;
	static TRI **tris_s = NULL, **tris_a = NULL;
	static int len_tris_s = 0, len_tris_a = 0;

	DEBUG_ENTER(append_buffer_surface1)
	
	if (DEBUG)
	{
	    (void) printf("dir = %d,nb = %d\n",dir,nb);
	    (void) printf("My surface\n");
	    print_surface(surf);
	    (void) printf("Buffer surface\n");
	    print_surface(adj_surf);
	}
	
	if (len_tris_s < surf->num_tri)
	{
	    len_tris_s = surf->num_tri;
	    if (tris_s != NULL)
		free(tris_s);
	    uni_array(&tris_s,len_tris_s,sizeof(TRI *));
	}
	if (len_tris_a < adj_surf->num_tri)
	{
	    len_tris_a = adj_surf->num_tri;
	    if (tris_a != NULL)
		free(tris_a);
	    uni_array(&tris_a,len_tris_a,sizeof(TRI *));
	}

	crx_coord = (nb == 0) ? grid->L[dir] : grid->U[dir];

	ns = na = 0;
	for (tri=first_tri(surf); !at_end_of_tri_list(tri,surf); tri=tri->next)
	{
	    if (tri_set_cross_line(tri,crx_coord,dir) == YES)
	    {
		tris_s[ns++] = tri;
	    }
	}

	for (tri = first_tri(adj_surf); !at_end_of_tri_list(tri,adj_surf); 
	     tri = tri->next)
	{
	    if (tri_set_cross_line(tri,crx_coord,dir) == YES)
	    {
		tris_a[na++] = tri;
	    }
	}

	/* Add matching points to the hashing table p_table */

	if (add_matching_pt_to_hash_table(tris_s,tris_a,ns,na,surf,
		      adj_surf,p_table,p_size))
	{
	    merge_gindex_of_tris(tris_a,na,grid,dir,nb,p_table,p_size);
	    new_adj_surf = copy_buffer_surface(adj_surf,p_table,p_size);
	    adj_surf = new_adj_surf;
	}
	else if (debugging("match_tris"))
	{
	    printf("dir = %d  nb = %d\n",dir,nb);
	    gview_plot_surface("surf_s",surf);
	    gview_plot_surface("surf_a",adj_surf);
	    gview_plot_tri_list("tris_s",tris_s,ns);
	    gview_plot_tri_list("tris_a",tris_a,na);
	    return NO;
	}

	synchronize_tris_at_subdomain_bdry(tris_a,tris_s,ns,p_table,p_size);
	
	na = 0;
	for (tri = first_tri(adj_surf); !at_end_of_tri_list(tri,adj_surf);
	     tri = tri->next)
	{
	    if (tri_set_cross_line(tri,crx_coord,dir) == YES)
	    {
	        tris_a[na++] = tri;
	    }
	}

	/*adjoin adj_surf tri list to surf tri list */
	last_tri(surf)->next = first_tri(adj_surf);
	first_tri(adj_surf)->prev = last_tri(surf);
	link_tri_list_to_surface(first_tri(surf),last_tri(adj_surf),surf);
	
	/* BEGIN curves in adj_surf should be added to surf. */
	for (c = adj_surf->pos_curves; c && *c; c++)
	    if (!delete_from_pointers(adj_surf,&(*c)->pos_surfaces))
	    {
	        printf("ERROR: in append_buffer_surface1, "
		       "adj_surf and pos_curves are not paired.\n");
		clean_up(ERROR);
	    }
	    else
	        install_curve_in_surface_bdry(surf,*c,POSITIVE_ORIENTATION);

	for (c = adj_surf->neg_curves; c && *c; c++)
	    if (!delete_from_pointers(adj_surf,&(*c)->neg_surfaces))
	    {
	        printf("ERROR: in append_buffer_surface1, "
		       "adj_surf and neg_curves are not paired.\n");
		clean_up(ERROR);
	    }
	    else
	        install_curve_in_surface_bdry(surf,*c,NEGATIVE_ORIENTATION);

	adj_surf->pos_curves = adj_surf->neg_curves = NULL;
	(void) delete_surface(adj_surf);
	adj_surf = NULL;

	/* average_btris*/
	if (!match_tris_at_subdomain_bdry(surf,adj_surf,tris_s,tris_a,ns,na))
	{
	    (void) printf("WARNING in append_buffer_surface1(), "
	                  "no match of tris at subdomain\n");
	    (void) printf("dir = %d, nb = %d\n",dir,nb);
	    DEBUG_LEAVE(append_buffer_surface1)
	    return NO;
	}

	DEBUG_LEAVE(append_buffer_surface1)
	return YES;
}		/*end append_buffer_surface1*/


LOCAL	void	synchronize_tris_at_subdomain_bdry(
	TRI    **tris_a,
	TRI    **tris_s,
	int    nt,
	P_LINK *p_table,
	int    p_size)
{
	POINT **ps, **pa, *p0, *p1, *p2;
	TRI   *ts, *ta;
	int   i, j, id, idp, idn;

	for (i = 0; i < nt; ++i)
	{
	    ts = tris_s[i];
	    ps = Point_of_tri(ts);
	    for (j = 0; j < nt; ++j)
	    {
		ta = tris_a[j];
		pa = Point_of_tri(ta);
		for (id = 0; id < 3; ++id)
		{
		    p0 = (POINT*) find_from_hash_table((POINTER)pa[id],
						       p_table,p_size);
		    if (p0 == ps[0])
		    {
		        idn = Next_m3(id);
		        p1 = (POINT*) find_from_hash_table((POINTER)pa[idn],
						           p_table,p_size);
		        if (p1 == ps[1])
		        {
		            idp = Prev_m3(id);
		            p2 = (POINT*) find_from_hash_table((POINTER)pa[idp],
						               p_table,p_size);
			    if (p2 == ps[2])
			    {
			        rotate_triangle(ts,(3-id)%3);
			        set_normal_of_tri(ts);
			        set_normal_of_tri(ta);
			        break;
			    }
		        }
		    }
		}
		if (id < 3)
		    break;
	    }
	    if(j == nt)
	    {
	        printf("WARNING, synchronize_tris_at_subdomain_bdry, "
		       "suitable triangle not found.\n");
	        print_tri(ts, ts->surf->interface);
	    }
	}
}		/*end synchronize_tris_at_subdomain_bdry*/

LOCAL	boolean	tri_cross_line(
	TRI		*tri,
	double		crx_coord,
	int		dir)
{
	double	min_coord, max_coord;
	double	crx_tol = tol1[dir];

	min_coord = max_coord = Coords(Point_of_tri(tri)[0])[dir];

	if (min_coord > Coords(Point_of_tri(tri)[1])[dir])
	    min_coord = Coords(Point_of_tri(tri)[1])[dir];
	if (min_coord > Coords(Point_of_tri(tri)[2])[dir])
	    min_coord = Coords(Point_of_tri(tri)[2])[dir];

	if (max_coord < Coords(Point_of_tri(tri)[1])[dir])
	    max_coord = Coords(Point_of_tri(tri)[1])[dir];
	if (max_coord < Coords(Point_of_tri(tri)[2])[dir])
	    max_coord = Coords(Point_of_tri(tri)[2])[dir];
	
	return (((min_coord - crx_coord) <= crx_tol) && 
	        ((crx_coord - max_coord) <= crx_tol)) ? YES : NO;
}		/*end tri_cross_line*/

/*NO    no bond side or no tri crx line
  YES   one tri crx line
*/

EXPORT boolean tri_set_cross_line(
	TRI		*tri,
	double		crx_coord,
	int		dir)
{
	int	i,j,n,l,num_bonds;
	TRI	**tris,*t;
	POINT	*p;
	BOND	*b,*bonds[2];
	BOND_TRI **btris;

    	if (tri_cross_line(tri,crx_coord,dir))
	    return YES;
	for (i = 0; i < 3; i++)
        {
            p = Point_of_tri(tri)[i];
            if (Boundary_point(p))
            {
		bonds_of_boundary_point(p,tri,bonds,&num_bonds);
                if (num_bonds == 0) continue;
		for (l = 0; l < num_bonds; ++l)
		{
		  b = bonds[l];
                  for (btris = Btris(b); btris && *btris; ++btris)
                  {
                    t = (*btris)->tri;
                    n = set_tri_list_around_point(p,t,&tris,t->surf->interface);
		    for (j = 0; j < n; ++j)
		    {
		    	if (tri_cross_line(tris[j],crx_coord,dir))
			    return YES;
		    }
		  }
		}
	    }
	}
	return NO;
}	/* end tri_bond_cross_line */


LOCAL boolean match_tris_at_subdomain_bdry(
	SURFACE		*surf,
	SURFACE		*adj_surf,
	TRI		**tri_s,
	TRI		**tri_a,
	int		ns,
	int		na)
{
	TRI		*ta, *ts;
	int		i, j;
	int		ums,uma;
	static boolean	*ms = NULL, *ma = NULL;
	static int      ms_len = 0, ma_len = 0;

	DEBUG_ENTER(match_tris_at_subdomain_bdry)
	if (DEBUG)
	{
	    (void) printf("Entered match_tris_at_subdomain_bdry()\n");
	    (void) printf("tri_a: na=%d\n",na);
	    for (i = 0; i < na; ++i)
		print_tri(tri_a[i],adj_surf->interface);
	    (void) printf("tri_s: ns=%d\n",ns);
	    for (i = 0; i < ns; ++i)
		print_tri(tri_s[i],surf->interface);
	}

	if (ms_len < ns)
	{
	    ms_len = ns;
	    if (ms != NULL)
		free(ms);
	    uni_array(&ms,ms_len,sizeof(boolean));
	}
	if (ma_len < na)
	{
	    ma_len = na;
	    if (ma != NULL)
		free(ma);
	    uni_array(&ma,ma_len,sizeof(boolean));
	}

	for (i = 0; i < ns; ++i)
	    ms[i] = NO;
	for (i = 0; i < na; ++i)
	    ma[i] = NO;

	for (i = 0; i < na; ++i)
	{
	    ta = tri_a[i];
	    
	    for (j = 0; j < ns; ++j)
	    {
	        ts = tri_s[j];
	        if (match_two_tris(ts,ta))
		{
		    ma[i] = ms[j] = YES;
		    
		    if(debugging("app_tri"))
		        printf("ns  %d \n", j);

	    	    average_btris(ts, surf, ta, adj_surf);
		    merge_two_tris(ts,ta,surf,adj_surf);
		    break;
		}
	    }

	    remove_from_debug("app_tri");
	}
	ums = uma = 0;
	for (i = 0; i < ns; ++i)
	    if (ms[i] == NO)
		++ums;
	for (i = 0; i < na; ++i)
	    if (ma[i] == NO)
		++uma;
	if (ums != 0 || uma != 0)
	{
	    return NO;
	    (void) printf("WARNING in match_tris_at_subdomain_bdry(), "
	                  "unmatched local tris\n"
	                  "na = %d ns = %d\n"
	                  "ums = %d uma = %d\n",na,ns,ums,uma);
	    for (i = 0; i < ns; ++i)
	    	if (ms[i] == NO)
		    print_tri(tri_s[i],surf->interface);
	    (void) printf("unmatched adj tris:\n");
	    for (i = 0; i < na; ++i)
	    	if (ma[i] == NO)
		    print_tri(tri_a[i],adj_surf->interface);
	    
	    DEBUG_LEAVE(match_tris_at_subdomain_bdry)
	    clean_up(ERROR);
	    return NO;
	}

	DEBUG_LEAVE(match_tris_at_subdomain_bdry)
	return YES;
}		/*end match_tris_at_subdomain_bdry*/

/*
*			match_two_tris():
*
*	Determines whether two triangles are the same triangle up to a rotation
*	of the vertex indices.  Returns YES if the tris are the same,
*	otherwise returns NO.
*/

LOCAL boolean match_two_tris(
	TRI		*tri,
	TRI		*atri)
{
	int		j;
	POINT		**p, **ap;

	p = Point_of_tri(tri);
	ap = Point_of_tri(atri);
	for (j = 0; j < 3; ++j)
	    if (p[0] == ap[j])
		break;
	if (j == 3)
	    return NO;
	return ((p[1]==ap[Next_m3(j)]) && (p[2]==ap[Prev_m3(j)])) ? YES : NO;
}		/*end match_two_tris*/

/*
*			merge_two_tris():
*
*	Merges two triangles ts and ta by replacing all linking
*	information to and from ta by linking with ts  while preserving
*	any nontrivial linking information in ts.
*/

EXPORT	void merge_two_tris(
	TRI	*ts,
	TRI	*ta,
	SURFACE	*s,
	SURFACE *as)
{
	TRI		*tri;
	BOND_TRI	*btri;
	int		i, j;

	DEBUG_ENTER(merge_two_tris)
	if (DEBUG)
	{
	    print_tri(ts,s->interface);
	    print_tri(ta,as->interface);
	}

	/* Set the links on the null side of ts by corresponding links to ta*/
	for (i = 0; i < 3; ++i)
	{
	    /* for grid based, this part is not necessary*/
            if(is_side_bdry(ta,i))
            {
	       if(Bond_on_side(ts,i) == NULL)
               {
                   if(Bond_on_side(ta,i) == NULL)
                           continue;
		       
	       	   printf("#merge_two_tris: bond on block face.\n");
		       /*set tri_neighbor btri*/
                   Bond_tri_on_side(ts,i) = btri = Bond_tri_on_side(ta,i);
                   btri->tri = ts;
                   btri->surface = s;
               }
	       continue;
            }

	    if (Tri_on_side(ts,i) == NULL)
	    {
	    	if (Tri_on_side(ta,i) == NULL)
		{
		    continue;
		}
		/* set tri_neighbor tri*/
	    	Tri_on_side(ts,i) = tri = Tri_on_side(ta,i);
		
		for (j = 0; j < 3; ++j)
	    	{
	    	    if (Tri_on_side(tri,j) == ta)
	    	    	Tri_on_side(tri,j) = ts;
	    	}
	    }
	}
	
	/*remove ta from s tri list */
	remove_tri_from_surface(ta,s,NO);
	if (DEBUG)
	{
	    (void) printf("after merge_two_tris\n");
	    print_tri(ts,s->interface);
	}
	DEBUG_LEAVE(merge_two_tris)
}		/*end merge_two_tris*/



/*
*			add_matching_pt_to_hash_table():
*
*	Creates a hashed list of matching points on intfc and adj_intfc,
*	where two points match if their positions are within a prescribed
*	floating point tolerance.  This is an extremely expensive function
*	and should be a candidate for an improved version.
*/

LOCAL boolean add_matching_pt_to_hash_table(
	TRI 	**tris_s,
	TRI	**tris_a,
	int	ns,
	int 	na,
	SURFACE	*ss,
	SURFACE	*sa,
	P_LINK	*p_table,
	int	p_size)
{
	boolean                  status;
	int 	                 i, tstn;
	POINT_LIST               *plists, *plista, *pls, *pla;
	static POINT_LIST_STORE  Pslist_store, Palist_store;

	plists = set_point_list(tris_s,ns,Hyper_surf(ss),&Pslist_store);
	plista = set_point_list(tris_a,na,Hyper_surf(sa),&Palist_store);

	tstn = 0;
	status = YES;
	for (pls = plists; pls != NULL; pls = pls->next)
	{
	    for (pla = plista; pla != NULL; pla = pla->next)
	    {
		for (i = 0; i < 3; ++i) /*Floating point TOLERANCE test*/
	            if (fabs(Coords(pla->p)[i]-Coords(pls->p)[i]) > tol1[i])
			break;
		if (i == 3)
		{
		    tstn++;
	            (void) add_to_hash_table((POINTER)pla->p,(POINTER)pls->p,
				             p_table,p_size);
		    set_full_average(NO);
		    (void) average_points(NO,pla->p,pla->hse,pla->hs,
				          pls->p,pls->hse,pls->hs);
		    set_full_average(YES);
		    
		    if (pla->prev == NULL)
			plista = pla->next;
		    else
			pla->prev->next = pla->next;
		    if (pla->next != NULL)
			pla->next->prev = pla->prev;
		    break;
		}
	    }
	    /*can not find a corresponding point in plista, in this case, 
	      we continue search just set status = NO.  */
	    if (pla == NULL) 
	        status = NO;
	}
	if (plista != NULL)
	    status = NO;

	if(status == NO)
	{
	    double	*pa, *ps;
	    double	len, min_len;

	    (void) printf("In add_matching_pt_to_hash_table(), status = NO\n");
	    for (pla = plista; pla != NULL; pla = pla->next)
	    {
	        min_len = HUGE_VAL;
		if (plists == NULL)
		{
		    (void) printf("plists is NULL\n");
		    continue;
		}
		for (pls = plists; pls != NULL; pls = pls->next)
	        {
		    len = distance_between_positions(Coords(pla->p),
				Coords(pls->p), 3);
		    if(len < min_len)
		    {
		        min_len = len;
			pa = Coords(pla->p);
			ps = Coords(pls->p);
		    }
		}
		
		(void) printf("Closest pair:\n");
		print_general_vector("pa", pa, 3, "\n");
		print_general_vector("ps", ps, 3, "\n");
		(void) printf("\n");
	    }
	}

	return (ns == na) ? status : NO;
}		/*end add_matching_pt_to_hash_table*/

boolean the_tri_rot(TRI *);

LOCAL	POINT_LIST	*set_point_list(
	TRI	         **tris,
	int              nt,
	HYPER_SURF       *hs,
	POINT_LIST_STORE *plist_store)
{
	POINT      **p;
	POINT_LIST *plist, *pl, Plist;
	TRI        *tri;
	int        i, j, max_np, tstnum;

	tstnum = 0;
	if (nt == 0)
	    return NULL;
	max_np = 3*nt + 1;
	if (max_np > plist_store->len)
	{
	    if (plist_store->pl != NULL)
		free(plist_store->pl);
	    plist_store->len = max_np;
	    uni_array(&plist_store->pl,sizeof(POINT_LIST),plist_store->len);
	}
	zero_scalar(plist_store->pl,sizeof(POINT_LIST)*plist_store->len);
	for (i = 0; i < nt; ++i)
	{
	    p = Point_of_tri(tris[i]);
	    for (j = 0; j < 3; ++j)
	        sorted(p[j]) = NO;
	}

	plist = plist_store->pl;
	pl = &Plist;
	for (i = 0; i < nt; ++i)
	{
	    tri = tris[i];
	    
	    p = Point_of_tri(tri);
	    for (j = 0; j < 3; ++j)
	    {
		if (!sorted(p[j]))
		{
		    pl->next = plist++;
		    pl->next->prev = pl;
		    pl = pl->next;
	            pl->p = p[j];
	            pl->hs = hs;
	            pl->hse = Hyper_surf_element(tri);
	            sorted(pl->p) = YES;
		    tstnum++; 
		}
	    }
	}
	Plist.next->prev = NULL;

	return Plist.next;
}		/*end set_point_list*/

/*copy surface as to s and replace the points in p_table in surface s.
  It will also put the newly generated points in p_table.
*/

EXPORT	SURFACE *copy_buffer_surface(
	SURFACE		*as,
	P_LINK		*p_table,
	int		p_size)
{
	SURFACE		*s;
	POINT		*p, *np, *ap;
	TRI		*tri, *atri;
	CURVE		**c, **pos_curves = NULL, **neg_curves = NULL;
	int		i;

	DEBUG_ENTER(copy_buffer_surface)

	/* First identify and copy nodes on as */

	for (c = as->pos_curves; c && *c; ++c)
	{
	    if (!add_to_pointers(matching_curve(*c,p_table,p_size),
				    &pos_curves))
	    {
	        screen("ERROR in copy_buffer_surface(), "
	               "add_to_pointers() failed\n");
	        clean_up(ERROR);
	    }
	}
	for (c = as->neg_curves; c && *c; ++c)
	{
	    if (!add_to_pointers(matching_curve(*c,p_table,p_size),
				    &neg_curves))
	    {
	        screen("ERROR in copy_buffer_surface(), "
	               "add_to_pointers() failed\n");
	        clean_up(ERROR);
	    }
	}

	/* set tri_array_numbers */
	for(i=0, tri=first_tri(as); !at_end_of_tri_list(tri,as); tri=tri->next)
	    Tri_index(tri) = i++;

	s = copy_surface(as,pos_curves,neg_curves,YES);

	/* null_tri_array_numbers */
	for(i=0, tri=first_tri(as); !at_end_of_tri_list(tri,as); tri=tri->next)
	    Tri_index(tri) = 0;

	for (tri = first_tri(s); !at_end_of_tri_list(tri,s); tri = tri->next)
	{
	    sorted(Point_of_tri(tri)[0]) = NO;
	    sorted(Point_of_tri(tri)[1]) = NO;
	    sorted(Point_of_tri(tri)[2]) = NO;
	}
	for (tri = first_tri(s), atri = first_tri(as);
	     !at_end_of_tri_list(tri,s); tri = tri->next, atri = atri->next)
	{
	    for (i = 0; i < 3; ++i)
	    {
	    	p = Point_of_tri(tri)[i];
	    	if (sorted(p) == NO)
	    	{
	    	    if (!vertex_on_bond(tri,i))  /* the if is useless */
		    {
			ap = Point_of_tri(atri)[i];
			np = (POINT*)find_from_hash_table((POINTER)ap,
					                  p_table,p_size);
			if (np == NULL)
			{
			    (void) add_to_hash_table((POINTER)ap,(POINTER)p,
						     p_table,p_size);
			}
			else if (p != np)
			{
			    p = np;
			    Point_of_tri(tri)[i] = p;
			}
		    }
		    sorted(p) = YES;
		}
	    }
	}
	DEBUG_LEAVE(copy_buffer_surface)
	return s;
}		/*end copy_buffer_surface*/

EXPORT	CURVE *matching_curve(
	CURVE		*ac,
	P_LINK		*p_table,
	int		p_size)
{
	BOND		*b, *ab;
	CURVE		**cc, *c;
	NODE		*ns, *ne;
	POINT		*ap, *p;

	c = NULL;
	ns = matching_node(ac->start,p_table,p_size);
	ne = matching_node(ac->end,p_table,p_size);
	for (cc = ns->out_curves; cc && *cc; ++cc)
	    if (curves_match(*cc,ac,p_table,p_size))
	    {
		/*
		if (Gindex(*cc) == 28)
		{
		    printf("In matching_curve() curves_match() found curve\n");
		    print_curve(*cc);
		}
		*/
		return *cc;
	    }
	
	c = copy_curve(ac,ns,ne);
	
	/*matching_node already put the starting and ending points in 
 	  the hash table do not deal the starting and the ending points here.
	*/
	for (b = c->first, ab = ac->first; b != c->last && ab != ac->last;
		b = b->next, ab = ab->next)
	{
	    ap = ab->end;
	    p = (POINT*) find_from_hash_table((POINTER)ap,p_table,p_size);
	    if (p != NULL)
	    {
	    	b->end = p;
		if (b->next)
	    	    b->next->start = p;
	    }
	    else
	    {
	    	(void) add_to_hash_table((POINTER)ap,(POINTER)b->end,
					 p_table,p_size);
	    }
	}
	
	test_curve_link(c);

	return c;
}		/*end matching_curve*/

EXPORT	boolean curves_match(
	CURVE		*c,
	CURVE		*ac,
	P_LINK		*p_table,
	int		p_size)
{
	BOND		*b, *ab;
	POINT		*p, *ap;

	ap = ac->start->posn;
	p = (POINT*)find_from_hash_table((POINTER)ap,p_table,p_size);
	if (p != c->start->posn)
	    return NO;
	for (b = c->first, ab = ac->first; b && ab; b = b->next, ab = ab->next)
	{
	    ap = ab->end;
	    p = (POINT*)find_from_hash_table((POINTER)ap,p_table,p_size);
	    if (p != b->end)
		return NO;
	}
	
	if(ab == NULL)
	    return YES;	 /* curve in orginal intfc is longer or equal */
	return NO;
}		/*end curves_match*/

EXPORT	NODE *matching_node(
	NODE		*an,
	P_LINK		*p_table,
	int		p_size)
{
	INTERFACE	*intfc = current_interface();
	NODE		*newn;
	POINT		*p, *ap = an->posn;

	p = (POINT*)find_from_hash_table((POINTER)ap,p_table,p_size);
	if (p == NULL)
	{
	    p = copy_point(ap);
	    (void) add_to_hash_table((POINTER)ap,(POINTER)p,p_table,p_size);
	}
	newn = node_of_point(p,intfc);
	if (newn == NULL)
	    newn = make_node(p);
	return newn;
}		/*end matching_node*/

EXPORT void open_surf_null_sides(
	SURFACE *surf,
	double *L,
	double *U,
	int dir,
	int nb)
{
	TRI *tri,*ntri;
	int i;
	POINT *p;

	ntri = NULL;
	for (tri=first_tri(surf); !at_end_of_tri_list(tri,surf); tri=ntri)
	{
	    ntri = tri->next;
	    /* bond with inside tri survives. */
	    if (tri_set_out_domain(tri,L,U,dir,nb))
	    {
	    	remove_out_domain_tri(tri,surf);
	    }
	}
}	/* end open_surf_null_sides */

EXPORT void open_null_bonds(
	INTERFACE *intfc,
	double *L,
	double *U,
	int dir,
	int nb)
{
	CURVE **c,**curves;
	BOND *b;
	NODE *node;
	int i,nc;

	nc = I_NumOfIntfcCurves(intfc);
	uni_array(&curves,3*nc,sizeof(CURVE*));
	set_floating_point_tolerance1(computational_grid(intfc)->h);

start_open_bonds:
	nc = 0;
	intfc_curve_loop(intfc,c)
	    curves[nc++] = *c;
	for (i = 0; i < nc; ++i)
	{
	    curve_bond_loop(curves[i],b)
	    {
		if (!bond_out_domain1(b,L,U,dir,nb)) continue;
		if (curves[i]->first == curves[i]->last)  /* only bond */
		{
		    delete_curve(curves[i]);
		    delete_node(curves[i]->start);
		    delete_node(curves[i]->end);
		    goto start_open_bonds;
		}
		if (is_closed_curve(curves[i]))
		{
		    I_MoveNodeToPoint(b->start,curves[i]);	
		    curves[i]->first = b->next;
		    curves[i]->first->prev = NULL;
		    node = make_node(b->end);
		    change_node_of_curve(curves[i],POSITIVE_ORIENTATION,node);
		    curves[i]->num_points--;
		}
		else if (b == curves[i]->first)
		{
		    curves[i]->first = b->next;
		    curves[i]->first->prev = NULL;
		    curves[i]->num_points--;
		    if (I_NumOfNodeCurves(curves[i]->start) == 1)
		    	curves[i]->start->posn = b->end;
		    else
		    {
		    	node = make_node(b->end);
			change_node_of_curve(curves[i],POSITIVE_ORIENTATION,
							node);
		    }
		}
		else if (b == curves[i]->last)
		{
		    curves[i]->last = b->prev;
		    curves[i]->last->next = NULL;
		    curves[i]->num_points--;
		    if (I_NumOfNodeCurves(curves[i]->end) == 1)
		    	curves[i]->end->posn = b->start;
		    else
		    {
		    	node = make_node(b->start);
			change_node_of_curve(curves[i],NEGATIVE_ORIENTATION,
							node);
		    }
		}
		else
		{
		    c = I_SplitCurve(b->start,curves[i]);
		    if (b == c[0]->last)
		    {
			c[0]->last = b->prev;
			c[0]->last->next = NULL;
		    	node = make_node(b->start);
		    	change_node_of_curve(c[0],NEGATIVE_ORIENTATION,node);
		    	c[0]->num_points--;
			
		    }
		    else if (b == c[1]->first)
		    {
			c[1]->first = b->next;
			c[1]->first->prev = NULL;
		    	node = make_node(b->end);
		    	change_node_of_curve(c[1],POSITIVE_ORIENTATION,node);
		    	c[1]->num_points--;
		    }
		    goto start_open_bonds;
		}
	    }		
	}
	free(curves);
}	/* end open_curve_bonds */

EXPORT void open_null_sides1(
	INTERFACE	*intfc,
	double		*L,
	double		*U,
	int		dir,
	int		nb)
{
	TRI		*tri;
	SURFACE 	**s;
	DEBUG_ENTER(open_null_sides1)

	if (DEBUG)
	{
	    char	      dname[1024];
	    static const char *xyz[] = { "x", "y", "z"};
	    static const char *sname[] = { "lower", "upper"};
	    static const char *strdir[3] = { "X", "Y", "Z" };
	    static const char *strnb[2] = { "LOWER", "UPPER" };
	    static int into[3][2];
	    (void) printf("Clipping interface in at %s %s side\n",
			  sname[nb],xyz[dir]);
	    (void) printf("U = ( %g %g %g )\n",U[0],U[1],U[2]);
	    (void) printf("L = ( %g %g %g )\n",L[0],L[1],L[2]);
	    (void) sprintf(dname,"fscatter/open_null_sides1/Into%d-%s_%s",
			   into[dir][nb],strdir[dir],strnb[nb]);
	    ++into[dir][nb];
	    summarize_interface(dname,"before",intfc,XY_PLANE,
				"open_null_sides1","before");
	}

	set_floating_point_tolerance1(computational_grid(intfc)->h);
	for (s = intfc->surfaces; s && *s; ++s)
	{
	    open_surf_null_sides(*s,L,U,dir,nb);
	}
	for (s = intfc->surfaces; s && *s; ++s)
	{
	    if (no_tris_on_surface(*s))
	    {
	    	(void) delete_surface(*s);
		--s;
	    }
	}
	DEBUG_LEAVE(open_null_sides1)
}		/*end open_null_sides1*/

EXPORT void remove_out_domain_tri(
	TRI		*tri,
	SURFACE		*s)
{
	TRI		*nbtri;
	int		i;

	for (i = 0; i < 3; ++i)
	{
	    if (!is_side_bdry(tri,i))
	    {
		nbtri = Tri_on_side(tri,i);
		if (nbtri != NULL)
		{
		    if (Tri_on_side01(nbtri) == tri)
		    	Tri_on_side01(nbtri) = NULL;
		    else if (Tri_on_side12(nbtri) == tri)
		    	Tri_on_side12(nbtri) = NULL;
		    else if (Tri_on_side20(nbtri) == tri)
		    	Tri_on_side20(nbtri) = NULL;
		}
	    }
	    else
	    {
		BOND_TRI *bt = Bond_tri_on_side(tri,i);
		if (bt != NULL)
		    (void) delete_from_pointers(bt,&Btris(bt->bond));
	    }
	}
	remove_tri_from_surface(tri,s,NO);
}		/*end remove_out_domain_tri*/


LOCAL boolean tri_out_domain(
	TRI		*tri,
	double		*L,
	double		*U,
	int		dir,
	int		nb)
{
	POINT **p;
	int   i;

	p = Point_of_tri(tri);
	if (nb == 0)
	{
	    for (i = 0; i < 3; ++i)
	    {
		if ((L[dir] - Coords(p[i])[dir]) <= tol1[dir])
	    	    return NO;
	    }
	}
	else
	{
	    for (i = 0; i < 3; ++i)
	    {
		if ((Coords(p[i])[dir] - U[dir]) <= tol1[dir])
	    	    return NO;
	    }
	}
	return YES;
}	/* end tri_out_domain */

LOCAL boolean tri_set_out_domain(
	TRI		*tri,
	double		*L,
	double		*U,
	int		dir,
	int		nb)
{
	int		i,j,n,l,num_bonds;
	POINT		*p;
	BOND		*b,*bonds[2];
	TRI		**tris,*t;
	BOND_TRI	**btris;

	if (!tri_out_domain(tri,L,U,dir,nb))
	    return NO;
	for (i = 0; i < 3; i++)
	{
	    p = Point_of_tri(tri)[i];
	    if (Boundary_point(p))
	    {
		bonds_of_boundary_point(p,tri,bonds,&num_bonds);
		if (num_bonds == 0) continue;
		for (l = 0; l < num_bonds; ++l)
		{
		  b = bonds[l];
		  for (btris = Btris(b); btris && *btris; ++btris)
		  {
		    t = (*btris)->tri;
		    n = set_tri_list_around_point(p,t,&tris,t->surf->interface);
		    for (j = 0; j < n; ++j)
		    {
			if (!tri_out_domain(tris[j],L,U,dir,nb))
			    return NO;
		    }
		  }
		}
	    }
	}
	return YES;
}	/* end tri_set_out_domain */


LOCAL void clip_intfc_at_grid_bdry1(
	INTERFACE	*intfc)
{
	INTERFACE	*cur_intfc = current_interface();
	RECT_GRID	*gr = computational_grid(intfc);
	int		dim = intfc->dim;
	int		dir, nb;
	double		L[MAXD],U[MAXD];

	DEBUG_ENTER(clip_intfc_at_grid_bdry1)
	strip_subdomain_bdry_curves(intfc);
	set_current_interface(intfc);
	for (dir = 0; dir < dim; ++dir)
	{
	    L[dir] = gr->L[dir];
	    U[dir] = gr->U[dir];
	    if (gr->lbuf[dir] == 0) L[dir] -= 0.5*gr->h[dir];
	    if (gr->ubuf[dir] == 0) U[dir] += 0.5*gr->h[dir];
	    
	    if(rect_boundary_type(intfc,dir,0) == OPEN_BOUNDARY)
		L[dir] = gr->VL[dir];
	    if(rect_boundary_type(intfc,dir,1) == OPEN_BOUNDARY)
		U[dir] = gr->VU[dir];

            /* do not cut the reflect part */
            if(rect_boundary_type(intfc,dir,0) == REFLECTION_BOUNDARY)
                L[dir] = gr->VL[dir];
            if(rect_boundary_type(intfc,dir,1) == REFLECTION_BOUNDARY)
                U[dir] = gr->VU[dir];
	}
	
	for (dir = 0; dir < dim; ++dir)
	{
	    for (nb = 0; nb < 2; ++nb)
	    	open_null_sides1(intfc,L,U,dir,nb);
	}
	for (dir = 0; dir < dim; ++dir)
	{
            for (nb = 0; nb < 2; ++nb)
            	open_null_bonds(intfc,L,U,dir,nb);
	}

	cut_out_curves_in_buffer(intfc);
	reset_intfc_num_points(intfc);
	set_current_interface(cur_intfc);
	DEBUG_LEAVE(clip_intfc_at_grid_bdry1)
}		/*end clip_intfc_at_grid_bdry1*/


LOCAL	boolean find_ending_null_side(TRI*,int,TRI**,int*);

LOCAL   boolean find_ending_null_side(
        TRI *start,
        int start_side,
        TRI **tri_start,
        int *side_start)
{
	POINT	*p, *p1, **recorded_pts;
	int	num_recorded_pts;
	TRI	*next_tri, *next_tri1;
	int	i, side, side1;
	CURVE	**c;

	uni_array(&recorded_pts, MAX_NULL_SIDE, sizeof(POINT*));

	next_tri = start;
	side = start_side;
        
	num_recorded_pts = 0;
        recorded_pts[0] = Point_of_tri(start)[side];

	while(1)
	{
	    num_recorded_pts++;
	    if(num_recorded_pts > MAX_NULL_SIDE)
	    {
		printf("ERROR find_ending_null_side"
		       "too many points on the null loop.\n");
		clean_up(ERROR);
	    }

	    p = Point_of_tri(next_tri)[Next_m3(side)];
	    recorded_pts[num_recorded_pts] = p;
	    
	    /* if a loop is found, return the ending tri and point. */
	    for(i=0; i<num_recorded_pts; i++)
	    {
		if (recorded_pts[i] == p)
		    break;
	    }
	    if(i<num_recorded_pts)
		break;

	    side1 = side;
	    p1 = p;
	    next_tri1 = find_following_null_tri(next_tri,&p1,
	                  &side1,COUNTER_CLOCK);
	    
	    /* hit another curve. */
	    if (next_tri1 == NULL)
	        break;
	    
	    next_tri = next_tri1;
	    side = side1;
	}

	*tri_start = next_tri;
	*side_start = side;
	
	free(recorded_pts);
	return YES;
}


LOCAL   boolean find_null_side_loop(
        INTERFACE *intfc,
        TRI *start,
        int start_side,
        TRI **tri_start,
        int *side_start)
{
        POINT *p, **recorded_pts, *ps;
        int i, j, k, num_recorded_pts,num_sides;
        TRI **null_tris, *next_tri = start;
        int *null_sides, side = start_side;
        CURVE **c;

        uni_array(&null_tris,MAX_NULL_SIDE,sizeof(TRI*));
        uni_array(&null_sides,MAX_NULL_SIDE,INT);
        uni_array(&recorded_pts,3*MAX_NULL_SIDE,sizeof(POINT*));

        num_sides = 0;
        num_recorded_pts = 0;

        recorded_pts[0] = Point_of_tri(start)[side];

        do
        {
            null_tris[num_sides] = next_tri;
            null_sides[num_sides] = side;
            ++num_sides;
            num_recorded_pts++;
            recorded_pts[num_recorded_pts] = p = Point_of_tri(next_tri)[Next_m3(side)];
            for (i = 0; i < num_recorded_pts; i++)
            {
                if (recorded_pts[i] == recorded_pts[num_recorded_pts])
                {
                    *tri_start = null_tris[i];
                    *side_start = null_sides[i];
                    free_these(3,null_tris,null_sides,recorded_pts);
                    return YES;
                }
	        ps = Point_of_tri(null_tris[i])[null_sides[i]];
		for (c = intfc->curves; c && *c; c++)
		{
		    if(is_closed_curve(*c))
		        continue;
		    if (ps == (*c)->first->start)
		    {
		        *tri_start = null_tris[i];
		        *side_start = null_sides[i];
		        free_these(3,null_tris,null_sides,recorded_pts);
		        return YES;
		    }
		}
	    }		
	    if (p == recorded_pts[0]) break;
	    next_tri = find_following_null_tri(next_tri,&p,
	                  &side,COUNTER_CLOCK);
	    if (next_tri == NULL)
	        break;
	}
	while (next_tri != start || side != start_side);

	*tri_start = start;
	*side_start = start_side;
	free_these(3,null_tris,null_sides,recorded_pts);
	return YES;
}       /*end find_null_side_loop*/

LOCAL void merge_point_pointers_at_subdomain_bdry(
	TRI             **tris_a,
	TRI		**tris_s,
	int		nt,
	P_LINK		*p_table,
	int		p_size)
{
	POINT		*p, *ap;
	int		i, j;

        for (i = 0; i < nt; ++i)
	{
	    sorted(Point_of_tri(tris_s[i])[0]) = NO; 
	    sorted(Point_of_tri(tris_s[i])[1]) = NO;
	    sorted(Point_of_tri(tris_s[i])[2]) = NO;
	    sorted(Point_of_tri(tris_a[i])[0]) = NO; 
	    sorted(Point_of_tri(tris_a[i])[1]) = NO;
	    sorted(Point_of_tri(tris_a[i])[2]) = NO;
	}
	for (i = 0; i < nt; ++i)
	{
	    for (j = 0; j < 3; ++j)
	    {
		ap = Point_of_tri(tris_a[i])[j];
		if (sorted(ap) == NO)
		{
	            p = (POINT*)find_from_hash_table((POINTER)ap,
						     p_table,p_size);
		    if(p != Point_of_tri(tris_a[i])[j])
		        printf("#merge_point  \n");
		    
		    Point_of_tri(tris_a[i])[j] = p;
		    sorted(p) = YES;
		}
	    }
	}
}	/*end merge_point_pointers_at_subdomain_bdry*/


LOCAL	void	copy_tri_state_to_btri(
	BOND_TRI	*btri,
	BOND		*b,
	ORIENTATION	orient,
	INTERFACE	*intfc)
{
	POINT		*p;
	Locstate	sl, sr;
	NODE		*n;
	CURVE		**c;
	BOND_TRI	**bt;
	boolean		found;
	size_t	 	sizest = size_of_state(intfc);

	if (orient == POSITIVE_ORIENTATION)
	{
	    p = b->start;
	    sl = left_start_btri_state(btri);
	    sr = right_start_btri_state(btri);
	}
	else
	{
	    p = b->end;
	    sl = left_end_btri_state(btri);
	    sr = right_end_btri_state(btri);
	}

	/*in install_subdomain_bdry_curve the new inserted point is always 
	  a node point.
	*/
	if ((n = node_of_point(p,intfc)) == NULL)
	{
	    ft_assign(sl,left_state(p),sizest);
	    ft_assign(sr,right_state(p),sizest);
	    return;
	}

	found = NO;
	/* deal with 3 comp curve */
	for(c = n->out_curves; c && *c; c++)
	{
	    if(*c == btri->curve)
	        continue;
	    found = YES;
	    for(bt=Btris((*c)->first); bt && *bt; bt++)
	        if((*bt)->surface == btri->surface)
		{
		    ft_assign(sl, left_start_btri_state(*bt), sizest);
		    ft_assign(sr, right_start_btri_state(*bt), sizest);
		    return;
		}
	}
	
	for(c = n->in_curves; c && *c; c++)
	{
	    if(*c == btri->curve)
	        continue;
	    found = YES;
	    for(bt=Btris((*c)->last); bt && *bt; bt++)
	        if((*bt)->surface == btri->surface)
		{
		    ft_assign(sl, left_end_btri_state(*bt), sizest);
		    ft_assign(sr, right_end_btri_state(*bt), sizest);
		    return;
		}
	}

	/* Now curve is the only curve in node n */
	if(!found)
	{
	    ft_assign(sl,left_state(p),sizest);
	    ft_assign(sr,right_state(p),sizest);
	    return;
	}

	printf("ERROR copy_tri_state_to_btri, no suitable btri found.\n");
	clean_up(ERROR);

}		/*end copy_tri_state_to_btri*/

LOCAL TRI *find_following_null_tri(
	TRI		*tri,
	POINT		**np,
	int		*side,
	ANGLE_DIRECTION dir)
{
	POINT *p;
	TRI   *newtri;
	int   nside;

	p = (dir == COUNTER_CLOCK) ? Point_of_tri(tri)[Next_m3(*side)] :
			             Point_of_tri(tri)[*side];

	newtri = tri;
	while (newtri != NULL)
	{
	    nside = Following_side_at_vertex(newtri,p,dir);
	    if (nside == ERROR)
	    {
	    	*side = -1;
	    	return NULL;
	    }
	    if (is_side_bdry(newtri,nside))
	    {
	    	if (Bond_tri_on_side(newtri,nside) == NULL)
	    	{
	    	    *side = nside;
	    	    *np = (dir == COUNTER_CLOCK) ?
	    			Point_of_tri(newtri)[Next_m3(*side)] :
	    			Point_of_tri(newtri)[*side];
	    	    return newtri;
	    	}
	    	return NULL;
	    }
	    else if (Tri_on_side(newtri,nside) == NULL)
	    {
	    	*side = nside;
	    	*np = (dir == COUNTER_CLOCK) ?
	    		Point_of_tri(newtri)[Next_m3(*side)] :
	    		Point_of_tri(newtri)[*side];
	    	return newtri;
	    }
	    newtri = Tri_on_side(newtri,nside);
	}
	return NULL;
}		/*end find_following_null_tri*/


EXPORT void strip_subdomain_bdry_curves(
	INTERFACE	*intfc)
{
	DEBUG_ENTER(strip_subdomain_bdry_curves)
	
	strip_bdry_curves(intfc, SUBDOMAIN_HSBDRY);

	DEBUG_LEAVE(strip_subdomain_bdry_curves)
}		/*end strip_subdomain_bdry_curves*/

EXPORT  void strip_bdry_curves(
	INTERFACE	*intfc,
	int		ctype)
{
	CURVE		**c;
	CURVE		**curves_to_delete = NULL;
	NODE		**nodes_to_delete = NULL;
	SURFACE		**s;
	NODE		**n;
	POINT		*p;

	for (c = intfc->curves; c && *c; ++c)
	{
	    if (hsbdry_type(*c) != ctype)
		continue;
	    for (s = (*c)->pos_surfaces; s && *s; ++s)
	    	strip_curve_from_surf(*c,*s,POSITIVE_ORIENTATION);
	    for (s = (*c)->neg_surfaces; s && *s; ++s)
	    	strip_curve_from_surf(*c,*s,NEGATIVE_ORIENTATION);
	    if (!add_to_pointers(*c,&curves_to_delete))
	    {
	    	screen("ERROR in strip_subdomain_bdry_curves(), "
	    	       "add_to_pointers() failed\n");
	    	clean_up(ERROR);
	    }
	}
	for (c = curves_to_delete; c && *c; ++c)
	    (void) delete_curve(*c);
	for (n = intfc->nodes; n && *n; ++n)
	{
	    if (!add_to_pointers(*n,&nodes_to_delete))
	    {
	    	screen("ERROR in strip_subdomain_bdry_curves(), "
	    	       "add_to_pointers() failed\n");
	    	clean_up(ERROR);
	    }
	}
	for (n = nodes_to_delete; n && *n; ++n)
	    (void) delete_node(*n);

}	/* end strip_bdry_curves */


/*it will make the Boundary_point for node to 0, it is not a problem
  (1) if the node contains only this curve, it is correct
  (2) if the node has several different curves. This case can happen when
      node has one 3-comp curve and 3 subdomain curves and on a subdomain 
      boundary. Then, in reconstruction part, the point is outside the 
      reconstruction domain in scatter_front part, the point will be 
      cutted because it is in a subdomain boundary.
*/
LOCAL void strip_curve_from_surf(
	CURVE		*curve,
	SURFACE		*surf,
	ORIENTATION	orient)
{
	BOND		*b;
	BOND_TRI	**btris;
	POINT		*p;
	TRI		*tri;
	int		i;
	size_t	 	sizest = size_of_state(surf->interface);
	
	DEBUG_ENTER(strip_curve_from_surf)
	for (b = curve->first; b != NULL; b = b->next)
	{
	    Boundary_point(b->start) = Boundary_point(b->end) = 0;
	    for (btris = Btris(b); btris && *btris; ++btris)
	    {
		tri = (*btris)->tri;
		for (i = 0; i < 3; ++i)
		{
		    if (Bond_tri_on_side(tri,i) == *btris)
		    {
			if(sizest > 0)
			{
			    p = (*btris)->bond->start;
			    ft_assign(left_state(p),
				left_start_btri_state(*btris),sizest);
			    ft_assign(right_state(p),
				right_start_btri_state(*btris),sizest);
			    
			    p = (*btris)->bond->end;
			    ft_assign(left_state(p),
				left_end_btri_state(*btris),sizest);
			    ft_assign(right_state(p),
				right_end_btri_state(*btris),sizest);
			}

		    	Bond_tri_on_side(tri,i) = NULL;
		    	set_side_bdry(Boundary_tri(tri),i,NO);
		    	break;
		    }
		}
		if (!delete_from_pointers(*btris,&Btris(b)))
		{
		    screen("ERROR in strip_curve_from_surf(), "
		           "delete_from_pointers() failed\n");
		    clean_up(ERROR);
		}
	    }
	}
	if (orient == POSITIVE_ORIENTATION)
	{
	    if (!(delete_from_pointers(curve,&surf->pos_curves) &&
		     delete_from_pointers(surf,&curve->pos_surfaces)))
	    {
		screen("ERROR in strip_curve_from_surf(), "
		       "delete_from_pointers() failed\n");
		clean_up(ERROR);
	    }
	}
	else
	{
	    if (!(delete_from_pointers(curve,&surf->neg_curves) &&
	    	     delete_from_pointers(surf,&curve->neg_surfaces)))
	    {
	    	screen("ERROR in strip_curve_from_surf(), "
		       "delete_from_pointers() failed\n");
		clean_up(ERROR);
	    }
	}
	DEBUG_LEAVE(strip_curve_from_surf)
}		/*end strip_curve_from_surf*/


EXPORT	INTERFACE  *cut_buf_interface1(
	INTERFACE	*intfc,
	int		dir,
	int		nb,
	int		*me,
	int		*him)
{
	INTERFACE	  *sav_intfc, *tmp_intfc, *buf_intfc;
	RECT_GRID	  *gr = computational_grid(intfc);
	RECT_GRID	  dual_gr;
	double		  L[3], U[3];
	boolean		  sav_copy = copy_intfc_states();
	char              dname[1024]; 
	static const char *strdir[3] = { "X", "Y", "Z" };
	static const char *strnb[2] = { "LOWER", "UPPER" };

	DEBUG_ENTER(cut_buf_interface1)

	if (DEBUG)
	{
	    static int into[3][2];
	    (void) sprintf(dname,"fscatter/cut_buf_interface1/Into%d-%s_%s",
			   into[dir][nb],strdir[dir],strnb[nb]);
	    (void) printf("cut_buf_interface1() at ENTRY\n"
			  "  dir = %s  nb = %s  me = ( %d, %d, %d )  "
			  "him = ( %d, %d, %d )\n",
			  strdir[dir],strnb[nb],me[0],me[1],me[2],
			  him[0],him[1],him[2]);
	    summarize_interface(dname,"intfc",intfc,XY_PLANE,
				"cut_buf_interface1","intfc at ENTER");
	    ++into[dir][nb];
	}

	set_copy_intfc_states(YES);
	sav_intfc = current_interface();

	set_size_of_intfc_state(size_of_state(intfc));
	tmp_intfc = copy_interface(intfc);
	set_dual_grid(&dual_gr,gr);

	if (nb == 0)
	{
	    L[dir] = gr->L[dir];
	    U[dir] = gr->L[dir]+(gr->L[dir]-dual_gr.VL[dir])+0.5*gr->h[dir];
	}
	else
	{
	    L[dir] = gr->U[dir]-(dual_gr.VU[dir]-gr->U[dir])-0.5*gr->h[dir];
	    U[dir] = gr->U[dir];
	}

	open_null_sides1(tmp_intfc,L,U,dir,(nb+1)%2);
	open_null_bonds(tmp_intfc,L,U,dir,(nb+1)%2);

	/*
	 * Split curves that have been disconnected from surfaces
	 * and delete those sections that lie entirely within the
	 * subdomain boundary.
	 */

	cut_out_curves_in_buffer(tmp_intfc);

	reset_intfc_num_points(tmp_intfc);

	if (me[dir] == him[dir])
	{
	    buf_intfc = tmp_intfc;
	}
	else
	{
	    set_size_of_intfc_state(size_of_state(intfc));
	    buf_intfc = copy_interface(tmp_intfc);
	    (void) delete_interface(tmp_intfc);
	}
	set_copy_intfc_states(sav_copy);
	set_current_interface(sav_intfc);

	if (DEBUG)
	{
	    static int outof[3][2];
	    (void) sprintf(dname,"fscatter/cut_buf_interface1/Outof%d-%s_%s",
			   outof[dir][nb],strdir[dir],strnb[nb]);
	    summarize_interface(dname,"intfc",intfc,XY_PLANE,
				"cut_buf_interface1","intfc at EXIT");
	    summarize_interface(dname,"buf_intfc",buf_intfc,XY_PLANE,
				"cut_buf_interface1","buf_intfc at EXIT");
	    ++outof[dir][nb];
	}
	DEBUG_LEAVE(cut_buf_interface1)
	return buf_intfc;
}		/*end cut_buf_interface1*/

EXPORT	void cut_out_curves_in_buffer(
	INTERFACE	*intfc)
{
	BOND  *bs, *be;
	CURVE **c;
	NODE  **n;
	NODE  **nodes_to_delete;

	for (c = intfc->curves; c && *c; ++c)
	{
	    for (bs = (*c)->first; bs != NULL; bs = bs->next)
	        if (Btris(bs) == NULL && hsbdry_type(*c) <
                        FIRST_PHYSICS_HSBDRY_TYPE)
		    break;
	    if (bs == NULL)
		continue;
	    for (be = bs; be->next != NULL; be = be->next)
		if (Btris(be->next) != NULL && hsbdry_type(*c) <
                        FIRST_PHYSICS_HSBDRY_TYPE)
		    break;
	    if ((bs == (*c)->first) && (be == (*c)->last))
	    {
	        /*Entire curve is in buffer region*/
		if (!delete_curve(*c))
		{
		    screen("ERROR in cut_out_curves_in_buffer(), "
		           "can't delete curve %llu\n",(long long unsigned int)curve_number(*c));
		    print_curve(*c);
		    print_interface(intfc);
		    clean_up(ERROR);
		}
		if (intfc->curves == NULL) /*No curves left*/
		    break;
		c = intfc->curves - 1;/*Start over on curve list*/
	    }
	    else if (be != (*c)->last)
	    {
	        /*without this, the closed curve will be cutted into 
		  two connected pieces it will make an inconsistent 
		  interface after copy_interface.
		*/
	        if(is_closed_curve(*c))
	            change_node_of_closed_curve(bs->start, *c);

		if (split_curve(be->end,be,*c,NO_COMP,NO_COMP,NO_COMP,NO_COMP)
		    == NULL)
		{
		    screen("ERROR in cut_out_curves_in_buffer(), can't split "
			   "start section of curve %llu\n",curve_number(*c));
		    print_curve(*c);
		    print_interface(intfc);
		    clean_up(ERROR);
		}
		c = intfc->curves - 1;/*Start over on curve list*/
	    }
	    else if (bs != (*c)->first)
	    {
	        if(is_closed_curve(*c))
	            change_node_of_closed_curve(be->end, *c);

		if (split_curve(bs->start,bs,*c,NO_COMP,NO_COMP,NO_COMP,NO_COMP)
		    == NULL)
		{
		    screen("ERROR in cut_out_curves_in_buffer(), can't split "
			   "end section of curve %llu\n",curve_number(*c));
		    print_curve(*c);
		    print_interface(intfc);
		    clean_up(ERROR);
		}
		c = intfc->curves - 1;/*Start over on curve list*/
	    }
	}
	nodes_to_delete = NULL;
	for (n = intfc->nodes; n && *n; ++n)
	{
	    if (!add_to_pointers(*n,&nodes_to_delete))
	    {
	    	screen("ERROR in cut_out_curves_in_buffer(), "
	    	       "add_to_pointers() failed\n");
	    	clean_up(ERROR);
	    }
	}
	for (n = nodes_to_delete; n && *n; ++n)
	    (void) delete_node(*n);
}		/*end cut_out_curves_in_buffer*/

EXPORT	void communicate_default_comp(
	Front *fr)
{
	INTERFACE    *intfc = fr->interf;
	PP_GRID      *pp_grid = fr->pp_grid;
	int          me[MAXD], him[MAXD];
	int          i,j,k,dim,myid, dst_id;
	int          *G;
	boolean	     status = NO;
	COMPONENT    comp_buf;
	RECT_GRID    *gr = computational_grid(intfc);
	int	     count,max_count;
	COMPONENT    save_default_comp = intfc->default_comp;
	int	     max_G;

	if (pp_numnodes() == 1) return;
	dim = intfc->dim;
	myid = pp_mynode();
        G = pp_grid->gmax;
        find_Cartesian_coordinates(myid,pp_grid,me);

	if (debugging("default_comp"))
	{
	    (void) printf("Entering communicate_default_comp()\n");
	    (void) printf("myid = %d\n",myid);
	    (void) printf("G = %d %d %d\n",G[0],G[1],G[2]);
	    (void) printf("me = %d %d %d\n",me[0],me[1],me[2]);
	    (void) printf("default_comp = %d\n",intfc->default_comp);
	}

	intfc->default_comp = NO_COMP;
	count = 0;
	max_G = G[0];
	for (i = 1; i < dim; ++i)
	    if (max_G < G[i]) max_G = G[i];
	max_count = 1;
	while (max_G != 1)
	{
	    max_count++;
	    max_G /= 2;
	}

	while (pp_min_status(status) == NO)
	{
	    if (debugging("default_comp"))
	    	(void) printf("Round %d\n",count);
	    
	    for (i = 0; i < dim; ++i)
	    {
	    	for (j = 0; j < 2; ++j)
		{
		    if(debugging("default_comp"))
		        printf("\n#comm buf %d %d\n", i, j);
		    pp_gsync();

		    /*send to adj proc */
		    for (k = 0; k < dim; ++k)
		    	him[k] = me[k];
		    if (rect_boundary_type(intfc,i,j) == SUBDOMAIN_BOUNDARY)
		    {
		    	him[i] = (me[i] + 2*j - 1 + G[i])%G[i];
			dst_id = domain_id(him,G,dim);
			comp_buf = buffer_component(intfc,i,j);
			if (comp_buf != NO_COMP)
	    		    intfc->default_comp = comp_buf;
			pp_send(0,&comp_buf,INT,dst_id);
			
			if (debugging("default_comp"))
			{
			    (void) printf("comp_buf = %d\n",comp_buf);
			    (void) printf("Send dst_id = %d\n",dst_id);
			    (void) printf("him = %d %d %d\n",him[0],
			    		him[1],him[2]);
			}
		    }
			
		    /*recv from adj proc */
		    for (k = 0; k < dim; ++k)
		    	him[k] = me[k];
		    if (rect_boundary_type(intfc,i,(j+1)%2) == 
		    			SUBDOMAIN_BOUNDARY)
		    {
		    	him[i] = (me[i] - 2*j + 1 + G[i])%G[i];
			dst_id = domain_id(him,G,dim);
			pp_recv(0,dst_id,&comp_buf,INT);
			if (comp_buf != NO_COMP)
			    intfc->default_comp = comp_buf;
			if (debugging("default_comp"))
			{
			    (void) printf("received comp_buf = %d\n",comp_buf);
			    (void) printf("Receive dst_id = %d\n",dst_id);
			    (void) printf("him = %d %d %d\n",him[0],
			    			him[1],him[2]);
			}
		    }
		}  /*for j: side */
	    }      /*for i: direction */

	    status = (intfc->default_comp == NO_COMP) ? NO : YES;
	    if (count++ > max_count)
	    {
		intfc->default_comp = save_default_comp;
	    }
	}
	if (debugging("default_comp"))
	    (void) printf("Leaving communicate_default_comp()\n");
}	/* end communicate_default_comp */

LOCAL COMPONENT buffer_component(
	INTERFACE *intfc,
	int dir,
	int nb)
{
	RECT_GRID *gr = computational_grid(intfc);
	double *L = gr->L;
	double *U = gr->U;
	double *h = gr->h;
	double coords[MAXD];
	int i, comp;

	if (intfc->surfaces == NULL) 
	    return intfc->default_comp;
	for (i = 0; i < 3; ++i)
	    coords[i] = 0.5*(L[i] + U[i]);
	coords[dir] = (nb == 0) ? (L[dir] - 0.5*h[dir]) : 
			(U[dir] + 0.5*h[dir]);

	comp = component(coords,intfc);
	if(comp == NO_COMP)
	    return intfc->default_comp;
	else
	    return comp;

}	/* end buffer_component */

EXPORT	void cut_surface(
	SURFACE *surf,
	boolean (*constr_func)(POINTER,double*), /* Constraint function */
        POINTER func_params,		     /* Constraint function params */
	boolean force_clip)
{
	INTERFACE *intfc = surf->interface;
	TRI *tri,*ntri;
	POINT *p;
	double *p1,*p2,*pc,pt[MAXD];
	int i,j,k,l,nside,dim = intfc->dim;
	double distance;
	double tol = 1e-5*MIN_SC_SEP(intfc);
	boolean p1_constr,p2_constr,pt_constr;
	TRI **tris;
	int num_tris;
	int count = 0;

	if (force_clip)
	{
	    set_current_interface(surf->interface);

top_of_loop:
	    count++;
	    for (tri=first_tri(surf); !at_end_of_tri_list(tri,surf); 
				tri=tri->next)
	    	for (i = 0; i < 3; ++i)
		    sorted(Point_of_tri(tri)[i]) = NO;

	    for (tri=first_tri(surf); !at_end_of_tri_list(tri,surf); 
				tri=tri->next)
	    {
	    	for (i = 0; i < 3; ++i)
		{
		    if (sorted(Point_of_tri(tri)[i]) &&
			sorted(Point_of_tri(tri)[(i+1)%3]))
			continue;
		    p1 = Coords(Point_of_tri(tri)[i]);
		    p2 = Coords(Point_of_tri(tri)[(i+1)%3]);
		    p1_constr = constr_func(func_params,p1);
		    p2_constr = constr_func(func_params,p2);
		    if ((p1_constr && !p2_constr) ||
		        (p2_constr && !p1_constr))
		    {
			pc = constr_position(p1,p2,constr_func,func_params);
                        distance = distance_between_positions(pc,p1,dim);
                        if (distance < tol)
                        {
                            for (j = 0; j < dim; ++j)
                                p1[j] = pc[j];
                            sorted(Point_of_tri(tri)[i]) = YES;
                            continue;
                        }
                        distance = distance_between_positions(pc,p2,dim);
                        if (distance < tol)
                        {
                            for (j = 0; j < dim; ++j)
                                p2[j] = pc[j];
                            sorted(Point_of_tri(tri)[(i+1)%3]) = YES;
                            continue;
                        }
                        p = Point(pc);
                        insert_point_in_tri_side(p,i,tri,surf);
                        sorted(p) = YES;  
		    }
		    else
		    {
			for (k = 1; k <= 10; ++k)
			{
			    for (j = 0; j < dim; ++j)
				pt[j] = p1[j] + k*(p2[j] - p1[j])/10.0;
		    	    pt_constr = constr_func(func_params,pt);
			    if (pt_constr == p1_constr) continue;

			    pc = constr_position(p1,pt,constr_func,func_params);
                            distance = distance_between_positions(pc,p1,dim);
                            if (distance < tol)
                            {
                            	for (j = 0; j < dim; ++j)
                                    p1[j] = pc[j];
				p = Point_of_tri(tri)[i];
                        	sorted(p) = YES;  
                            }
			    else
			    {
                        	p = Point(pc);
                        	insert_point_in_tri_side(p,i,tri,surf);
                        	sorted(p) = YES;  
			    }
			    ntri = NULL;
			    num_tris = set_tri_list_around_point(p,tri,
							&tris,intfc);
			    for (l = 0; l < num_tris; ++l)
			    {
				for (j = 0; j < 3; ++j)
				{
				    if ((Point_of_tri(tris[l])[j] == p &&
					Coords(Point_of_tri(tris[l])[(j+1)%3]) 
						== p2) ||
					(Point_of_tri(tris[l])[(j+1)%3] == p &&
					Coords(Point_of_tri(tris[l])[j]) == p2))
				    {
					ntri = tris[l];
					nside = j;
					break;
				    }
				}
			    }
			    if (ntri == NULL) 
			    {
				(void) printf("In cut_surface(): ");
				(void) printf("cannot find tri of p2\n");
				clean_up(ERROR);
			    }
			    pc = constr_position(p2,pt,constr_func,func_params);
                            distance = distance_between_positions(pc,p2,dim);
                            if (distance < tol)
                            {
                            	for (j = 0; j < dim; ++j)
                                    p2[j] = pc[j];
				sorted(p) = YES;
                            }
			    else
			    {
                        	p = Point(pc);
                        	insert_point_in_tri_side(p,nside,ntri,surf);
                        	sorted(p) = YES;  
			    }
			    break;
			}
		    }
		}
	    }
	    if (count < 4) goto top_of_loop;
	}
	for (tri=first_tri(surf); !at_end_of_tri_list(tri,surf); tri=ntri)
	{
	    ntri = tri->next;
	    for (i = 0; i < 3; ++i)
	    {
		pc[i] = (Coords(Point_of_tri(tri)[0])[i] +
                         Coords(Point_of_tri(tri)[1])[i] +
                         Coords(Point_of_tri(tri)[2])[i])/3.0;
	    }
	    if (!constr_func(func_params,pc))
	    {
		remove_out_domain_tri(tri,surf);
	    }
	}
}	/* end cut_surface */

LOCAL 	double *constr_position(
	double *p1,
	double *p2,
	boolean (*constr_func)(POINTER,double*), /* Constraint function */
        POINTER func_params)		     /* Constraint function params */
{
	static double p_in[MAXD];
	double p_out[MAXD],pc[MAXD];
	int i,j,N = 20;		/* accurate to 2^-20*(p1-p2) */
	
	if (constr_func(func_params,p1) && !constr_func(func_params,p2))
	{
	    for (i = 0; i < MAXD; ++i)
	    {
		p_in[i]  = p1[i];
		p_out[i] = p2[i];
	    }
	}
	else if (constr_func(func_params,p2) && !constr_func(func_params,p1))
	{
	    for (i = 0; i < MAXD; ++i)
	    {
		p_in[i]  = p2[i];
		p_out[i] = p1[i];
	    }
	}
	else
	{
	    screen("ERROR: In constr_position(), same side!\n");
	    clean_up(ERROR);
	}
	for (i = 0; i < N; ++i)
	{
	    for (j = 0; j < MAXD; ++j)
		pc[j] = 0.5*(p_in[j] + p_out[j]);
	    if (constr_func(func_params,pc))
		for (j = 0; j < MAXD; ++j)
		    p_in[j] = pc[j];
	    else
		for (j = 0; j < MAXD; ++j)
		    p_out[j] = pc[j];
	}
	return p_in;
}	/* end constr_position */
	
EXPORT void install_hsbdry_on_surface(
	SURFACE	*surf,
	int hsbdry_type)
{
	BOND	 *b;
	BOND_TRI *btri;
	CURVE	 *curve, **c;
	NODE	 *ns, *ne, **n, **n1;
	POINT	 *p, *ps, *pe;
	TRI	 *tri_start, *start_tri;
	Locstate sl, sr;
	boolean  	 sav_intrp, dup_nodes;
	int	 side_start, start_side, ntris, i;
	INTERFACE *intfc = surf->interface;
	int	 dim = intfc->dim;
	TRI	 **null_tris = NULL;
	int	 num_tri_allocated = 0;

	DEBUG_ENTER(install_hsbdry_on_surface)
	
	sav_intrp = interpolate_intfc_states(intfc);
	interpolate_intfc_states(intfc) = NO;
	
	incremental_alloc_tris(&null_tris,&num_tri_allocated);
	while (null_side_on_surface(surf,&start_tri,&start_side))
	{
	    find_ending_null_side(start_tri,start_side,&tri_start,&side_start);
		
	    ps = Point_of_tri(tri_start)[side_start];
	    ns = make_node(ps);
	    pe = Point_of_tri(tri_start)[Next_m3(side_start)];
	    ne = make_node(pe);
	    curve = make_curve(NO_COMP,NO_COMP,ns,ne);

	    hsbdry_type(curve) = hsbdry_type;
	    install_curve_in_surface_bdry(surf,curve,POSITIVE_ORIENTATION);
	    	
	    null_tris[0] = tri_start;
	    ntris = 1;

	    b = curve->first;
	    while((tri_start = find_following_null_tri(tri_start,&p,
					&side_start,CLOCKWISE)) != NULL)
	    {
		if (insert_point_in_bond(ns->posn,b,curve) !=
			FUNCTION_SUCCEEDED)
		{
		    screen("ERROR in install_hsbdry_on_surface(), "
		               "insert_point_in_bond() failed\n");
		    clean_up(ERROR);
		}
		    
		ns->posn = b->start = p;
		set_bond_length(b,dim);
	
		null_tris[ntris] = tri_start;
		ntris++;
		if(ntris >= num_tri_allocated)
		{
		    incremental_alloc_tris(&null_tris,&num_tri_allocated);
		}
	    
		if (ns->posn == ne->posn) /* Closed loop formed */
		{
		    change_node_of_curve(curve,POSITIVE_ORIENTATION,ne);
		    if(!delete_node(ns))
		    {
			printf("ERROR install_hsbdry_on_surface "
			    	   "can not delete ns node for a loop.\n");
			print_node(ns);
			clean_up(ERROR);
		    }
		    break;
		}
	    }
		
	    for(i=0, b=curve->first; b!=NULL; b=b->next, i++)
	    {
		btri = link_tri_to_bond(NULL,null_tris[ntris-1-i],surf,b,curve);
		copy_tri_state_to_btri(btri,b,NEGATIVE_ORIENTATION,intfc);
		copy_tri_state_to_btri(btri,b,POSITIVE_ORIENTATION,intfc);
	    }
	}   /* null_side_on_surface. */

	dup_nodes = YES;
	while(dup_nodes)
	{
	    dup_nodes = NO;
	    for(n=intfc->nodes; n && *n; n++)
	    {
		for(n1=n+1; n1 && *n1; n1++)
		{
		    if((*n)->posn != (*n1)->posn)
			continue;
		    dup_nodes = YES;

		    for(c=(*n1)->in_curves; c && *c; c++)
			change_node_of_curve(*c,NEGATIVE_ORIENTATION,*n);
		    for(c=(*n1)->out_curves; c && *c; c++)
			change_node_of_curve(*c,POSITIVE_ORIENTATION,*n);
		    
		    if(!delete_node(*n1))
		    {
			printf("ERROR install_hsbdry_on_surface "
		    	       "can not delete node.\n");
			print_node(*n1);
			clean_up(ERROR);
		    }
		    break;
		}
		if(dup_nodes)
		    break;
	    }
	}
	free_these(1,null_tris);

	if (debugging("consistency"))
	{
	    if (!consistent_interface(intfc))
	    {
		screen("ERROR in install_hsbdry_on_surface(), "
		       "intfc is inconsistent\n");
		clean_up(ERROR);
	    }
	}
	interpolate_intfc_states(intfc) = sav_intrp;
	DEBUG_LEAVE(install_hsbdry_on_surface)
}		/*end install_hsbdry_on_surface*/

EXPORT 	boolean surfaces_matched(
	SURFACE *as,
	SURFACE *s)
{
	if (Gindex(as) != -1 && Gindex(s) != -1)
	    return (Gindex(as) == Gindex(s)) ? YES : NO;
	else if ((Gindex(as) == -1 && Gindex(s) != -1) ||
		 (Gindex(as) != -1 && Gindex(s) == -1))
	{
	    (void) printf("In surfaces_matched(): "
			  "Gindex incompletely assigned\n");
	    clean_up(ERROR);	
	}
	/* The following is for cases without Gindex */
	if (wave_type(s) != wave_type(as))
	    return NO;
	else if ((is_bdry(as) && !is_bdry(s)) || (is_bdry(s) && !is_bdry(as)))
	    return NO;
	else if (is_bdry(as) && is_bdry(s)) 
	    return (Hyper_surf_index(as) == Hyper_surf_index(s)) ? YES : NO;
	else if (negative_component(as) == negative_component(s) &&
	    	 positive_component(as) == positive_component(s))
	{
	    if (wave_type(s) == ICE_PARTICLE_BOUNDARY ||
		wave_type(as) == ICE_PARTICLE_BOUNDARY)
	    {
		return (body_index(s) == body_index(as)) ? YES : NO;
	    }
	    else if (wave_type(s) == ELASTIC_BOUNDARY && 
		     wave_type(as) == ELASTIC_BOUNDARY)
		return (body_index(s) == body_index(as)) ? YES : NO;
	    else
	    	return YES;
	}
	else
	    return NO;
}	/* end surfaces_matched */

int compare_points(const void *, const void *);

int compare_points(const void *a, const void *b)
{
        POINT   **a1=(POINT**)a, **b1=(POINT**)b;
        POINT   *p1 = *a1, *p2 = *b1;
        double   *crds1, *crds2;
        int     i;

        crds1 = Coords(p1);
        crds2 = Coords(p2);
        for(i=0; i<3; i++)
        {
            if(fabs(crds1[i]-crds2[i]) < tol1[i])
                continue;
            else if(crds1[i] < crds2[i])
                return -1;
            else
                return 1;
        }
        return 0;
}

LOCAL	void incremental_alloc_tris(
	TRI ***tris,
	int *num_allocated_tris)
{
	TRI **new_tris;
	int i,old_num_allocated_tris = *num_allocated_tris;

	*num_allocated_tris += MAX_SUBDOMAIN_TRIS;
	uni_array(&new_tris,*num_allocated_tris,sizeof(TRI*));
	for (i = 0; i < old_num_allocated_tris; ++i)
	    new_tris[i] = (*tris)[i];
	if (old_num_allocated_tris != 0 || *tris != NULL)
	    free_these(1,*tris);
	*tris = new_tris;
}	/* end incremental_alloc_tri */

LOCAL	void merge_point_pair_gindex(
	POINT *ps,		/* local point */
	POINT *pa,		/* adjacent point */
	RECT_GRID *grid,	/* computational grid */
        int dir,
        int nb)
{
	double *L = grid->L;
	double *U = grid->U;
	if (Gindex(ps) == -1 && Gindex(pa) == -1)
	    return;		/* global index not assigned */
	if (nb == 0)
	{
	    if (Coords(pa)[dir] < L[dir])
		Gindex(ps) = Gindex(pa);
	    else
		Gindex(pa) = Gindex(ps);
	}
	else
	{
	    if (Coords(pa)[dir] >= U[dir])
		Gindex(ps) = Gindex(pa);
	    else
		Gindex(pa) = Gindex(ps);
	}
}	/* end merge_point_pair_gindex */

LOCAL	void merge_gindex_of_tris(
	TRI **tris_a,
	int na,
	RECT_GRID *grid,
	int dir,
	int nb,
	P_LINK		*p_table,
	int		p_size)
{
	int i,j;
	POINT *ps,*pa;
	TRI *t;
	double *L = grid->L;
	double *U = grid->U;

	for (i = 0; i < na; ++i)
	{
	    t = tris_a[i];
	    for (j = 0; j < 3; ++j)
	    {
		pa = Point_of_tri(t)[j];
		ps = (POINT*)find_from_hash_table((POINTER)pa,p_table,p_size);
		merge_point_pair_gindex(ps,pa,grid,dir,nb);
	    }
	}
}	/* end merge_gindex_of_tris */

LOCAL boolean bond_out_domain1(
	BOND		*bond,
	double		*L,
	double		*U,
	int		dir,
	int		nb)
{
	POINT *ps,*pe;

	if (Btris(bond) != NULL) return NO;	/* not a null bond */

	ps = bond->start;	pe = bond->end;
	if (nb == 0)
	{
	    if ((L[dir] - Coords(ps)[dir]) <= tol1[dir])
	        return NO;
	    if ((L[dir] - Coords(pe)[dir]) <= tol1[dir])
	        return NO;
	}
	else
	{
	    if ((Coords(ps)[dir] - U[dir]) <= tol1[dir])
	        return NO;
	    if ((Coords(pe)[dir] - U[dir]) <= tol1[dir])
	        return NO;
	}
	return YES;
}	/* end bond_out_domain1 */

#define		MAX_NUM_OTHER_CURVES	500

LOCAL boolean append_other_curves1(
	INTERFACE *intfc,
	INTERFACE *adj_intfc,
	RECT_GRID *grid,
	int dir,
	int nb,
	P_LINK *p_table,
	int p_size)
{
	CURVE **cc;
	CURVE *c[MAX_NUM_OTHER_CURVES];
	CURVE *ac[MAX_NUM_OTHER_CURVES];
	CURVE *mc[MAX_NUM_OTHER_CURVES];
	int i,j,num_c,num_ac;
	BOND *b,*ba;
	boolean bond_matched;
	POINT *p;

	num_c = num_ac = 0;
	for (cc = intfc->curves; cc && *cc; ++cc)
	{
	    if (I_NumOfCurveSurfaces(*cc) != 0)
	    	continue;
	    if (num_c >= MAX_NUM_OTHER_CURVES)
	    {
		printf("In append_other_curves1(): num_c = %d\n",num_c);
		printf("MAX_NUM_OTHER_CURVES too small!\n");
		clean_up(ERROR);
	    }
	    c[num_c++] = *cc;
	}
	for (cc = adj_intfc->curves; cc && *cc; ++cc)
	{
	    if (I_NumOfCurveSurfaces(*cc) != 0)
	    	continue;
	    if (num_ac >= MAX_NUM_OTHER_CURVES)
	    {
		printf("In append_other_curves1(): num_ac = %d\n",num_ac);
		printf("MAX_NUM_OTHER_CURVES too small!\n");
		clean_up(ERROR);
	    }
	    ac[num_ac++] = *cc;
	}
	for (i = 0; i < num_ac; ++i)
	{
	    mc[i] = NULL;
	    for (j = 0; j < num_c; ++j)
	    {
		if (c[j] == NULL)	/* already matched */
		    continue;
	    	for (ba = ac[i]->first; ba != NULL; ba = ba->next)
	    	{
		    for (b = c[j]->first; b != NULL; b = b->next)
		    {
			bond_matched = NO;
		    	if (bond_match1(b,ba))
		    	{
			    bond_matched = YES;
			    p = (POINT*)find_from_hash_table((POINTER)ba->start,
                                                p_table,p_size);
			    if (p == NULL)
				(void) add_to_hash_table((POINTER)ba->start,
						(POINTER)b->start,
                                                p_table,p_size);
			    else if (p != b->start)
			    {
				printf("Bond start not from hashing table!\n");
				clean_up(ERROR);
			    }
			    p = (POINT*)find_from_hash_table((POINTER)ba->end,
                                                p_table,p_size);
			    if (p == NULL)
				(void) add_to_hash_table((POINTER)ba->end,
						(POINTER)b->end,
                                                p_table,p_size);
			    else if (p != b->end)
			    {
				printf("Bond end not from hashing table!\n");
				clean_up(ERROR);
			    }
			    mc[i] = c[j];
		    	}
			if (bond_matched) break;
		    }
	    	}	    
		if (mc[i] != NULL) 
		{
		    c[j] = NULL;
		    break;
		}
	    }
	}
}	/* end append_other_curves1 */

LOCAL	boolean bond_match1(
	BOND *b,
	BOND *ba)
{
	int i;
	for (i = 0; i < 3; ++i)
	{
	    if (fabs(Coords(b->start)[i] - Coords(ba->start)[i]) > tol1[i])
		return NO;
	    if (fabs(Coords(b->end)[i] - Coords(ba->end)[i]) > tol1[i])
		return NO;
	}
	return YES;
}	/* end bond_match1 */

LOCAL	void merge_gindex_of_curve(
	CURVE *ac,
	RECT_GRID *grid,
	int dir,
	int nb,
	P_LINK *p_table,
	int p_size)
{
	POINT *ps,*pa;
	BOND *b;

	pa = ac->first->start;
	ps = (POINT*)find_from_hash_table((POINTER)pa,p_table,p_size);
	if (ps == NULL)
	{
	    (void) printf("ERROR in merge_gindex_of_curve()\n");
	    (void) printf("Cannot find ps from hash table!\n");
	    clean_up(ERROR);
	}
	merge_point_pair_gindex(ps,pa,grid,dir,nb);

	curve_bond_loop(ac,b)
	{
	    pa = b->end;
	    ps = (POINT*)find_from_hash_table((POINTER)pa,p_table,p_size);
	    if (ps == NULL)
	    {
	    	(void) printf("ERROR in merge_gindex_of_curve()\n");
	    	(void) printf("Cannot find ps from hash table!\n");
	    	clean_up(ERROR);
	    }
	    merge_point_pair_gindex(ps,pa,grid,dir,nb);
	}
}	/* end merge_gindex_of_curve */
