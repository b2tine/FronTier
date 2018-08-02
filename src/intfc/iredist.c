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
*				iredist.c:
*
*       Copyright 1999 by The University at Stony Brook, All rights reserved.
*
*	Containing function of redistribute interface.
*
*/


#include <intfc/int.h>

#define DEBUG_STRING "i_redistribute"

LOCAL	int rb_cell_index(double*,int,RECT_GRID*);
LOCAL  	void equi_redist_curve_seg_o1(CURVE*,BOND*,BOND*,int,double,double,
                                        RECT_GRID*);

enum _SPQ_FLAG {
        SHORTEST = 0,
        MIDDLE = 1,
        LONGEST  = 2
};
typedef enum _SPQ_FLAG SPQ_FLAG;

/* 3D redistribute surface functions */
LOCAL   int     find_scaled_extrem_edge(TRI*,RECT_GRID*,SPQ_FLAG);
LOCAL   boolean    is_tri_in_queue(TRI*,POINTER_Q*);
LOCAL   boolean    delete_min_side_of_tri(TRI*,int,SURFACE*,POINTER_Q**,
                                        INTERFACE*,boolean*);
LOCAL   void    sort_pointer_queue(POINTER_Q*,INTERFACE*,SPQ_FLAG);
LOCAL	int	remove_tris_and_seal(TRI**,TRI**,int,SURFACE*,POINTER_Q**,
					INTERFACE*);
LOCAL 	void 	exchange_queues(POINTER_Q*,POINTER_Q*);
LOCAL 	boolean 	compare_pointers(POINTER_Q*,POINTER_Q*);
LOCAL  	boolean redistribute_surf_o1(SURFACE*,RECT_GRID*,
			SCALED_REDIST_PARAMS);
LOCAL	void reset_neighbors_of_adjacent_sides(TRI*,int);
/*	Functions dealing with tetra tris */
LOCAL 	boolean is_tetra_point(POINT*,TRI*,INTERFACE*);
LOCAL 	boolean remove_interior_tetra_pt(POINT*,TRI*,SURFACE*,POINTER_Q**,
			INTERFACE *);
LOCAL 	void set_neighbors_from_tetra_tris(TRI*,TRI**,int*,boolean*);
LOCAL	void dequeue_tris_around_bond_point(POINT*,BOND*,INTERFACE*,
			POINTER_Q**);
LOCAL 	boolean will_form_foldback(POINT**,TRI*,TRI*,SURFACE*,TRI***,int*);
/*TMP*/ 
LOCAL	boolean tri_in_intfc(TRI*,INTERFACE*);

EXPORT	void equi_redist_curve_seg(
	CURVE		*c,
	BOND		*bs,
	BOND		*be,
	int		nbds,
	double		seg_len,
	double		space,
	RECT_GRID	*rgr)
{
	switch (c->redist_order)
	{
	case 1:
	    equi_redist_curve_seg_o1(c,bs,be,nbds,seg_len,space,rgr);
	    return;
	case 2:
	    equi_redist_curve_seg_o2(c,bs,be,nbds,seg_len,space,rgr);
	    return;
	}
}	/* end equi_redist_curve_seg */

/*
*			equi_redist_curve_seg_o1():
*
*	Redistributes a curve segment from bs->start to be->end by
*	deleting all points between, and inserting N equally spaced points.
*	NOTE that the resultant point spacing is equal only if measured
*	with respect to the arc length of the ORIGINAL, input curve
*	segment and NOT with respect to the final redistributed curve
*	segment.
*	If nbds <= 0, the routine computes N such that the spacing between
*	points on the redistributed curve is <= space.
*	If nbds > 0, the routine uses N == nbds.
*	This is the first order version.
*/

LOCAL	void equi_redist_curve_seg_o1(
	CURVE		*c,
	BOND		*bs,
	BOND		*be,
	int		nbds,
	double		seg_len,
	double		space,
	RECT_GRID	*rgr)
{
	BOND		*b, *bstart, *bend;
	double		b_len, sc_len, offset, s, oms;
	double		coords[MAXD];
	double		*h = rgr->h;
	int		dim = rgr->dim;
	int		new_nbds;

	DEBUG_ENTER(equi_redist_curve_seg)
	if (debugging("high_order_redist"))
	    printf("Entering equi_redist_curve_seg_o1()\n");
	if (nbds <= 0)
	{
	    new_nbds = (int)(ceil((double)seg_len/(double)space));
	    b_len = seg_len/(double)new_nbds;
	}
	else
	{
	    new_nbds = nbds;
	    b_len = seg_len/(double)new_nbds;
	}
	if (new_nbds <= 1)
	{
	    if (is_closed_curve(c) && bs == c->first && be == c->last)
	    {
	    	new_nbds = c->num_points - 1;
	    	if (new_nbds > 1)
	            equi_redist_curve_seg(c,bs,be,new_nbds,
					  seg_len,space,rgr);
	    }
	    else
	        replace_curve_seg_by_bond(c,bs,be);

	    DEBUG_LEAVE(equi_redist_curve_seg)
	    return;
	}

	offset = b_len;
	bstart = bs;		bend = be->next;
	while (bstart != bend)
	{
	    b = bstart;
	    while ((sc_len = scaled_bond_length(b,h,dim)) < offset)
	    {
	    	if (b->next == bend)
	    	{
	    	    replace_curve_seg_by_bond(c,bstart,b);
	    	    goto leave;
	    	}
	    	offset -= sc_len;
	    	b = b->next;
	    }
	    if ((b->next != bend) ||
	        (sc_len >= offset + MIN_SC_SEP(c->interface)))
	    {
	    	s = offset/sc_len;
		interpolate_coords(coords,Coords(b->start),Coords(b->end),s,2);
	    	if (insert_point_in_bond(Point(coords),b,c) !=
		    FUNCTION_SUCCEEDED)
	    	{
	    	    screen("ERROR in equi_redist_curve_seg(), "
	    	           "insert_point_in_bond failed\n");
	    	    clean_up(ERROR);
	    	}
	    }
	    replace_curve_seg_by_bond(c,bstart,b);
	    bstart = bstart->next;
	    offset = b_len;
	}

leave:
	DEBUG_LEAVE(equi_redist_curve_seg)
	return;
}		/*end equi_redist_curve_seg_o1*/


EXPORT	void	rect_bdry_curve_redist(
	CURVE		*c,
	ORIENTATION	orient,
	RECT_GRID	*gr,
	double		*tol)
{
	BOND		*b, *bb;
	double		tolx = tol[0], toly = tol[1];
	double		*ps, *pe;
	double		x, y, coords[MAXD];
	int		ixs, ixe, iys, iye, ix, iy;
	int		ibegin, ifinish;
	ORIENTATION	opor = Opposite_orient(orient);
	int		di;

	DEBUG_ENTER(rect_bdry_curve_redist)
	ps = Coords(c->start->posn);		pe = Coords(c->end->posn);

	if (fabs(ps[1]-pe[1]) < fabs(ps[0]-pe[0]))
	{

	    /* horizontal boundary */

	    ixs = rb_cell_index(ps,0,gr);
	    ixe = rb_cell_index(pe,0,gr);

	    if (ixe == ixs)
	    {
	    	while (c->first != c->last)
		    i_delete_point_adjacent_to_node(c,POSITIVE_ORIENTATION);
	    	DEBUG_LEAVE(rect_bdry_curve_redist)
	    	return;
	    }
	    if (orient == POSITIVE_ORIENTATION)
	    {
	    	di = (ps[0] < pe[0]) ? 1 : -1;
	    	ibegin = ixs;
	    	ifinish = ixe;
	    }
	    else
	    {
	    	di = (ps[0] < pe[0]) ? -1 : 1;
	    	ibegin = ixe;
	    	ifinish = ixs;
	    }

			/* Check positions at nodes */
		
	    x = cell_center(ibegin,0,gr);
	    if ((fabs(x - Coords(Node_of(c,orient)->posn)[0]) < tolx) ||
	    				(!Between(x,ps[0],pe[0])))
	    	ibegin += di;

	    x = cell_center(ifinish,0,gr);
	    if ((fabs(x - Coords(Node_of(c,opor)->posn)[0]) < tolx) ||
	    				(!Between(x,ps[0],pe[0])))
	    	ifinish -= di;

	    if (di*(ibegin - ifinish) > 0)
	    {
	    	while (c->first != c->last)
		    i_delete_point_adjacent_to_node(c,POSITIVE_ORIENTATION);
	    	DEBUG_LEAVE(rect_bdry_curve_redist)
	    	return;
	    }

	    coords[1] = y = Coords(c->start->posn)[1];
	    b = Bond_at_node(c,orient);
	    for (ix = ibegin; di*(ifinish - ix) >= 0; ix += di)
	    {
	        x = cell_center(ix,0,gr);
	        for (bb = b; bb != NULL; bb = Following_bond(bb,orient))
	        {
	    	    if (fabs(x-Coords(Point_of_bond(bb,opor))[0]) < tolx)
	    	    {
	    		Coords(Point_of_bond(bb,opor))[0] = x;
	    		break;
	    	    }
	    	    else if (Between(x,Coords(bb->start)[0],Coords(bb->end)[0]))
	    	    {
	    	        coords[0] = x;
	    	        if (insert_point_in_bond(Point(coords),bb,c) !=
			    FUNCTION_SUCCEEDED)
	    	        {
	    	            screen("ERROR in rect_bdry_curve_redist(), "
	    	                   "insert_point_in_bond failed\n");
	    	            clean_up(ERROR);
	    	        }
	    	        if (orient == NEGATIVE_ORIENTATION) 
	    	        {
	    		    if (bb == b)
				b = bb->next;
	    		    bb = bb->next;
	    	        }
	    	        break;
	            }
	        }
		if (bb == NULL)
		{
		    bb = Bond_at_node(c,opor);
		    if (bb == b) break;
		}
		if (orient == POSITIVE_ORIENTATION)
		{
		    replace_curve_seg_by_bond(c,b,bb);
		    b = b->next;
		}
		else
		{
		    replace_curve_seg_by_bond(c,bb,b);
		    b = bb->prev;
		}
		if (b == NULL) break;
	    }
	    if (orient == POSITIVE_ORIENTATION)
	    	replace_curve_seg_by_bond(c,b,Bond_at_node(c,opor));
	    else
	    	replace_curve_seg_by_bond(c,Bond_at_node(c,opor),b);
	}
	else if (fabs(ps[0]-pe[0]) < fabs(ps[1]-pe[1]))
	{	
	    /* vertical boundary */

	    iys = rb_cell_index(ps,1,gr);
	    iye = rb_cell_index(pe,1,gr);
	    if (iye == iys)
	    {
	    	while (c->first != c->last)
		    i_delete_point_adjacent_to_node(c,POSITIVE_ORIENTATION);
	    	DEBUG_LEAVE(rect_bdry_curve_redist)
	    	return;
	    }
	    if (orient == POSITIVE_ORIENTATION)
	    {
	    	di = (ps[1] < pe[1]) ? 1 : -1;
	    	ibegin = iys;
	    	ifinish = iye;
	    }
	    else
	    {
	    	di = (ps[1] < pe[1]) ? -1 : 1;
	    	ibegin = iye;
	    	ifinish = iys;
	    }

	    /* Check positions at nodes */

	    y = cell_center(ibegin,1,gr);
	    if ((fabs(y - Coords(Node_of(c,orient)->posn)[1]) < toly) ||
	    				(!Between(y,ps[1],pe[1])))
	    	ibegin += di;

	    y = cell_center(ifinish,1,gr);
	    if ((fabs(y - Coords(Node_of(c,opor)->posn)[1]) < toly) ||
	    				(!Between(y,ps[1],pe[1])))
	        ifinish -= di;

	    if (di*(ibegin - ifinish) > 0)
	    {
	    	while (c->first != c->last)
		    i_delete_point_adjacent_to_node(c,POSITIVE_ORIENTATION);
	    	DEBUG_LEAVE(rect_bdry_curve_redist)
	    	return;
	    }

	    coords[0] = x = Coords(c->start->posn)[0];
	    b = Bond_at_node(c,orient);
	    for (iy = ibegin; di*(ifinish - iy) >= 0; iy += di)
	    {
	        y = cell_center(iy,1,gr);
	        for (bb = b; bb != NULL; bb = Following_bond(bb,orient))
	        {
	    	    if (fabs(y-Coords(Point_of_bond(bb,opor))[1]) < toly)
	    	    {
	    		Coords(Point_of_bond(bb,opor))[1] = y;
	    		break;
	    	    }
	    	    else if (Between(y,Coords(bb->start)[1],Coords(bb->end)[1]))
	    	    {
	    	        coords[1] = y;

	    	        if (insert_point_in_bond(Point(coords),bb,c) !=
			    FUNCTION_SUCCEEDED)
	    	        {
	    	            screen("ERROR in rect_bdry_curve_redist(), "
	    	                   "insert_point_in_bond failed\n");
	    	            clean_up(ERROR);
	    	        }
	    	        if (orient == NEGATIVE_ORIENTATION)
	    	        {
	    	            if (bb == b) b = bb->next;
	    	            bb = bb->next;
	    	        }
	    	        break;
	    	    }
	        }  
		if (bb == NULL)
		{
		    bb = Bond_at_node(c,opor);
		    if (bb == b) break;
		}
		if (orient == POSITIVE_ORIENTATION)
		{
		    replace_curve_seg_by_bond(c,b,bb);
		    b = b->next;
		}
		else
		{
		    replace_curve_seg_by_bond(c,bb,b);
		    b = bb->prev;
		}
		if (b == NULL)
		    break;
	    }
	    if (orient == POSITIVE_ORIENTATION)
	    	replace_curve_seg_by_bond(c,b,Bond_at_node(c,opor));
	    else
	    	replace_curve_seg_by_bond(c,Bond_at_node(c,opor),b);
	}
	DEBUG_LEAVE(rect_bdry_curve_redist)
}		/*end rect_bdry_curve_redist*/


EXPORT boolean i_delete_point_adjacent_to_node(
	CURVE           *c,
	ORIENTATION     orient)
{
	BOND		*b = Bond_at_node(c,orient);
	boolean		stat;

	if (orient == POSITIVE_ORIENTATION) 
	{
		if (b->next == NULL) return FUNCTION_FAILED;
		if (debugging("delete_point_adjacent_to_node"))
		{
			(void) printf("deleting start of bond - ");
			print_bond(b->next);
		}

		stat = delete_start_of_bond(b->next,c);
	}
	else 
	{
		if (b->prev == NULL) return FUNCTION_FAILED;
		if (debugging("delete_point_adjacent_to_node"))
		{
			(void) printf("deleting start of bond - ");
			print_bond(b);
		}

		stat = delete_start_of_bond(b,c);
	}

	return stat;
}	/* end i_delete_point_adjacent_to_node */


LOCAL	int rb_cell_index(
	double		*p,
	int		i,
	RECT_GRID	*gr)
{
	int		indx;
	int		*gmax = gr->gmax;
	int		*lbuf = gr->lbuf, *ubuf = gr->ubuf;

	indx = cell_index(p[i],i,gr);
	if (indx < -lbuf[i])
	    indx = -lbuf[i];
	if (indx >= gmax[i] + ubuf[i])
	    indx = gmax[i] + ubuf[i] - 1;
	return indx;
}		/*end rb_cell_index*/

/*
*			replace_curve_seg_by_bond():
*
*	Replaces the curve segment from the start of bond bs to the end
*	of bond be by a single bond. The resultant single bond has the
*	same pointer value as bs.
*/

EXPORT   void replace_curve_seg_by_bond(
	CURVE		*c,
	BOND		*bs,
	BOND		*be)
{
	BOND		*bnext;

	if ((bs == NULL) || (be == NULL))
	    return;
	if (bs == be)
	    return;
	bnext = be->next;
	while (bs->next != bnext)
		(void) delete_start_of_bond(bs->next,c);
}		/*end replace_curve_seg_by_bond*/

EXPORT  void rect_bdry_redist2d(
	INTERFACE	*intfc,
	RECT_GRID	*rgr,
	int		iperm)
{
	CURVE		**c;
	double		tol[MAXD];
	int		i;
	ORIENTATION	orient;

	DEBUG_ENTER(rect_bdry_redist2d)
	for (i = 0; i < rgr->dim; ++i)
	    tol[i] = MIN_SC_SEP(intfc) * rgr->h[i];
	orient = (iperm % 2) ? POSITIVE_ORIENTATION : NEGATIVE_ORIENTATION;
	for (c = intfc->curves; c && *c; ++c) 
	{
	    if (!is_bdry(*c))
	        continue;

	    rect_bdry_curve_redist(*c,orient,rgr,tol);
	}
	DEBUG_LEAVE(rect_bdry_redist2d)
}		/*end rect_bdry_redist2d*/


/*
*			closed_curve_node_redistribute():
*
*	A node on a closed curve never gets redistributed. With the
*	assumption that a node on a closed curve is more of a convenience
*	than an important point, this routine randomly changes the
*	position of the node on a closed curve each time redistribute
*	is called.
*/

EXPORT	boolean closed_curve_node_redistribute(
	INTERFACE	*intfc,
	boolean		status)
{
	BOND		*b;
	CURVE		**c;

        if (debugging("noCCNR"))
	    return status;
	for (c = intfc->curves; *c ; ++c)
	{
	    if (is_bdry(*c))
		continue;
	    if (!is_closed_curve(*c))
		continue;

			/* pick a random bond - any bond */

	    b = random_bond_on_curve(*c);
	    if (!move_closed_loop_node(*c,b))
		status = NO;
	}
	return status;
}		/*end closed_curve_node_redistribute*/

LOCAL	const int	Num_pqs_in_block = 1000;

struct _TRI_SURF {
	TRI	*tri;
	SURFACE *surf;
	CURVE   *c01, *c12, *c20;
	double   sqr_norm, dist;
	int     side;
};
typedef struct _TRI_SURF	TRI_SURF; 

#define Tri_surf(p)			((TRI_SURF *) (p)->pointer)
#define PQ_for_tri(tri)			((POINTER_Q *) Tri_workspace(tri))
#define tri_surface_from_queue(tri)	(Tri_surf(PQ_for_tri(tri)))
#define Bond_of_q(pq)		        ((BOND *)(pq)->pointer)
#define Tri_of_q(pq)		        (Tri_surf(pq)->tri)


LOCAL POINTER_Q *dequeue(
	TRI       *tri,
	POINTER_Q *pq)
{
	if (PQ_for_tri(tri))
	{
	    if (head_of_pointer_queue(PQ_for_tri(tri))
				    == head_of_pointer_queue(pq))
	    {
	        pq = delete_from_pointer_queue(PQ_for_tri(tri));
	        Tri_workspace(tri) = NULL;
	    }
	}
	return pq;
}		/*end dequeue*/

LOCAL POINTER_Q *alloc_and_add_to_queue(
	TRI       *t,
	SURFACE   *s,
	POINTER_Q *pq,
	int       nside)
{
	TRI_SURF    *ts;
	const double *nor;
	int	    i, j;

	pq = add_to_pointer_queue(NULL,pq);
	Tri_workspace(t) = (POINTER) pq;
	ts = tri_surface_from_queue(t);
	if (ts == NULL)
	{
	    screen("ERROR in alloc_and_add_to_queue(), "
	           "tri_surface_from_queue() returns NULL\n");
	    clean_up(ERROR);
	}
	ts->tri = t;
	nor = Tri_normal(t);
	ts->sqr_norm = Dot3d(nor,nor);

	/* a distant for further check, should be shift invariant */
	ts->dist = 0.0;
	for(i=0; i<3; i++)
	    for(j=0; j<3; j++)
	        ts->dist += Coords(Point_of_tri(t)[i])[j];
	
	ts->surf = s;
	ts->c01 = ts->c12 = ts->c20 = NULL;
	ts->side = nside;
	return pq;
}		/*end alloc_and_add_to_queue*/

LOCAL  void  tecplot_tri_queue(
	const char	*msg,
	FILE		*file,
	POINTER_Q 	*p_q)
{
	POINTER_Q 	*q;
	TRI_SURF 	*t_surf;
	POINT		*p;
	TRI 		*tri;
	int		k, i, cnt = 0;
 
	double	pt[3] = { 0.9733961893900565,     -0.301512040795133,     -15.64756169180028  };
	double	pt1[3] = { -1.026603810609944,     -0.301512040795133,     -15.64756169180028 };
		
 
	if (p_q == NULL) 
	{
	    (void) printf("tecplot_tri_queue NULL POINTER_Q %s\n", msg);
	    return;
        }

	q = head_of_pointer_queue(p_q);
	while (q != tail_of_pointer_queue(p_q))
        {
	    t_surf = Tri_surf(q);
	    tri = t_surf->tri;
	    q = q->next;
	    cnt++;
	}
	(void) fprintf(file, "ZONE T=\"%s\" N=%d E=%d\nF=FEPOINT, ET=TRIANGLE\n",
			     msg, 3*cnt, cnt);

	k = 0;
	q = head_of_pointer_queue(p_q);
	while (q != tail_of_pointer_queue(p_q))
        {
	    t_surf = Tri_surf(q);
	    tri = t_surf->tri;
	    for (i = 0; i < 3; i++)
	    {
		p = Point_of_tri(tri)[i];
		fprintf(file,"%-9g %-9g %-9g\n",Coords(p)[0],
				 Coords(p)[1],Coords(p)[2]);
		
		if(distance_between_positions(Coords(Point_of_tri(tri)[i]), pt, 3)<0.05  || 
		   distance_between_positions(Coords(Point_of_tri(tri)[i]), pt1, 3)<0.05 )
		{
		        /*printf("#de tri sort  k = %d\n", k); */
			/*printf("sqr_norm = %24.16e, dist = %24.16e\n", t_surf->sqr_norm, t_surf->dist); */
			/*print_tri(tri, t_surf->surf->interface); */
		}
	    }
	    k++;
	    q = q->next;
	}
	
	for(i=0; i<cnt; i++)
	{
	    fprintf(file, "%d %d %d\n", 3*i+1, 3*i+2, 3*i+3);
	}
}

LOCAL  void  tri_queue_test(
	const char	*msg,
	POINTER_Q 	*p_q)
{
	POINTER_Q 	*q;
	TRI_SURF 	*t_surf;
	POINT		*p;
	TRI 		*tri;
	int		k, i;
	double		tst_pt[3] = {  0.0,     0.8985,                 1.47066};
	double		tst_pt1[3] = { 1.0,     0.8985,                 1.47066};
	double		tol = 1.0/40.0;
 
	if (p_q == NULL) 
	{
	    (void) printf("tri_queue_test NULL POINTER_Q %s\n", msg);
	    return;
        }
	
	k = 0;
	q = head_of_pointer_queue(p_q);
	while (q != tail_of_pointer_queue(p_q))
        {
	    t_surf = Tri_surf(q);
	    tri = t_surf->tri;

	    for (i = 0; i < 3; i++)
	    {
		p = Point_of_tri(tri)[i];
		
		if(distance_between_positions(Coords(Point_of_tri(tri)[i]), tst_pt, 3)<tol  || 
		   distance_between_positions(Coords(Point_of_tri(tri)[i]), tst_pt1, 3)<tol )
		{
		    printf("#de tri sort  k = %d\n", k);
		    printf("sqr_norm = %24.15e, dist = %24.15e\n", t_surf->sqr_norm, t_surf->dist);
		    print_tri(tri, t_surf->surf->interface);
		    
		    break;
		}
	    }

	    k++;
	    q = q->next;
	}
}

EXPORT  boolean redistribute_surf(
	SURFACE		*s,
	RECT_GRID 	*gr,
        SCALED_REDIST_PARAMS scaled_redist_params)
{
	switch (s->redist_order)
	{
	case 1:
	    return redistribute_surf_o1(s,gr,scaled_redist_params);
	case 2:
	    return redistribute_surf_o2(s,gr,scaled_redist_params);
	default:
	    screen("Redistribute order %d not implemented!\n",s->redist_order);
	    clean_up(ERROR);
	}
}	/* end redistribute_surface */

/*	
*	First order surface redistribution function
*/

LOCAL  boolean redistribute_surf_o1(
	SURFACE		*s,
	RECT_GRID 	*gr,
	SCALED_REDIST_PARAMS scaled_redist_params)
{
	INTERFACE *intfc;
	POINT	  *midp;
	POINTER_Q *insert_queue, *delete_queue;
	TRI	  *tri, *oppt;
	boolean   nothing_done = YES;
	int       nside, nf, nl, ns, nt, i;
	int	  dim;
	FILE	  *db_file;
	double	  coords[3];
	boolean	  status,change_side;

	if (debugging("trace"))
	    (void) printf("Entering redistribute_surf_o1()\n");

	set_pointer_queue_opts(PQ_BLOCK_SIZE,Num_pqs_in_block,PQ_ALLOC_TYPE,
		               "vmalloc",PQ_ALLOC_SIZE_FOR_POINTERS,
			       sizeof(TRI_SURF),0);

	if (debugging("high_order_redist"))
	{
	    (void) printf("Entering redistribute_surf_o1()\n");
	    (void) printf("redist_order = %d\n",s->redist_order);
	}
	intfc = s->interface;
	dim = intfc->dim;
	
	/*set the tolerance for tri_status */
	
	insert_queue = delete_queue = NULL;
	nt = nf = nl = ns = 0;
	
	for (tri=first_tri(s); !at_end_of_tri_list(tri,s); tri=tri->next)
	{
	    nt++;
	    switch (tri_scaled_status(tri,gr,scaled_redist_params))
	    {
	    case BAD_ANGLE:	
		++nf;
		break;
	    case LARGE:
		++nl;
		insert_queue = alloc_and_add_to_queue(tri,s,insert_queue,nside);
		break;
	    case SMALL:
		++ns;
		delete_queue = alloc_and_add_to_queue(tri,s,delete_queue,nside);
		break;
	    case GOOD_ANGLE:
	    default:
		Tri_workspace(tri) = NULL;
		break;
	    }
	}
	
	if(insert_queue==NULL && delete_queue==NULL)
	{
	    if (debugging("trace"))
	    	(void) printf("Leaving redistribute_surf_o1()\n");
	    return nothing_done;
	}

	sort_pointer_queue(insert_queue,intfc,LONGEST);
	while (insert_queue)
	{
	    insert_queue = head_of_pointer_queue(insert_queue);
	    tri = Tri_of_q(insert_queue);
	    
	    nside = find_scaled_extrem_edge(tri,gr,LONGEST);
	    insert_queue = dequeue(tri,insert_queue);
	    
	    if (!is_side_bdry(tri,nside))
	    {
	    	oppt = Tri_on_side(tri,nside);	    
	    	if(is_tri_in_queue(oppt,insert_queue))
		    insert_queue = dequeue(oppt,insert_queue);
	    	if(is_tri_in_queue(oppt,delete_queue))
		    delete_queue = dequeue(oppt,delete_queue);
	    }
	    

	    /* find and make tri side mid point */
	    for(i = 0; i < dim; ++i)
		coords[i] = 0.5*(Coords(Point_of_tri(tri)[nside])[i] +
		                 Coords(Point_of_tri(tri)[Next_m3(nside)])[i]);
	    midp = Point(coords);

	    if (!insert_point_in_tri_side(midp,nside,tri,s))
	    {
		printf("WARNING redistribute_surf_o1, "
		       "insert_point_in_tri_side fails.\n");
	    }
	    nothing_done = NO;
	}
	    
	sort_pointer_queue(delete_queue,intfc,SHORTEST);

	while (delete_queue)
	{
	    delete_queue = head_of_pointer_queue(delete_queue);
	    tri = Tri_of_q(delete_queue);

	    nside = find_scaled_extrem_edge(tri,gr,SHORTEST);
		
	    delete_queue = dequeue(tri, delete_queue);

	    status = delete_min_side_of_tri(tri,nside,s,&delete_queue,intfc,
					&change_side);
	    if (!status && change_side)
	    {
		if (debugging("foldback"))
		    (void) printf("WARNING: delete shortest side failed\n");
	    	nside = find_scaled_extrem_edge(tri,gr,MIDDLE);
	    	status = delete_min_side_of_tri(tri,nside,s,&delete_queue,
					intfc,&change_side);
		if (!status && debugging("foldback"))
		    (void) printf("WARNING: delete middle side also failed\n");
	    }
	    if (!status && change_side)
	    {
	    	nside = find_scaled_extrem_edge(tri,gr,LONGEST);
	    	status = delete_min_side_of_tri(tri,nside,s,&delete_queue,
					intfc,&change_side);
		if (!status && debugging("foldback"))
		    (void) printf("WARNING: delete longest side also failed\n");
	    }
	    if (!status)
	        delete_queue = dequeue(tri,delete_queue);
	    if (!status && debugging("foldback"))
	    {
		printf("WARNING, redistribute_surf_o1, "
		       "delete_min_side_of_tri fails.\n");
	    }
	    nothing_done = NO;
	}
	    
	interface_reconstructed(intfc) = NO;
	
	if (debugging("trace"))
	    (void) printf("Leaving redistribute_surf_o1()\n");
	return nothing_done;
}		/*end redistribute_surf_o1*/

LOCAL  int find_scaled_extrem_edge(
	TRI		*tri,
	RECT_GRID	*grid,
	SPQ_FLAG	to_find)
{
	const double* const *s;
	double	h0 = grid->h[0], h1 = grid->h[1], h2 = grid->h[2];
	double	s00, s01, s02;
	double	s10, s11, s12;
	double	s20, s21, s22;
	double	len[3],len_tmp;
	int  side[3],side_tmp;
	int i,j;

	s = side_vector(tri);
	s00 = s[0][0]/h0; s01 = s[0][1]/h1; s02 = s[0][2]/h2;
	s10 = s[1][0]/h0; s11 = s[1][1]/h1; s12 = s[1][2]/h2;
	s20 = s[2][0]/h0; s21 = s[2][1]/h1; s22 = s[2][2]/h2;
	len[0] = QDot3d(s0,s0); len[1] = QDot3d(s1,s1); len[2] = QDot3d(s2,s2);
	for (i = 0; i < 3; ++i) side[i] = i;
	for (i = 0; i < 3; ++i)
	for (j = i+1; j < 3; ++j)
	{
	    if (len[i] > len[j])
	    {
		len_tmp = len[i]; 	side_tmp = side[i];
		len[i] = len[j];	side[i] = side[j];
		len[j] = len_tmp;	side[j] = side_tmp;
	    }
	}

	switch (to_find)
	{
	case SHORTEST:
	    return side[0];
	case MIDDLE:
	    return side[1];
	case LONGEST:
	    return side[2];
	default:
	    return -1;
	}
}		/*end find_scaled_extrem_edge*/


/*
*				tri_status():
*
*	Determines whether redistribution of a triangle is needed by
*	comparing the triangles normalized area with two tolerances.
*	In addition, this routine also checks the squared edge lengths and
*	determines their aspect ratios.
*	This routine has the following return values:
*
*	if (norm_area > max_sqr_area) return LARGE;
*	if (norm_area < min_sqr_area) return SMALL;
*	if (aspect ratio < aspect_tol) return BAD_ANGLE;
*	return GOOD_ANGLE;
*
*	The aspect ratio of a triangle is defined as A/(l0^2+l1^2+l2^2)
*	where A is the area of the triangle in the grid scaled metric
*	and the li are the lengths of the sides of the triangle in
*	grid scaled metric.
*/

EXPORT	TRI_STATUS tri_scaled_status(
	TRI *tri,
	RECT_GRID *gr,
	SCALED_REDIST_PARAMS scaled_redist_params)
{
	double	*p0 = Coords(Point_of_tri(tri)[0]);
	double	*p1 = Coords(Point_of_tri(tri)[1]);
	double	*p2 = Coords(Point_of_tri(tri)[2]);
	double	s00, s01, s02;
	double	s10, s11, s12;
	double	s20, s21, s22;
	double  N0, N1, N2;
	double	h0 = gr->h[0], h1 = gr->h[1], h2 = gr->h[2];
	double	a_ratio,len[3],tri_area;
	double  len_max,len_min;
	double	len2[3],cos_angle[3],angle[3],min_angle;
	int 	i;

	s00 = (p1[0]-p0[0])/h0; s01 = (p1[1]-p0[1])/h1; s02 = (p1[2]-p0[2])/h2;
	s10 = (p2[0]-p1[0])/h0; s11 = (p2[1]-p1[1])/h1; s12 = (p2[2]-p1[2])/h2;
	s20 = (p0[0]-p2[0])/h0; s21 = (p0[1]-p2[1])/h1; s22 = (p0[2]-p2[2])/h2;
	len2[0] = QDot3d(s0,s0);
	len2[1] = QDot3d(s1,s1);
	len2[2] = QDot3d(s2,s2);
	for (i = 0; i < 3; ++i)
	{
	    len[i] = sqrt(len2[i]);
	    if (len[i] > scaled_redist_params.max_scaled_side_length)
	    	return LARGE;
	}

	QCross3d(s0,s2,N);
	tri_area = 0.5*sqrt(QDot3d(N,N));
	if (tri_area < scaled_redist_params.min_scaled_tri_area)
	{
	    return SMALL;
	}
	else if (tri_area > scaled_redist_params.max_scaled_tri_area)
	{
	    return LARGE;
	}
			/* Check aspect ratio	*/
	len_max = -HUGE;
	len_min =  HUGE;
	min_angle = HUGE;
	for (i = 0; i < 3; ++i)
	{
	    cos_angle[i] = (len2[i] + len2[(i+1)%3] - len2[(i+2)%3])
			/2.0/len[i]/len[(i+1)%3];
	    angle[i] = acos(cos_angle[i]);
	    if (len_min > len[i]) len_min = len[i];
	    if (len_max < len[i]) len_max = len[i];
	    if (min_angle > angle[i]) min_angle = angle[i];
	}
	if (min_angle < PI/12.0)
	{
	    return BAD_ANGLE;
	}
	a_ratio = len_max/len_min;
	if (a_ratio > scaled_redist_params.aspect_tol)
	{
	    return BAD_ANGLE;
	}
	return GOOD_ANGLE;
}		/*end tri_status*/


/*TMP*/
LOCAL	boolean tri_in_intfc(
	TRI *tri,
	INTERFACE *intfc)
{
	SURFACE **s;
	TRI *t;
	for (s = intfc->surfaces; s && *s; ++s)
	{
	    for (t = first_tri(*s); !at_end_of_tri_list(t,*s); t = t->next)
		if (t == tri) return YES;
	}
	return NO;
}	

LOCAL boolean delete_min_side_of_tri(
	TRI	  *tri,
	int	  side,
	SURFACE	  *s,
	POINTER_Q **pq,
	INTERFACE *intfc,
	boolean	  *change_side)
{
	TRI	*nbtri, *t, *nbt, **tmp_tris;
	static TRI **new_tris, **in_tris, ***tris;
	POINT	*p[4], *pt, *pmid, *plist[2][100];
	int	i, j, k, nt, np[2], nside, ntris[2];
	boolean	rm_flag;
	static	int	cnt = 0;
	FILE	*file;
	char	fname[100];
	int pr_side,nx_side;
	boolean status;

	DEBUG_ENTER(delete_min_side_of_tri)
	/*
	printf("tri = %d  tri_in_intfc = %d\n",tri,tri_in_intfc(tri,intfc));
	*/
	*change_side = NO;
	if (is_side_bdry(tri,side))
	{
	    BOND *b = Bond_on_side(tri,side);
	    BOND_TRI **btris = Btris(b);
	    CURVE *curve;
	    if (!btris || !*btris)
	    {
		(void) printf("Warning: in delete_min_side_of_tri():\n");
		(void) printf("boundary case, btris or *btris is NULL\n");
		return YES;
	    }
	    curve = (*btris)->curve;
	    pr_side = Prev_m3(side);
	    nx_side = Next_m3(side);
	    if (is_side_bdry(tri,pr_side))
	    {
	    	BOND *bp = Bond_on_side(tri,pr_side);
		if (bp == b->next)
		{
		    dequeue_tris_around_bond_point(bp->start,bp,intfc,pq);
		    status = delete_start_of_bond(bp,curve);
		}
		else if (b->prev != NULL)
		{
		    dequeue_tris_around_bond_point(b->start,b,intfc,pq);
		    status = delete_start_of_bond(b,curve);
		}
		if (!status)
		    *change_side = YES;
		return status;
	    }
	    else if (is_side_bdry(tri,nx_side))
	    {
	    	BOND *bn = Bond_on_side(tri,nx_side);
		if (bn == b->next)
		{
		    dequeue_tris_around_bond_point(bn->start,bn,intfc,pq);
		    status = delete_start_of_bond(bn,curve);
		}
		else if (b->prev != NULL)
		{
		    dequeue_tris_around_bond_point(b->start,b,intfc,pq);
		    status = delete_start_of_bond(b,curve);
		}
		if (!status)
		    *change_side = YES;
		return status;
	    }
	    else
	    {
		if (b->prev != NULL)
		{
		    dequeue_tris_around_bond_point(b->start,b,intfc,pq);
		    status = delete_start_of_bond(b,curve);
		}
		else if (b->next != NULL)
		{
		    dequeue_tris_around_bond_point(b->end,b,intfc,pq);
		    status = delete_start_of_bond(b->next,curve);
		}
		if (!status)
		    *change_side = YES;
	    	return status;
	    }
	}
	if (new_tris == NULL)
	{
	    stat_vector(&new_tris,500,sizeof(TRI*));
	    stat_vector(&in_tris,200,sizeof(TRI*));
	    stat_matrix(&tris,2,200,sizeof(TRI*));
	}

	p[0] = Point_of_tri(tri)[side];
	p[1] = Point_of_tri(tri)[Next_m3(side)];
	if (Boundary_point(p[0]) && Boundary_point(p[1]))
	    return YES;

	nbtri = Tri_on_side(tri,side);
	for (nside = 0; nside < 3; nside++)
	{
	    if (Tri_on_side(nbtri,nside) == tri)
		break;
	}

	p[2] = Point_of_tri(tri)[Prev_m3(side)];
	if (is_tetra_point(p[2],tri,intfc))
	{
	    if (debugging("tetra_pt"))
	    	(void) printf("p[2] is tetra point\n");
	    return remove_interior_tetra_pt(p[2],tri,s,pq,intfc);
	}
	p[3] = Point_of_tri(nbtri)[Prev_m3(nside)];
	if (is_tetra_point(p[3],nbtri,intfc))
	{
	    if (debugging("tetra_pt"))
	    	(void) printf("p[3] is tetra point\n");
	    return remove_interior_tetra_pt(p[3],nbtri,s,pq,intfc);
	}

	for(k=0; k<2; k++)
	{
	    ntris[k] = set_tri_list_around_point(p[k],tri,&tmp_tris,intfc);
	    for (i = 0; i < ntris[k]; i++)
	    {
		tris[k][i] = tmp_tris[i];
	    }
	}
	if (ntris[0] + ntris[1] > 12) /* form too small angle */
	{
	    *change_side = YES;
	    return NO;
	}

	for(k=0; k<2; k++)
	{
	    for (i = 0; i < ntris[k]; i++)
	        *pq = dequeue(tris[k][i],*pq);
	    np[k] = 0;
	    /*finding bounding points except the 4 common points. */
	    for(i=0; i<ntris[k]; i++)
	    {
		int l;
		t = tris[k][i];
		for (l = 0; l < 3; ++l)
		{
		    boolean not_in_lists = YES;
		    pt = Point_of_tri(t)[l];
		    for (j = 0; j < 4; ++j)
		    {
		    	if (pt == p[j])
			{
			    not_in_lists = NO;
			    break;
			}
		    }
		    for (j = 0; j < np[k]; ++j)
		    {
		    	if (pt == plist[k][j])
			{
			    not_in_lists = NO;
			    break;
			}
		    }
		    if (not_in_lists == YES)
		    {
			plist[k][np[k]] = pt;
			np[k]++;
		    }
		}
	    }
	}

	/*check if there are duplicate points in the bounding tris. */
	rm_flag = NO;
	if(np[0] > 0 && np[1] > 0)
	{
	    /*the general case, test duplicate points */
	    for(i=0; i<np[0]; i++)
		for(j=0; j<np[1]; j++)
		    if(plist[0][i] == plist[1][j])
			rm_flag = YES;
	}
	else if(np[0] == 0 && np[1] == 0)
	{
	    /*the tetrahedron case */
	    rm_flag = YES;
	}

	if(rm_flag)
	{
	    nt = 0;
	    for(k=0; k<2; k++)
	    for(k=0; k<2; k++)
		nt = merge_tris_set(in_tris, nt, tris[k], ntris[k]);

	    nt = remove_tris_and_seal(new_tris, in_tris, nt, s, pq, intfc);
	
	    DEBUG_LEAVE(delete_min_side_of_tri)
	    return   nt==-1 ? NO : YES;
	}

	if (will_form_foldback(p,tri,nbtri,s,tris,ntris))
	{
	    *change_side = YES;
	    return NO;
	}
		

	/*collapse two tris. */
	if (Boundary_point(p[0]))
	    pmid = p[0];
	else if (Boundary_point(p[1]))
	    pmid = p[1];
	else
	    pmid = average_points(YES,p[0],Hyper_surf_element(tri),
		    Hyper_surf(s),p[1],Hyper_surf_element(tri),Hyper_surf(s));
	
	/*change the point for the surrounding tris. */
	for (i = 0; i < 2; ++i)
	{
	    if (p[i] == pmid) continue;
	    for (j = 0; j < ntris[i]; ++j)
	    {
		t = tris[i][j];
		k = Vertex_of_point(t,p[i]);
		Point_of_tri(t)[k] = pmid;
		if ((t != tri) && (t != nbtri))
		{
		    set_normal_of_tri(t);
		}
	    }
	}

	/*change tri neighbor for tri. */
	reset_neighbors_of_adjacent_sides(tri,side);
	/*change tri neighbor for nbtri. */
	reset_neighbors_of_adjacent_sides(nbtri,nside);

	remove_tri_from_surface(tri,s,YES);
	remove_tri_from_surface(nbtri,s,YES);
	
	DEBUG_LEAVE(delete_min_side_of_tri)
	return YES;
}	/* end delete_min_side_of_tri */


/*ARGSUSED*/
LOCAL void sort_pointer_queue(
	POINTER_Q	*pq,
	INTERFACE	*intfc,
	SPQ_FLAG	flag)
{
	POINTER_Q	*pq1,*pq2;
	
	if (pq == NULL)
	    return;
	pq1 = head_of_pointer_queue(pq);

	pq1 = head_of_pointer_queue(pq);
	while (pq1 != tail_of_pointer_queue(pq))
	{
	    pq2 = pq1->next;
	    if (flag == SHORTEST)
	    {
		if (compare_pointers(pq1, pq2))
	    	    exchange_queues(pq1,pq2);
	        
		while (pq2 != tail_of_pointer_queue(pq))
	        {
	    	    pq2 = pq2->next;
		    if (compare_pointers(pq1, pq2))
		    {
	    	        exchange_queues(pq1,pq2);
	    	    }
	        }
	    }
	    else if (flag == LONGEST)
	    {
		if (!compare_pointers(pq1, pq2))
	    	    exchange_queues(pq1,pq2);
	        
		while (pq2 != tail_of_pointer_queue(pq))
	        {
	    	    pq2 = pq2->next;
		    if (!compare_pointers(pq1, pq2))
		    {
	    	        exchange_queues(pq1,pq2);
	    	    }
	        }
	    }
	    pq1 = pq1->next;
	}
}		/*end sort_pointer_queue*/

/*COND: there should be no boundary point in the bound of tris */
LOCAL	int	remove_tris_and_seal(
	TRI		**new_tris,
	TRI		**tris,
	int		nt,
	SURFACE		*s,
	POINTER_Q	**pq,
	INTERFACE	*intfc)
{
	TRI	*out_tris[500], **new_out_tris;
	int	i, num_out_tris, num_new_tris;

	DEBUG_ENTER(remove_tris_and_seal);
	
	num_out_tris = bound_tris_set(out_tris, tris, nt);

	/*tris are already dequeue in delete_min_side_of_tri */
	for(i=0; i<nt; i++)
	    remove_tri_from_surface(tris[i], s, NO);
	    
	if(num_out_tris == 0)
	{
	    DEBUG_LEAVE(remove_tris_and_seal);
	    return 0;
	}

	sep_common_point_from_loop(out_tris, num_out_tris, NULL, NULL, intfc);
	
	/*new null side tris can be added into out_tris */
	num_out_tris = sep_common_edge_from_tris(&new_out_tris, 
				out_tris, num_out_tris, intfc);

	/*since a smooth_null_loop is applied above, the positions of  */
	/*3 vertics of a tri is changed, all the bound tris should be  */
	/*removed from the que. */
	
	for(i=0; i<num_out_tris; i++)
	    *pq = dequeue(new_out_tris[i], *pq);

	num_new_tris = 0;
	nt = seal_all_loops_wo_constraint(new_tris, &num_new_tris, 
			new_out_tris, num_out_tris, 0, NO);
	nt = merge_tris_set(new_tris, num_new_tris, new_out_tris, nt);

	DEBUG_LEAVE(remove_tris_and_seal);
	return nt;
}	/* end remove_tris_and_seal */

LOCAL	void	exchange_queues(
	POINTER_Q *pq1,
	POINTER_Q *pq2)
{
	TRI_SURF *ts1, *ts2, T;

	ts1 = Tri_surf(pq1);
	ts2 = Tri_surf(pq2);
	T = *ts1;
	*ts1 = *ts2;
	*ts2 = T;
	Tri_workspace(ts1->tri) = (POINTER) pq1;
	Tri_workspace(ts2->tri) = (POINTER) pq2;
}		/*end exchange_queues*/

LOCAL 	boolean 	compare_pointers(
	POINTER_Q	*pq1,
	POINTER_Q	*pq2)
{
	double	ave_norm, norm1, norm2;
	double	tol = 1.0e-8;

	norm1 = Tri_surf(pq1)->sqr_norm;
	norm2 = Tri_surf(pq2)->sqr_norm;
	ave_norm = (norm1 + norm2)*0.5;

	/*two tris have very similar area, compare the postions. */
	if( fabs(norm1 - norm2) < ave_norm*tol ) 
	    return   Tri_surf(pq1)->dist > Tri_surf(pq2)->dist;
	else
	    return   norm1 > norm2;
}

LOCAL boolean is_tri_in_queue(
	TRI		*tri,
	POINTER_Q		*pq)
{
	POINTER_Q		*tri_q;

	for (tri_q = head_of_pointer_queue(pq); tri_q; tri_q = tri_q->next)
	{
	    if (Tri_of_q(tri_q) == tri)
		return YES;
	}
	return NO;
}		/*end is_tri_in_queue*/


EXPORT boolean redistribute_curve(
	CURVE *curve,
	RECT_GRID *gr,
	SCALED_REDIST_PARAMS scaled_redist_params)
{
	BOND **ordered_bonds,*b;
	double *length_bonds,len;
	int i,j,nb;
	double len_ubound,len_lbound;
	boolean nothing_done = YES;
	double *h = gr->h;
	int dim = gr->dim;

	len_ubound = scaled_redist_params.max_scaled_bond_length;
	len_lbound = scaled_redist_params.min_scaled_bond_length;

	if (debugging("redist_curve"))
	{
	    (void) printf("\nEntering redistribute_curve()\n");
	    (void) printf("len_ubound = %f  len_lbound = %f\n",
				len_ubound,len_lbound);
	}

	uni_array(&ordered_bonds,curve->num_points,sizeof(BOND*));
	uni_array(&length_bonds,curve->num_points,sizeof(double));

	nb = 0;
	for (b = curve->first; b != NULL; b = b->next)
	{
	    length_bonds[nb] = scaled_bond_length(b,h,dim);
	    ordered_bonds[nb++] = b;
	}

	for (i = 0; i < nb-1; ++i)
	{
	    for (j = i+1; j < nb; ++j)
	    {
		if (length_bonds[j] < length_bonds[i])
		{
		    b = ordered_bonds[i];
		    ordered_bonds[i] = ordered_bonds[j];
		    ordered_bonds[j] = b;
		    len = length_bonds[i];
		    length_bonds[i] = length_bonds[j];
		    length_bonds[j] = len;
		}
	    }
	}
	if (debugging("redist_curve"))
	{
	    for (i = 0; i < nb; ++i)
	    	(void) printf("ordered bond_length [%d] = %f\n",
					i,length_bonds[i]);
	}
	for (i = 0; i < nb; ++i)
	{
	    double forward_len,backward_len;
	    b = ordered_bonds[i];
	    if (b == NULL) continue;
	    forward_len = (b != curve->last) ? scaled_bond_length(b,h,dim) + 
			scaled_bond_length(b->next,h,dim) : HUGE;
	    backward_len = (b != curve->first) ? scaled_bond_length(b,h,dim) + 
			scaled_bond_length(b->prev,h,dim) : HUGE;
	    if (forward_len < backward_len && 
		scaled_bond_length(b,h,dim) < len_lbound)
	    {
		nothing_done = NO;
		for (j = 0; j < nb; ++j)
		{
		    if (ordered_bonds[j] == b || ordered_bonds[j] == b->next)
			ordered_bonds[j] = NULL;
		}
		delete_start_of_bond(b->next,curve);
	    }
	    else if (backward_len < forward_len && 
		     scaled_bond_length(b,h,dim) < len_lbound)
	    {
		nothing_done = NO;
		for (j = 0; j < nb; ++j)
		{
		    if (ordered_bonds[j] == b || ordered_bonds[j] == b->prev)
			ordered_bonds[j] = NULL;
		}
		delete_start_of_bond(b,curve);
	    }
	    else if (scaled_bond_length(b,h,dim) > len_ubound)
	    {
		double newp[MAXD];
		nothing_done = NO;
		for (j = 0; j < nb; ++j)
		{
		    if (ordered_bonds[j] == b || ordered_bonds[j] == b->prev)
			ordered_bonds[j] = NULL;
		}
		interpolate_coords(newp,Coords(b->start),Coords(b->end),0.5,3);
		insert_point_in_bond(Point(newp),b,curve);
	    }
	}
	free_these(2,ordered_bonds,length_bonds);
	if (debugging("redist_curve"))
	{
	    (void) printf("Leaving redistribute_curve()\n");
	    (void) printf("Checking consistency of interface\n");
	    consistent_interface(curve->interface);
	    (void) printf("Check complete\n");
	}
	return nothing_done;
}	/* end redistribute_curve */

LOCAL	void reset_neighbors_of_adjacent_sides(
	TRI *tri,
	int side)
{
	int prev_side = Prev_m3(side);
	int next_side = Next_m3(side);
	TRI *next_tri,*prev_tri;
	BOND *b;
	BOND_TRI **btris;
	int i;

	if (is_side_bdry(tri,prev_side) && is_side_bdry(tri,next_side))
	{
	    (void) printf("ERROR: In reset_neighbors_of_adjacent_sides():\n");
	    (void) printf("both adjacent sides are boundaries\n");
	    clean_up(ERROR);
	}
	else if (is_side_bdry(tri,prev_side))
	{
	    int bside;
	    b = Bond_on_side(tri,prev_side);
	    next_tri = Tri_on_side(tri,next_side);
	    for (i = 0; i < 3; ++i)
	    {
		if (Tri_on_side(next_tri,i) == tri)
		{
		    set_side_bdry(Boundary_tri(next_tri),i,YES);
		    bside = i;
	    	    for (btris = Btris(b); btris && *btris; ++btris)
	    	    {
			if ((*btris)->tri == tri)
			{
		    	    (*btris)->tri = next_tri;
		    	    Bond_tri_on_side(next_tri,i) = *btris;
			}
	    	    }
		    break;
		}
	    }
	    b = Bond_on_side(next_tri,bside);
	}
	else if (is_side_bdry(tri,next_side))
	{
	    int bside;
	    b = Bond_on_side(tri,next_side);
	    prev_tri = Tri_on_side(tri,prev_side);
	    for (i = 0; i < 3; ++i)
		if (Tri_on_side(prev_tri,i) == tri)
		{
		    set_side_bdry(Boundary_tri(prev_tri),i,YES);
		    bside = i;
	    	    for (btris = Btris(b); btris && *btris; ++btris)
	    	    {
			if ((*btris)->tri == tri)
			{
		    	    (*btris)->tri = prev_tri;
		    	    Bond_tri_on_side(prev_tri,i) = *btris;
			}
	    	    }
		    break;
		}
	    b = Bond_on_side(prev_tri,bside);
	}
	else
	{
	    prev_tri = Tri_on_side(tri,prev_side);
	    next_tri = Tri_on_side(tri,next_side);
	    for (i = 0; i < 3; ++i)
		if (Tri_on_side(prev_tri,i) == tri)
		    Tri_on_side(prev_tri,i) = next_tri;
	    for (i = 0; i < 3; ++i)
		if (Tri_on_side(next_tri,i) == tri)
		    Tri_on_side(next_tri,i) = prev_tri;
	}
}	/* end reset_neighbors_of_adjacent_sides */

LOCAL	void dequeue_tris_around_bond_point(
	POINT *p,
	BOND *b,
	INTERFACE *intfc,
	POINTER_Q **pq)
{
	int i,nt;
	BOND_TRI **btris;
	TRI *tri,**tris;
	for (btris = Btris(b); btris && *btris; ++btris)
	{
	    tri = (*btris)->tri;
	    nt = set_tri_list_around_point(p,tri,&tris,intfc);
	    for (i = 0; i < nt; ++i)
	    {
	    	*pq = dequeue(tris[i],*pq);
	    }
	}
}	/* end dequeue_tris_around_bond_point */
	
LOCAL boolean is_tetra_point(
	POINT *p,
	TRI *tri,
	INTERFACE *intfc)
{
	int nt;
	TRI **tris;
	if (Boundary_point(p)) return NO;
	nt = set_tri_list_around_point(p,tri,&tris,intfc);
	return (nt == 3) ? YES : NO;
}	/* end is_tetra_point */

LOCAL boolean remove_interior_tetra_pt(
	POINT *p,
	TRI *tri,
	SURFACE *s,
	POINTER_Q **pq,
        INTERFACE *intfc)
{
	int i,j,nt;
	TRI *new_t,*t,**tris;
	int side[3];
	boolean is_bdry[3];
	POINT *pts[3];
	boolean status;

	if (debugging("tetra_pt"))
	{
	    (void) printf("Entering remove_interior_tetra_pt()\n");
	    status = consistent_interface(intfc);
	    (void) printf("Consisteny status = %d\n",status);
	}
	nt = set_tri_list_around_point(p,tri,&tris,intfc);
	if (nt != 3 || Boundary_point(p))
	{
	    (void) printf("Is not a tetra point, should not have entered"
			  "this function remove_interior_tetra_pt()\n");
	}
	for (i = 0; i < nt; ++i)
	{
	    t = tris[i];
	    for (j = 0; j < 3; ++j)
	    {
		if (Point_of_tri(t)[j] == p)
		{
		    side[i] = Next_m3(j);
		    pts[i] = Point_of_tri(t)[side[i]];
		    if (is_side_bdry(t,side[i]))
			is_bdry[i] = YES;
		    else
			is_bdry[i] = NO;
		}
	    }
	}
	new_t = make_tri(pts[0],pts[1],pts[2],NULL,NULL,NULL,0);
	insert_tri_at_tail_of_list(new_t,s);
	set_neighbors_from_tetra_tris(new_t,tris,side,is_bdry);
	for (i = 0; i < nt; ++i)
	{
	    remove_tri_from_surface(tris[i],s,YES);
	    *pq = dequeue(tris[i],*pq);
	}
	if (debugging("tetra_pt"))
	{
	    (void) printf("Leaving remove_interior_tetra_pt()\n");
	    status = consistent_interface(intfc);
	    (void) printf("Consisteny status = %d\n",status);
	}
	return YES;
}	/* end remove_interior_tetra_pt */

LOCAL void set_neighbors_from_tetra_tris(
	TRI *new_tri,
	TRI **old_tris,
	int *side,
	boolean *is_bdry)
{
	int i,j;
	TRI *t,*nbt;
	BOND *b;
	BOND_TRI *btri,**bt;

	for (i = 0; i < 3; ++i)
	{
	    t = old_tris[i];
	    if (is_bdry[i])
	    {
		b = Bond_on_side(t,side[i]);
		for (bt = Btris(b); bt && *bt; ++bt)
		{
		    btri = *bt;
		    if (btri->tri == t) break;
		}
		link_tri_to_bond(btri,new_tri,btri->surface,b,btri->curve);
		
	    }
	    else
	    {
		nbt = Tri_on_side(t,side[i]);
		Tri_on_side(new_tri,i) = nbt;
		for (j = 0; j < 3; ++j)
		{
		    if (Tri_on_side(nbt,j) == t)
		    {
			Tri_on_side(nbt,j) = new_tri;
			break;
		    }
		}
	    }
	}
}	/* end set_neighbors_from_tetra_tris */

LOCAL boolean will_form_foldback(
	POINT **p,
	TRI *tri,
	TRI *nbtri,
	SURFACE *s,
	TRI ***tris,
	int *ntris)
{
	int i,j,k;
	TRI *t;
	double old_nor[MAXD],new_nor[MAXD];
	double v1[MAXD],v2[MAXD],ptmp[MAXD],*p1,*p2;
	double old_length,new_length,angle;

	/* First attempt */
	if (Boundary_point(p[0]))
	{
	    for (i = 0; i < 3; ++i)
		ptmp[i] = Coords(p[0])[i];
	}
	else if (Boundary_point(p[1]))
	{
	    for (i = 0; i < 3; ++i)
		ptmp[i] = Coords(p[1])[i];
	}
	else
	{
	    for (i = 0; i < 3; ++i)
		ptmp[i] = 0.5*(Coords(p[0])[i] + Coords(p[1])[i]);
	}

	for (i = 0; i < 2; ++i)
	{
	    for (j = 0; j < ntris[i]; ++j)
	    {
		t = tris[i][j];
		if ((t != tri) && (t != nbtri))
		{
		    for (k = 0; k < 3; ++k)
		    	old_nor[k] = Tri_normal(t)[k];
		    old_length = Mag3d(old_nor);
		    for (k = 0; k < 3; ++k)
		    	old_nor[k] /= old_length;
		    k = Vertex_of_point(t,p[i]);
		    p1 = Coords(Point_of_tri(t)[Next_m3(k)]);
		    p2 = Coords(Point_of_tri(t)[Prev_m3(k)]);
		    for (k = 0; k < 3; ++k)
		    {
			v1[k] = p1[k] - ptmp[k];
			v2[k] = p2[k] - ptmp[k];
		    }
		    Cross3d(v1,v2,new_nor);
		    new_length = Mag3d(new_nor);
		    for (k = 0; k < 3; ++k)
			new_nor[k] /= new_length;
		    if (new_length < 1.0e-8*old_length)
                        return YES;
		    angle = acos(Dot3d(old_nor,new_nor));
		    if (angle > PI/2.0)
		    {
		    	if (debugging("foldback"))
		    	{
		    	    (void) printf("old_nor = %f %f %f\n",old_nor[0],
					old_nor[1],old_nor[2]);
		    	    (void) printf("new_nor = %f %f %f\n",new_nor[0],
					new_nor[1],new_nor[2]);
			    (void) printf("angle = %f\n",angle);
		    	}
			return YES;
		    }
		}
	    }
	}
	return NO;
}	/* end will_form_foldback */

EXPORT void interpolate_coords(
	double *ans,
	double *coords1,
	double *coords2,
	double t,
	int dim)
{
	int i;
	for (i = 0; i < dim; ++i)
	    ans[i] = (coords1[i] == coords2[i]) ? coords1[i] :
			(1.0 - t)*coords1[i] + t*coords2[i];
}	/* end interpolate_coords */
