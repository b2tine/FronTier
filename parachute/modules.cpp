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

static void initSingleModule(Front*);
static void initMultiModule(Front*,int);
static void MergeTwoIntfc(INTERFACE*,INTERFACE*);
static void CopyNodeInfo(INTERFACE*,INTERFACE*);
static void modifyCanopySet(FILE*,Front*,SURFACE*);
static boolean curve_of_boundary_hs(CURVE*);
static void setCanopyBodyIndex(Front*);

static void setCanopyBodyIndex(Front *front)
{
	int count = 1;
	SURFACE **s;

	intfc_surface_loop(front->interf, s)
	{
	    if (wave_type(*s) == MOVABLE_BODY_BOUNDARY)
		count++;
	}
	intfc_surface_loop(front->interf, s)
	{
	    if (wave_type(*s) == ELASTIC_BOUNDARY)
		body_index(*s) = count++;
	}
}	/* end setCanopyBodyIndex */

static void divideAtGoreBdry(INTERFACE *intfc)
{
	SURFACE **s,*surf;
	int num_pieces;

	if (debugging("divide_gore_surf"))
	    (void) printf("Entering divideAtGoreBdry()\n");

start_loop:
	intfc_surface_loop(intfc,s)
	{
	    if (wave_type(*s) != ELASTIC_BOUNDARY) continue;
	    surf = *s;
	    num_pieces = color_bounded_surf_pieces(surf);
	    if (num_pieces == 1) continue;
	    if (debugging("divide_gore_surf"))
	    	(void) printf("divide into num_pieces: %d\n",num_pieces);
	    split_surf_with_colors(surf, num_pieces);
	    goto start_loop;
	}
	if (debugging("divide_gore_surf"))
	{
	    int count = 0;
            CURVE **c;
            intfc_surface_loop(intfc, s)
            {
            	if (wave_type(*s) == ELASTIC_BOUNDARY) 
		{
                    int pc = 0, nc = 0;
                    surf_pos_curve_loop(*s, c)
                    	pc++;
                    surf_neg_curve_loop(*s, c)
                    	nc++;
                    (void) printf("for surf %d: pc = %d nc = %d total = %d\n",
                        count, pc, nc, pc+nc);
                    count++;
            	}
            }
	    (void) printf("%d elastic surfaces\n", count);
	    if (!consistent_interface(intfc))
             	(void) printf("intfc is not consistent\n");
            else
             	(void) printf("consistent interface\n");
            (void) printf("Leaving divideAtGoreBdry()\n");
	}
}	/* end divideAtGoreBdry */

extern void initParachuteModules(Front *front)
{
	int i,num_canopy;
	FILE *infile = fopen(InName(front),"r");
	SURFACE *surf; //not being used?
	static RG_PARAMS rgb_params;

	if (debugging("trace"))
	    (void) printf("Entering initParachuteModules()\n");

	if (debugging("set_module"))
	    gview_plot_interface("module-step-1",front->interf);

	initRigidBody(front);
	rgb_init(front,rgb_params);

	CursorAfterString(infile,"Enter number of canopy surfaces:");
        fscanf(infile,"%d",&num_canopy);
        (void) printf("%d\n",num_canopy);
	fclose(infile);

	if (num_canopy == 1)
	    initSingleModule(front);
	else
	    initMultiModule(front,num_canopy);

	divideAtGoreBdry(front->interf);
	setCanopyBodyIndex(front);

	if (debugging("trace"))
	{
	    gview_plot_interface("ginit",front->interf);
	    (void) printf("Leaving initParachuteModules()\n");
	}
}	/* end initParachuteModules */

extern void initParachuteDefault(
	Front *front)
{
	AF_PARAMS *af_params = (AF_PARAMS*)front->extra2;
	FILE *infile = fopen(InName(front),"r");
        char string[100];
	af_params->is_parachute_system = YES;
	af_params->num_opt_round = 20;
        af_params->spring_model = MODEL1;
	af_params->gore_len_fac = 1.0;
	if (CursorAfterStringOpt(infile,
            "Enter yes to attach gores to canopy:"))
        {
            fscanf(infile,"%s",string);
            if (string[0] == 'y' || string[0] == 'Y')
                af_params->attach_gores = YES;
        }
        fclose(infile);
}	/* end initParachuteDefault */

static void initSingleModule(
        Front *front)
{
	FILE *infile = fopen(InName(front),"r");
	SURFACE *surf;

	CgalCanopySurface(infile,front,&surf);
	fclose(infile);

	InstallNewLoadNode(front,1);
}	/* end initSingleModule */

static void initMultiModule(
	Front *front,
	int num_canopy)
{
	INTERFACE *intfc = front->interf;
	SURFACE **surfs;
	FILE *infile = fopen(InName(front),"r");
	double center[MAXD];
	double phi,theta;
	char string[100];
	int i;
        INTERFACE *cur_intfc;
        cur_intfc = current_interface();

	FT_VectorMemoryAlloc((POINTER*)&surfs,num_canopy,sizeof(SURFACE*));

	for (i = 0; i < num_canopy; ++i)
	{
	    sprintf(string,"For module %d", i+1);
	    CursorAfterString(infile, string);
	    printf("\n");
	    CgalCanopySurface(infile,front,&surfs[i]);
	    modifyCanopySet(infile,front,surfs[i]);
	}
	fclose(infile);

	InstallNewLoadNode(front,num_canopy);
	set_current_interface(cur_intfc);
	FT_FreeThese(1,surfs);
}	/* end initMultiModule */

static void CopyNodeInfo(
	INTERFACE *intfc,
	INTERFACE *newintfc)
{
	AF_NODE_EXTRA *extra, *tmp;
        INTERFACE *cur_intfc;
	NODE **n,**newn;

        cur_intfc = current_interface();
	set_current_interface(newintfc);
	intfc_node_loop(intfc,n)
	{
	    intfc_node_loop(newintfc,newn)
	    {
		if (fabs(Coords((*n)->posn)[0]-Coords((*newn)->posn)[0]) < 1e-6
		 && fabs(Coords((*n)->posn)[1]-Coords((*newn)->posn)[1]) < 1e-6
		 && fabs(Coords((*n)->posn)[2]-Coords((*newn)->posn)[2]) < 1e-6)
		{
		    tmp = (AF_NODE_EXTRA*)((*n)->extra);
		    if (tmp != NULL)
		    {
            	        FT_ScalarMemoryAlloc((POINTER*)&extra,
				sizeof(AF_NODE_EXTRA));
            	    	extra->af_node_type = tmp->af_node_type;
            	        (*newn)->extra = (POINTER)extra;
            	        (*newn)->size_of_extra = sizeof(AF_NODE_EXTRA);
		    } 
		    else
			(*newn)->extra = NULL;
		    break;
		}
	    }
	}
	set_current_interface(cur_intfc);
}	/* end CopyNodeInfo */


static void MergeTwoIntfc(
	INTERFACE *intfc,
	INTERFACE *tmp_intfc)
{
	P_LINK    *p_table;
	int p_size;
	SURFACE   **tmp_s,*news;
	CURVE **tmp_c,*newc;
	NODE **tmp_n,*newn;
	AF_NODE_EXTRA *extra, *tmp_ex;
        INTERFACE *cur_intfc;

        cur_intfc = current_interface();
	set_current_interface(intfc);

	p_size = 4*(tmp_intfc->num_points) + 1;
        uni_array(&p_table,p_size,sizeof(P_LINK));
        reset_hash_table(p_table,p_size);

	for (tmp_s = tmp_intfc->surfaces; tmp_s && *tmp_s; ++tmp_s)
	{
	    if (is_bdry_hs(Hyper_surf(*tmp_s)))
		continue;
	    news = copy_buffer_surface(*tmp_s,p_table,p_size);
	    Hyper_surf_index(news) = Hyper_surf_index((*tmp_s));
	    
	}
        intfc_curve_loop(tmp_intfc,tmp_c)
        {
            if (curve_of_boundary_hs(*tmp_c))
                continue;
	    newc = matching_curve(*tmp_c,p_table,p_size);
	    hsbdry_type(newc) = hsbdry_type(*tmp_c);
	}
	intfc_node_loop(tmp_intfc,tmp_n)
	{
	    tmp_ex = (AF_NODE_EXTRA*)((*tmp_n)->extra);
	    if (NULL == tmp_ex)
		continue;
	    newn = matching_node((*tmp_n),p_table,p_size);
	    FT_ScalarMemoryAlloc((POINTER*)&extra,sizeof(AF_NODE_EXTRA));
	    extra->af_node_type = tmp_ex->af_node_type;
	    newn->extra = (POINTER)extra;	    
	    newn->size_of_extra = sizeof(AF_NODE_EXTRA);	    
	}
	free(p_table);
	set_current_interface(cur_intfc);
}	/* end MergeTwoIntfc */

static boolean curve_of_boundary_hs(
        CURVE *c)
{
        SURFACE **s;
        curve_pos_surf_loop(c,s)
        {
            if (Boundary_hs(Hyper_surf(*s)))
                return YES;
        }
        curve_neg_surf_loop(c,s)
        {
            if (Boundary_hs(Hyper_surf(*s)))
                return YES;
        }
        return NO;
}       /* end curve_of_boundary_hs */

static void modifyCanopySet(
	FILE *infile,
	Front *front,
	SURFACE *canopy)
{
	char string[200];
	double displacement[MAXD];	// Translation displacement
	double center[MAXD];		// Center of rotation
	double phi,theta;		// Spherical angles of rotation
	int i;
	int nc,nn;			// nc, nn: number of curves and nodes;
	CURVE **c,*curves[500];
	NODE **n,*nodes[500];
	TRI *tri;
	BOND *b;
	POINT *p;

        if (CursorAfterStringOpt(infile,
            "Entering yes to modify initialization:"))
        {
            fscanf(infile,"%s",string);
            (void) printf("%s\n",string);
            if (string[0] != 'y' && string[0] != 'Y')
                return;
        }
	else
	    return;

        CursorAfterString(infile,
                "Enter yes for rotation of canopy:");
        fscanf(infile,"%s",string);
        (void) printf("%s\n",string);
        if (string[0] == 'y' || string[0] == 'Y')
        {
            CursorAfterString(infile,"Enter center of rotation:");
            fscanf(infile,"%lf %lf %lf",center,center+1,center+2);
            (void) printf("%f %f %f\n",center[0],center[1],center[2]);
            CursorAfterString(infile,"Enter azimuthal and polar angles:");
            fscanf(infile,"%lf %lf",&phi,&theta);
            (void) printf("%f %f\n",phi,theta);
            theta *= PI/180.0;
            phi *= PI/180.0;
        }

	/* Assemble curves and nodes */
	nc = nn = 0;
	surf_pos_curve_loop(canopy,c)
	{
	    if (!pointer_in_list(*c,nc,(POINTER*)curves))
	    {
		curves[nc++] = *c;
		if (!pointer_in_list((*c)->start,nn,(POINTER*)nodes))
		    nodes[nn++] = (*c)->start;
		if (!pointer_in_list((*c)->end,nn,(POINTER*)nodes))
		    nodes[nn++] = (*c)->end;
	    }
	}
	surf_neg_curve_loop(canopy,c)
	{
	    if (!pointer_in_list(*c,nc,(POINTER*)curves))
	    {
		curves[nc++] = *c;
		if (!pointer_in_list((*c)->start,nn,(POINTER*)nodes))
		    nodes[nn++] = (*c)->start;
		if (!pointer_in_list((*c)->end,nn,(POINTER*)nodes))
		    nodes[nn++] = (*c)->end;
	    }
	}
	for (i = 0; i < nn; ++i)
	{
	    node_in_curve_loop(nodes[i],c)
	    {
		if (!pointer_in_list(*c,nc,(POINTER*)curves))
		{
		    curves[nc++] = *c;
		    if (!pointer_in_list((*c)->start,nn,(POINTER*)nodes))
		    	nodes[nn++] = (*c)->start;
		    if (!pointer_in_list((*c)->end,nn,(POINTER*)nodes))
		    	nodes[nn++] = (*c)->end;
		}
	    }
	    node_out_curve_loop(nodes[i],c)
	    {
		if (!pointer_in_list(*c,nc,(POINTER*)curves))
		{
		    curves[nc++] = *c;
		    if (!pointer_in_list((*c)->start,nn,(POINTER*)nodes))
		    	nodes[nn++] = (*c)->start;
		    if (!pointer_in_list((*c)->end,nn,(POINTER*)nodes))
		    	nodes[nn++] = (*c)->end;
		}
	    }
	}
	/* Systematically rotate all points */
	I_SphericalRotateInteriorSurfPoints(canopy,center,phi,theta);
	for (i = 0; i < nc; ++i)
	    I_SphericalRotateInteriorCurvePoints(curves[i],center,phi,theta);
	for (i = 0; i < nn; ++i)
	    I_SphericalRotatePoint(nodes[i]->posn,center,phi,theta,NO);
}	/* end modifyCanopySet */
