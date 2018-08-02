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
typedef struct {
        int dim;
        double cen[MAXD];
	double ang[MAXD];
        double rad;
} ARC_PARAMS;

typedef struct {
        /*equation for star shaped function*/
        /*x = a[0]*cos(t)+x0*/
        /*y = b[1]*sin(t)+y0*/
        double a[2];       
        double ang[2];
        double cen[MAXD];
} AIRBAG_PARAMS;

static void initSingleModule(Front*);
static void initMultiModule(Front*,int);
static void MergeTwoIntfc(INTERFACE*,INTERFACE*);
static void CopyNodeInfo(INTERFACE*,INTERFACE*);
static void modifyCanopySet(FILE*,Front*,SURFACE*);
static boolean curve_of_boundary_hs(CURVE*);
static boolean line_seg_func(POINTER,double,double*);
static boolean parab_seg_func(POINTER,double,double*);
static boolean arc_func(POINTER,double,double*);
static void InstallLoadNode2d(Front*,double*,boolean);
static void initAirBag(Front*,FILE*);
static void initParachute(Front*,FILE*);
static void airbag_shaped_func(POINTER,double,double*);


extern void init2DModules(Front *front)
{
	int i,num_canopy;
	FILE *infile = fopen(InName(front),"r");
	SURFACE *surf;
	RG_PARAMS rgb_params;

	if (debugging("trace"))
	    (void) printf("Entering init2DModules()\n");

	if (debugging("set_module"))
	    gview_plot_interface("module-step-1",front->interf);

	num_canopy = 0;
	if (CursorAfterStringOpt(infile,"Enter number of canopy curve:"))
	{
            fscanf(infile,"%d",&num_canopy);
	    fclose(infile);
            (void) printf("%d\n",num_canopy);
	}

	if (num_canopy == 1)
	    initSingleModule(front);
	else if (num_canopy > 1)
	    initMultiModule(front,num_canopy);

	initRigidBody(front);
	rgb_init(front,rgb_params);

	if (debugging("trace"))
	    (void) printf("Leaving init2DModules()\n");
}	/* end init2DModules */

extern void initParachuteDefault(
	Front *front)
{
	AF_PARAMS *af_params = (AF_PARAMS*)front->extra2;
	af_params->is_parachute_system = YES;
	af_params->num_opt_round = 20;
        af_params->spring_model = MODEL1;
	af_params->gore_len_fac = 1.0;
}	/* end initParachuteDefault */

static void initSingleModule(
        Front *front)
{
	FILE *infile = fopen(InName(front),"r");
	SURFACE *surf;
	char string[200];

	(void) printf("Available type of 2d problem include:\n");
	(void) printf("\tParachute (P)\n");
	(void) printf("\tAirbag (A)\n");

	CursorAfterString(infile,"Enter type of parachute:");
	fscanf(infile,"%s",string);
	(void) printf("%s\n",string);
	switch (string[0])
	{
	case 'p':
	case 'P':
	    initParachute(front,infile);    
	    break;
	case 'A':
	case 'a':
	    initAirBag(front,infile);
	    break;
	default:
	    (void) printf("Unknow type of parachute %s!\n",string);
	    clean_up(ERROR);
	}
	fclose(infile);
	if (debugging("init_module"))
	{
	    gview_plot_interface("init_module",front->interf);
	    CURVE** c;
	    int i = 0;
	    char fname[200];
            intfc_curve_loop(front->interf,c)
	    {
			printf("curve[%d]: [%f %f] --> [%f %f]\n",
                        i, Coords(((*c)->start)->posn)[0],
                           Coords(((*c)->start)->posn)[1],
                           Coords(((*c)->end)->posn)[0],
                           Coords(((*c)->end)->posn)[1]);
			sprintf(fname,"curve-%d",i++);
		 gview_plot_curve(*c,"curves",fname,SURFACE_COLOR(i%8),2);
            }
	}
}	/* end initSingleModule */

static void initParachute(Front* front,FILE* infile)
{
        NODE* nload;
        int w_type;
        char string[200];
        double len,dist,load_pos[MAXD];
        int i, num, k, dim = FT_Dimension();
        int neg_comp,pos_comp;
        CURVE *curve,*string_curve;
        LINE_SEG_PARAMS line_params;
        boolean fixed = NO;
	AF_NODE_EXTRA* extra;
	LEVEL_FUNC_PACK level_func_pack;

        w_type = ELASTIC_BOUNDARY;	
	CursorAfterString(infile,"Enter start coordinate:");
	fscanf(infile,"%lf %lf",line_params.coords_start,
				line_params.coords_start+1);
	(void) printf("%f %f\n",line_params.coords_start[0],
				line_params.coords_start[1]);
	CursorAfterString(infile,"Enter end coordinate:");
	fscanf(infile,"%lf %lf",line_params.coords_end,
				line_params.coords_end+1);
	(void) printf("%f %f\n",line_params.coords_end[0],
				line_params.coords_end[1]);
	CursorAfterString(infile,"Enter load coordinate:");
        fscanf(infile,"%lf %lf",&load_pos[0],&load_pos[1]);
        (void) printf("%f %f\n",load_pos[0],load_pos[1]);
	CursorAfterString(infile,"Enter shape of parachute canopy:");
        fscanf(infile,"%s",string);
        (void) printf("%s\n",string);
	neg_comp = LIQUID_COMP2;
	pos_comp = LIQUID_COMP2;
	line_params.dim = 2;
	/*make 2d canopy*/
	switch (string[0])
	{
	    case 'L':
	    case 'l':
	        curve = FT_MakeParametricCurve(front,neg_comp,pos_comp,w_type,
				line_seg_func,(POINTER)&line_params,1,NO);
		break;
	    case 'P':
	    case 'p':
		curve = FT_MakeParametricCurve(front,neg_comp,pos_comp,w_type,
                                parab_seg_func,(POINTER)&line_params,1,NO);
		break;
	    default:
		printf("Warning: unknown parachute shape %s\n",string);
		printf("Avaiable parachute shapes are: LINE, PARABOLIC\n");
	}
    	node_type(curve->start) = node_type(curve->end) = MONO_STRING_NODE;
    	FT_ScalarMemoryAlloc((POINTER*)&extra,sizeof(AF_NODE_EXTRA));
    	extra->af_node_type = STRING_NODE;
   	curve->start->extra = (POINTER)extra;
    	curve->start->size_of_extra = sizeof(AF_NODE_EXTRA);
    	curve->end->extra = (POINTER)extra;
    	curve->end->size_of_extra = sizeof(AF_NODE_EXTRA);
    	//hsbdry_type(curve) = MONO_COMP_HSBDRY;

	/*make curve between load node and 2d canopy*/
	CursorAfterString(infile,"Enter yes to fix load node:");
        fscanf(infile,"%s",string);
        (void) printf("%s\n",string);
	if (string[0] == 'Y' || string[0] == 'y')
	    fixed = YES;
	else 
	    fixed = NO;
	InstallLoadNode2d(front,load_pos,fixed);
   	set_equilibrium_mesh(front); /*very important*/
}

static void initAirBag(Front* front,FILE* infile)
{
	CURVE *curve;
	CURVE **c;
	CURVE **cc, **delete_curves;
	INTERFACE* intfc;
	COMPONENT neg_comp, pos_comp;
	AIRBAG_PARAMS airbag_params;
	int w_type = ELASTIC_BOUNDARY;
        AF_NODE_EXTRA* extra;
	int *gmax = front->rect_grid->gmax;
	neg_comp = LIQUID_COMP2; pos_comp = SOLID_COMP;
	double **node_coords;
	RECT_GRID* rgr = front->rect_grid;
	int np, num_nodes = 2;

	/*delete boundary curves and reset*/
	delete_curves = NULL;
        for (cc = front->interf->curves; cc && *cc; ++cc)
            add_to_pointers(*cc,&delete_curves);
        for (cc = delete_curves; cc && *cc; ++cc)
            (void) delete_curve(*cc);

	/*read params*/
	CursorAfterString(infile,"Enter airbag center:");
        fscanf(infile,"%lf %lf",airbag_params.cen,
                                airbag_params.cen+1);
	(void) printf("%f %f\n",airbag_params.cen[0],
                                airbag_params.cen[1]);
	CursorAfterString(infile,"Enter airbag radius:");
            fscanf(infile,"%lf",&airbag_params.a[0]);
	(void) printf("%f\n",airbag_params.a[0]);
	airbag_params.a[1] = airbag_params.a[0]/4.0; 
	airbag_params.ang[1] = -4.0/9.0*M_PI;
	airbag_params.ang[0] = 13.0/9.0*M_PI;
	airbag_params.cen[1] = rgr->L[1]-airbag_params.a[1]
			     * sin(airbag_params.ang[1]);
	LEVEL_FUNC_PACK level_func_pack;
	level_func_pack.func_params = NULL;
	level_func_pack.func = NULL;
	level_func_pack.closed_curve = NO;
	level_func_pack.neg_component = LIQUID_COMP2;
        level_func_pack.pos_component = LIQUID_COMP2;
	level_func_pack.wave_type = ELASTIC_BOUNDARY;

	level_func_pack.num_points = np = (int)gmax[0];
	FT_MatrixMemoryAlloc((POINTER*)&level_func_pack.point_array,np,
                                2,sizeof(double));
	double t = airbag_params.ang[0];
	double dt = (airbag_params.ang[1]-airbag_params.ang[0])/(np-1);
	for (int i = 0; i < np; i++)
	{
	    airbag_shaped_func((POINTER)&airbag_params,t,
				level_func_pack.point_array[i]);
	    t += dt; 
	}
	level_func_pack.point_array[np-1][1] = front->rect_grid->L[1];
	level_func_pack.point_array[0][1] = front->rect_grid->L[1];
	FT_InitIntfc(front,&level_func_pack);
	set_equilibrium_mesh(front);
	intfc_curve_loop(front->interf,c)
	{
	    if(wave_type(*c) == ELASTIC_BOUNDARY)
	    {
		node_type((*c)->start) = node_type((*c)->end) = FIXED_NODE;
    		FT_ScalarMemoryAlloc((POINTER*)&extra,sizeof(AF_NODE_EXTRA));
    		extra->af_node_type = PRESET_NODE;
   		(*c)->start->extra = (POINTER)extra;
    		(*c)->start->size_of_extra = sizeof(AF_NODE_EXTRA);
    		(*c)->end->extra = (POINTER)extra;
    		(*c)->end->size_of_extra = sizeof(AF_NODE_EXTRA);
	    }
	}
}	/* end initAirBag */

static void airbag_shaped_func(
	POINTER params,
	double t,
	double *coords)
{
	/*star shaped level function, pos: outter, negative: inner*/
	AIRBAG_PARAMS *airbag_params = (AIRBAG_PARAMS*)params;
	double *a = airbag_params->a;
	double *cen = airbag_params->cen;
	double theta;
	theta = t;
	coords[0] = a[0] * cos(theta) + cen[0];
	coords[1] = a[1] * sin(theta) + cen[1];
	return;
}

static void InstallLoadNode2d(
	Front *front,
	double *load_pos,
	boolean fixed)
{
	INTERFACE *intfc = front->interf;
        NODE **n, *nload;
        CURVE **string_curves;
        AF_NODE_EXTRA *extra;
        BOND *bond;
        double center[MAXD],dir[MAXD],coords[MAXD];
        double spacing,*h = front->rect_grid->h;
        int i,j,k,nb;
        INTERFACE *cur_intfc;
	int dim = FT_Dimension();

	/*make load node*/
	FT_VectorMemoryAlloc((POINTER*)&string_curves,2,
                                sizeof(CURVE*));
        nload = make_node(Point(load_pos));
        FT_ScalarMemoryAlloc((POINTER*)&extra,sizeof(AF_NODE_EXTRA));
	if (fixed)
	{
            extra->af_node_type = PRESET_NODE;
	    node_type(nload) = FIXED_NODE;
	}
	else
	{
            extra->af_node_type = LOAD_NODE;
	    node_type(nload) = MONO_STRING_NODE;
	}
        nload->extra = (POINTER)extra;
        nload->size_of_extra = sizeof(AF_NODE_EXTRA);

	i = 0;
	/*link intfc and load node*/
	intfc_node_loop(intfc,n)
        {
            extra = (AF_NODE_EXTRA*)((*n)->extra);
            if (extra == NULL)
                continue;
            if (extra->af_node_type == STRING_NODE)
            {
                string_curves[i] = make_curve(LIQUID_COMP2,LIQUID_COMP2,
						(*n),nload);
		wave_type(string_curves[i]) = ELASTIC_STRING;
                spacing = separation((*n)->posn,nload->posn,dim);
                for (j = 0; j < dim; ++j)
                    dir[j] = (Coords(nload->posn)[j] -
                        Coords((*n)->posn)[j])/spacing;
                nb = (int)(spacing/h[0]);
                spacing /= (double)nb;
                bond = string_curves[i]->first;
                for (j = 1; j < nb; ++j)
                {
                    for (k = 0; k < 3; ++k)
                        coords[k] = Coords((*n)->posn)[k] +
                                        j*dir[k]*spacing;
                    insert_point_in_bond(Point(coords),bond,string_curves[i]);
                    bond->length0 = spacing;
                    bond = bond->next;
                }
                bond->length0 = spacing;
                i++;
            }
        }
 
}


static void initMultiModule(
	Front *front,
	int num_canopy)
{
	INTERFACE *intfc = front->interf;
	SURFACE **surfs;
	FILE *infile = fopen(InName(front),"r");
	double center[MAXD];
	double phi,theta;
	char string[10];
	int i;
        INTERFACE *cur_intfc;
        cur_intfc = current_interface();

	/*
	FT_VectorMemoryAlloc((POINTER*)&surfs,num_canopy,sizeof(SURFACE*));

	for (i = 0; i < num_canopy; ++i)
	{
	    CgalCanopySurface(infile,front,&surfs[i]);
	    modifyCanopySet(infile,front,surfs[i]);
	}

	InstallNewLoadNode(front,num_canopy);
	set_current_interface(cur_intfc);
	FT_FreeThese(1,surfs);
	*/
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
            fscanf(infile,"%lf %lf",center,center+1);
            (void) printf("%f %f\n",center[0],center[1]);
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

extern void initRigidBody(
	Front *front)
{
	FILE *infile = fopen(InName(front),"r");
	char string[100];
	double cen[MAXD];
	double radius,radii[MAXD];
	double len,dist;
	int w_type, x_type;
	int i, num, k, dim = FT_Dimension();
	int neg_comp,pos_comp;
	CURVE *curve;
	LINE_SEG_PARAMS line_params;
	ARC_PARAMS arc_params;

	if (CursorAfterStringOpt(infile,"Enter yes to add rigid body:"))
	{
	    fscanf(infile,"%s",string);
	    (void) printf("%s\n",string);
	    if (string[0] != 'y' && string[0] != 'Y')
		return;
	}
	else
	    return;

	w_type = MOVABLE_BODY_BOUNDARY;
	(void) printf("Rigid body can be fixed (F) or Movable (M)\n");
	(void) printf("The default is Movable (M)\n");
	if (CursorAfterStringOpt(infile,"Type yes if the rigid body is fixed:"))
	{
	    fscanf(infile,"%s",string);
	    (void) printf("%s\n",string);
	    if (string[0] == 'y' || string[0] == 'Y')
		w_type = NEUMANN_BOUNDARY;
	}

	(void) printf("Available type of rigid body include:\n");
	(void) printf("\tCircle        (c)\n");
	(void) printf("\tLine segment  (l)\n");
	(void) printf("\tLine segments (n)\n");
	(void) printf("\tCircular arc  (a)\n");
	(void) printf("\tTwo arcs      (t)\n");
	(void) printf("\tTwo segments  (s)\n");

	CursorAfterString(infile,"Enter type of rigid body:");
	fscanf(infile,"%s",string);
	(void) printf("%s\n",string);
	switch (string[0])
	{
	case 'c':
	case 'C':
	    CursorAfterString(infile,"Enter center of the circle:");
	    fscanf(infile,"%lf %lf",cen,cen+1);
	    (void) printf("%f %f\n",cen[0],cen[1]);
	    CursorAfterString(infile,"Enter radius of the circle:");
	    fscanf(infile,"%lf",&radius);
	    (void) printf("%f\n",radius);
	    for (i = 0; i < dim; ++i) radii[i] = radius;
	    neg_comp = SOLID_COMP;
	    pos_comp = LIQUID_COMP2;
	    FT_MakeEllipticCurve(front,cen,radii,neg_comp,pos_comp,w_type,
				2.0,&curve);
	    node_type(curve->start) = CLOSED_NODE;
	    break;
	case 'L':
	case 'l':
	    CursorAfterString(infile,"Enter start coordinate:");
	    fscanf(infile,"%lf %lf",line_params.coords_start,
				line_params.coords_start+1);
	    (void) printf("%f %f\n",line_params.coords_start[0],
				line_params.coords_start[1]);
	    CursorAfterString(infile,"Enter end coordinate:");
	    fscanf(infile,"%lf %lf",line_params.coords_end,
				line_params.coords_end+1);
	    (void) printf("%f %f\n",line_params.coords_end[0],
				line_params.coords_end[1]);
	    neg_comp = LIQUID_COMP2;
	    pos_comp = LIQUID_COMP2;
	    line_params.dim = 2;
	    curve = FT_MakeParametricCurve(front,neg_comp,pos_comp,w_type,
				line_seg_func,(POINTER)&line_params,2,NO);
	    node_type(curve->start) = node_type(curve->end) = FIXED_NODE;
	    break;
	case 'N':
        case 'n':
            CursorAfterString(infile,"Enter number of segments:");
            fscanf(infile,"%d",&num);
            (void) printf("%d\n",num);
            CursorAfterString(infile,"Enter length of segments:");
            fscanf(infile,"%lf",&len);
            (void) printf("%f\n",len);
            CursorAfterString(infile,"Enter distance between segments:");
            fscanf(infile,"%lf",&dist);
	    (void) printf("%f\n",dist);
            CursorAfterString(infile,
				"Enter start coordinates for first segment:");
            fscanf(infile,"%lf %lf",line_params.coords_start,
                                line_params.coords_start+1);
            (void) printf("%f %f\n",line_params.coords_start[0],
                                line_params.coords_start[1]);
	    for (k = 0; k < num; k++)
	    {
	    	line_params.coords_end[0] =  line_params.coords_start[0] + len;
	    	line_params.coords_end[1] = line_params.coords_start[1];
            	neg_comp = LIQUID_COMP2;
            	pos_comp = LIQUID_COMP2;
            	line_params.dim = 2;
            	curve = FT_MakeParametricCurve(front,neg_comp,pos_comp,w_type,
                                line_seg_func,(POINTER)&line_params,2,NO);
            	node_type(curve->start) = node_type(curve->end) = FIXED_NODE;
	    	line_params.coords_start[0] = line_params.coords_end[0] + dist;            
            }
	    break;
	case 'A':
        case 'a':
            CursorAfterString(infile,"Enter center coordinates:");
            fscanf(infile,"%lf %lf",arc_params.cen,
                                arc_params.cen+1);
            (void) printf("%f %f\n",arc_params.cen[0],
                                arc_params.cen[1]);
            CursorAfterString(infile,"Enter radius:");
            fscanf(infile,"%lf",&(arc_params.rad));
            (void) printf("%f\n",arc_params.rad);
	    CursorAfterString(infile,"Enter start end angles in degrees:");
            fscanf(infile,"%lf %lf",arc_params.ang,
                                arc_params.ang+1);
            (void) printf("%f %f\n",arc_params.ang[0],
                                arc_params.ang[1]);

            neg_comp = LIQUID_COMP2;
            pos_comp = LIQUID_COMP2;
            arc_params.dim = 3;
            curve = FT_MakeParametricCurve(front,neg_comp,pos_comp,w_type,
                                arc_func,(POINTER)&arc_params,3,NO);
            node_type(curve->start) = node_type(curve->end) = FIXED_NODE;
            break;
	case 'S':
        case 's':
            CursorAfterString(infile,
				"Enter start coordinates of first segment:");
            fscanf(infile,"%lf %lf",line_params.coords_start,
				line_params.coords_start+1);
            (void) printf("%f %f\n",line_params.coords_start[0],
                                line_params.coords_start[1]);
            CursorAfterString(infile,"Enter end coordinates of first segment:");
            fscanf(infile,"%lf %lf",line_params.coords_end,
                                line_params.coords_end+1);
            (void) printf("%f %f\n",line_params.coords_end[0],
                                line_params.coords_end[1]);
            neg_comp = LIQUID_COMP2;
            pos_comp = LIQUID_COMP2;
            line_params.dim = 2;
            curve = FT_MakeParametricCurve(front,neg_comp,pos_comp,w_type,
                                line_seg_func,(POINTER)&line_params,2,NO);
            node_type(curve->start) = node_type(curve->end) = FIXED_NODE;

            CursorAfterString(infile,
				"Enter start coordinates of second segment:");
            fscanf(infile,"%lf %lf",line_params.coords_start,
                                line_params.coords_start+1);
            (void) printf("%f %f\n",line_params.coords_start[0],
                                line_params.coords_start[1]);
            CursorAfterString(infile,
				"Enter end coordinates of second segment:");
            fscanf(infile,"%lf %lf",line_params.coords_end,
                                line_params.coords_end+1);
            (void) printf("%f %f\n",line_params.coords_end[0],
                                line_params.coords_end[1]);
            neg_comp = LIQUID_COMP2;
            pos_comp = LIQUID_COMP2;
            line_params.dim = 2;
            curve = FT_MakeParametricCurve(front,neg_comp,pos_comp,w_type,
                                line_seg_func,(POINTER)&line_params,2,NO);
            node_type(curve->start) = node_type(curve->end) = FIXED_NODE;
            break;
	case 'T':
        case 't':
            CursorAfterString(infile,"Enter center coordinates for first arc:");
            fscanf(infile,"%lf %lf",arc_params.cen,
                                arc_params.cen+1);
            (void) printf("%f %f\n",arc_params.cen[0],
                                arc_params.cen[1]);
            CursorAfterString(infile,"Enter radius for first arc:");
            fscanf(infile,"%lf",&(arc_params.rad));
            (void) printf("%f\n",arc_params.rad);
            CursorAfterString(infile,"Enter start end angles for first arc:");
            fscanf(infile,"%lf %lf",arc_params.ang,
                                arc_params.ang+1);
            (void) printf("%f %f\n",arc_params.ang[0],
                                arc_params.ang[1]);
	    
	    neg_comp = LIQUID_COMP2;
            pos_comp = LIQUID_COMP2;
            arc_params.dim = 3;
            curve = FT_MakeParametricCurve(front,neg_comp,pos_comp,w_type,
                                arc_func,(POINTER)&arc_params,3,NO);
            node_type(curve->start) = node_type(curve->end) = FIXED_NODE;
		
	    CursorAfterString(infile,
				"Enter center coordinates for second arc:");
            fscanf(infile,"%lf %lf",arc_params.cen,
                                arc_params.cen+1);
            (void) printf("%f %f\n",arc_params.cen[0],
                                arc_params.cen[1]);
	    CursorAfterString(infile,"Enter radius for second arc:");
            fscanf(infile,"%lf",&(arc_params.rad));
            (void) printf("%f\n",arc_params.rad);
            CursorAfterString(infile,"Enter start end angles for second arc:");
            fscanf(infile,"%lf %lf",arc_params.ang,
                                arc_params.ang+1);
            (void) printf("%f %f\n",arc_params.ang[0],
                                arc_params.ang[1]);

            neg_comp = LIQUID_COMP2;
            pos_comp = LIQUID_COMP2;
            arc_params.dim = 3;
            curve = FT_MakeParametricCurve(front,neg_comp,pos_comp,w_type,
                                arc_func,(POINTER)&arc_params,3,NO);
            node_type(curve->start) = node_type(curve->end) = FIXED_NODE;
            break;

	default:
	    (void) printf("Unknow type of rigid body!\n");
	    clean_up(ERROR);
	}

	fclose(infile);
}	/* end initRigidBody */

static boolean line_seg_func(
	POINTER params,
	double t,
	double *coords)
{
	LINE_SEG_PARAMS *l_params = (LINE_SEG_PARAMS*)params;
	double *coords_start = l_params->coords_start;
	double *coords_end = l_params->coords_end;
	int i,dim = l_params->dim;
	for (i = 0; i < dim; ++i)
	{
	    if (coords_end[i] == coords_start[i])
		coords[i] = coords_start[i];
	    else
	    	coords[i] = coords_start[i] + t*(coords_end[i] - 
				coords_start[i]);
	}
	return YES;
}	/* end line_seg_func */


static boolean arc_func(
        POINTER params,
        double t,
        double *coords)
{
        ARC_PARAMS *l_params = (ARC_PARAMS*)params;
        double *cen = l_params->cen;
        double ang[2];
        double rad = l_params->rad;
        int i,dim = l_params->dim;
	
	ang[0] = PI*l_params->ang[0]/180;
	ang[1] = PI*l_params->ang[1]/180;
        coords[0] = cen[0] + rad*cos(ang[1] - ang[1]*t*(ang[1]-ang[0])/ang[1]);
        coords[1] = cen[1] + rad*sin(ang[1] - ang[1]*t*(ang[1]-ang[0])/ang[1]);
        return YES;
}       /* end arc_func */

extern void initInnerBoundary(
	Front *front,
	LEVEL_FUNC_PACK *level_func_pack)
{
	char string[100];
	FILE *infile = fopen(InName(front),"r");
	double cen[MAXD],R;
	static CIRCLE_PARAMS circle_params;

	if (!CursorAfterStringOpt(infile,"Enter yes to add inner boundary: "))
	    return;
	fscanf(infile,"%s",string);
	(void) printf("%s\n",string);
	if (string[0] != 'y' && string[0] != 'Y')
	    return;

	level_func_pack->pos_component = SOLID_COMP;
	level_func_pack->neg_component = LIQUID_COMP2;
	(void) printf("Available inner boundary types are\n");
	(void) printf("\tCircle (c)\n");
	CursorAfterString(infile,"Entering inner boundary type: ");
	fscanf(infile,"%s",string);
	(void) printf("%s\n",string);
	switch(string[0])
	{
	case 'c':
	case 'C':
	    circle_params.dim = 2;
	    CursorAfterString(infile,"Entering center of the circle: ");
	    fscanf(infile,"%lf %lf",circle_params.cen,circle_params.cen+1);
	    (void) printf("%f %f\n",circle_params.cen[0],circle_params.cen[1]);
	    CursorAfterString(infile,"Entering radius of the circle: ");
	    fscanf(infile,"%lf",&circle_params.R);
	    (void) printf("%f\n",circle_params.R);
	    level_func_pack->func_params = (POINTER)&circle_params;
	    level_func_pack->func = level_circle_func;
	    level_func_pack->wave_type = NEUMANN_BOUNDARY;
	    break;
	default:
	    (void) printf("Unknown inner boundary type\n");
	    clean_up(ERROR);
	}
}	/* end initInnerBoundary */


static boolean parab_seg_func(
	POINTER params,
	double t,
	double *coords)
{
	LINE_SEG_PARAMS *l_params = (LINE_SEG_PARAMS*)params;
	double *coords_start = l_params->coords_start;
	double *coords_end = l_params->coords_end;
	double coords_top[2], coords_cen[2], nor[2], vec[2];
	double ratio = 2.2, mag_vec = 0.0;
	double A[3][3], b[3];
	static double a[3];
	static boolean first = YES;
	if (first)
	{
	    first = NO;
	    /*rotate and translate two points to cartesian coords*/
	    coords_cen[0] = (coords_start[0] + coords_end[0])/2.0;
	    coords_cen[1] = (coords_start[1] + coords_end[1])/2.0;
	    vec[0] =  coords_end[0]-coords_start[0];
	    vec[1] =  coords_end[1]-coords_start[1];
	    nor[0] = -vec[1]; nor[1] = vec[0];
	    mag_vec = Mag2d(nor);
	    nor[0] /= mag_vec; nor[1] /= mag_vec;
	    mag_vec = Mag2d(vec)*ratio;
	    coords_top[0] = coords_cen[0] + mag_vec * nor[0];
	    coords_top[1] = coords_cen[1] + mag_vec * nor[1];
	
	    for (int j = 0; j < 3; j++)
	    {
	        A[j][0] = pow(coords_start[0],2-j); 
	        A[j][1] = pow(coords_end[0],2-j); 
	        A[j][2] = pow(coords_top[0],2-j); 
	    }
	    b[0] = coords_start[1]; 
	    b[1] = coords_end[1];
	    b[2] = coords_top[1];
	    a[0] = Det3d(b,A[1],A[2])/Det3d(A[0],A[1],A[2]);
	    a[1] = Det3d(A[0],b,A[2])/Det3d(A[0],A[1],A[2]);
	    a[2] = Det3d(A[0],A[1],b)/Det3d(A[0],A[1],A[2]);
	}
	coords[0] = coords_start[0] + t*(coords_end[0] - 
			 coords_start[0]);
	coords[1] = a[0]*coords[0]*coords[0]+a[1]*coords[0]+a[2];
	
	return YES;
}	/* end parab_seg_func */

