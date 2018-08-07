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

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Delaunay_mesher_2.h>
#include <CGAL/Delaunay_mesh_face_base_2.h>
#include <CGAL/Delaunay_mesh_size_criteria_2.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>

#include <iostream>
#include <fstream>
#include <iFluid.h>
#include <airfoil.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Triangulation_vertex_base_with_info_2<unsigned, K> Vb;
typedef CGAL::Delaunay_mesh_face_base_2<K> Fb;
typedef CGAL::Triangulation_data_structure_2<Vb, Fb> Tds;
typedef CGAL::Exact_predicates_tag Itag;
typedef CGAL::Constrained_Delaunay_triangulation_2<K, Tds,Itag> CDT;
typedef CGAL::Delaunay_mesh_size_criteria_2<CDT> Criteria;

typedef CDT::Vertex_handle Vertex_handle;
typedef CDT::Point Cgal_Point;

static void CgalFlatSurface(FILE*,Front*,SURFACE**);
static void CgalParabolicSurface(FILE*,Front*,SURFACE**);
static void CgalCross(FILE*,Front*,SURFACE**);
static void CgalCircle(FILE*,Front*,SURFACE**);
static void CgalRectangular(FILE*,Front*,SURFACE**);
static void CgalEllipse(FILE*,Front*,SURFACE**);
static void GenerateCgalSurf(Front*,SURFACE**,CDT*,int*,double);
static void InstallGore(Front*,SURFACE*,int,double*,double*);
static void InstallInCurve(Front*,SURFACE*,double,double,double*,int);
static void SplitCirBdry(Front*,SURFACE*,int,double*,ORIENTATION);
static void linkCurveTriBond(CURVE*,SURFACE*);
static bool ptinbox(double *c, double *l, double *u);
static bool ptoutcircle(double*,double*,double);
static void foldSurface(FILE*,SURFACE*);
static void findStringNodePoints(SURFACE*,double*,POINT**,int,CURVE**);
static void installString(FILE*,Front*,SURFACE*,CURVE*,POINT**,int);
static void resetStringNodePoints(SURFACE*,POINT**,int*,CURVE**);
static void setCurveZeroLength(CURVE*,double);
static void setSurfZeroMesh(SURFACE*);
static void resetGoreBdryZerolength(SURFACE*);
static void setMonoCompBdryZeroLength(SURFACE*);
static boolean sewSurface(FILE*,SURFACE*);

static void connectStringtoRGB(Front*,SURFACE*,NODE**,int);
static void linkSurfaceTriPoint(INTERFACE*,SURFACE*);
static void findPointsonRGB(Front*,SURFACE*,std::vector<POINT*>&);
static void vector_extreme(std::vector<POINT*>&,int,std::vector<POINT*>&);
static void checkReducedTri(SURFACE*);

extern void CgalCanopySurface(
	FILE *infile,
	Front *front,
	SURFACE **surf)
{
	char string[10];
        (void) printf("Available canopy surface types are:\n");
        (void) printf("\tFLAT (F)\n");
        (void) printf("\tPARABOLIC (P)\n");
        (void) printf("\tELLIPTIC (E)\n");
        CursorAfterString(infile,"Enter canopy surface type:");
        fscanf(infile,"%s",string);
        (void) printf("%s\n",string);
        switch (string[0])
	{
	    case 'F':
	    case 'f':
		CgalFlatSurface(infile,front,surf);
		break;
	    case 'P':
	    case 'p':
		CgalParabolicSurface(infile,front,surf);
            break;
	    default:
		(void) printf("Unknown canopy surface type\n");
		clean_up(ERROR);
	}
}	/* end CgalCanopySurface */

static void CgalParabolicSurface(
	FILE *infile,
	Front *front,
	SURFACE **surf)
{
	double cen[3], rad[2];
	TRI *tri;
	POINT *p;
	int i;

	CgalCircle(infile,front,surf);

        CursorAfterString(infile,"Enter vertex coordinate of the paraboloid:");
        fscanf(infile,"%lf %lf %lf",&cen[0],&cen[1],&cen[2]);
        (void) printf("%f %f %f\n",cen[0],cen[1],cen[2]);
        CursorAfterString(infile,"Enter coefficients of the paraboloid:");
        fscanf(infile,"%lf %lf",&rad[0],&rad[1]);
        (void) printf("%f %f\n",rad[0],rad[1]);

	unsort_surf_point(*surf);
	surf_tri_loop(*surf,tri)
	{
	    for (i = 0; i < 3; i++)
	    {
		tri->side_length0[i] = -1.0;
		p = Point_of_tri(tri)[i];
		if (sorted(p) == NO)
		{
		    Coords(p)[2] = cen[2] - rad[0]*sqr(Coords(p)[0]-cen[0])
					  - rad[1]*sqr(Coords(p)[1]-cen[1]); 
		    sorted(p) = YES;
		}
	    }
	}
	setSurfZeroMesh(*surf);
	resetGoreBdryZerolength(*surf);
}	/* end CgalParabolicSurface */

static void CgalFlatSurface(
	FILE *infile,
	Front *front,
	SURFACE **surf)
{
        char string[200];

        (void) printf("Available types of canopy boundaries are:\n");
        (void) printf("\tCircular (C)\n");
        (void) printf("\tRectangular (R)\n");
        (void) printf("\tElliptic (E)\n");
        (void) printf("\tCross (X)\n");
        (void) printf("\tWing (W)\n");
        CursorAfterString(infile,"Enter type of canopy boundary:");
        fscanf(infile,"%s",string);
        (void) printf("%s\n",string);
        switch (string[0])
        {
        case 'X':
        case 'x':
            CgalCross(infile,front,surf);
            break;
        case 'C':
        case 'c':
            CgalCircle(infile,front,surf);
            break;
        case 'E':
        case 'e':
            CgalEllipse(infile,front,surf);
            break;
        case 'R':
        case 'r':
            CgalRectangular(infile,front,surf);
            break;
        }
}	/* end CgalFlatSurface */

static void CgalRectangular(
	FILE *infile,
	Front *front,
	SURFACE **surf)
{
	int *flag;
	int i;
	int num_strings = 4;
	double height;
	double Lower[2];
	double Upper[2];
	double cri_dx = 0.6*computational_grid(front->interf)->h[0];
	std::list<Cgal_Point> list_of_seeds;
	CDT cdt;
	CDT::Finite_faces_iterator fit;
        Vertex_handle *v_out;

        CursorAfterString(infile,"Enter the height of the plane:");
        fscanf(infile,"%lf",&height);
        (void) printf("%f\n",height);
        CursorAfterString(infile,"Enter lower bounds of the rectangle:");
        fscanf(infile,"%lf %lf",Lower, Lower+1);
        (void) printf("%f %f\n",Lower[0], Lower[1]);
        CursorAfterString(infile,"Enter upper bounds of the rectangle:");
        fscanf(infile,"%lf %lf",Upper, Upper+1);
        (void) printf("%f %f\n",Upper[0], Upper[1]);

	v_out = new Vertex_handle[num_strings];

	v_out[0] = cdt.insert(Cgal_Point(Lower[0], Lower[1]));
	v_out[1] = cdt.insert(Cgal_Point(Upper[0], Lower[1]));
	v_out[2] = cdt.insert(Cgal_Point(Upper[0], Upper[1]));
	v_out[3] = cdt.insert(Cgal_Point(Lower[0], Upper[1]));

	for (i = 0; i < num_strings-1; i++)
		cdt.insert_constraint(v_out[i],v_out[i+1]);
	cdt.insert_constraint(v_out[0],v_out[num_strings-1]);

	CGAL::refine_Delaunay_mesh_2(cdt, list_of_seeds.begin(), 
			list_of_seeds.end(),Criteria(0.3, cri_dx));

	flag = new int[cdt.number_of_faces()];

	i = 0;
        for (fit = cdt.finite_faces_begin(); fit != cdt.finite_faces_end();
		++fit)
	    flag[i++] = 1;	

	GenerateCgalSurf(front,surf,&cdt,flag,height);
        wave_type(*surf) = ELASTIC_BOUNDARY;
        FT_InstallSurfEdge(*surf,MONO_COMP_HSBDRY);
	setSurfZeroMesh(*surf);
	setMonoCompBdryZeroLength(*surf);
	if (consistent_interface(front->interf) == NO)
	    clean_up(ERROR);
}	/* end CgalRectangular */

static void CgalCircle(
	FILE *infile,
	Front *front,
	SURFACE **surf)
{
	double height;
	double CirCenter[2];
	double CirR[2];
	int num_strings;
	CDT cdt;
	CDT::Finite_faces_iterator fit;
        Vertex_handle *v_out, *v_in;
        int num_out_vtx,num_in_vtx;
	POINT **string_node_pts;
        double *out_nodes_coords,*in_nodes_coords;
	double *out_vtx_coords,*in_vtx_coords;
	double ang_out, ang_in;
	int out_vtx_oneside = 15, in_vtx_oneside = 2;
	char gore_bool[10],vent_bool[10], string_bool[10];
	std::list<Cgal_Point> list_of_seeds;
	double cri_dx = 0.6*computational_grid(front->interf)->h[0];
	AF_PARAMS *af_params = (AF_PARAMS*)front->extra2;
	int i;
	CURVE *cbdry;

        CursorAfterString(infile,"Enter the height of the plane:");
        fscanf(infile,"%lf",&height);
        (void) printf("%f\n",height);
        CursorAfterString(infile,"Enter circle center:");
        fscanf(infile,"%lf %lf",&CirCenter[0],&CirCenter[1]);
        (void) printf("%f %f\n",CirCenter[0],CirCenter[1]);
        CursorAfterString(infile,"Enter circle radius:");
        fscanf(infile,"%lf",&CirR[0]);
        (void) printf("%f\n",CirR[0]);

	CirR[1] = 0;
	CursorAfterStringOpt(infile,"Enter yes to attach gores to canopy:");
        fscanf(infile,"%s",gore_bool);
        (void) printf("%s\n",gore_bool);
        if (gore_bool[0]=='y' || gore_bool[0]=='Y')
        {
	    CirR[1] = 0.1 * CirR[0];
	    af_params->attach_gores = YES;
	}
	else
	    af_params->attach_gores = NO;
	CursorAfterStringOpt(infile,"Enter yes to cut a vent on canopy:");
	fscanf(infile,"%s",vent_bool);
	(void) printf("%s\n",vent_bool);
	if (vent_bool[0]=='y' || vent_bool[0]=='Y')
        {
            CursorAfterString(infile,"Enter radius of the vent:");
	    fscanf(infile,"%lf",&CirR[1]);
	    (void) printf("%f\n",CirR[1]);
        }

	num_strings = 28;   //default
	CursorAfterStringOpt(infile,"Enter yes to attach strings to canopy:");
	fscanf(infile,"%s",string_bool);
	(void) printf("%s\n",string_bool);
	if (string_bool[0]=='y' || string_bool[0]=='Y')
	{
	    CursorAfterString(infile,"Enter number of chords:");
	    fscanf(infile,"%d",&num_strings);
	    (void) printf("%d\n",num_strings);
	}
	FT_VectorMemoryAlloc((POINTER*)&string_node_pts,num_strings,
                                sizeof(POINT*));

	double offset = PI/num_strings;
	num_out_vtx = num_strings * out_vtx_oneside;
	num_in_vtx = num_strings * in_vtx_oneside;
	ang_out = 2*PI/num_out_vtx;
	ang_in = 2*PI/num_in_vtx;

	v_out = new Vertex_handle[num_out_vtx];
	v_in = new Vertex_handle[num_in_vtx];
	out_nodes_coords = new double[num_strings*2];
	in_nodes_coords = new double[num_strings*2];

	for (i = 0; i < num_out_vtx; i++)
	{
	    v_out[i] = cdt.insert(Cgal_Point(
				CirCenter[0]+CirR[0]*cos(i*ang_out+offset),
				CirCenter[1]+CirR[0]*sin(i*ang_out+offset)));
	    if (0 == i%out_vtx_oneside)
	    {
		out_nodes_coords[i/out_vtx_oneside] = CirCenter[0]+
				CirR[0]*cos(i*ang_out+offset);
		out_nodes_coords[i/out_vtx_oneside+num_strings] = 
				CirCenter[1]+CirR[0]*sin(i*ang_out+offset);
	    }
	}

	for (i = 0; i < num_out_vtx-1; i++)
		cdt.insert_constraint(v_out[i],v_out[i+1]);
	cdt.insert_constraint(v_out[0],v_out[num_out_vtx-1]);


	if (CirR[1] != 0)
	{
	    for (i = 0; i < num_in_vtx; i++)
	    {
	        v_in[i] = cdt.insert(Cgal_Point(CirCenter[0]+
				CirR[1]*cos(i*ang_in+offset),
				CirCenter[1]+CirR[1]*sin(i*ang_in+offset)));
	        if (0 == i%in_vtx_oneside)
	        {
		    in_nodes_coords[i/in_vtx_oneside] = CirCenter[0]+
				CirR[1]*cos(i*ang_in+offset);
		    in_nodes_coords[i/in_vtx_oneside+num_strings] = 
				CirCenter[1]+CirR[1]*sin(i*ang_in+offset);
	        }
	    }
	    for (i = 0; i < num_in_vtx-1; i++)
		cdt.insert_constraint(v_in[i],v_in[i+1]);
	    cdt.insert_constraint(v_in[0],v_in[num_in_vtx-1]);

	    if (gore_bool[0]=='y'|| gore_bool[0]=='Y')
	    {
	    	for (i = 0; i < num_strings; i++)
		    cdt.insert_constraint(v_out[i*out_vtx_oneside],
				v_in[i*in_vtx_oneside]);
	    }
	    if (vent_bool[0]=='y'|| vent_bool[0]=='Y')
	    {
		list_of_seeds.push_back(Cgal_Point(CirCenter[0], CirCenter[1]));
	    }
	}
	
	CGAL::refine_Delaunay_mesh_2(cdt, list_of_seeds.begin(), 
			list_of_seeds.end(),Criteria(0.3, cri_dx));

	int *flag;
	flag = new int[cdt.number_of_faces()];
	double tri_center[2];

        i=0;
	if (vent_bool[0]=='y'|| vent_bool[0]=='Y')
	{
            for (fit = cdt.finite_faces_begin(); fit != cdt.finite_faces_end();
				 ++fit)
            {
                tri_center[0] = (fit->vertex(0)->point()[0] + 
				fit->vertex(1)->point()[0]
                    + fit->vertex(2)->point()[0]) / 3.0;
                tri_center[1] = (fit->vertex(0)->point()[1] + 
				fit->vertex(1)->point()[1]
                    + fit->vertex(2)->point()[1]) / 3.0;

                if (ptoutcircle(tri_center, CirCenter, CirR[1]))
                    flag[i] = 1;
                else
                    flag[i] = 0;
                i++;
            }
	}
	else
	{
            for (fit = cdt.finite_faces_begin(); fit != cdt.finite_faces_end();
				++fit)
	    	flag[i++] = 1;	
	}

	GenerateCgalSurf(front,surf,&cdt,flag,height);
	checkReducedTri(*surf);
    wave_type(*surf) = ELASTIC_BOUNDARY;
    FT_InstallSurfEdge(*surf,MONO_COMP_HSBDRY);
	setMonoCompBdryZeroLength(*surf);

    if (string_bool[0] == 'y' || string_bool[0] == 'Y')
	{
	    findStringNodePoints(*surf,out_nodes_coords,string_node_pts,
                                num_strings,&cbdry);
	    installString(infile,front,*surf,cbdry,string_node_pts,num_strings);
	}

    if (gore_bool[0]=='y'|| gore_bool[0]=='Y')
        {
            if (vent_bool[0] !='y' && vent_bool[0] !='Y')
		InstallInCurve(front,*surf,in_nodes_coords[0],
				in_nodes_coords[num_strings],
				CirCenter,num_in_vtx);

	    SplitCirBdry(front,*surf,num_strings,in_nodes_coords,
			POSITIVE_ORIENTATION);
	    if (string_bool[0] !='y' && string_bool[0] !='Y')
	    {
	    	SplitCirBdry(front,*surf,num_strings,out_nodes_coords,
			NEGATIVE_ORIENTATION);
	    }
	    InstallGore(front,*surf,num_strings,out_nodes_coords,
				in_nodes_coords);
	}
	setSurfZeroMesh(*surf);
	setMonoCompBdryZeroLength(*surf);
	FT_FreeThese(1,string_node_pts);
	if (consistent_interface(front->interf) == NO)
	    clean_up(ERROR);
}	/* end CgalCircle */

static void InstallInCurve(
        Front *front,
        SURFACE *surf,
	double x,
	double y,
	double *CirCenter,
	int num_vtx)
{
	double r = sqrt(sqr(x-CirCenter[0]) + sqr(y-CirCenter[1]));
	double cri = 1.1*r*(1-cos(PI/num_vtx));
	POINT *p,*ptmp,*p_end,*p_pre;
	TRI *tri,**tris;
	CURVE *vent_curve;
	BOND *bond;
	double tri_center[3];
	int i,j;
	NODE *node;
	AF_NODE_EXTRA *extra;
	AF_PARAMS *af_params = (AF_PARAMS*)front->extra2;
	double len_fac = af_params->gore_len_fac;

	p = NULL;
	surf_tri_loop(surf,tri)
	{
	    for (i = 0; i < 3; i++)
	    {
		ptmp = Point_of_tri(tri)[i];
		if (fabs(Coords(ptmp)[0] - x) < 1e-6 && 
		    fabs(Coords(ptmp)[1] - y) < 1e-6)
		{
		    p = ptmp;
		    break;
		}
	    }
	    if (p != NULL)
		break;
	}
	node = make_node(p);
        FT_ScalarMemoryAlloc((POINTER*)&extra,sizeof(AF_NODE_EXTRA));
        extra->af_node_type = GORE_NODE;
        node->extra = (POINTER)extra;
        node->size_of_extra = sizeof(AF_NODE_EXTRA);
	vent_curve = make_curve(0,0,node,node);
	install_curve_in_surface_bdry(surf,vent_curve,POSITIVE_ORIENTATION);
	hsbdry_type(vent_curve) = GORE_HSBDRY;
	bond = vent_curve->first;

	p_end = p;
	do
	{
	    tris = p->tris;
	    for (i = 0; i < p->num_tris; i++)
	    {
		POINT **pts = Point_of_tri(tris[i]);
		for (j = 0; j < 2; j++)
		    tri_center[j] = (Coords(pts[0])[j] + Coords(pts[1])[j]
				+ Coords(pts[2])[j])/3.0;

		if (sqr(tri_center[0]-CirCenter[0]) +
		    sqr(tri_center[1]-CirCenter[1]) > r*r)
		    continue;

                for (j = 0; j < 3; j++)
                    if (pts[j] == p)
                        break;
		p_pre = pts[(j+2)%3];

		if (p_pre == p_end)
		{
		    p = p_pre;
		    break;
		}
		
		if (fabs(sqrt(sqr(Coords(p_pre)[0]-CirCenter[0])+
			sqr(Coords(p_pre)[1]-CirCenter[1]))-r) < cri )
		{
                    insert_point_in_bond(p_pre,bond,vent_curve);
                    bond = bond->next;
                    p = p_pre;
		    break;
		}
	    }
	} while(p != p_end);

	linkCurveTriBond(vent_curve,surf);
	setCurveZeroLength(vent_curve,len_fac);
}	/* end InstallInCurve */

static void linkCurveTriBond(
	CURVE *c,
	SURFACE *surf)
{
	TRI **tris;
	TRI **pos_tris,**neg_tris;
	BOND *b;
	int num_pts;
	int i,j,k;

	num_pts = c->num_points;
        FT_VectorMemoryAlloc((POINTER*)&pos_tris,num_pts-1,sizeof(TRI*));
        FT_VectorMemoryAlloc((POINTER*)&neg_tris,num_pts-1,sizeof(TRI*));

	for (i = 0, b = c->first; b!= NULL; b = b->next, i++)
	{
	    tris = b->start->tris;
	    pos_tris[i] = neg_tris[i] = NULL;
	    for (j = 0; j < b->start->num_tris; j++)
	    {
		for (k = 0; k < 3; ++k)
		{
		    if (b->start == Point_of_tri(tris[j])[k] &&
		        b->end == Point_of_tri(tris[j])[(k+1)%3])
			pos_tris[i] = tris[j];
		    if (b->end == Point_of_tri(tris[j])[k] &&
		        b->start == Point_of_tri(tris[j])[(k+1)%3])
			neg_tris[i] = tris[j];
		}
	    }
            if (pos_tris[i] == NULL)
            {
                printf("pos_tris[%d] not found\n",i);
            }
            if (neg_tris[i] == NULL)
            {
                printf("neg_tris[%d] not found\n",i);
            }
	}
	for (i = 0, b = c->first; b!= NULL; b = b->next, i++)
	{
	    link_tri_to_bond(NULL,pos_tris[i],surf,b,c);
	    link_tri_to_bond(NULL,neg_tris[i],surf,b,c);
	}
	FT_FreeThese(2,pos_tris,neg_tris);
}	/* end linkCurveTriBond */

static void InstallGore(
        Front *front,
        SURFACE *surf,
        int num_gores,
        double *sta_coords,
        double *end_coords)
{
	INTERFACE *intfc = front->interf;
	CURVE **gore_curves;
	BOND *bond;
	POINT *sta_p;
	NODE *sta_node = NULL, *end_node = NULL, **n;
	TRI **tris;
	TRI *tri;
	double dir[3],dir_tmp[3];
	POINT *p_tmp[2];
	int i,j,k,l;
	AF_PARAMS *af_params = (AF_PARAMS*)front->extra2;
	double len_fac = af_params->gore_len_fac;

        FT_VectorMemoryAlloc((POINTER*)&gore_curves,num_gores,sizeof(CURVE*));

	for (i = 0; i < num_gores; i++)
	{
	    sta_node = NULL;
	    end_node = NULL;
	    intfc_node_loop(intfc,n)
	    {
		if (Coords((*n)->posn)[0] == sta_coords[i] &&
		    Coords((*n)->posn)[1] == sta_coords[i+num_gores])
		    sta_node = *n;
		else if (Coords((*n)->posn)[0] == end_coords[i] &&
			 Coords((*n)->posn)[1] == end_coords[i+num_gores])
		    end_node = *n;
	    }

	    gore_curves[i] = make_curve(0,0,sta_node,end_node);
	    install_curve_in_surface_bdry(surf,gore_curves[i],
				POSITIVE_ORIENTATION);
	    hsbdry_type(gore_curves[i]) = GORE_HSBDRY;
	    direction_vector(Coords(sta_node->posn),
				Coords(end_node->posn),dir,3);
	    bond = gore_curves[i]->first;
	    sta_p = sta_node->posn;
	    while (sta_p != end_node->posn)
	    {
		tris = sta_p->tris;
		for (j = 0; j < sta_p->num_tris; j++)
		{
		    for (k = 0; k < 3; k++)
			if (Point_of_tri(tris[j])[k] == sta_p)
			    break;
	    
		    p_tmp[0] = Point_of_tri(tris[j])[(k+1)%3];
		    p_tmp[1] = Point_of_tri(tris[j])[(k+2)%3];
		    for (l = 0; l < 2; l++)
		    {
		    	direction_vector(Coords(sta_p),Coords(p_tmp[l]),
					dir_tmp,3);
		    	for (k = 0; k < 3; k++)
			{
			    if (fabs(dir[k] - dir_tmp[k]) > 1e-6)
			    	break;
			}
			if (p_tmp[l] == end_node->posn)
			{
			    sta_p = p_tmp[l];
			    break;
			}
		    	if (3 == k)
		    	{
			    insert_point_in_bond(p_tmp[l],bond,gore_curves[i]);
			    bond = bond->next;
			    sta_p = p_tmp[l];
			    break;
		    	}
		    }
		    if (3 == k || p_tmp[l] == end_node->posn)
			break; 
		}
	    }
	    linkCurveTriBond(gore_curves[i],surf);
	    setCurveZeroLength(gore_curves[i],len_fac);
	}
}	/* end InstallGore */

static void SplitCirBdry(
        Front *front,
        SURFACE *surf,
        int num_strings,
        double *nodes_coords,
	ORIENTATION ORIEN)
{
	INTERFACE *intfc = front->interf;
	CURVE *vent_bdry,**c;
	BOND *bond, *pre_bond;
	NODE **vent_nodes;
        AF_NODE_EXTRA *extra;
        boolean node_moved;
	int i;

        FT_VectorMemoryAlloc((POINTER*)&vent_nodes,num_strings,
                                sizeof(NODE*));

	vent_bdry = NULL;
	surf_pos_curve_loop(surf,c)
        {
            if (vent_bdry != NULL)
                break;
            curve_bond_loop(*c,bond)
            if (Coords(bond->start)[0] == nodes_coords[0] &&
                Coords(bond->start)[1] == nodes_coords[num_strings])
            {
                vent_bdry = *c;
                break;
            }
        }
        surf_neg_curve_loop(surf,c)
        {
            if (vent_bdry != NULL)
                break;
            curve_bond_loop(*c,bond)
            if (Coords(bond->start)[0] == nodes_coords[0] &&
                Coords(bond->start)[1] == nodes_coords[num_strings])
            {
                vent_bdry = *c;
                break;
            }
        }

        node_moved = NO;
        for (i = 0; i < num_strings; ++i)
        {
            if (!node_moved)
            {
                node_moved = YES;
                curve_bond_loop(vent_bdry,bond)
                {
                    if (Coords(bond->start)[0] == nodes_coords[i] &&
                        Coords(bond->start)[1] == nodes_coords[i+num_strings])
                    {
                        move_closed_loop_node(vent_bdry,bond);
                        vent_nodes[i] = I_NodeOfPoint(intfc,bond->start);
                        FT_ScalarMemoryAlloc((POINTER*)&extra,
					sizeof(AF_NODE_EXTRA));
                        extra->af_node_type = GORE_NODE;
                        vent_nodes[i]->extra = (POINTER)extra;
                        vent_nodes[i]->size_of_extra = sizeof(AF_NODE_EXTRA);
                        break;
                    }
                }
                continue;
            }
            curve_bond_loop(vent_bdry,bond)
            {
                if (Coords(bond->start)[0] == nodes_coords[num_strings-i] &&
                    Coords(bond->start)[1] == nodes_coords[2*num_strings-i])
                {
		    pre_bond = bond->prev;
                    split_curve(bond->start,bond,vent_bdry,0,0,0,0);
                    vent_nodes[i] = I_NodeOfPoint(intfc,bond->start);
                    FT_ScalarMemoryAlloc((POINTER*)&extra,
					sizeof(AF_NODE_EXTRA));
                    extra->af_node_type = GORE_NODE;
                    vent_nodes[i]->extra = (POINTER)extra;
                    vent_nodes[i]->size_of_extra = sizeof(AF_NODE_EXTRA);
                    break;
                }
            }
	    if (ORIEN == POSITIVE_ORIENTATION)
                vent_bdry = I_CurveOfPoint(intfc,bond->start,&bond);
	    else if (ORIEN == NEGATIVE_ORIENTATION)
                vent_bdry = I_CurveOfPoint(intfc,pre_bond->start,&bond);
        }
}	/* end SplitInCirBdry */
            
static void CgalEllipse(
	FILE *infile,
	Front *front,
	SURFACE **surf)
{
	double height;
	double EliCenter[2];
	double EliR[2];
	double xrange[2];
	char gore_bool[10];
        int num_strings;
        int num_out_vtx;
        int num_gore_oneside;
	POINT **string_node_pts;
	double startx, distance;
        double *out_nodes_coords,*out_vtx_coords;
        CDT cdt;
        CDT::Finite_faces_iterator fit;
	AF_PARAMS *af_params = (AF_PARAMS*)front->extra2;
	CURVE *cbdry;

	Vertex_handle *v_out;
	double cri_dx = 0.6*computational_grid(front->interf)->h[0];
        double width=cri_dx;
        int i,j;

	CursorAfterString(infile,"Enter the height of the plane:");
        fscanf(infile,"%lf",&height);
        (void) printf("%f\n",height);
	CursorAfterString(infile,"Enter ellipse center:");
        fscanf(infile,"%lf %lf",&EliCenter[0],&EliCenter[1]);
        (void) printf("%f %f\n",EliCenter[0],EliCenter[1]);
        CursorAfterString(infile,"Enter ellipse radius:");
        fscanf(infile,"%lf %lf",&EliR[0],&EliR[1]);
        (void) printf("%f %f\n",EliR[0],EliR[1]);
        CursorAfterString(infile,"Enter x range of ellipse:");
        fscanf(infile,"%lf %lf",&xrange[0],&xrange[1]);
        (void) printf("%f %f\n",xrange[0],xrange[1]);
	CursorAfterStringOpt(infile,"Enter yes to attach gores to canopy:");
	fscanf(infile,"%s",gore_bool);
        (void) printf("%s\n",gore_bool);		
	if (gore_bool[0] == 'y' || gore_bool[0] == 'Y')
	    af_params->attach_gores = YES;
	else
	    af_params->attach_gores = NO;
	CursorAfterStringOpt(infile,"Enter number of vertical gores:");
	fscanf(infile,"%d",&num_gore_oneside);
        (void) printf("%d\n",num_gore_oneside);		
	CursorAfterString(infile,"Enter start x-coordinate of gore:");
        fscanf(infile,"%lf",&startx);
        (void) printf("%f\n",startx);
	CursorAfterString(infile,"Enter distance between gores:");
        fscanf(infile,"%lf",&distance);
        (void) printf("%f\n",distance);

	num_strings = 2*(num_gore_oneside+2);
        out_nodes_coords = new double[num_strings*2];
	FT_VectorMemoryAlloc((POINTER*)&string_node_pts,num_strings,
                                sizeof(POINT*));

	out_nodes_coords[0]=out_nodes_coords[num_strings-1]=xrange[0];
	out_nodes_coords[num_gore_oneside+1]=
			out_nodes_coords[num_gore_oneside+2]=xrange[1];
        for (i=1; i<num_gore_oneside+1; i++)
            out_nodes_coords[i]=out_nodes_coords[num_strings-1-i]=
			startx+(i-1)*distance;

	for (i=0; i<num_gore_oneside+2; i++)
	{
	    out_nodes_coords[i+num_strings]=EliCenter[1] - 
			sqrt(1-sqr(out_nodes_coords[i]-EliCenter[0])
			/sqr(EliR[0])) * EliR[1];
	    out_nodes_coords[2*num_strings-1-i]=2*EliCenter[1]-
			out_nodes_coords[i+num_strings];
	}

	int *N_sub;
	int gl_idx=0;
	N_sub = new int[num_gore_oneside+1];    
	num_out_vtx = 0;
	for (i=0; i<num_gore_oneside+1; i++)
	{
	    N_sub[i] = (out_nodes_coords[i+1] - out_nodes_coords[i]) / width;
    	    num_out_vtx+= N_sub[i];
	}
	num_out_vtx=(num_out_vtx+1)*2;

        out_vtx_coords = new double[num_out_vtx*2];
        v_out = new Vertex_handle[num_out_vtx];

	for(j=0; j<num_gore_oneside+1; j++)
	{
	    width = (out_nodes_coords[j+1]-out_nodes_coords[j])/N_sub[j];

    	    for(i=0; i<N_sub[j]; i++)
    	    {
		out_vtx_coords[gl_idx]=out_vtx_coords[num_out_vtx-1-gl_idx]=
			out_nodes_coords[j]+i*width;
		out_vtx_coords[gl_idx+num_out_vtx]=EliCenter[1] - 
			sqrt(1-sqr(out_vtx_coords[gl_idx]-EliCenter[0])
			/sqr(EliR[0])) * EliR[1];
		out_vtx_coords[2*num_out_vtx-1-gl_idx]=2*EliCenter[1]-
			out_vtx_coords[gl_idx+num_out_vtx];
      		gl_idx++;
	    }
	}
	out_vtx_coords[num_out_vtx/2-1]=out_nodes_coords[num_gore_oneside+1];
	out_vtx_coords[num_out_vtx*3/2-1]=out_nodes_coords[num_gore_oneside+
			1+num_strings];
	out_vtx_coords[num_out_vtx/2]=out_nodes_coords[num_gore_oneside+2];
	out_vtx_coords[num_out_vtx*3/2]=out_nodes_coords[num_gore_oneside+
			2+num_strings];

	for (i=0; i<num_out_vtx; i++)
            v_out[i]=cdt.insert(Cgal_Point(out_vtx_coords[i],
			out_vtx_coords[i+num_out_vtx]));
        for (i=0; i<num_out_vtx-1; i++)
            cdt.insert_constraint(v_out[i],v_out[i+1]);
        cdt.insert_constraint(v_out[0],v_out[num_out_vtx-1]);

	gl_idx = 0;
	if (gore_bool[0]=='y'|| gore_bool[0]=='Y')
	{
	    for (i = 0; i < num_gore_oneside; i++)
	    {
		gl_idx += N_sub[i];
		cdt.insert_constraint(v_out[gl_idx],
			v_out[num_out_vtx-1-gl_idx]);
	    }
	}

	CGAL::refine_Delaunay_mesh_2(cdt, Criteria(0.125, cri_dx)); 

	int *flag;
	flag = new int[cdt.number_of_faces()];

	i = 0;
        for (fit = cdt.finite_faces_begin(); fit != cdt.finite_faces_end(); 
			++fit)
	    flag[i++] = 1;	

	GenerateCgalSurf(front,surf,&cdt,flag,height);
        wave_type(*surf) = ELASTIC_BOUNDARY;
        FT_InstallSurfEdge(*surf,MONO_COMP_HSBDRY);
	/* curve the surface */
	double coeff, cen[3];
	TRI *tri;
	CursorAfterString(infile,"Enter vertex coordinate of the paraboloid:");
	fscanf(infile,"%lf %lf %lf\n",&cen[0],&cen[1],&cen[2]);
	(void) printf("%f %f %f\n",cen[0],cen[1],cen[2]);
        CursorAfterString(infile,"Enter coefficient of the paraboloid:");
        fscanf(infile,"%lf",&coeff);
        (void) printf("%f\n",coeff);
        unsort_surf_point(*surf);
        surf_tri_loop(*surf,tri)
        {   
            for (i = 0; i < 3; i++)
            {   
                tri->side_length0[i] = -1.0;
                POINT* p = Point_of_tri(tri)[i];
                if (sorted(p) == NO)
                {   
                    Coords(p)[2] = cen[2] - coeff*sqr(Coords(p)[0]-cen[0]);
                    sorted(p) = YES;
                }
            }
        }
	setSurfZeroMesh(*surf);
	/* end curve the surface */
	setMonoCompBdryZeroLength(*surf);

	findStringNodePoints(*surf,out_nodes_coords,string_node_pts,
                                num_strings,&cbdry);
	installString(infile,front,*surf,cbdry,string_node_pts,num_strings);

        if (gore_bool[0]=='y'|| gore_bool[0]=='Y')
        {
	    double *gore_sta_nodes,*gore_end_nodes;

	    gore_sta_nodes = new double[num_gore_oneside*2];
	    gore_end_nodes = new double[num_gore_oneside*2];
	    for (i = 0; i < num_gore_oneside; i++)
	    {
	        gore_sta_nodes[i] = out_nodes_coords[i+1];
	        gore_sta_nodes[i+num_gore_oneside] = 
			out_nodes_coords[num_strings+i+1];
	        gore_end_nodes[i] = out_nodes_coords[num_strings-2-i];
	        gore_end_nodes[i+num_gore_oneside] = 
			out_nodes_coords[2*num_strings-2-i];
	    }
	    InstallGore(front,*surf,num_gore_oneside,gore_sta_nodes,
			gore_end_nodes);
	}
	setMonoCompBdryZeroLength(*surf);
}	/* end CgalEllipse */

static void CgalCross(
	FILE *infile,
	Front *front,
	SURFACE **surf)
{
	double height;
	double L[2][2],U[2][2];
	char gore_bool[10];
	int num_strings;
	int num_out_vtx;
	int num_gore_oneside;
	POINT **string_node_pts;
	double *out_nodes_coords,*out_vtx_coords;
	CDT cdt;
	CDT::Finite_faces_iterator fit;
	AF_PARAMS *af_params = (AF_PARAMS*)front->extra2;
	CURVE *cbdry;

	Vertex_handle *v_out;
	double width;
	double cri_dx = 0.6*computational_grid(front->interf)->h[0];
	int i,j;

        CursorAfterString(infile,"Enter the height of the plane:");
        fscanf(infile,"%lf",&height);
        (void) printf("%f\n",height);
	(void) printf("Input the two crossing rectangles\n");
        CursorAfterString(infile,"Enter lower bounds of first rectangle:");
        fscanf(infile,"%lf %lf",&L[0][0],&L[0][1]);
        (void) printf("%f %f\n",L[0][0],L[0][1]);
        CursorAfterString(infile,"Enter upper bounds of first rectangle:");
        fscanf(infile,"%lf %lf",&U[0][0],&U[0][1]);
        (void) printf("%f %f\n",U[0][0],U[0][1]);
        CursorAfterString(infile,"Enter lower bounds of second rectangle:");
        fscanf(infile,"%lf %lf",&L[1][0],&L[1][1]);
        (void) printf("%f %f\n",L[1][0],L[1][1]);
        CursorAfterString(infile,"Enter upper bounds of second rectangle:");
        fscanf(infile,"%lf %lf",&U[1][0],&U[1][1]);
        (void) printf("%f %f\n",U[1][0],U[1][1]);
	CursorAfterStringOpt(infile,"Enter yes to attach gores to canopy:");
	fscanf(infile,"%s",gore_bool);
        (void) printf("%s\n",gore_bool);		
	if (gore_bool[0] == 'y' || gore_bool[0] == 'Y')
	    af_params->attach_gores = YES;
	else
	    af_params->attach_gores = NO;
	CursorAfterStringOpt(infile,"Enter number of chords on one side:");
        fscanf(infile,"%d",&num_gore_oneside);
        (void) printf("%d\n",num_gore_oneside);


	num_strings = 4*num_gore_oneside;
	num_out_vtx = num_strings + 4;
	out_nodes_coords = new double[num_strings*2];
	out_vtx_coords = new double[num_out_vtx*2];
	v_out = new Vertex_handle[num_out_vtx];
	FT_VectorMemoryAlloc((POINTER*)&string_node_pts,num_strings,
				sizeof(POINT*));

	width=(U[0][0]-L[0][0])/(num_gore_oneside-1);

	for (i=0; i<num_gore_oneside; i++)
	{
	    out_nodes_coords[i]=L[0][0]+i*width;
	    out_nodes_coords[i+num_strings]=L[0][1];
	    out_nodes_coords[i+num_gore_oneside]=U[1][0];
	    out_nodes_coords[i+num_gore_oneside+num_strings]=U[1][1]-
			(num_gore_oneside-i-1)*width;
	    out_nodes_coords[i+2*num_gore_oneside]=U[0][0]-i*width;
            out_nodes_coords[i+2*num_gore_oneside+num_strings]=U[0][1];
	    out_nodes_coords[i+3*num_gore_oneside]=L[1][0];
            out_nodes_coords[i+3*num_gore_oneside+num_strings]=L[1][1]+
			(num_gore_oneside-i-1)*width;
	}
	for (i=0; i<4; i++)
	{
	    for(j=0; j<num_gore_oneside; j++)
	    {
		out_vtx_coords[i*(num_gore_oneside+1)+j]=
			out_nodes_coords[i*num_gore_oneside+j];	
		out_vtx_coords[i*(num_gore_oneside+1)+j+num_out_vtx]=
		    out_nodes_coords[i*num_gore_oneside+j+num_strings];
	    }
	}
	out_vtx_coords[num_gore_oneside]=U[0][0];
	out_vtx_coords[num_gore_oneside+num_out_vtx]=L[1][1];
	out_vtx_coords[2*num_gore_oneside+1]=U[0][0];
        out_vtx_coords[2*num_gore_oneside+1+num_out_vtx]=U[1][1];
	out_vtx_coords[3*num_gore_oneside+2]=L[0][0];
        out_vtx_coords[3*num_gore_oneside+2+num_out_vtx]=U[1][1];
	out_vtx_coords[4*num_gore_oneside+3]=L[0][0];
        out_vtx_coords[4*num_gore_oneside+3+num_out_vtx]=L[1][1];

	for (i=0; i<num_out_vtx; i++)
	    v_out[i]=cdt.insert(Cgal_Point(out_vtx_coords[i],
			out_vtx_coords[i+num_out_vtx]));
	for (i=0; i<num_out_vtx-1; i++)
	    cdt.insert_constraint(v_out[i],v_out[i+1]);
	cdt.insert_constraint(v_out[0],v_out[num_out_vtx-1]);

	CGAL::refine_Delaunay_mesh_2(cdt, Criteria(0.125, cri_dx)); 

	int *flag;
	flag = new int[cdt.number_of_faces()];
	double tri_center[2];

	i=0;
	for (fit = cdt.finite_faces_begin(); fit != cdt.finite_faces_end(); 
			++fit)
        {
            tri_center[0] = (fit->vertex(0)->point()[0] + 
			fit->vertex(1)->point()[0] 
		+ fit->vertex(2)->point()[0]) / 3.0;
            tri_center[1] = (fit->vertex(0)->point()[1] + 
			fit->vertex(1)->point()[1] 
		+ fit->vertex(2)->point()[1]) / 3.0;

            if (ptinbox(tri_center, L[0], U[0]) || ptinbox(tri_center, L[1], 
			U[1]))
                flag[i] = 1;
            else
                flag[i] = 0;
            i++;
        }

	GenerateCgalSurf(front,surf,&cdt,flag,height);
        wave_type(*surf) = ELASTIC_BOUNDARY;
        FT_InstallSurfEdge(*surf,MONO_COMP_HSBDRY);
	setMonoCompBdryZeroLength(*surf);
	findStringNodePoints(*surf,out_nodes_coords,string_node_pts,
				num_strings,&cbdry);
	foldSurface(infile,*surf);
	if (sewSurface(infile,*surf))
	{
	    resetStringNodePoints(*surf,string_node_pts,&num_strings,&cbdry);
	}
	setSurfZeroMesh(*surf);
	setMonoCompBdryZeroLength(*surf);
	installString(infile,front,*surf,cbdry,string_node_pts,num_strings);
	FT_FreeThese(1,string_node_pts);
}	/* end CgalCross */

static void foldSurface(
	FILE *infile,
	SURFACE *surf)
{
	double dir[MAXD],axis[MAXD],angle;
	char string[200];
	SIDE side;
	int i,num_foldings;
	boolean first;

	if (CursorAfterStringOpt(infile,
            "Entering yes to fold surface:"))
        {
            fscanf(infile,"%s",string);
            (void) printf("%s\n",string);
            if (string[0] != 'y' && string[0] != 'Y')
                return;
        }
        else
            return;

	CursorAfterString(infile,"Enter number of folding steps: ");
        fscanf(infile,"%d",&num_foldings);
        (void) printf("%d\n",num_foldings);
	first = YES;
	for (i = 0; i < num_foldings; ++i)
	{
	    sprintf(string,"Folding step %d",i+1);
	    CursorAfterString(infile,string);
	    (void) printf("\n");
	    CursorAfterString(infile,"Enter direction of folding: ");
            fscanf(infile,"%lf %lf %lf",dir,dir+1,dir+2);
            (void) printf("%f %f %f\n",dir[0],dir[1],dir[2]);
	    CursorAfterString(infile,"Enter axis point of folding: ");
            fscanf(infile,"%lf %lf %lf",axis,axis+1,axis+2);
            (void) printf("%f %f %f\n",axis[0],axis[1],axis[2]);
	    CursorAfterString(infile,"Enter angle of folding: ");
            fscanf(infile,"%lf",&angle);
            (void) printf("%f\n",angle);
	    angle *= PI/180.0;
	    CursorAfterString(infile,"Enter folding side: ");
            fscanf(infile,"%s",string);
            (void) printf("%s\n",string);
	    if (string[0] == 'p' || string[0] == 'P')
	    	side = POSITIVE_SIDE;
	    else if (string[0] == 'n' || string[0] == 'N')
	    	side = NEGATIVE_SIDE;
	    I_FoldSurface(surf,dir,axis,angle,side,first);
	    first = NO;
	}
}	/* end foldSurface */

static boolean sewSurface(
	FILE *infile,
	SURFACE *surf)
{
	double crds_start[MAXD],crds_end[MAXD];
	char string[200];
	int i,num_stitches;

	if (debugging("sewing"))
	    (void) printf("Entering  sewSurface()\n");

	if (CursorAfterStringOpt(infile,
            "Entering yes to sew surface:"))
        {
            fscanf(infile,"%s",string);
            (void) printf("%s\n",string);
            if (string[0] != 'y' && string[0] != 'Y')
                return NO;
        }
        else
            return NO;

        (void) printf("Stitches must be along existing curves\n");
	CursorAfterString(infile,"Enter number of sewing stitches: ");
        fscanf(infile,"%d",&num_stitches);
        (void) printf("%d\n",num_stitches);
	for (i = 0; i < num_stitches; ++i)
	{
	    sprintf(string,"For stitche %d",i+1);
	    CursorAfterString(infile,string);
	    (void) printf("\n");
	    CursorAfterString(infile,"Enter start point: ");
            fscanf(infile,"%lf %lf %lf",crds_start,crds_start+1,crds_start+2);
            (void) printf("%f %f %f\n",crds_start[0],crds_start[1],
					crds_start[2]);
	    CursorAfterString(infile,"Enter end point: ");
            fscanf(infile,"%lf %lf %lf",crds_end,crds_end+1,crds_end+2);
            (void) printf("%f %f %f\n",crds_end[0],crds_end[1],crds_end[2]);
	    I_SewSurface(surf,crds_start,crds_end);
	}
	if (debugging("sewing"))
	{
	    (void) printf("Leaving sewSurface()\n");
	}
	return YES;
}	/* end sewSurface */

static void GenerateCgalSurf(
	Front *front,
	SURFACE **surf,
	CDT *cdt,
	int *flag,
	double height)
{
	COMPONENT amb_comp = front->interf->default_comp;
	COMPONENT  neg_comp, pos_comp;
	INTERFACE *intfc;
	SURFACE *newsurf;
	int num_vtx,num_tris,num_point_tris,num_ptris;
	int* index;
	double *vertex;
	int i,j,k;
	POINT **points,*p;
	TRI **tris,**ptris;
	INTERFACE *sav_intfc;

	CDT::Finite_vertices_iterator vit;
	CDT::Finite_faces_iterator fit;

	neg_comp = amb_comp;
	pos_comp = amb_comp;
	
	intfc = front->interf;
	sav_intfc = current_interface();
	set_current_interface(intfc);

	newsurf = make_surface(neg_comp,pos_comp,NULL,NULL);
	num_vtx = cdt->number_of_vertices();
        uni_array(&vertex,3*num_vtx,FLOAT);
        uni_array(&index,num_vtx,INT);

        for (i = 0, vit = cdt->finite_vertices_begin();
		vit != cdt->finite_vertices_end(); ++vit)
        {
	    vit->info() = index[i] = i;
	    vertex[3*i] = vit->point()[0];
	    vertex[3*i+1] = vit->point()[1];
	    vertex[3*i+2] = height;
	    i++;
        }
	uni_array(&points,num_vtx,sizeof(POINT*));
	for ( i = 0; i < num_vtx; ++i)
	{
	    points[i] = Point(vertex+3*i);
	    points[i]->num_tris = 0;
	} 

	num_tris = 0;
	for (i = 0, fit = cdt->finite_faces_begin(); 
		fit != cdt->finite_faces_end(); ++fit)
	{
	    if (!flag[i++])
		continue;
	    num_tris++;
	}
	num_point_tris = num_tris * 3;
	uni_array(&tris,num_tris,sizeof(TRI*));
	for (i = 0, j = 0, fit = cdt->finite_faces_begin(); 
		fit != cdt->finite_faces_end(); ++fit)
	{
	    if (!flag[i++])
		continue;

	    int i1,i2,i3;

	    i1 = index[fit->vertex(0)->info()];
	    i2 = index[fit->vertex(1)->info()];
	    i3 = index[fit->vertex(2)->info()];
	    tris[j] = make_tri(points[i1],points[i2],points[i3],NULL,NULL,
				NULL,NO);

	    tris[j]->surf = newsurf;
            points[i1]->num_tris++;
            points[i2]->num_tris++;
            points[i3]->num_tris++;
	    j++;
	}
	intfc->point_tri_store = (TRI**)store(num_point_tris*sizeof(TRI*));
	ptris = intfc->point_tri_store;
	for (i = 0; i < num_vtx; ++i)
        {
            points[i]->tris = ptris;
            ptris += points[i]->num_tris;
            points[i]->num_tris = 0;
        }
        for (i = 0; i < num_tris; ++i)
        {
            if (i != 0)
            {
                tris[i]->prev = tris[i-1];
                tris[i-1]->next = tris[i];
            }
            for (j = 0; j < 3; ++j)
            {
                p = Point_of_tri(tris[i])[j];
                p->tris[p->num_tris++] = tris[i];
            }
        }
        for (i = 0; i < num_vtx; ++i)
        {
	    int m,l;
            ptris = points[i]->tris;
            num_ptris = points[i]->num_tris;
            for (j = 0; j < num_ptris; ++j)
            for (k = 0; k < j; ++k)
            {
                TRI *tri1 = ptris[j];
                TRI *tri2 = ptris[k];
                for (m = 0; m < 3; ++m)
                for (l = 0; l < 3; ++l)
                {
                    if ((Point_of_tri(tri1)[m] == Point_of_tri(tri2)[(l+1)%3] 
			&& 
			Point_of_tri(tri1)[(m+1)%3] == Point_of_tri(tri2)[l]))
                    {
                        Tri_on_side(tri1,m) = tri2;
                        Tri_on_side(tri2,l) = tri1;
                    }
                }
            }
        }
        newsurf->num_tri = num_tris;

        first_tri(newsurf) = tris[0];
        last_tri(newsurf) = tris[num_tris-1];
        last_tri(newsurf)->next = tail_of_tri_list(newsurf);
        first_tri(newsurf)->prev = head_of_tri_list(newsurf);
        reset_intfc_num_points(newsurf->interface);

	*surf = newsurf;	
	setSurfZeroMesh(newsurf);
	set_current_interface(sav_intfc);
}	/* end GenerateCgalSurf */

bool ptinbox(
	double *c,
	double *l,
	double *u)
{
	if (c[0] > l[0] && c[0] < u[0] && c[1] > l[1] && c[1] < u[1])
	    return 1;
	else
	    return 0;
}	/* end ptinbox */

bool ptoutcircle(
	double *p, 
	double *c, 
	double r)
{
	if ((p[0]-c[0])*(p[0]-c[0]) + (p[1]-c[1])*(p[1]-c[1]) > r*r)
	    return 1;
	else
	    return 0;
}	/* end ptoutcircle */

static void findStringNodePoints(
	SURFACE *surf,
	double *nodes_coords,
	POINT **string_node_pts,
	int num_strings,
	CURVE **cbdry)
{
	CURVE *canopy_bdry = NULL;
	CURVE **c;
	BOND *bond;
	int i;
	surf_pos_curve_loop(surf,c)
	{
	    if (canopy_bdry != NULL)
		break;
	    curve_bond_loop(*c,bond)
	    if (Coords(bond->start)[0] == nodes_coords[0] &&
            	Coords(bond->start)[1] == nodes_coords[num_strings])
	    {
	    	canopy_bdry = *c; 
	    	break;
	    }
	}
	surf_neg_curve_loop(surf,c)
	{
	    if (canopy_bdry != NULL)
		break;
	    curve_bond_loop(*c,bond)
	    if (Coords(bond->start)[0] == nodes_coords[0] &&
            	Coords(bond->start)[1] == nodes_coords[num_strings])
	    {
	    	canopy_bdry = *c; 
	    	break;
	    }
	}

	for (i = 0; i < num_strings; ++i)
	{
	    curve_bond_loop(canopy_bdry,bond)
	    {
		if (Coords(bond->start)[0] == nodes_coords[i] &&
		    Coords(bond->start)[1] == nodes_coords[i+num_strings])
		{
		    string_node_pts[i] = bond->start;
		}
	    }
	    if (Coords(canopy_bdry->end->posn)[0] == nodes_coords[i] &&
		Coords(canopy_bdry->end->posn)[1] == nodes_coords[i+num_strings])
	    {
		string_node_pts[i] = canopy_bdry->end->posn;
	    }
	}
	*cbdry = canopy_bdry;
}	/* end findStringNodePoints */

static void installString(
	FILE *infile,
	Front *front,
	SURFACE *surf,
	CURVE *canopy_bdry,
	POINT **string_node_pts,
	int num_strings)
{
	INTERFACE *intfc = front->interf;
	char string[100];
	double cload[MAXD],coords[MAXD],dir[MAXD];
	NODE *nload,**string_nodes;
	CURVE **string_curves,**c;
	BOND *bond;
	AF_NODE_EXTRA *extra;
	boolean node_moved;
	double spacing,*h = computational_grid(intfc)->h;
	int nb;
	int i, j,k;
	double length,len_fac;
	AF_PARAMS *af_params = (AF_PARAMS*)front->extra2;

        CursorAfterString(infile,"Enter initial position of load:");
        fscanf(infile,"%lf %lf %lf",&cload[0],&cload[1],&cload[2]);
        (void) printf("%f %f %f\n",cload[0],cload[1],cload[2]);

	nload = make_node(Point(cload));
	FT_ScalarMemoryAlloc((POINTER*)&extra,sizeof(AF_NODE_EXTRA));
        extra->af_node_type = LOAD_NODE;
	if (CursorAfterStringOpt(infile,"Enter yes to fix the load node:"))
	{
	    char s[100];
	    fscanf(infile,"%s",s);
            (void) printf("%s\n",s);
	    if (s[0] == 'Y' || s[0] == 'y')
		extra->af_node_type = PRESET_NODE;
	}
        nload->extra = (POINTER)extra;
        nload->size_of_extra = sizeof(AF_NODE_EXTRA);

	FT_VectorMemoryAlloc((POINTER*)&string_nodes,num_strings,
                                sizeof(NODE*));

	node_moved = NO;
	for (i = 0; i < num_strings; ++i)
	{
	    if (!node_moved)
	    {
		node_moved = YES;
		curve_bond_loop(canopy_bdry,bond)
		{
		    if (bond->start == string_node_pts[i])
		    {
			move_closed_loop_node(canopy_bdry,bond);
			string_nodes[i] = I_NodeOfPoint(intfc,bond->start);
			FT_ScalarMemoryAlloc((POINTER*)&extra,
					sizeof(AF_NODE_EXTRA));
                	extra->af_node_type = STRING_NODE;
                	string_nodes[i]->extra = (POINTER)extra;
                	string_nodes[i]->size_of_extra = sizeof(AF_NODE_EXTRA);
			break;
		    }
		}
		continue;		
	    }
	    curve_bond_loop(canopy_bdry,bond)
	    {
		if (bond->start == string_node_pts[i])
		{
		    split_curve(bond->start,bond,canopy_bdry,0,0,0,0);
		    string_nodes[i] = I_NodeOfPoint(intfc,bond->start);
                    FT_ScalarMemoryAlloc((POINTER*)&extra,
					sizeof(AF_NODE_EXTRA));
                    extra->af_node_type = STRING_NODE;
                    string_nodes[i]->extra = (POINTER)extra;
                    string_nodes[i]->size_of_extra = sizeof(AF_NODE_EXTRA);
                    break;
		}
	    }
	    canopy_bdry = I_CurveOfPoint(intfc,bond->start,&bond);
	}

	if (CursorAfterStringOpt(infile,
			"Enter yes to install the strings to RGB:"))
	{
	    fscanf(infile,"%s",string);
	    (void) printf("%s\n",string);
	    if (string[0] == 'y' || string[0] == 'Y')
	    {
		int rg_index;
		SURFACE **rg_surf;
		CursorAfterString(infile, 
				"Enter the body index of the target RGB:");
		fscanf(infile,"%d",&rg_index);
		(void) printf("%d\n",rg_index);
		intfc_surface_loop(intfc, rg_surf)
		{
		    if ((wave_type(*rg_surf) == NEUMANN_BOUNDARY) || 
			(wave_type(*rg_surf) == MOVABLE_BODY_BOUNDARY))
		    {
			if (body_index(*rg_surf) == rg_index)
			    break;
		    }
		}
		connectStringtoRGB(front,*rg_surf,string_nodes,num_strings);
		delete_node(nload);
		return;
	    }
	}

	/* make the all initial springs at their equilibruim length */
        FT_VectorMemoryAlloc((POINTER*)&string_curves,num_strings,
                                sizeof(CURVE*));
	for (i = 0; i < num_strings; ++i)
	{
	    string_curves[i] = make_curve(0,0,string_nodes[i],nload);
	    hsbdry_type(string_curves[i]) = STRING_HSBDRY;
	    spacing = separation(string_nodes[i]->posn,nload->posn,3);
	    for (j = 0; j < 3; ++j)
		dir[j] = (Coords(nload->posn)[j] -
			Coords(string_nodes[i]->posn)[j])/spacing;
	    nb = (int)spacing/(0.40*h[0]);
	    spacing /= (double)nb;
	    bond = string_curves[i]->first;
	    for (j = 1; j < nb; ++j)
	    {
		for (k = 0; k < 3; ++k)
		    coords[k] = Coords(string_nodes[i]->posn)[k] +
					j*dir[k]*spacing;
		insert_point_in_bond(Point(coords),bond,string_curves[i]);
		bond->length0 = spacing;
		bond = bond->next;
	    }
	    bond->length0 = spacing;
	    af_params->string_curves.push_back(string_curves[i]);
	}
	FT_FreeThese(1,string_curves);
}	/* end installString */

static void resetStringNodePoints(
	SURFACE *surf,
	POINT **string_node_pts,
	int *num_strings,
	CURVE **cbdry)
{
	CURVE *canopy_bdry = NULL;
	CURVE **c;
	BOND *bond;
	int i,j,nv;
	boolean node_found;

	nv = *num_strings;
	if (debugging("sewing"))
	{
	    (void) printf("Entering resetStringNodePoints()\n");
	    (void) printf("nv = %d\n",nv);
	}
	surf_pos_curve_loop(surf,c)
	{
	    if (canopy_bdry != NULL)
		break;
	    curve_bond_loop(*c,bond)
	    {
		for (i = 0; i < nv; ++i)
		{
		    if (bond->start == string_node_pts[i])
		    {
	    		canopy_bdry = *c; 
	    		break;
		    }	
		    if (canopy_bdry != NULL)
			break;
		}
	    }
	}
	surf_neg_curve_loop(surf,c)
	{
	    if (canopy_bdry != NULL)
		break;
	    curve_bond_loop(*c,bond)
	    {
		for (i = 0; i < nv; ++i)
		{
		    if (bond->start == string_node_pts[i])
		    {
	    		canopy_bdry = *c; 
	    		break;
		    }	
		    if (canopy_bdry != NULL)
			break;
		}
	    }
	}

	for (i = 0; i < nv; ++i)
	{
	    node_found = NO;
	    curve_bond_loop(canopy_bdry,bond)
	    {
		if (bond->start == string_node_pts[i])
		    node_found = YES;
	    }
	    if (canopy_bdry->last->end == string_node_pts[i])
		node_found = YES;
	    if (!node_found)
	    {
		for (j = i; j < nv-1; ++j)
		    string_node_pts[j] = string_node_pts[j+1];
		nv--;
	    }
	}
	*num_strings = nv;
	*cbdry = canopy_bdry;
	if (debugging("sewing"))
	{
	    (void) printf("Leaving resetStringNodePoints()\n");
	    (void) printf("nv = %d\n",nv);
	}
}	/* end resetStringNodePoints */

static void setSurfZeroMesh(
	SURFACE *surf)
{
	TRI *t;
	int i,j;
	double total_num_sides;
	double max_len,min_len,ave_len,len;

	if (debugging("zero_mesh"))
	    (void) printf("Entering setSurfZeroMesh()\n");

	ave_len = 0.0;
	max_len = 0.0;
	min_len = HUGE;
	total_num_sides = 0.0;

	surf_tri_loop(surf,t)
	{
	    for (i = 0; i < 3; ++i)
	    {
		if (t->side_length0[i] != -1.0) continue;
		t->side_length0[i] = separation(Point_of_tri(t)[i],
			Point_of_tri(t)[(i+1)%3],3);
		for (j = 0; j < 3; ++j)
		{
		    t->side_dir0[i][j] = 
				(Coords(Point_of_tri(t)[(i+1)%3])[j] -
				 Coords(Point_of_tri(t)[i])[j])/
				 t->side_length0[i];
		}
		if (max_len < t->side_length0[i]) 
		    max_len = t->side_length0[i];
		if (min_len > t->side_length0[i])
		    min_len = t->side_length0[i];
		ave_len += t->side_length0[i];
		total_num_sides += 1.0;
	    }
	}
	never_redistribute(Hyper_surf(surf)) = YES;

	if (debugging("zero_mesh"))
	{
	    (void) printf("Leaving setSurfZeroMesh()\n");
	    if (total_num_sides == 0.0)
	    {
		(void) printf("Nothing done.\n");
	    }
	    else
	    {
	    	(void) printf("Equilibrium length:\n");
	    	(void) printf("total_num_sides = %d\n",(int)total_num_sides);
	    	(void) printf("min_len = %16.12f\n",min_len);
	    	(void) printf("max_len = %16.12f\n",max_len);
	    	(void) printf("ave_len = %16.12f\n",ave_len/total_num_sides);
	    }
	}
}	/* end setSurfZeroMesh */

static void setCurveZeroLength(
	CURVE *curve,
	double len_fac)
{
	int i,j;
	double max_len,min_len,ave_len,len;
	double total_num_sides;
	BOND *b;

	if (debugging("zero_mesh"))
	    (void) printf("Entering setCurveZeroLength()\n");

	ave_len = 0.0;
        max_len = 0.0;
        min_len = HUGE;
        total_num_sides = 0.0;

	curve_bond_loop(curve,b)
	{
	    if (b->length0 != -1.0) continue;
	    set_bond_length(b,3);
	    b->length0 = len_fac*bond_length(b);
	    for (i = 0; i < 3; ++i)
		b->dir0[i] = (Coords(b->end)[i] - Coords(b->start)[i])/
					b->length0;	
	    if (max_len < b->length0) max_len = b->length0;
	    if (min_len > b->length0) min_len = b->length0;
	    ave_len += b->length0;
	    total_num_sides += 1.0;
	}
	if (debugging("zero_mesh"))
	{
	    (void) printf("Leaving setCurveZeroLength()\n");
	    if (total_num_sides == 0.0)
	    {
		(void) printf("Nothing done.\n");
	    }
	    else
	    {
	    	(void) printf("Equilibrium length:\n");
	    	(void) printf("total_num_sides = %d\n",(int)total_num_sides);
	    	(void) printf("min_len = %16.12f\n",min_len);
	    	(void) printf("max_len = %16.12f\n",max_len);
	    	(void) printf("ave_len = %16.12f\n",ave_len/total_num_sides);
	    }
	}
}	/* end setCurveZeroLength */

static void setMonoCompBdryZeroLength(
	SURFACE *surf)
{
	CURVE **c;

	surf_pos_curve_loop(surf,c)
	{
	    if (hsbdry_type(*c) != MONO_COMP_HSBDRY)
		continue;
	    setCurveZeroLength(*c,1.0);
	}
	surf_neg_curve_loop(surf,c)
	{
	    if (hsbdry_type(*c) != MONO_COMP_HSBDRY)
		continue;
	    setCurveZeroLength(*c,1.0);
	}
}	/* end setMonoCompBdryZeroLength */

static void resetGoreBdryZerolength(
        SURFACE *surf)
{
	CURVE **c;
	BOND *b;
	surf_pos_curve_loop(surf,c)
        {
            if (hsbdry_type(*c) != GORE_HSBDRY)
                continue;
	    curve_bond_loop(*c,b)
	    {
		b->length0 = -1.0;
	    }
            setCurveZeroLength(*c,1.0);
        }
        surf_neg_curve_loop(surf,c)
        {
            if (hsbdry_type(*c) != GORE_HSBDRY)
                continue;
	    curve_bond_loop(*c,b)
	    {
		b->length0 = -1.0;
	    }
            setCurveZeroLength(*c,1.0);
        }
}       /* end resetGoreBdryZerolength */

static void connectStringtoRGB(
        Front *front,
        SURFACE *rg_surf,
	NODE **string_nodes,
	int num_strings)
{
	int i, j, k;
	int num;
	INTERFACE *intfc = front->interf;
	std::vector<POINT*> target;
	std::vector<CURVE*> delete_curves;
	std::vector<POINT*>::iterator it;
	NODE **rg_string_nodes, **n, *nload;
	CURVE **c, **string_curves, **rg_curves;
	BOND *b;
	AF_NODE_EXTRA *extra;
	double dist, dist1;
	double spacing,*h = computational_grid(intfc)->h;
	double coords[MAXD], dir[MAXD];
	int nb;
	AF_PARAMS *af_params = (AF_PARAMS*)front->extra2;

	if (debugging("trace"))
	    printf("Entering connectStringtoRGB() \n");
	linkSurfaceTriPoint(intfc, rg_surf);

	/* find and make rg_string_node */
	findPointsonRGB(front, rg_surf, target);
	num = target.size();
	FT_VectorMemoryAlloc((POINTER*)&rg_string_nodes, num, sizeof(NODE*));
	for (i = 0; i < num; ++i)
	{
	    rg_string_nodes[i] = make_node(target[i]);
	    FT_ScalarMemoryAlloc((POINTER*)&extra, sizeof(AF_NODE_EXTRA));
	    extra->af_node_type = RG_STRING_NODE;
	    rg_string_nodes[i]->extra = (POINTER)extra;
	    rg_string_nodes[i]->size_of_extra = sizeof(AF_NODE_EXTRA);
	}
	/* install a one-bond curve at each rg_string_node */
	FT_VectorMemoryAlloc((POINTER*)&rg_curves,num,sizeof(CURVE*));
	for (i = 0; i < num; ++i)
	{
	    TRI *t = target[i]->tris[0];
            int vertex = Vertex_of_point(t,target[i]);
	    for (int j = 0; j < target[i]->num_tris; ++j)
	    {
		t = target[i]->tris[j];
		vertex = Vertex_of_point(t,target[i]);
		if (Coords(Point_of_tri(t)[Next_m3(vertex)])[2] < 
		    Coords(target[i])[2] && 
		    Coords(Point_of_tri(t)[Prev_m3(vertex)])[2] <
		    Coords(target[i])[2])
		    break;
	    }
	    NODE *rg_node = make_node(Point_of_tri(t)[Next_m3(vertex)]);
            rg_curves[i] = make_curve(0,0,rg_string_nodes[i], rg_node);
            install_curve_in_surface_bdry(rg_surf,rg_curves[i],
					POSITIVE_ORIENTATION);
	    hsbdry_type(rg_curves[i]) = PASSIVE_HSBDRY;
	    linkCurveTriBond(rg_curves[i],rg_surf);
	}

	/* distinguish single parachute system from multi-parachute system */
	boolean multi_para = NO;
	int string_curve_onenode = 1; //test
	if (num_strings == 1)
	{
	    num_strings = num;
	    multi_para = YES;
	    string_curve_onenode = 10;
	}
	FT_VectorMemoryAlloc((POINTER*)&string_curves,num_strings,
						sizeof(CURVE*));
	for (k = 0; k < num_strings; ++k)
	{
	    NODE *start, *end;
	    if (multi_para == NO)
	    {
		start = string_nodes[k];
		end = rg_string_nodes[0];
		dist = distance_between_positions(Coords(start->posn),
					Coords(target[0]),3);
		for (i = 1; i < num; ++i)
		{
		    dist1 = distance_between_positions(Coords(start->posn),
					Coords(target[i]),3);
		    if (dist1 < dist)
		    {
			end = rg_string_nodes[i];
			dist = dist1;
		    }
		}
	    }
	    else
	    {
		start = *string_nodes;
		end = rg_string_nodes[k];
	    }
	    for (int l = 0; l < string_curve_onenode; ++l)
	    {
		string_curves[k] = make_curve(0,0,start,end);
		hsbdry_type(string_curves[k]) = STRING_HSBDRY;
		spacing = separation(start->posn,end->posn,3);
		for (j = 0; j < 3; ++j)
		    dir[j] = (Coords(end->posn)[j] - Coords(start->posn)[j])
							/ spacing;
		nb = (int)spacing/(0.40 * h[0]);
		spacing /= (double)nb;
		b = string_curves[k]->first;
		for (i = 1; i < nb; ++i)
		{
		    for (j = 0; j < 3; ++j)
			coords[j] = Coords(start->posn)[j] + i*dir[j]*spacing;
		    insert_point_in_bond(Point(coords),b,string_curves[k]);
		    b->length0 = spacing;
		    b = b->next;
		}
		b->length0 = spacing;
		af_params->string_curves.push_back(string_curves[k]);
	    }
	}
	FT_FreeThese(1,string_curves);

	if (debugging("trace"))
	    printf("Leaving connectStringtoRGB() \n");
}	/* end connectStringtoRGB */

static void findPointsonRGB(
	Front *front,
	SURFACE *rg_surf,
	std::vector<POINT*> & target)
{
	int i;
	TRI* tri;
	POINT *pt;
	std::vector<POINT*> candidate;
	std::vector<POINT*>::iterator it;
	INTERFACE *intfc = front->interf;
	double *h = computational_grid(intfc)->h;
	boolean find_x = YES, find_y = YES;
	FILE *infile = fopen(InName(front),"r");
	char string[100];

	if (debugging("trace"))
	    printf("Entering findPointsonRGB() \n");
	sprintf(string, "Enter yes to skip rg_nodes in x direction for %d:", 
				body_index(rg_surf));
	if (CursorAfterStringOpt(infile,string))
	{
	    fscanf(infile,"%s",string);
	    (void) printf("%s\n",string);
	    if (string[0] == 'y' || string[0] == 'Y')
		find_x = NO;
	}
	sprintf(string, "Enter yes to skip rg_nodes in y direction for %d:", 
				body_index(rg_surf));
	if (CursorAfterStringOpt(infile,string))
        {
            fscanf(infile,"%s",string);
            (void) printf("%s\n",string);
            if (string[0] == 'y' || string[0] == 'Y')
                find_y = NO;
        }
	fclose(infile);
	/* find points with the max z coordinate */
	surf_tri_loop(rg_surf, tri)
	{
	    for (i = 0; i < 3; ++i)
	    {
		pt = Point_of_tri(tri)[i];
		if (candidate.empty())
		    candidate.push_back(pt);
		else
		{
		    POINT* tmp = candidate.front();
		    if (fabs(Coords(pt)[2] - Coords(tmp)[2]) < 1.0e-6)
		    {
			if (std::find(candidate.begin(),candidate.end(),pt) 
				== candidate.end())
			    candidate.push_back(pt);
		    }
		    else if (Coords(pt)[2] > Coords(tmp)[2])
		    {
			candidate.clear();
			candidate.push_back(pt);
		    }
		}
	    }
	}
	/* for one-point cases, repeat finding with a lower z coordinate */
	if (candidate.size() == 1)
	{
	    double z_coord = Coords(candidate.front())[2] - 0.5*h[2];
	    candidate.clear();
	    surf_tri_loop(rg_surf, tri)
	    {
		for (i = 0; i < 3; ++i)
		{
		    pt = Point_of_tri(tri)[i];
		    if ( (fabs(Coords(pt)[2] - z_coord) < 0.1*h[2]) && 
			 (std::find(candidate.begin(),candidate.end(),pt)
					== candidate.end()) )
			candidate.push_back(pt);
		}
	    }
	}
	if (candidate.size() <= 4)
	{
	    target = candidate;
	    if (debugging("trace"))
		printf("Leaving findPointsonRGB() \n");
	    return;
	}

	std::vector<POINT*> x_max, x_min;
	x_max.push_back(candidate.front());
	x_min.push_back(candidate.front());
	std::vector<POINT*> y_max, y_min;
	y_max.push_back(candidate.front());
	y_min.push_back(candidate.front());
	POINT* tmp;
	for (it = candidate.begin() + 1; it != candidate.end(); ++it)
	{
	    /* find points with the max and min x coordinate */
	    tmp = x_max.front();
	    if (fabs(Coords(*it)[0] - Coords(tmp)[0]) < 1.0e-6)
		x_max.push_back(*it);
	    else if (Coords(*it)[0] > Coords(tmp)[0])
	    {
		x_max.clear();
		x_max.push_back(*it);
	    }
	    tmp = x_min.front();
	    if (fabs(Coords(*it)[0] - Coords(tmp)[0]) < 1.0e-6)
		x_min.push_back(*it);
	    else if (Coords(*it)[0] < Coords(tmp)[0])
	    {
		x_min.clear();
		x_min.push_back(*it);
	    }
	    /* find points with the max and min y coordinate */
	    tmp = y_max.front();
	    if (fabs(Coords(*it)[1] - Coords(tmp)[1]) < 1.0e-6)
		y_max.push_back(*it);
	    else if (Coords(*it)[1] > Coords(tmp)[1])
	    {
		y_max.clear();
		y_max.push_back(*it);
	    }
	    tmp = y_min.front();
	    if (fabs(Coords(*it)[1] - Coords(tmp)[1]) < 1.0e-6)
		y_min.push_back(*it);
	    else if (Coords(*it)[1] < Coords(tmp)[1])
	    {
		y_min.clear();
		y_min.push_back(*it);
	    }
	}
	candidate.clear();
	if (find_x == YES)
	{
	    if (x_max.size() == 1)
		candidate.push_back(x_max.front());
	    else
		vector_extreme(x_max, 1, candidate);
	    if (x_min.size() == 1)
		candidate.push_back(x_min.front());
	    else
		vector_extreme(x_min, 1, candidate);
	}
	if (find_y == YES)
	{
	    if (y_max.size() == 1)
		candidate.push_back(y_max.front());
	    else
		vector_extreme(y_max, 0, candidate);
	    if (y_min.size() == 1)
		candidate.push_back(y_min.front());
	    else
		vector_extreme(y_min, 0, candidate);
	}
	target = candidate;

	if (debugging("rigid_body"))
	{
	    for (i = 0; i < target.size(); ++i)
	    	printf("%d, %f %f %f\n", i, Coords(target[i])[0], 
			Coords(target[i])[1], Coords(target[i])[2]);
	}
	if (debugging("trace"))
	    printf("Leaving findPointsonRGB() \n");
}	/* end findPointsonRGB */

static void vector_extreme(
	std::vector<POINT*> &v,
	int dir,
	std::vector<POINT*> &candidate)
{
	std::vector<POINT*>::iterator it;
	POINT *p_max = v.front();
	POINT *p_min = v.front();
	for (it = v.begin() + 1; it != v.end(); ++it)
	{
	    if (Coords(*it)[dir] > Coords(p_max)[dir])
		p_max = *it;
	    if (Coords(*it)[dir] < Coords(p_min)[dir])
		p_min = *it;
	}
	if (std::find(candidate.begin(),candidate.end(),p_max)==candidate.end())
	    candidate.push_back(p_max);
	if (std::find(candidate.begin(),candidate.end(),p_min)==candidate.end())
	    candidate.push_back(p_min);
}	/* end vector_extreme */

static void linkSurfaceTriPoint(
        INTERFACE *intfc,
        SURFACE *surf)
{
        TRI *tri,**ptris;
        int i,j,num_point_tris;
        POINT **pts;
        std::vector<POINT*> pts_on_rgb;

        if (debugging("trace"))
            printf("Entering linkSurfaceTriPoint() \n");
        num_point_tris = 0;
        surf_tri_loop(surf,tri)
	{
	    pts = Point_of_tri(tri);
	    for (j = 0; j < 3; j++)
		pts[j]->num_tris = 0;
	}
        surf_tri_loop(surf,tri)
        {
            pts = Point_of_tri(tri);
            for (j = 0; j < 3; j++)
            {
                if (pts[j]->num_tris == 0)
                    pts_on_rgb.push_back(pts[j]);
                pts[j]->num_tris++;
            }
            num_point_tris += 3;
        }
        intfc->point_tri_store_rgb = (TRI**)store(num_point_tris*sizeof(TRI*));
        ptris = intfc->point_tri_store_rgb;
        for (i = 0; i < pts_on_rgb.size(); ++i)
        {
            pts_on_rgb[i]->tris = ptris;
            ptris += pts_on_rgb[i]->num_tris;
            pts_on_rgb[i]->num_tris = 0;
        }
        surf_tri_loop(surf,tri)
        {
            pts = Point_of_tri(tri);
            for (j = 0; j < 3; j++)
            {
                pts[j]->tris[pts[j]->num_tris++] = tri;
            }
        }
        if (debugging("trace"))
            printf("Leaving linkSurfaceTriPoint() \n");
        return;
}       /* end linkSurfaceTriPoint */

extern void InstallNewLoadNode(
	Front *front,
	int num_canopy)
{
	INTERFACE *intfc = front->interf;
	FILE *infile = fopen(InName(front),"r");
	NODE **n, *sec_nload, *nload;
	CURVE **string_curves;
	AF_NODE_EXTRA *extra;
	BOND *bond;
	double connection[MAXD],newload[MAXD],dir[MAXD],coords[MAXD];
	double spacing,*h = front->rect_grid->h;
	int i,j,k,nb;
 	INTERFACE *cur_intfc;
	AF_PARAMS *af_params = (AF_PARAMS*)front->extra2;
	int string_curve_onenode = 10; // test

	if (CursorAfterStringOpt(infile,"Enter new load position:"))
	{
	    fscanf(infile,"%lf %lf %lf",newload,newload+1,newload+2);
	    (void) printf("%f %f %f\n",newload[0],newload[1],newload[2]);
	}
	else
	{
	    fclose(infile);
	    return;
	}

	CursorAfterString(infile,"Enter connection position:");
	fscanf(infile,"%lf %lf %lf",connection,connection+1,connection+2);
	(void) printf("%f %f %f\n",connection[0],connection[1],connection[2]);

	cur_intfc = current_interface();
	set_current_interface(intfc);
	FT_VectorMemoryAlloc((POINTER*)&string_curves,num_canopy+1,
							sizeof(CURVE*));
	sec_nload = make_node(Point(connection));
	FT_ScalarMemoryAlloc((POINTER*)&extra,sizeof(AF_NODE_EXTRA));
	extra->af_node_type = SEC_LOAD_NODE;
	sec_nload->extra = (POINTER)extra;
	sec_nload->size_of_extra = sizeof(AF_NODE_EXTRA);

	i = 0;
	intfc_node_loop(intfc,n)
	{
	    extra = (AF_NODE_EXTRA*)((*n)->extra);
	    if (extra == NULL) continue;
	    if (extra->af_node_type != LOAD_NODE) continue;
	    extra->af_node_type = THR_LOAD_NODE;
	    for (int l = 0; l < string_curve_onenode; ++l)
	    {
		string_curves[i] = make_curve(0,0,(*n),sec_nload);
		hsbdry_type(string_curves[i]) = STRING_HSBDRY;
		spacing = separation((*n)->posn,sec_nload->posn,3);
		for (j = 0; j < 3; ++j)
		    dir[j] = (Coords(sec_nload->posn)[j] - 
					Coords((*n)->posn)[j])/spacing;
		nb = (int)spacing/(0.40*h[0]);
		spacing /= (double)nb;
		bond = string_curves[i]->first;
		for (j = 1; j < nb; ++j)
		{
		    for (k = 0; k < 3; ++k)
			coords[k] = Coords((*n)->posn)[k] + j*dir[k]*spacing;
		    insert_point_in_bond(Point(coords),bond,string_curves[i]);
		    bond->length0 = spacing;
		    bond = bond->next;
		}
		bond->length0 = spacing;
		af_params->string_curves.push_back(string_curves[i]);
	    }
	    i++;
	}

	nload = make_node(Point(newload));
	FT_ScalarMemoryAlloc((POINTER*)&extra,sizeof(AF_NODE_EXTRA));
	extra->af_node_type = LOAD_NODE;
	if(CursorAfterStringOpt(infile,"Enter yes to fix new load node:"))
        {
            char string[256];
            fscanf(infile,"%s",string);
            (void) printf("%s\n",string);
            if (string[0] == 'Y' || string[0] == 'y')
                extra->af_node_type = PRESET_NODE;
        }
	nload->extra = (POINTER)extra;
	nload->size_of_extra = sizeof(AF_NODE_EXTRA);

	if (CursorAfterStringOpt(infile,
			"Enter yes to install the multi-parachute to RGB:"))
	{
	    char string[100];
	    fscanf(infile,"%s",string);
	    (void) printf("%s\n",string);
	    if (string[0] == 'y' || string[0] == 'Y')
	    {
		int rg_index;
		SURFACE **rg_surf;
		CursorAfterString(infile,
				"Enter the body index of the target RGB:");
		fscanf(infile,"%d",&rg_index);
		(void) printf("%d\n",rg_index);
		intfc_surface_loop(intfc, rg_surf)
		{
		    if ((wave_type(*rg_surf) == NEUMANN_BOUNDARY) ||
			(wave_type(*rg_surf) == MOVABLE_BODY_BOUNDARY))
		    {
			if (body_index(*rg_surf) == rg_index)
			    break;
		    }
		}
		connectStringtoRGB(front,*rg_surf,&sec_nload,1);
		delete_node(nload);
		set_current_interface(cur_intfc);
		fclose(infile);
		FT_FreeThese(1,string_curves);
		return;
	    }
	}

	for (int l = 0; l < string_curve_onenode; ++l)
	{
	    string_curves[i] = make_curve(0,0,sec_nload,nload);
	    hsbdry_type(string_curves[i]) = STRING_HSBDRY;
	    spacing = separation(sec_nload->posn,nload->posn,3);
	    for (j = 0; j < 3; ++j)
		dir[j] = (Coords(nload->posn)[j] -
				Coords(sec_nload->posn)[j])/spacing;
	    nb = (int)spacing/(0.40*h[0]);
	    spacing /= (double)nb;
	    bond = string_curves[i]->first;
	    for (j = 1; j < nb; ++j)
	    {
		for (k = 0; k < 3; ++k)
		    coords[k] = Coords(sec_nload->posn)[k] + j*dir[k]*spacing;
		insert_point_in_bond(Point(coords),bond,string_curves[i]);
		bond->length0 = spacing;
		bond = bond->next;
	    }
	    bond->length0 = spacing;
	    af_params->string_curves.push_back(string_curves[i]);
	}
	set_current_interface(cur_intfc);
	fclose(infile);
	FT_FreeThese(1,string_curves);
}       /* end InstallNewLoadNode */

static void checkReducedTri(SURFACE* s)
{
	TRI *t;
	double tol = 1.0e-10;

	surf_tri_loop(s, t)
	{
	    double u[3], v[3];
	    for (int i = 0; i < 3; ++i)
	    {
		u[i] = Coords(Point_of_tri(t)[1])[i] - 
			Coords(Point_of_tri(t)[0])[i];
		v[i] = Coords(Point_of_tri(t)[2])[i] - 
			Coords(Point_of_tri(t)[0])[i];
	    }
	    double u_mag = Mag3d(u), v_mag = Mag3d(v);
	    if (u_mag < tol || v_mag < tol)
	    {
		printf("Error: reduced tri exists after initialization.\n");
		print_tri(t, s->interface);
		clean_up(ERROR);
	    }
	    for (int i = 0; i < 3; ++i)
	    {
		u[i] /= u_mag;
		v[i] /= v_mag;
	    }
	    double tmp[3];
	    Cross3d(u, v, tmp);
	    if (Mag3d(tmp) < tol)
	    {
		printf("Error: reduced tri exists after initialization.\n");
		print_tri(t, s->interface);
		clean_up(ERROR);
	    }
	}
}	/* end checkReducedTri */
