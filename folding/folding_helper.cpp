#include <FronTier.h>
#include "spring_solver.h"
#include <map>
#include "../iFluid/ifluid_state.h"
#include "folding.h"

//helper functions to convert FronTier data structures
//to spring mesh
static void unsortIntfc(INTERFACE*);
static void unsortSurface(SURFACE*);
static void unsortCurve(CURVE*);

static void assemblePointsFromSurf(SURFACE*, std::vector<SpringVertex*>&,
				std::map<POINT*, size_t>&);
static void assemblePointsFromCurve(CURVE*, std::vector<SpringVertex*>& pts,
				std::map<POINT*, size_t>&);
static void assemblePointsFromNode(NODE*, std::vector<SpringVertex*>& pts,
				std::map<POINT*, size_t>&);

static void setConnection(POINT* p1, POINT* p2, double length0, std::vector<SpringVertex*>& pts, std::map<POINT*, size_t>&);

static  SpringVertex* FT_Point2SpringVertex(POINT*);

static bool isElasticCurve(CURVE* c) {
    return (hsbdry_type(c) == STRING_HSBDRY ||
            hsbdry_type(c) == MONO_COMP_HSBDRY ||
            hsbdry_type(c) == GORE_HSBDRY);
}

static SpringVertex* FT_Point2SpringVertex(POINT* p) {
    SpringVertex* spv = new SpringVertex();
    STATE* state = (STATE*)left_state(p);
    for (int i = 0; i < 3; ++i)
        spv->ext_accel[i] = state->fluid_accel[i] +
                   state->other_accel[i];
    spv->x = Coords(p);
    spv->v = state->vel;
    spv->org_vtx = (void*)p;
    return spv;
}

extern void FT_Intfc2SpringMesh(INTERFACE* intfc, std::vector<SpringVertex*>&pts) {
    //assemble tris list from input intfc
    //this function should be called before
    //spring interior dynamics computed

    std::map<POINT*, size_t> ht;
    //clear hash table and point list
    ht.clear();
    pts.clear();
    unsortIntfc(intfc);

    //clone points from surface
    SURFACE** s;
    intfc_surface_loop(intfc,s) {
        if (wave_type(*s) == ELASTIC_BOUNDARY)
            assemblePointsFromSurf(*s, pts, ht);
    }
    printf("#%lu num of points in surf\n",pts.size());

    //clone points from curve
    CURVE** c;
    intfc_curve_loop(intfc,c) {
        //do not assemble strings
        if (hsbdry_type(*c) == STRING_HSBDRY) continue;
        if (isElasticCurve(*c))
            assemblePointsFromCurve(*c, pts, ht);
    }
    printf("#%lu num of points in curve\n",pts.size());

    //clone points from node
    NODE** n;
    intfc_node_loop(intfc,n) {
        assemblePointsFromNode(*n, pts, ht);
    }
    printf("#%lu num of points in nodes\n",pts.size());
}

static void assemblePointsFromSurf(SURFACE* s, 
			std::vector<SpringVertex*>& pts,
			std::map<POINT*, size_t>& ht) {
        TRI* tris[20];
        TRI* tri;
        surf_tri_loop(s,tri)
        for (int i = 0; i < 3; ++i)
        {
            POINT* p = Point_of_tri(tri)[i];
            if(Boundary_point(p) || sorted(p)) continue;
            sorted(p) = YES;

            int nt = 0;
            PointAndFirstRingTris(p,Hyper_surf_element(tri),
                            Hyper_surf(s),&nt,tris);
            for (int j = 0; j < nt; ++j)
            for (int l = 0; l < 3; ++l)
            if (Point_of_tri(tris[j])[l] == p)
            {
                POINT* p_nb = Point_of_tri(tris[j])[(l+1)%3];
                setConnection(p,p_nb,tris[j]->side_length0[l], pts, ht);
            }
        }
}

static void assemblePointsFromCurve(CURVE *curve, 
				std::vector<SpringVertex*>& pts,
				std::map<POINT*, size_t>& ht) {
 	//assemble from curve
	for (BOND* b = curve->first; b != curve->last; b = b->next)
        for (int i = 0; i < 2; ++i)
        {
            POINT* p = b->end;
            POINT* p_nb = (i == 0) ? b->start : b->next->end;
            double dist = (i == 0) ? b->length0 : b->next->length0;
	    setConnection(p,p_nb,dist, pts, ht);
        }

	//connect with btris 
	TRI** tris;
	for (BOND* b = curve->first; b != curve->last; b = b->next)
	{
	    POINT* p = b->end;
	    for (BOND_TRI** btris = Btris(b); btris && *btris; ++btris)
            {
                int nt = I_FirstRingTrisAroundPoint(p,(*btris)->tri,&tris);
		for (int j = 0; j < nt; ++j)
                {
                    for (int side = 0; side < 3; ++side)
                    {
                        if (p == Point_of_tri(tris[j])[side])
                        {
                            if (is_side_bdry(tris[j],side))
                                continue;
                            POINT* p_nb = Point_of_tri(tris[j])[(side+1)%3];
			    setConnection(p,p_nb,tris[j]->side_length0[side], pts, ht);
                        }
                    }
                }
	    }
	}
}

static void assemblePointsFromNode(NODE* n, 
			std::vector<SpringVertex*>& pts,
			std::map<POINT*, size_t>& ht) {
    for (CURVE** c = n->in_curves; c && *c; ++c)
    {
        if (!isElasticCurve(*c)) continue;
        //do not fold strings
        if (hsbdry_type(*c) == STRING_HSBDRY) continue;
        BOND* b = (*c)->last;
        POINT* p_nb = b->start;
        setConnection(n->posn,p_nb,b->length0, pts, ht);
    }

    for (CURVE** c = n->out_curves; c && *c; ++c)
    {
        if (!isElasticCurve(*c)) continue;
        //do not fold strings
        if (hsbdry_type(*c) == STRING_HSBDRY) continue;
        BOND* b = (*c)->first;
        POINT* p_nb = b->end;
        setConnection(n->posn,p_nb,b->length0, pts, ht);
    }
}

static void setConnection(POINT* p1, POINT* p2, double length0, 
			std::vector<SpringVertex*>& pts,
			std::map<POINT*, size_t>& ht) {
    POINT* points[2] = {p1,p2};
    for (int i = 0; i < 2; ++i) {
        if (ht.find(points[i]) == ht.end()) {
            //not inserted yet
            ht[points[i]] = pts.size();
            pts.push_back(FT_Point2SpringVertex(points[i]));
        }
    }
    if (length0 < 0)
        length0 = separation(p1,p2,3);
    size_t ind = ht[p1];
    size_t ind_nb = ht[p2];
    pts[ind]->addNeighbor(ind_nb,length0);
}

static void unsortIntfc(INTERFACE* intfc) {
    SURFACE** s;
    intfc_surface_loop(intfc,s)
        unsortSurface(*s);
    CURVE** c;
    intfc_curve_loop(intfc,c)
        unsortCurve(*c);
}

static void unsortSurface(SURFACE* s) {
    TRI* t;
    surf_tri_loop(s,t) {
        for (int i = 0; i < 3; ++i) {
            POINT* p = Point_of_tri(t)[i];
            sorted(p) = NO;
        }
    }
}

static void unsortCurve(CURVE* c) {
    BOND* b;
    curve_bond_loop(c,b) {
        sorted(b->start) = NO;
        sorted(b->end) = NO;
    }
}
