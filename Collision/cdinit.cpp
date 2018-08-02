#include <FronTier.h>
#include <fstream>
#include <sstream>
#include <string>
#include <map>
#include "collid.h"
static CURVE* make3dCurve(Front&,double[][3],int);
static void initSurface(Front&,char*);
static void initCurves(Front&,char*); 

void initTestModule(Front &front, char* inName)
{
	initSurface(front,inName);
	initCurves(front,inName);
	gview_plot_interface("init",front.interf);
}

std::map<std::string,int> hashMap({
	{"FIRST_PHYSICS_WAVE_TYPE",FIRST_PHYSICS_WAVE_TYPE},
	{"MOVABLE_BODY_BOUNDARY",MOVABLE_BODY_BOUNDARY},
	{"NEUMANN_BOUNDARY",NEUMANN_BOUNDARY},
	{"DIRICHLET_BOUNDARY",DIRICHLET_BOUNDARY},
	{"STRING_HSBDRY",STRING_HSBDRY}});

static void initSurface(Front &front,char* inName) {
	SURFACE *surf;
	std::ifstream fileStream(inName);
	std::string str, objName;
	int count = 0;
	while (!fileStream.eof()) {
	    getline(fileStream,str);
	    if (str.find("Sphere") != std::string::npos) {
		std::stringstream ss(str);
		double center[3], vel[3], radius[3];
		std::string wave_type;
		ss >> objName;
		std::cout << "#"<<count++<<" "<<objName << std::endl;
		for (int i = 0; i < 3; ++i)
		    ss >> center[i];
		std::cout << "    Center: "<< center[0] <<" " << center[1] 
			  << " " << center[2]<< std::endl;

		for (int i = 0; i < 3; ++i)
		    ss >> radius[i];
		std::cout << "    Radius: "<< radius[0]<<" "<<radius[1]
			  <<" "<<radius[2]<<std::endl;

		for (int i = 0; i < 3; ++i) 
		    ss >> vel[i];
		std::cout << "    Vel: "<< vel[0] << " " << vel[1] 
			  << " " << vel[2] << std::endl;

		ss >> wave_type;
		if (hashMap.find(wave_type) == hashMap.end())
		    std::cout << "Unknown supported type: "<<wave_type<<std::endl;
		std::cout << "    Wave type: " << wave_type << std::endl; 
		FT_MakeEllipticSurf(&front,center,radius,
                        1,2,    // negative and positive components
                        hashMap[wave_type],
                        1,      // refinement level
                        &surf);
	        initSurfaceState(surf,vel);
	    }
	    if (str.find("Cuboid") != std::string::npos) {
		std::stringstream ss(str);
                double center[3], vel[3], edge[3];
                std::string wave_type;
                ss >> objName;
                std::cout << "#"<<count++<<" "<<objName << std::endl;
                for (int i = 0; i < 3; ++i)
                    ss >> center[i];
                std::cout << "    Center: "<< center[0] <<" " << center[1] << " " 
			  << center[2]<< std::endl;

                for (int i = 0; i < 3; ++i)
                    ss >> edge[i];
                std::cout << "    Edge: "<< edge[0]<<" "<<edge[1]<<" "
			  << edge[2]<< std::endl;

                for (int i = 0; i < 3; ++i)
                    ss >> vel[i];
                std::cout << "    Vel: "<< vel[0] << " " << vel[1] << " " << vel[2]
                << std::endl;

                ss >> wave_type;
                if (hashMap.find(wave_type) == hashMap.end())
                    std::cout << "Unknown supported type: "<<wave_type<<std::endl;
                std::cout << "    Wave type: " << wave_type << std::endl; 
		FT_MakeCuboidSurf(&front,center,edge,1,2,
                        hashMap[wave_type],
                        1,
                        &surf);
	        initSurfaceState(surf,vel);
	    }
	}
	fileStream.close();
}

static void initCurves(Front &front,char* inName) {
	CURVE* curve;
	std::ifstream fileStream(inName);
        std::string str;
        int count = 0;
        while (!fileStream.eof()) {
            getline(fileStream,str);
            if (str.find("Curve") != std::string::npos) {
                std::stringstream ss(str);
		std::string objName, hsbdry_type;
		double endPoint[2][3], vel[3];
		ss >> objName;
		std::cout <<"#"<< count++ << " " << objName <<std::endl;

		for (int i = 0; i < 3; ++i)
		    ss >> endPoint[0][i];
		std::cout << "    Start: "<<endPoint[0][0]<<" "<<endPoint[0][1]<<" "
			  <<endPoint[0][2] << std::endl;
		
		for (int i = 0; i < 3; ++i)
		    ss >> endPoint[1][i];
		std::cout << "    End: "<< endPoint[1][0] << " " << endPoint[1][1] << " "
			  << endPoint[1][2] << std::endl;

		for (int i = 0; i < 3; ++i)
		    ss >> vel[i];
		std::cout << "    Vel: " << vel[0] << " " << vel[1] << " "
			  << vel[2]<< std::endl;
		
		ss >> hsbdry_type;
		std::cout << "    Wave type: " << hsbdry_type << std::endl;
		if (hashMap.find(hsbdry_type) == hashMap.end())
		    std::cout << "Unknown hypersurface boundary type" << std::endl;
		else
		{
		    curve = make3dCurve(front,endPoint,hashMap[hsbdry_type]);
        	    initCurveState(curve,vel);
		}
	    }
	    
	}
	fileStream.close();
}

static CURVE* make3dCurve(Front& front,double pt[][3],int hsb_type){
	NODE      *string_nodes[2];
        CURVE     *curve;
        INTERFACE *intfc = front.interf;
        for (int i = 0; i < 2; ++i)
            string_nodes[i] = make_node(Point(pt[i]));
        curve = make_curve(0,0,string_nodes[0],string_nodes[1]);
        hsbdry_type(curve) = hsb_type;
        double spacing = separation(string_nodes[0]->posn,
                                    string_nodes[1]->posn,3);
        double dir[3];
        for (int j = 0; j < 3; ++j)
            dir[j] = (Coords(string_nodes[1]->posn)[j] -
                      Coords(string_nodes[0]->posn)[j])/spacing;
        double *h = computational_grid(intfc)->h;
        int nb = (int)(spacing/(0.25*h[0]));
        spacing /= (double)nb;  
        BOND* b = curve->first; 
        double coords[3];       
        for (int j = 1; j < nb; ++j)
        {
            for (int k = 0; k < 3; ++k)
                coords[k] = Coords(string_nodes[0]->posn)[k] +
                                   j*dir[k]*spacing;
            insert_point_in_bond(Point(coords),b,curve);
	    b->length0 = spacing;
            b = b->next;
        }
	b->length0 = spacing;
	return curve;
}

void initSurfaceState(
        SURFACE* surf,
        const double* vel)
{
        TRI* tri;
        POINT* p;
        STATE* sl, *sr;
	if (!surf) return; 
        surf_tri_loop(surf,tri)
        {
            for (int i = 0; i < 3; ++i)
            {
                p = Point_of_tri(tri)[i];
                sl = (STATE*)left_state(p);
                sr = (STATE*)right_state(p);
                for (int j = 0; j < 3; ++j)
                {
                    sl->vel[j] = sr->vel[j] = vel[j];
                    sl->collsnImpulse[j] = sr->collsnImpulse[j] = 0.0;
                    sl->collsnImpulse_RG[j] = sr->collsnImpulse_RG[j] = 0.0;
                    sl->friction[j] = sr->friction[j] = 0.0;
                    sl->collsn_num = sr->collsn_num = 0;
                    sl->collsn_num_RG = sr->collsn_num_RG = 0;
                    sl->x_old[j] = Coords(p)[j];
		    sl->is_fixed = false;
		    sl->is_movableRG = false;
		    if (wave_type(p->hs) == NEUMANN_BOUNDARY)
			sl->is_fixed = true;
		    if (wave_type(p->hs) == MOVABLE_BODY_BOUNDARY)
			sl->is_movableRG = true;
                }
            }
        }
}

void initCurveState(CURVE* curve, const double* vel)
{
	BOND* b;
	STATE* sl;
	POINT* p[2];
	if (!curve) return;
	curve_bond_loop(curve,b)
	{
	    p[0] = b->start;
	    p[1] = b->end;
	    for (int i = 0; i < 2; ++i)
	    {
		sl = (STATE*)left_state(p[i]);
		for (int j = 0; j < 3; ++j)
                {
                    sl->vel[j] = vel[j];
                    sl->collsnImpulse[j] = 0.0;
                    sl->collsnImpulse_RG[j] = 0.0;
                    sl->friction[j] = 0.0;
                    sl->collsn_num = 0;
                    sl->collsn_num_RG = 0;
                    sl->x_old[j] = Coords(p[i])[j];
                }	
	    }
	}
}
