#include "folding.h"
#include "spring_solver.h"
#include <fstream>
#include <sstream>
#include <map>
#include <iomanip>
#include <cstdlib>

static void setCollisionFreePoints3d(INTERFACE*, Drag*);
static bool string_double(const std::string& s); 

bool findAndLocate(std::ifstream&, const char*);
bool getString(std::ifstream&, const char*);

double Folder::m_thickness = 0.001;
double Folder::max_dt = 0.01;

bool string_double(const std::string& s) {
    int numD = 0;

    for (int i = 0; i < s.size(); i++) {
         if (i == 0 && s[i] == '-') continue;
         if (s[i] < '0' || s[i] > '9')
             if (s[i] == '.' && numD < 1) {
                 numD++;
                 continue;
             }
             else return false;
        else continue;
    }
    return true;

}

void Folder::addDragsFromFile(std::string fname) {
    std::ifstream ifs(fname.c_str());
    if (!ifs.is_open()) {
        std::cout << "File \"" << fname <<
        "\" doesn't exist" << std::endl;
        return;
    }

    Drag::Info info;
    std::set<std::string> foldset;

    for (std::vector<Drag*>::iterator it = Drag::prototypes.begin();
                it != Drag::prototypes.end(); it++)
         foldset.insert((*it) -> id());

    std::string s;
    info.clear();
    while (ifs >> s)
    {
	//three cases: drag name, numeric, garbage
	if (find(foldset.begin(), foldset.end(), s) != foldset.end())
        {
	    // new drag found, create drag using info
	    // begin to collect information for new drag
	    if (info.data().size() != 0)
             {
                addDrag(Drag::dragFactory(info));
             }
             info.clear();
             info.id() = s;
        }
        else
        {
	    //check if it is a double
             if (!string_double(s)) continue;

	     std::stringstream ss(s);
	     double temp; 
	     ss >> temp;  

             info.data().push_back(temp);
        }
    }
    //last drag
    if (info.data().size() != 0)
        addDrag(Drag::dragFactory(info));
}

void Folder::setSpringParameters(double k, double lambda, double m) {
    spring_params.k = k;
    spring_params.lambda = lambda;
    spring_params.m = m;
}

void Folder::setParaFromFile(const char* s)
{
    std::ifstream fin(s);

    if (!fin.is_open())
    {
	std::cerr << "Can't open file " << s << std::endl; 
	clean_up(ERROR);
    }
    // we may use the default value
    // so it's ok not find these strings
    if (findAndLocate(fin, "fabric spring constant:"))
    {
	fin >> spring_params.k; 
    	std::cout << spring_params.k << std::endl; 
    }
    else
	std::cout << "use default value!\n"; 
    if (findAndLocate(fin, "fabric friction constant:"))
    {
	fin >> spring_params.lambda;
    	std::cout << spring_params.lambda << std::endl;
    }
    else
        std::cout << "use default value!\n";
    if (findAndLocate(fin, "fabric point mass:"))
    {
	fin >> spring_params.m;
        std::cout << spring_params.m << std::endl;
    }
    else
        std::cout << "use default value!\n";
    if (findAndLocate(fin, "frame step size:"))
    {
	fin >> max_dt;
        std::cout << max_dt << std::endl;
    }
    else
        std::cout << "use default value!\n";

    fin.close(); 
}

void Folder3d::doFolding() {
    SpringSolver* sp_solver = SpringSolver::createSpringSolver(
						getOdeScheme()); 
    CollisionSolver* cd_solver = new CollisionSolver3d();
    
    //configure collision solver
    cd_solver->assembleFromInterface(m_intfc,getFrameStepSize());
    cd_solver->setSpringConstant(getSpringParams().k);
    cd_solver->setFrictionConstant(getSpringParams().lambda);
    cd_solver->setPointMass(getSpringParams().m);
    cd_solver->setFabricThickness(getThickness());
  
    //configure spring solver
    FT_Intfc2SpringMesh(m_intfc, sp_solver->getSpringMesh());
    if (getSpringParams().k == 0 || getSpringParams().m == 0)
	throw std::invalid_argument("tensile stiffness and mass cannot zero");
    sp_solver->setParameters(getSpringParams());
    
    //default bend coefficient and bend damping coefficient is 
    //0.01 and 0.0
    BendingForce* btemp = new BendingForce(m_intfc); 
    std::string temp(getInputFile()); 

    btemp->getParaFromFile(temp.c_str()); 
    sp_solver->ext_forces.push_back(btemp);
    
    Drag::setTolerance(m_intfc->table->rect_grid.h[0]*0.5);
    Drag::setThickness(0.001);

    for (std::vector<Drag*>::iterator it = drags.begin();
	 it != drags.end(); ++it) 
    {
	if (*it == NULL) continue;
	printf("folding base on %lu drag: %s ...\n",
		it-drags.begin(),(*it)->id().c_str());
	doFolding(*it,sp_solver,cd_solver);
	sp_solver->resetVelocity();
    }
    straightenStrings();
}

double Folder3d::computePotentialEnergy()
{
    SURFACE** s;
    double E = 0;
    intfc_surface_loop(m_intfc, s)
    {
	if (wave_type(*s) != ELASTIC_BOUNDARY)
	    continue;
	TRI* t;
	surf_tri_loop(*s, t)
	{
	    for (int i = 0; i < 3; ++i)
	    {
		double len = separation(Point_of_tri(t)[i],
				Point_of_tri(t)[(i+1)%3], 3); 
		E += 0.5*spring_params.k*pow(len - t->side_length0[i], 2.0);	
	    }
	}
    } 
    return E*0.5;
}

void Folder3d::doFolding(
     Drag* drag, 
     SpringSolver* sp_solver,
     CollisionSolver* cd_solver) 
{
    double min_dt = sp_solver->getTimeStepSize();
    double max_dt = std::max(min_dt,getFrameStepSize());
    
    sp_solver->setTimeStepSize(max_dt);
    sp_solver->setDrag(drag);

    setCollisionFreePoints3d(m_intfc,drag);

    //call spring solver and collision solver
    //based on adaptive dt
    // t0 is the time control for one drag plan
    double t0 = 0; 
    static double t = 0;
    double dt = max_dt;
    printf("drag = %s, m_t = %f\n", drag->id().c_str(), drag->m_t);
    movie->recordMovieFrame();
    while (t0 < drag->m_t - MACH_EPS) {
	printf("--------------------------------\n");
	printf("dt = %f, t = %f, total time = %f\n",
		dt,t0,drag->m_t);
	printf("--------------------------------\n");
	if (t0+dt > drag->m_t) 
	    dt = drag->m_t-t0;
	
	cd_solver->recordOriginPosition();
	cd_solver->setTimeStepSize(dt);

    	sp_solver->solve(dt);
	
	recordData(t,movie->out_dir);

	//cd_solver->resolveCollision();

	t += dt;
        t0 += dt; 
	if (movie->isMovieTime(t))
	    movie->recordMovieFrame();
    }
    movie->recordMovieFrame();
}

// friend function to Folder 
// used to check force for the unregistered point whose other two points on 
// the same triangle are registered points
void Folder3d::check_force(SpringSolver* sp_solver)
{
    std::vector<SpringVertex*> pts = sp_solver->getSpringMesh(); 
    std::map<POINT*, SpringVertex*> pToVer; 
    int i;  

    for (size_t j = 0; j < pts.size(); j++)
    {
	 pToVer[(POINT*)(pts[j]->org_vtx)] = pts[j]; 
    }

    SURFACE ** surf; 
    std::ofstream fout1("forceregis.txt"); 
    std::ofstream fout2("force.txt");

    if (!fout1.is_open() || !fout2.is_open())
    {
	std::cerr << "Open file error!\n"; 
	exit(EXIT_FAILURE);
    }
    
    fout1 << std::fixed << std::setprecision(14); 
    fout2 << std::fixed << std::setprecision(14);

    intfc_surface_loop(m_intfc, surf)
    { 
	if (wave_type(*surf) != ELASTIC_BOUNDARY) continue; 
	if (is_bdry(*surf)) continue; 

        TRI *tri;
	     
	surf_tri_loop(*surf, tri)
	    for (i = 0; i < 3; i++)
		 sorted(Point_of_tri(tri)[i]) = NO; 
	surf_tri_loop(*surf, tri)
	{
	    for (i = 0; i < 3; i++)
	    {
		 POINT *p = Point_of_tri(tri)[i];
                 SpringVertex *vp = pToVer[p];
 		 
		 if (sorted(p) || vp->isRegistered()) continue;
                 fout2 << "Point is: " << std::setw(20) << Coords(p)[0]
                        << std::setw(20) << Coords(p)[1] << std::setw(20)
                        << Coords(p)[2] << std::endl;
                     fout2 << "Force is: " << std::setw(20) << p->force[0]
                        << std::setw(20) << p->force[1] << std::setw(20)
                        << p->force[2] << std::endl;
		 
		 int count = 0;
		 std::vector<size_t> indexNb = vp->index_nb;  
		 
		 for (size_t j = 0; j < indexNb.size(); j++)
		 {
		      if (pts[indexNb[j]]->isRegistered())
			  count++; 
		 }
		 if (count == 2)
		 {
		     fout1 << "Point is: " << std::setw(20) << Coords(p)[0] 
			<< std::setw(20) << Coords(p)[1] << std::setw(20) 
			<< Coords(p)[2] << std::endl; 
		     fout1 << "Force is: " << std::setw(20) << p->force[0]
			<< std::setw(20) << p->force[1] << std::setw(20) 
			<< p->force[2] << std::endl; 
		 }
	    }
	}
    }
    fout1.close(); 
    fout2.close();
}

void Folder3d::appendDataToFile(double x, double y, std::string fname)
{
    //if new file, clear it and add to list
    std::ofstream ofs;
    if (dataFileSet.find(fname) == dataFileSet.end())
    {
	dataFileSet.insert(fname);
	ofs.open(fname.c_str(), std::ofstream::trunc);
	ofs.close();
    }
    ofs.open(fname.c_str(), std::ofstream::app);
    ofs << x << " " << y << std::endl;
    ofs.close();
}

void Folder3d::recordData(double t, std::string out_dir)
{
    appendDataToFile(t, computePotentialEnergy(), out_dir+"/"+"energy");
}

void Folder3d::setupMovie(std::string dname, std::string oname, 
			double dt) {
    movie = new Movie();
    movie->mv_dt = dt;
    std::string pathgv = oname + '/' + dname + "/gview";
    std::string pathvtk = oname + '/' + dname + "/vtk"; 
    movie->mv_gv_dir = pathgv;
    movie->mv_vtk_dir = pathvtk; 
    movie->out_dir = oname;
    movie->doMakeMovie = true;
    movie->mv_intfc = this->m_intfc;
    std::string sys_cmd = "mkdir -p " + pathgv;
    system(sys_cmd.c_str());
    sys_cmd = "mkdir -p " + pathvtk; 
    system(sys_cmd.c_str());
}

void Movie::recordMovieFrame() {
    std::stringstream ss; 
    std::string count; 
    
    ss << std::setw(6) << std::setfill('0') << mv_count; 
    count = ss.str(); 

    std::string fname = mv_gv_dir + "/intfc-" + count;
    gview_plot_interface(fname.c_str(),mv_intfc);

    // keep marking number to have fixed number of digits
    fname = mv_vtk_dir + "/vtk-" + count;
    std::string sys_cmd = "mkdir -p " + fname; 
    system(sys_cmd.c_str());
    vtk_interface_plot(fname.c_str(),mv_intfc, YES,
                                        mv_dt * mv_count ,mv_count);
    mv_count++; 
}

bool Movie::isMovieTime(double t) {
    if (doMakeMovie)
        return fabs(t-mv_count*mv_dt) < 1e-10;
    else
	return false;
}

static void setCollisionFreePoints3d(INTERFACE* intfc, Drag* drag)
{
        POINT *p;
        HYPER_SURF *hs;
        HYPER_SURF_ELEMENT *hse;
        assert (intfc->dim == 3);

        next_point(intfc,NULL,NULL,NULL);
        while(next_point(intfc,&p,&hse,&hs)){
            STATE* sl = (STATE*)left_state(p);
            sl->is_fixed = false;
            if (wave_type(hs) == NEUMANN_BOUNDARY ||
                wave_type(hs) == MOVABLE_BODY_BOUNDARY)
            {
                sl->is_fixed = true;
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
}       /* setCollisionFreePoints3d() */

Folder3d::Folder3d(INTERFACE* intfc, SURFACE* s) : m_intfc(intfc)
{}

Folder3d::~Folder3d() {
    if (movie) delete movie;
}

static int numConnectCurves(NODE* n) {
    int ans = 0;
    for (CURVE** c = n->in_curves; c && *c; ++c)
	ans++;

    for (CURVE** c = n->out_curves; c && *c; ++c)
	ans++;
    return ans;
}

void Folder3d::straightenStrings() {
    std::cout << "Straighten the string" << std::endl;
    //straighten the string curves
    CURVE** c;
    double z_avg = 0;
    int c_count = 1;
    intfc_curve_loop(m_intfc,c) {
        if (hsbdry_type(*c) != STRING_HSBDRY) continue;
	//count how many bonds in a curve
        int count = 0;
        BOND* b;
	curve_bond_loop(*c,b)
	    count++;

	//adjust load node
	POINT* p1, *p2;
	p1 = (*c)->start->posn;
	p2 = (*c)->end->posn;
	if (numConnectCurves((*c)->start) < numConnectCurves((*c)->end))
	    std::swap(p1,p2); //make sure p1 is load node;

	double l = (*c)->first->length0 * count;
	double new_z  =  Coords(p2)[2] - sqrt(sqr(l) 
		       - sqr(Coords(p1)[0]-Coords(p2)[0]) 
		       - sqr(Coords(p1)[1]-Coords(p2)[1]));
	z_avg = z_avg + (new_z - z_avg)/(c_count++);
	Coords(p1)[2] = z_avg;
    }

    intfc_curve_loop(m_intfc,c) {
        if (hsbdry_type(*c) != STRING_HSBDRY) continue;
	//count how many bonds in a curve
        int count = 0;
        BOND* b;
	curve_bond_loop(*c,b)
	    count++;

        double dir[3] = {0};
        for (size_t i = 0; i < 3; ++i)
            dir[i] = Coords((*c)->end->posn)[i] -
                     Coords((*c)->start->posn)[i];
        //normalize dir
        double len = Mag3d(dir);
        for (size_t i = 0; i < 3; ++i)
            dir[i] /= len;

        //move bonds
        double len0 = len/count;
        count = 1;
        for (b = (*c)->first; b != (*c)->last; b = b->next) {
            POINT* p = b->end;
            for (size_t i = 0; i < 3; ++i)
                Coords(p)[i] = Coords((*c)->start->posn)[i]
                             + (count)*len0 * dir[i];
            count++;
        }
   }
}
