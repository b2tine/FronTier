#include <FronTier.h>
#include <vector>
#include "collid.h"

//test module for 3d surface
//proximity and collision detection
char *in_name,*restart_state_name,*restart_name,*out_name;
boolean RestartRun;
int RestartStep;
static void propagation_driver(Front*);
static void dummySpringSolver(Front*);
static void collision_point_propagate(Front*,POINTER,POINT*,
        			      POINT *newp,HYPER_SURF_ELEMENT *,
        			      HYPER_SURF*,double,double*);
static void collision_curve_propagate(Front*,POINTER,CURVE*,CURVE*,double);

int main(int argc, char** argv)
{
	static Front front;
        static F_BASIC_DATA f_basic;
        static LEVEL_FUNC_PACK level_func_pack;
	f_basic.dim = 3;
        FT_Init(argc,argv,&f_basic);
        f_basic.size_of_intfc_state = sizeof(STATE);

        /* Initialize basic computational data */
        in_name                 = f_basic.in_name;
        restart_state_name      = f_basic.restart_state_name;
        out_name                = f_basic.out_name;
        restart_name            = f_basic.restart_name;
        RestartRun              = f_basic.RestartRun;
        RestartStep             = f_basic.RestartStep;
	//initialize interface and velocity
	FT_ReadSpaceDomain(in_name,&f_basic);
        FT_StartUp(&front,&f_basic);
        FT_InitDebug(in_name);

	level_func_pack.pos_component = 2;
        FT_InitIntfc(&front,&level_func_pack);

	FT_ReadTimeControl(in_name,&front);
	
	//Custom function for test
	initTestModule(front,in_name);

	front.vfunc = NULL;
	//PointPropagationFunction(&front) = fourth_order_point_propagate;
	PointPropagationFunction(&front) = collision_point_propagate;
	front.curve_propagate = collision_curve_propagate;
	char dname[256];
	sprintf(dname,"%s/intfc",OutName(&front));
	geomview_interface_plot(dname,front.interf,front.rect_grid);

	propagation_driver(&front);
	clean_up(0);
}

static  void propagation_driver(
        Front *front)
{
        double CFL;
	CollisionSolver *collision_solver = new CollisionSolver3d();

        CFL = Time_step_factor(front);
	Tracking_algorithm(front) = STRUCTURE_TRACKING;

	Frequency_of_redistribution(front,GENERAL_WAVE) = 100000;
        printf("CFL = %f\n",CFL);
        printf("Frequency_of_redistribution(front,GENERAL_WAVE) = %d\n",
                Frequency_of_redistribution(front,GENERAL_WAVE));

        if (!RestartRun)
        {
            FT_RedistMesh(front);
            FT_ResetTime(front);

            // Always output the initial interface.
            FT_Save(front);
            FT_Draw(front);

            // This is a virtual propagation to get maximum front 
            // speed to determine the first time step.

            FT_Propagate(front);
            FT_SetTimeStep(front);
            FT_SetOutputCounter(front);
        }
        else
        {
            FT_SetOutputCounter(front);
        }

        FT_TimeControlFilter(front);

        FT_PrintTimeStamp(front);
        for (;;)
        {
            /* Propagating interface for time step dt */


            FT_Propagate(front);

	    dummySpringSolver(front);
	    //collision detect and handling
	    collision_solver->assembleFromInterface(front->interf,front->dt);
	    collision_solver->setFrictionConstant(0.0);
	    collision_solver->resolveCollision();

            FT_AddTimeStepToCounter(front);

            //Next time step determined by maximum speed of previous
            //step, assuming the propagation is hyperbolic and
            //is not dependent on second order derivatives of
            //the interface such as curvature, and etc.

            FT_SetTimeStep(front);

            /* Output section */

            FT_TimeControlFilter(front);
            FT_PrintTimeStamp(front);
            if (FT_IsSaveTime(front))
                FT_Save(front);
            if (FT_IsDrawTime(front))
                FT_Draw(front);

            if (FT_TimeLimitReached(front))
                    break;

        }
	delete collision_solver;
}       /* end propagation_driver */

static void collision_curve_propagate(
	Front* front,
	POINTER wave,
	CURVE* oldc,
	CURVE* newc,
	double dt)
{
	BOND *oldb,*newb;
        POINT *oldp,*newp;
	int dim = 3;
	oldp = oldc->start->posn;
        newp = newc->start->posn;
	ft_assign(left_state(newp),left_state(oldp),front->sizest);
        ft_assign(right_state(newp),right_state(oldp),front->sizest);
	STATE* newsl = (STATE*)left_state(newp);
	STATE* oldsl = (STATE*)left_state(oldp);
	for (int i = 0; i < dim; ++i)
        {
            newsl->vel[i] = oldsl->vel[i];
            Coords(newp)[i] = Coords(oldp)[i] + dt*oldsl->vel[i];
            newsl->collsnImpulse[i] = 0.0;
            newsl->x_old[i] = Coords(oldp)[i];
        }

	oldp = oldc->end->posn;
        newp = newc->end->posn;
	newsl = (STATE*)left_state(newp);
	oldsl = (STATE*)left_state(oldp);
        ft_assign(left_state(newp),left_state(oldp),front->sizest);
        ft_assign(right_state(newp),right_state(oldp),front->sizest);
	for (int i = 0; i < dim; ++i)
        {
            newsl->vel[i] = oldsl->vel[i];
            Coords(newp)[i] = Coords(oldp)[i] + dt*oldsl->vel[i];
            newsl->collsnImpulse[i] = 0.0;
            newsl->x_old[i] = Coords(oldp)[i];
        }


	for (oldb = oldc->first, newb = newc->first; oldb != oldc->last;
                oldb = oldb->next, newb = newb->next)
        {
            oldp = oldb->end;
            newp = newb->end;
	    newsl = (STATE*)left_state(newp);
	    oldsl = (STATE*)left_state(oldp);
            ft_assign(left_state(newp),left_state(oldp),front->sizest);
            ft_assign(right_state(newp),right_state(oldp),front->sizest);
	    for (int i = 0; i < dim; ++i)
            {
                newsl->vel[i] = oldsl->vel[i];
                Coords(newp)[i] = Coords(oldp)[i] + dt*oldsl->vel[i];
                newsl->collsnImpulse[i] = 0.0;
                newsl->x_old[i] = Coords(oldp)[i];
            }
        }
}

static void collision_point_propagate(
        Front *front,
        POINTER wave,
        POINT *oldp,
        POINT *newp,
        HYPER_SURF_ELEMENT *oldhse,
        HYPER_SURF         *oldhs,
        double              dt,
        double              *V)
{
        double vel[MAXD],s;
        int i, dim = front->rect_grid->dim;

	STATE *newsl;
        STATE *sl;
	sl = (STATE*)left_state(oldp);
        newsl = (STATE*)left_state(newp);
        ft_assign(left_state(newp),left_state(oldp),front->sizest);

	for (i = 0; i < dim; ++i)
	    vel[i] = sl->vel[i];

	for (i = 0; i < dim; ++i)
	{
	    newsl->vel[i] = vel[i];
            Coords(newp)[i] = Coords(oldp)[i];
	    newsl->collsnImpulse[i] = 0.0;
	    newsl->x_old[i] = Coords(oldp)[i];
	}
	newsl->collsn_num = 0;
        s = mag_vector(V,dim);
        set_max_front_speed(dim,s,NULL,Coords(newp),front);
}       /* fourth_order_point_propagate */

static void dummySpringSolver(Front* front) {
	INTERFACE *intfc = front->interf;
	POINT *p;
	HYPER_SURF              *hs;
        HYPER_SURF_ELEMENT      *hse;
	(void) next_point(intfc,NULL,NULL,NULL);
	while (next_point(intfc,&p,&hse,&hs)) 
        {
            if(Boundary_point(p))
                continue;
	    STATE* sl = (STATE*)left_state(p);
	    for (int i = 0; i < FT_Dimension(); ++i)
	    	Coords(p)[i] = sl->x_old[i] + front->dt*(sl->vel[i]);
        }
	CURVE**c;
	BOND* bond;
	intfc_curve_loop(intfc,c) {
	    for ((bond) = (*c)->first; (bond) != (*c)->last; 
		 (bond) = (bond)->next){
		p = bond->end;
		STATE* sl = (STATE*)left_state(p);
		for (int i = 0; i < FT_Dimension(); ++i)	
		    Coords(p)[i] = sl->x_old[i] + front->dt*(sl->vel[i]);
	    }
	}
	NODE** n;
	intfc_node_loop(intfc,n) {
	     p = (*n)->posn;
	     STATE* sl = (STATE*)left_state(p);
                for (int i = 0; i < FT_Dimension(); ++i)
                    Coords(p)[i] = sl->x_old[i] + front->dt*(sl->vel[i]);
	}
}
