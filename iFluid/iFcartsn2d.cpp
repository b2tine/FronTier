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

/*******************************************************************
 * 	               iFcartsn2d.cpp	
 *******************************************************************/
#include "iFluid.h"
#include "solver.h"

static double (*getStateVel[3])(POINTER) = {getStateXvel,getStateYvel,
                                        getStateZvel};


void Incompress_Solver_Smooth_2D_Cartesian::computeAdvection(void)
{
	int i,j,index;
	static HYPERB_SOLVER hyperb_solver(*front);
	static double *rho;
	double **vel = field->vel;

	if (rho == NULL)
	{
	    int size = (top_gmax[0]+1)*(top_gmax[1]+1);
	    FT_VectorMemoryAlloc((POINTER*)&rho,size,sizeof(double));
	}
	for (j = 0; j <= top_gmax[1]; j++)
	for (i = 0; i <= top_gmax[0]; i++)
	{
	    index = d_index2d(i,j,top_gmax);
	    rho[index] = field->rho[index];
	}
	if (debugging("field_var"))
	{
	    for (j = jmin; j <= jmax; j++)
	    for (i = imin; i <= imax; i++)
	    {
	    	index = d_index2d(i,j,top_gmax);
		field->old_var[0][index] = vel[0][index];
		field->old_var[1][index] = vel[1][index];
	    }
	}
	hyperb_solver.rho = rho;

	hyperb_solver.obst_comp = SOLID_COMP;
	switch (iFparams->adv_order)
	{
	case 1:
	    hyperb_solver.order = 1;
	    hyperb_solver.numericalFlux = upwind_flux;
	    break;
	case 4:
	    hyperb_solver.order = 4;
	    hyperb_solver.numericalFlux = weno5_flux;
	    break;
	default:
	    (void) printf("Advection order %d not implemented!\n",
					iFparams->adv_order);
	    clean_up(ERROR);
	}
	hyperb_solver.dt = m_dt;
	hyperb_solver.var = field->vel;
	hyperb_solver.soln = field->vel;
	hyperb_solver.soln_comp1 = LIQUID_COMP1;
	hyperb_solver.soln_comp2 = LIQUID_COMP2;
	hyperb_solver.getStateVel[0] = getStateXvel;
	hyperb_solver.getStateVel[1] = getStateYvel;
	hyperb_solver.rho1 = iFparams->rho1;
	hyperb_solver.rho2 = iFparams->rho2;
	hyperb_solver.findStateAtCrossing = findStateAtCrossing;
	hyperb_solver.solveRungeKutta();
	if (debugging("field_var"))
	{
	    (void) printf("\nIn computeAdvection(), \n");
	    (void) printf("one step increment for v[0]:\n");
	    computeVarIncrement(field->old_var[0],vel[0],NO);
	    (void) printf("one step increment for v[1]:\n");
	    computeVarIncrement(field->old_var[1],vel[1],NO);
	    (void) printf("\n");
	}
}

void Incompress_Solver_Smooth_2D_Cartesian::computeProjection(void)
{
	setIndexMap();
	switch (iFparams->num_scheme.ellip_method)
	{
	case SIMPLE_ELLIP:
	case DUAL_ELLIP:
	    computeProjectionSimple();
	    return;
	case DOUBLE_ELLIP:
	    computeProjectionDouble();
	    return;
	case CIM_ELLIP:
	    computeProjectionCim();
	    return;
	}
}	/* end computeProjection */

void Incompress_Solver_Smooth_2D_Cartesian::computeProjectionCim(void)
{
	static CIM_ELLIPTIC_SOLVER elliptic_solver(*front);
	int index;
	int i,j,l,icoords[MAXD];
	double **vel = field->vel;
	double *phi = field->phi;
	double *div_U = field->div_U;
	double sum_div;
	double value;

	sum_div = 0.0;
	max_value = 0.0;

	/* Compute velocity divergence */
	for (j = jmin; j <= jmax; j++)
        for (i = imin; i <= imax; i++)
	{
	    icoords[0] = i;
	    icoords[1] = j;
	    index  = d_index(icoords,top_gmax,dim);
	    source[index] = computeFieldPointDiv(icoords,vel);
	    diff_coeff[index] = 1.0/field->rho[index];
	}
	FT_ParallelExchGridArrayBuffer(source,front,NULL);
	FT_ParallelExchGridArrayBuffer(diff_coeff,front,NULL);
	for (j = 0; j <= top_gmax[1]; j++)
        for (i = 0; i <= top_gmax[0]; i++)
	{
	    index  = d_index2d(i,j,top_gmax);
	    div_U[index] = source[index];
	    source[index] /= -accum_dt;
	    array[index] = phi[index];
	}

	if(debugging("step_size"))
	{
	    for (j = jmin; j <= jmax; j++)
	    for (i = imin; i <= imax; i++)
	    {
		index = d_index2d(i,j,top_gmax);
	        value = fabs(div_U[index]);
		sum_div = sum_div + div_U[index];
	    }
	    pp_global_sum(&sum_div,1);
	    (void) printf("\nThe summation of divergence of U is %.16g\n",
					sum_div);
	    pp_global_max(&max_value,1);
	    (void) printf("\nThe max value of divergence of U is %.16g\n",
					max_value);
	    max_value = 0.0;
	}
	elliptic_solver.w_type = FIRST_PHYSICS_WAVE_TYPE;
	elliptic_solver.neg_comp = LIQUID_COMP2;
	elliptic_solver.pos_comp = LIQUID_COMP1;
	elliptic_solver.solutionJump = p_jump;
        elliptic_solver.gradJumpDotN = grad_p_jump_n;
        elliptic_solver.gradJumpDotT = grad_p_jump_t;

        elliptic_solver.diff_coeff[0] = 1.0/m_rho[0];
        elliptic_solver.diff_coeff[1] = 1.0/m_rho[1];
        elliptic_solver.kappa = diff_coeff;
        elliptic_solver.source = source;
        elliptic_solver.ilower = ilower;
        elliptic_solver.iupper = iupper;
        elliptic_solver.soln = array;
        elliptic_solver.ij_to_I = ij_to_I;
        elliptic_solver.I_to_ij = I_to_ij;
	elliptic_solver.getStateVar = getStatePhi;
	elliptic_solver.findStateAtCrossing = ifluid_find_state_at_crossing;
	elliptic_solver.size = iupper - ilower;
	elliptic_solver.set_solver_domain();
	elliptic_solver.solve(array);

	for (j = 0; j <= top_gmax[1]; j++)
	for (i = 0; i <= top_gmax[0]; i++)
	{
	    index  = d_index2d(i,j,top_gmax);
	    phi[index] = array[index];
	}
}	/* end computeProjectionCim */

void Incompress_Solver_Smooth_2D_Cartesian::computeProjectionDouble(void)
{
	iFparams->total_div_cancellation = YES;
	computeProjectionSimple(); 
	return;
}	/* end computeProjectionDouble */

void Incompress_Solver_Smooth_2D_Cartesian::computeProjectionSimple(void)
{
	static ELLIPTIC_SOLVER elliptic_solver(*front);
	int index;
	int i,j,l,icoords[MAXD];
	double **vel = field->vel;
	double *phi = field->phi;
	double *div_U = field->div_U;
	double sum_div;
	double value;
	int num_colors;

	sum_div = 0.0;
	max_value = 0.0;
	for (l = 0; l < dim; ++l)
	{
	    vmin[l] = HUGE;
	    vmax[l] = -HUGE;
	}

	if (debugging("field_var"))
	{
	    for (j = jmin; j <= jmax; j++)
	    for (i = imin; i <= imax; i++)
	    {
	    	index  = d_index2d(i,j,top_gmax);
		field->old_var[0][index] = field->phi[index];
	    }
	}
	/* Compute velocity divergence */
	for (j = jmin; j <= jmax; j++)
        for (i = imin; i <= imax; i++)
	{
	    icoords[0] = i;
	    icoords[1] = j;
	    index  = d_index(icoords,top_gmax,dim);
	    if (!ifluid_comp(top_comp[index]))
		continue;
	    source[index] = computeFieldPointDiv(icoords,vel);
	    diff_coeff[index] = 1.0/field->rho[index];
	    if (debugging("check_div"))
	    {
		for (l = 0; l < dim; ++l)
		{
		    if (vmin[l] > field->vel[l][index])
			vmin[l] = field->vel[l][index];
		    if (vmax[l] < field->vel[l][index])
			vmax[l] = field->vel[l][index];
		}
	    }
	}
	if (debugging("field_var"))
	{
	    (void) printf("\nCheck one step increment of div_U:\n");
	    computeVarIncrement(field->div_U,source,NO);
	    (void) printf("\n");
	}
	FT_ParallelExchGridArrayBuffer(source,front,NULL);
	FT_ParallelExchGridArrayBuffer(diff_coeff,front,NULL);
	for (j = 0; j <= top_gmax[1]; j++)
        for (i = 0; i <= top_gmax[0]; i++)
	{
	    index  = d_index2d(i,j,top_gmax);
	    if (!ifluid_comp(top_comp[index]))
		continue;
	    source[index] = (source[index])/accum_dt;
            /*Compute pressure jump due to porosity*/
            icoords[0] = i; icoords[1] = j;
            source[index] += computeFieldPointPressureJump(icoords,
                             iFparams->porous_coeff[0],
                             iFparams->porous_coeff[1]);
            /*end of computing pressure jump*/
	    array[index] = phi[index];
	}

	if(debugging("step_size"))
	{
	    for (j = jmin; j <= jmax; j++)
	    for (i = imin; i <= imax; i++)
	    {
		index = d_index2d(i,j,top_gmax);
	        if (!ifluid_comp(top_comp[index]))
		    continue;
	        value = fabs(div_U[index]);
		sum_div = sum_div + div_U[index];
		if (max_value < value)
		    max_value = value;
	    }
	    pp_global_sum(&sum_div,1);
	    (void) printf("\nThe summation of divergence of U is %.16g\n",
					sum_div);
	    pp_global_max(&max_value,1);
	    (void) printf("\nThe max value of divergence of U is %.16g\n",
					max_value);
	    max_value = 0.0;
	}
	if (debugging("check_div"))
        {
	    checkVelocityDiv("Before computeProjection()");
        }
        elliptic_solver.D = diff_coeff;
        elliptic_solver.source = source;
        elliptic_solver.soln = array;
	elliptic_solver.set_solver_domain();
	elliptic_solver.getStateVar = getStatePhi;
	elliptic_solver.findStateAtCrossing = findStateAtCrossing;
	elliptic_solver.skip_neumann_solver = skip_neumann_solver;
	num_colors = drawColorMap();
	paintAllGridPoint(NOT_SOLVED);
	for (i = 1; i < num_colors; ++i)
	{
	    paintToSolveGridPoint2(i);
	    setGlobalIndex();
            setIndexMap();
            elliptic_solver.ij_to_I = ij_to_I;
            elliptic_solver.ilower = ilower;
            elliptic_solver.iupper = iupper;
	    if (iFparams->total_div_cancellation)
	    	elliptic_solver.dsolve(array);
	    else
	    	elliptic_solver.solve(array);
	    paintSolvedGridPoint();
	}
	FT_ParallelExchGridArrayBuffer(array,front,NULL);

	for (j = 0; j <= top_gmax[1]; j++)
	for (i = 0; i <= top_gmax[0]; i++)
	{
	    index  = d_index2d(i,j,top_gmax);
	    phi[index] = array[index];
	}
	if (debugging("field_var"))
	{
	    printf("\nIn computeProjectionSimple()\n");
	    printf("Check one step increment of phi:\n");
	    computeVarIncrement(field->old_var[0],phi,NO);
	    printf("\n");
	    
	}
}	/* end computeProjectionSimple */

void Incompress_Solver_Smooth_2D_Cartesian::computeNewVelocity(void)
{
	int i,j,k,l,index;
	double grad_phi[2],rho;
	COMPONENT comp;
	int icoords[MAXD];
	double **vel = field->vel;
	double *phi = field->phi;

	if (iFparams->num_scheme.ellip_method == DUAL_ELLIP)
	{
	    computeNewVelocityDual();
	    return;
	}
	for (j = 0; j <= top_gmax[1]; j++)
	for (i = 0; i <= top_gmax[0]; i++)
	{
	    index  = d_index2d(i,j,top_gmax);
	    array[index] = phi[index];
	}

	for (j = jmin; j <= jmax; j++)
        for (i = imin; i <= imax; i++)
	{
	    index = d_index2d(i,j,top_gmax);
	    comp = top_comp[index];
	    if (!ifluid_comp(comp))
	    {
		vel[0][index] = 0.0;
		vel[1][index] = 0.0;
		continue;
	    }
	    rho = field->rho[index];
	    icoords[0] = i;
	    icoords[1] = j;
	    if (iFparams->with_porosity)
                computeFieldPointGradJump(icoords,phi,grad_phi);
            else
                computeFieldPointGrad(icoords,phi,grad_phi);
	    vel[0][index] -= accum_dt/rho*grad_phi[0];
	    vel[1][index] -= accum_dt/rho*grad_phi[1];
	}
	FT_ParallelExchGridVectorArrayBuffer(vel,front);
	extractFlowThroughVelocity();
	computeVelDivergence();

	if (debugging("check_div"))
        {
	    checkVelocityDiv("After computeNewVelocity()");
        }
}	/* end computeNewVelocity */

void Incompress_Solver_Smooth_2D_Cartesian::computeNewVelocityDual(void)
{
	int i,j,k,l,m,n,index,dual_index;
	double grad_phi[2],rho;
        COMPONENT comp;
        double v_ave;
        int icoords[MAXD];
	int icu[MAXD],icl[MAXD];
	int index_l,index_u;
        double **vel = field->vel;
        double *d_phi = field->d_phi;
	static double **V;
	int symmetry[MAXD];

	computeProjectionDual();

	if (V == NULL)
            FT_MatrixMemoryAlloc((POINTER*)&V,2,(top_gmax[0]+1)*(top_gmax[1]+1),
                			sizeof(double));
        for (j = 0; j <= top_gmax[1]; ++j)
        for (i = 0; i <= top_gmax[0]; ++i)
	{
            index = d_index2d(i,j,top_gmax);
	    for (l = 0; l < dim; ++l)
		V[l][index] = vel[l][index];
	}

        for (j = jmin; j <= jmax; ++j)
        for (i = imin; i <= imax; ++i)
        {
	    icoords[0] = i;
	    icoords[1] = j;
            index = d_index(icoords,top_gmax,dim);
            comp = top_comp[index];

	    rho = field->rho[index];
	    computeDualFieldPointGrad(icoords,d_phi,grad_phi);
	    if (!ifluid_comp(comp))
	    {
	    	int index1,index2,index3,index4;
	    	index1 = d_index2d(i-1,j,top_gmax);
	    	index2 = d_index2d(i+1,j,top_gmax);
	    	index3 = d_index2d(i,j-1,top_gmax);
	    	index4 = d_index2d(i,j+1,top_gmax);
		if (!ifluid_comp(top_comp[index1])&&
		    !ifluid_comp(top_comp[index2])&&
		    !ifluid_comp(top_comp[index3])&&
		    !ifluid_comp(top_comp[index4]))
		{
		    vel[0][index] = 0.0;
                    vel[1][index] = 0.0;
		}
		else if (!(ifluid_comp(top_comp[index1]))&&
                        !(ifluid_comp(top_comp[index2])))
		{
		    if (ifluid_comp(top_comp[index3]))
		    	vel[0][index] = (V[0][index]+V[0][index3])/2.0 -
				accum_dt/rho*grad_phi[0];
		    else
			vel[0][index] = (V[0][index]+V[0][index4])/2.0 -
                                accum_dt/rho*grad_phi[0];

		    vel[1][index] = (V[1][index1]+V[1][index2] +
                        2*V[1][index])/4.0 - accum_dt/rho*grad_phi[1]; 
		}
		else if (!(ifluid_comp(top_comp[index3]))&&
                        !(ifluid_comp(top_comp[index4])))
		{
		    vel[0][index] = 0.0;
		    vel[0][index] = (V[0][index3]+V[0][index4]+
                        2*V[0][index])/4.0 - accum_dt/rho*grad_phi[0];
		    if (ifluid_comp(top_comp[index2]))
			vel[1][index] = (V[1][index]+V[1][index2])/2.0 -
				accum_dt/rho*grad_phi[1];
		    else
			vel[1][index] = (V[1][index]+V[1][index1])/2.0 -
                                accum_dt/rho*grad_phi[1];
		} else
		{
		    vel[0][index] = (V[0][index3]+V[0][index4]+
                        2*V[0][index])/4.0 - accum_dt/rho*grad_phi[0];
		    vel[1][index] = (V[1][index1]+V[1][index2] +
                        2*V[1][index])/4.0 - accum_dt/rho*grad_phi[1];
		}
		continue;
	    }

	    double denom;
	    for (l = 0; l < dim; ++l)
	    {
		v_ave = 0.0;
                for (n = 1; n < dim; ++n)
                {
                    denom = 0.0;
                    for (m = 0; m < dim; ++m)
                        icl[m] = icu[m] = icoords[m];
                    icl[(l+n)%dim] -= 1;        icu[(l+n)%dim] += 1;
                    index_l = d_index(icl,top_gmax,dim);
                    index_u = d_index(icu,top_gmax,dim);
                    v_ave += 0.5*(V[l][index] + V[l][index_l]);
                    v_ave += 0.5*(V[l][index] + V[l][index_u]);
                    denom += 2.0;
                }
		v_ave /= denom;
		vel[l][index] = v_ave - accum_dt/rho*grad_phi[l];
	    }
	}

	FT_ParallelExchGridVectorArrayBuffer(vel,front);

        if (debugging("check_div"))
	    checkVelocityDiv("Before extractFlowThroughVelocity()");
	extractFlowThroughVelocity();
        if (debugging("check_div"))
            checkVelocityDiv("After extractFlowThroughVelocity()");
}	/* end computeNewVelocityDual */

void Incompress_Solver_Smooth_2D_Cartesian::computeSourceTerm(
	double *coords, 
	double *source) 
{
        int i,j,k,l;
	double x,y,a,f,Force;
	double t = front->time;
	double phi[MAXD+1];
	x = coords[0];
	y = coords[1];
	UNIFORM_PARAMS uniform_params;
	short unsigned int xsubi[3]; 

	if(iFparams->if_buoyancy)
	{
	    int ic[MAXD],index;
            rect_in_which(coords,ic,top_grid);
            index = d_index(ic,top_gmax,dim);
            for (i = 0; i < dim; ++i)
                source[i] = field->ext_accel[i][index];
	}
	else
	{
            for (i = 0; i < dim; ++i)
                source[i] = iFparams->gravity[i];
	}
}	/* end computeSourceTerm */

void Incompress_Solver_Smooth_2D_Cartesian::solve(double dt)
{
	static boolean first = YES;
	double **vel = field->vel;
	if (first)
	{
	    accum_dt = 0.0;
	    first = NO;
	}
	m_dt = dt;
	int j,icoords[MAXD];
	double div;

	start_clock("solve");
	setDomain();

	setComponent();
	if (iFparams->num_scheme.ellip_method == DUAL_ELLIP)
	    updateComponent();
	if (debugging("trace"))
	    printf("Passed setComponent()\n");

	paintAllGridPoint(TO_SOLVE);
	setGlobalIndex();
	if (debugging("trace"))
	    printf("Passed setGlobalIndex()\n");
	start_clock("setSmoothedProperties");
	setSmoothedProperties();
	stop_clock("setSmoothedProperties");
	if (debugging("trace"))
	    printf("Passed setSmoothedProperties()\n");
	
	// 1) solve for intermediate velocity
	start_clock("computeAdvection");
	computeAdvection();
	if (debugging("check_div"))
            checkVelocityDiv("After computeAdvection()");
	if (debugging("check_div") || debugging("step_size"))
	{
	    computeMaxSpeed();
	    (void) printf("max_speed after   computeAdvection(): %20.14f ",
				max_speed);
	    (void) printf("occured at (%d, %d)\n",icrds_max[0],icrds_max[1]);
	}
	stop_clock("computeAdvection");

	if (debugging("sample_velocity"))
	    sampleVelocity();
	start_clock("computeDiffusion");
	computeDiffusion();
	if (debugging("check_div"))
            checkVelocityDiv("After computeDiffusion()");
	if (debugging("check_div") || debugging("step_size"))
	{
	    computeMaxSpeed();
	    (void) printf("max_speed after   computeDiffusion(): %20.14f ",
				max_speed);
	    (void) printf("occured at (%d, %d)\n",icrds_max[0],icrds_max[1]);
	}
	stop_clock("computeDiffusion");

	if (debugging("sample_velocity"))
	    sampleVelocity();

	// 2) projection step
	accum_dt += m_dt;
	if (debugging("step_size"))
	    (void) printf("min_dt = %f  accum_dt = %f\n",min_dt,accum_dt);
	if (accum_dt >= min_dt)
	{
	    start_clock("computeProjection");
	    computeProjection();
	    stop_clock("computeProjection");

	    start_clock("computePressure");
	    computePressure();
	    stop_clock("computePressure");

	    start_clock("computeNewVelocity");
	    computeNewVelocity();
	    stop_clock("computeNewVelocity");
	    accum_dt = 0.0;
	}
	computeMaxSpeed();

	if (debugging("sample_velocity"))
	    sampleVelocity();
	if (debugging("step_size"))
	{
	    (void) printf("max_speed after computeNewVelocity(): %20.14f ",
				max_speed);
	    (void) printf("occured at (%d, %d)\n",icrds_max[0],icrds_max[1]);
	}
	start_clock("copyMeshStates");
	copyMeshStates();
	stop_clock("copyMeshStates");

	/* The following is for climate modelling */
        if(iFparams->if_ref_pres == YES)
            setReferencePressure();

	setAdvectionDt();
	stop_clock("solve");
}	/* end solve */


double Incompress_Solver_Smooth_2D_Cartesian::getVorticity(int i, int j)
{
        int icoords[MAXD],icnb[MAXD];
        int index,index_nb;
        COMPONENT comp;
        double div,u_edge[2][2];
        double crx_coords[MAXD];
        POINTER intfc_state;
        HYPER_SURF *hs;
        int status;
        GRID_DIRECTION dir[2][2] = {{WEST,EAST},{SOUTH,NORTH}};
        int k,idir,nb;
        double u0;
        double dx,dy;
        double vorticity;
        double **vel = field->vel;
        int dim = 2;

        icoords[0] = i;
        icoords[1] = j;

        index = d_index(icoords,top_gmax,dim);
        comp = top_comp[index];

        for (idir = 0; idir < dim; idir++)
        {
            u0 = vel[(idir+1)%dim][index];
            for (k = 0; k < dim; ++k)
                icnb[k] = icoords[k];
            for (nb = 0; nb < 2; nb++)
            {
                dx = top_h[0];
                dy = top_h[1];
                icnb[idir] = (nb == 0) ? icoords[idir] - 1 : icoords[idir] + 1;
                index_nb = d_index(icnb,top_gmax,dim);
                status = (*findStateAtCrossing)(front,icoords,dir[idir][nb],
                                comp,&intfc_state,&hs,crx_coords);
                if (status == NO_PDE_BOUNDARY)
                    u_edge[idir][nb] = vel[(idir+1)%dim][index_nb];
                else if (status ==CONST_P_PDE_BOUNDARY)
                    u_edge[idir][nb] = u0;
                else if (status ==CONST_V_PDE_BOUNDARY &&
                        wave_type(hs) == DIRICHLET_BOUNDARY)
                    u_edge[idir][nb] = getStateVel[(idir+1)%dim](intfc_state);
                else
                    u_edge[idir][nb] = u0;
            }
        }

        vorticity = (u_edge[0][1]-u_edge[0][0])/2/dx -
                        (u_edge[1][1]-u_edge[1][0])/2/dy;

	return vorticity;
}	/* end getVorticity */

void Incompress_Solver_Smooth_2D_Cartesian::copyMeshStates(void)
{
	int i,j,index;
	double **vel = field->vel;
	double *pres = field->pres;
	double *vort = field->vort;
	int symmetry[MAXD];
	
	for (i = imin; i <= imax; ++i)
	for (j = jmin; j <= jmax; ++j)
	{
	    index  = d_index2d(i,j,top_gmax);
	    if (ifluid_comp(top_comp[index]))
		vort[index] = getVorticity(i,j);
	    else
		vort[index] = 0.0;
	}
	symmetry[0] = symmetry[1] = ODD;
	FT_ParallelExchGridArrayBuffer(vort,front,symmetry);
}	/* end copyMeshStates */

void Incompress_Solver_Smooth_2D_Cartesian::
	computeDiffusion(void)
{
	return computeDiffusionCN();
}

void Incompress_Solver_Smooth_2D_Cartesian::
	computeDiffusionCN(void)
{
        COMPONENT comp;
        int index,index_nb[4],size;
        int I,I_nb[4];
        int i,j,k,l,nb,icoords[MAXD];
        double coords[MAXD], crx_coords[MAXD];
        double coeff[4],mu[4],mu0,rho,rhs,U_nb[4];
        double *x;
        GRID_DIRECTION dir[4] = {WEST,EAST,SOUTH,NORTH};
        POINTER intfc_state;
        HYPER_SURF *hs;
        PetscInt num_iter;
        double rel_residual;
        double aII;
        double source[MAXD];
        double **vel = field->vel;
        double **f_surf = field->f_surf;
        INTERFACE *grid_intfc = front->grid_intfc;
	int status;

        if (debugging("trace"))
            (void) printf("Entering Incompress_Solver_Smooth_2D_Cartesian::"
                        "computeDiffusionCN()\n");
	start_clock("computeDiffusionCN");

        setIndexMap();

        size = iupper - ilower;
        FT_VectorMemoryAlloc((POINTER*)&x,size,sizeof(double));

	if (debugging("field_var"))
	{
	    for (j = jmin; j <= jmax; j++)
	    for (i = imin; i <= imax; i++)
	    {
	    	index = d_index2d(i,j,top_gmax);
		field->old_var[0][index] = vel[0][index];
		field->old_var[1][index] = vel[1][index];
	    }
	}
	for (l = 0; l < dim; ++l)
	{
            PETSc solver;
            solver.Create(ilower, iupper-1, 5, 5);
	    solver.Reset_A();
	    solver.Reset_b();
	    solver.Reset_x();

            for (j = jmin; j <= jmax; j++)
            for (i = imin; i <= imax; i++)
            {
            	I  = ij_to_I[i][j];
            	if (I == -1) continue;

            	index  = d_index2d(i,j,top_gmax);
            	index_nb[0] = d_index2d(i-1,j,top_gmax);
            	index_nb[1] = d_index2d(i+1,j,top_gmax);
            	index_nb[2] = d_index2d(i,j-1,top_gmax);
            	index_nb[3] = d_index2d(i,j+1,top_gmax);
		icoords[0] = i;
		icoords[1] = j;
		comp = top_comp[index];

            	I_nb[0] = ij_to_I[i-1][j]; // left or west
            	I_nb[1] = ij_to_I[i+1][j]; // right or east
            	I_nb[2] = ij_to_I[i][j-1]; // down or south
            	I_nb[3] = ij_to_I[i][j+1]; // up or north


		mu0 = field->mu[index];
		rho = field->rho[index];

            	for (nb = 0; nb < 4; nb++)
            	{
                    if ((*findStateAtCrossing)(front,icoords,dir[nb],comp,
					&intfc_state,&hs,crx_coords) &&
                                wave_type(hs) != FIRST_PHYSICS_WAVE_TYPE)
                    {
			if (wave_type(hs) == DIRICHLET_BOUNDARY &&
                            boundary_state_function(hs) &&
                            strcmp(boundary_state_function_name(hs),
                            "flowThroughBoundaryState") == 0)
                        {
                            U_nb[nb] = vel[l][index];
                        }
                        else
			{
			    U_nb[nb] = getStateVel[l](intfc_state);
			}
			
			if (wave_type(hs) == DIRICHLET_BOUNDARY ||
			    neumann_type_bdry(wave_type(hs)))
                            mu[nb] = mu0;
			else
			    mu[nb] = 1.0/2*(mu0 + field->mu[index_nb[nb]]);
                    }
                    else
		    {
                        U_nb[nb] = vel[l][index_nb[nb]];
			mu[nb] = 1.0/2*(mu0 + field->mu[index_nb[nb]]);
		    }
            	}

            	coeff[0] = 0.5*m_dt/rho*mu[0]/(top_h[0]*top_h[0]);
            	coeff[1] = 0.5*m_dt/rho*mu[1]/(top_h[0]*top_h[0]);
            	coeff[2] = 0.5*m_dt/rho*mu[2]/(top_h[1]*top_h[1]);
            	coeff[3] = 0.5*m_dt/rho*mu[3]/(top_h[1]*top_h[1]);

            	getRectangleCenter(index, coords);
            	computeSourceTerm(coords, source);

        	//first equation  decoupled, some terms may be lost
		aII = 1+coeff[0]+coeff[1]+coeff[2]+coeff[3];
            	rhs = (1-coeff[0]-coeff[1]-coeff[2]-coeff[3])
				*vel[l][index];

            	for (nb = 0; nb < 4; nb++)
            	{
		    status = (*findStateAtCrossing)(front,icoords,dir[nb],comp,
                                &intfc_state,&hs,crx_coords);
	    	    if (status == NO_PDE_BOUNDARY)
		    {
            	    	solver.Set_A(I,I_nb[nb],-coeff[nb]);
            	    	rhs += coeff[nb]*U_nb[nb];
		    }
		    else
		    {
                        if (status == CONST_P_PDE_BOUNDARY)
                        {
                            aII -= coeff[nb];
                            rhs += coeff[nb]*U_nb[nb];
                        }
                        else
                            rhs += 2.0*coeff[nb]*U_nb[nb];
                    }
		}
            	rhs += m_dt*source[l];
		rhs += m_dt*f_surf[l][index];
	    	//rhs -= m_dt*grad_q[l][index]/rho;

		solver.Set_A(I,I,aII);
            	solver.Set_b(I, rhs);
            }

            solver.SetMaxIter(40000);
            solver.SetTol(1e-14);

	    start_clock("Befor Petsc solve");
            solver.Solve();
            solver.GetNumIterations(&num_iter);
            solver.GetFinalRelativeResidualNorm(&rel_residual);

	    stop_clock("After Petsc solve");

            // get back the solution
            solver.Get_x(x);

            if (debugging("PETSc"))
                (void) printf("Incompress_Solver_Smooth_2D_Cartesian::"
			"computeDiffusion: "
       			"num_iter = %d, rel_residual = %g. \n", 
                        num_iter,rel_residual);

            for (j = jmin; j <= jmax; j++)
            for (i = imin; i <= imax; i++)
            {
                I = ij_to_I[i][j];
                index = d_index2d(i,j,top_gmax);
                if (I >= 0)
                    vel[l][index] = x[I-ilower];
                else
                    vel[l][index] = 0.0;
            }
        }
	FT_ParallelExchGridVectorArrayBuffer(vel,front);
	stop_clock("computeDiffusionCN");
	if (debugging("field_var"))
	{
	    (void) printf("\nIn computeDiffusionCN(), \n");
	    (void) printf("one step increment for v[0]:\n");
	    computeVarIncrement(field->old_var[0],vel[0],NO);
	    (void) printf("one step increment for v[1]:\n");
	    computeVarIncrement(field->old_var[1],vel[1],NO);
	    (void) printf("\n");
	}

        FT_FreeThese(1,x);
        if (debugging("trace"))
            (void) printf("Leaving Incompress_Solver_Smooth_3D_Cartesian::"
                        "computeDiffusionCN()\n");
}	/* end computeDiffusionCN */

void Incompress_Solver_Smooth_2D_Cartesian::computePressurePmI(void)
{
        int i,j,index;
	double *pres = field->pres;
	double *phi = field->phi;
	double *q = field->q;

	for (j = 0; j <= top_gmax[1]; j++)
        for (i = 0; i <= top_gmax[0]; i++)
	{
            index = d_index2d(i,j,top_gmax);
            pres[index] += phi[index];
	    q[index] = pres[index];
	}
}        /* end computePressurePmI2d */


void Incompress_Solver_Smooth_2D_Cartesian::computePressurePmII(void)
{
        int i,j,index;
        double mu0;
	double *pres = field->pres;
	double *phi = field->phi;
	double *q = field->q;
	double *div_U = field->div_U;

	for (j = 0; j <= top_gmax[1]; j++)
        for (i = 0; i <= top_gmax[0]; i++)
	{
            index = d_index2d(i,j,top_gmax);
            mu0 = 0.5*field->mu[index];
	    pres[index] += phi[index] - accum_dt*mu0*div_U[index];
	    q[index] = pres[index];
	}
}        /* end computePressurePmII2d */

void Incompress_Solver_Smooth_2D_Cartesian::computePressurePmIII(void)
{
        int i,j,index;
        double mu0;
	double *pres = field->pres;
	double *phi = field->phi;
	double *q = field->q;
	double *div_U = field->div_U;

	if (debugging("field_var"))
	{
	    for (j = jmin; j <= jmax; j++)
            for (i = imin; i <= imax; i++)
	    {
            	index = d_index2d(i,j,top_gmax);
		field->old_var[0][index] = pres[index];
	    }
	}
	for (j = 0; j <= top_gmax[1]; j++)
        for (i = 0; i <= top_gmax[0]; i++)
	{
            index = d_index2d(i,j,top_gmax);
            mu0 = 0.5*field->mu[index];
            pres[index] = phi[index] -
                        	accum_dt*mu0*div_U[index];
	    q[index] = 0.0;
	}
	if (debugging("field_var"))
	{
	    (void) printf("\nCheck one step increment of Pressure:\n");
	    computeVarIncrement(field->old_var[0],pres,NO);
	    (void) printf("\n");
	}
}        /* end computePressurePmIII */

void Incompress_Solver_Smooth_2D_Cartesian::computePressure(void)
{
	switch (iFparams->num_scheme.projc_method)
	{
	case BELL_COLELLA:
	    computePressurePmI();
	    break;
	case KIM_MOIN:
	    computePressurePmII();
	    break;
	case SIMPLE:
	case PEROT_BOTELLA:
	    computePressurePmIII();
	    break;
	case ERROR_PROJC_SCHEME:
	default:
	    (void) printf("Unknown computePressure scheme!\n");
	    clean_up(ERROR);
	}
	computeGradientQ();
}	/* end computePressure */

void Incompress_Solver_Smooth_2D_Cartesian::computeGradientQ(void)
{
	int i,j,l,index;
	double **grad_q = field->grad_q;
	int icoords[MAXD];
	double *phi = field->phi;
	double point_grad_q[MAXD];

	if (debugging("field_var"))
	{
	    int comp;
	    for (j = jmin; j <= jmax; j++)
            for (i = imin; i <= imax; i++)
	    {
	    	index = d_index2d(i,j,top_gmax);
	    	comp = top_comp[index];
	    	if (!ifluid_comp(comp))
		    continue;
	    	icoords[0] = i;
	    	icoords[1] = j;
		field->old_var[0][index] = field->grad_q[0][index];
		field->old_var[1][index] = field->grad_q[1][index];
	    }
	}
	for (j = 0; j < top_gmax[1]; ++j)
	for (i = 0; i < top_gmax[0]; ++i)
	{
	    index = d_index2d(i,j,top_gmax);
	    array[index] = phi[index];
	}
	for (j = jmin; j <= jmax; j++)
	for (i = imin; i <= imax; i++)
	{
	    index = d_index2d(i,j,top_gmax);
	    icoords[0] = i;
	    icoords[1] = j;
	    computeFieldPointGrad(icoords,array,point_grad_q);
	    for (l = 0; l < dim; ++l)
	    	grad_q[l][index] = point_grad_q[l];
	}
	FT_ParallelExchGridVectorArrayBuffer(grad_q,front);
	if (debugging("field_var"))
	{
	    (void) printf("\nIn computeGradientQ(), \n");
	    (void) printf("one step increment for grad_phi[0]:\n");
	    computeVarIncrement(field->old_var[0],field->grad_q[0],NO);
	    (void) printf("one step increment for grad_phi[1]:\n");
	    computeVarIncrement(field->old_var[1],field->grad_q[1],NO);
	    printf("\n");
	}
}	/* end computeGradientQ2d */


void Incompress_Solver_Smooth_2D_Cartesian::surfaceTension(
	double *coords,
	HYPER_SURF_ELEMENT *hse,
	HYPER_SURF *hs,
	double *force,
	double sigma)
{
	int i,k,nb;
	BOND *pb,*bonds[100];
	double kappa,nor[MAXD];
	double kappa0,nor0[MAXD];
	double kappa1,nor1[MAXD];
	double len,delta;
	HYPER_SURF_ELEMENT *phse;
	double p[MAXD];
	

	BondAndNeighbors(hse,hs,&nb,bonds,3);

	for (i = 0; i < dim; ++i) force[i] = 0.0;
	for (k = 0; k < nb; ++k)
	{
	    pb = bonds[k];
	    for (i = 0; i < dim; ++i) 
		p[i] = 0.5*(Coords(pb->start)[i] + Coords(pb->end)[i]);
	    delta = smoothedDeltaFunction(coords,p);
	    if (delta == 0.0) continue;

	    len = bond_length(pb);
	    phse = Hyper_surf_element(pb);
	    GetFrontNormal(pb->start,phse,hs,nor0,front);
	    GetFrontNormal(pb->end,phse,hs,nor1,front);
	    for (i = 0; i < dim; ++i) nor[i] = 0.5*(nor0[i] + nor1[i]);
	    GetFrontCurvature(pb->start,phse,hs,&kappa0,front);
	    GetFrontCurvature(pb->end,phse,hs,&kappa1,front);
	    kappa = 0.5*(kappa0 + kappa1);
	    for (i = 0; i < dim; ++i) 
	    {
		force[i] += delta*sigma*len*kappa*nor[i];
	    }
	}
}	/* end surfaceTension2d */

void Incompress_Solver_Smooth_2D_Cartesian::setInitialCondition()
{
	int i;
	COMPONENT comp;
	double coords[MAXD];
	int size = (int)cell_center.size();

	FT_MakeGridIntfc(front);
	if (iFparams->num_scheme.ellip_method == DUAL_ELLIP)
	    FT_MakeCompGridIntfc(front);
	setDomain();

        m_rho[0] = iFparams->rho1;
        m_rho[1] = iFparams->rho2;
        m_mu[0] = iFparams->mu1;
        m_mu[1] = iFparams->mu2;
	m_comp[0] = iFparams->m_comp1;
	m_comp[1] = iFparams->m_comp2;
	m_smoothing_radius = iFparams->smoothing_radius;
	m_sigma = iFparams->surf_tension;
	mu_min = rho_min = HUGE;
	for (i = 0; i < 2; ++i)
	{
	    if (ifluid_comp(m_comp[i]))
	    {
        	mu_min = std::min(mu_min,m_mu[i]);
        	rho_min = std::min(rho_min,m_rho[i]);
	    }
	}

	// Initialize state at cell_center
        for (i = 0; i < size; i++)
        {
            getRectangleCenter(i, coords);
	    //cell_center[i].m_state.setZero();
	    comp = top_comp[i];
	    if (getInitialState != NULL)
	    	(*getInitialState)(comp,coords,field,i,dim,iFparams);
        }
	computeGradientQ();
        copyMeshStates();
	setAdvectionDt();
}       /* end setInitialCondition */

void Incompress_Solver_Smooth_2D_Cartesian::setParallelVelocity()
{
        FILE *infile;
        int i,j,id,k,l,index,G_index;
        char fname[100];
        COMPONENT comp;
        double coords[MAXD];
        int size = (int)cell_center.size();
	int myid = pp_mynode();
	int numprocs = pp_numnodes();

        int G_icoords[MAXD],pp_icoords[MAXD],icoords[MAXD];
        int local_gmax[MAXD], global_gmax[MAXD];
        int G_size, L_size;
	PP_GRID *pp_grid = front->pp_grid;
	double *local_L = pp_grid->Zoom_grid.L;
	double *local_U = pp_grid->Zoom_grid.U;
	double *GU_buff,*GV_buff, *U_buff, *V_buff;
        for (i = 0; i < dim; i++)
        {
            global_gmax[i] = pp_grid->Global_grid.gmax[i]-1;
            local_gmax[i] = pp_grid->Zoom_grid.gmax[i]-1;
        }
        FT_MakeGridIntfc(front);
        setDomain();
        G_size = 1;
        L_size = 1;
	for (i = 0; i < dim; i++)
	{
	    G_size = G_size * (global_gmax[i]+1);
	    L_size = L_size * (top_gmax[i]+1);
	}
	uni_array(&U_buff,L_size,sizeof(double));
	uni_array(&V_buff,L_size,sizeof(double));
	if (myid == 0)
	{
	    uni_array(&GU_buff,G_size,sizeof(double));
	    uni_array(&GV_buff,G_size,sizeof(double));
	 
	    if (setInitialVelocity != NULL)
                (*setInitialVelocity)(comp,pp_grid->Global_grid.gmax,
				   GU_buff,GV_buff,NULL,
				    &(pp_grid->Global_grid),iFparams);
	    for (id = 0; id < numprocs; id++)
	    {            
		find_Cartesian_coordinates(id,pp_grid,pp_icoords);
		for (j = jmin; j <= jmax; ++j)
            	for (i = imin; i <= imax; ++i)
            	{
	            icoords[0] = i;
	            icoords[1] = j;
	            G_icoords[0] = pp_icoords[0]*(local_gmax[0]+1)+icoords[0]-imin;
	            G_icoords[1] = pp_icoords[1]*(local_gmax[1]+1)+icoords[1]-jmin;
	            G_index = d_index(G_icoords,global_gmax,dim);
	            index = d_index(icoords,top_gmax,dim);
	            U_buff[index] = GU_buff[G_index];
	            V_buff[index] = GV_buff[G_index];
	        }
		if (id == 0)
		{	
		    for (i = 0; i < L_size; i++)
		    {
			field->vel[0][i] = U_buff[i];
			field->vel[1][i] = V_buff[i];
		    }
		}
		else
		{
	            pp_send(1,(POINTER)(U_buff),sizeof(double)*L_size,id);
	            pp_send(2,(POINTER)(V_buff),sizeof(double)*L_size,id);
		}
	    }
	    FT_FreeThese(2,GU_buff,GV_buff);
	}
	else
	{
	    pp_recv(1,0,(POINTER)(U_buff),sizeof(double)*L_size);
	    pp_recv(2,0,(POINTER)(V_buff),sizeof(double)*L_size);
	    for (i = 0; i < L_size; i++)
	    {
	        field->vel[0][i] = U_buff[i];
		field->vel[1][i] = V_buff[i];
	    }
	}


        m_rho[0] = iFparams->rho1;
        m_rho[1] = iFparams->rho2;
        m_mu[0] = iFparams->mu1;
        m_mu[1] = iFparams->mu2;
        m_comp[0] = iFparams->m_comp1;
        m_comp[1] = iFparams->m_comp2;
        m_smoothing_radius = iFparams->smoothing_radius;
        m_sigma = iFparams->surf_tension;
        mu_min = rho_min = HUGE;
        for (i = 0; i < 2; ++i)
        {
            if (ifluid_comp(m_comp[i]))
            {
                mu_min = std::min(mu_min,m_mu[i]);
                rho_min = std::min(rho_min,m_rho[i]);
            }
        }
	FT_FreeThese(2,U_buff,V_buff);

	computeGradientQ();
        copyMeshStates();
        setAdvectionDt();
}

void Incompress_Solver_Smooth_2D_Cartesian::computeProjectionDual(void)
{
	static DUAL_ELLIPTIC_SOLVER dual_elliptic_solver(*front);
	int index;
	int i,j,l,icoords[MAXD];
	double **vel = field->vel;
	double *phi = field->phi;
	double *div_U = field->div_U;
	double *d_phi = field->d_phi;
	double sum_div;
	double value;
	setDualGlobalIndex();
	setDualIndexMap();

	updateComponent();

	sum_div = 0.0;
	max_value = 0.0;

	/* Compute velocity divergence */
	for (j = jmin; j <= jmax; j++)
        for (i = imin; i <= imax; i++)
	{
	    index  = d_index2d(i,j,top_gmax);
	    diff_coeff[index] = 1.0/field->rho[index];
	}
	FT_ParallelExchGridArrayBuffer(diff_coeff,front,NULL);
	for (j = cjmin; j <= cjmax; j++)
        for (i = cimin; i <= cimax; i++)
	{
	    icoords[0] = i;
	    icoords[1] = j;
	    index  = d_index(icoords,ctop_gmax,dim);
	    if (!ifluid_comp(ctop_comp[index]))
		continue;
	    source[index] = computeDualFieldPointDiv(icoords,vel);
	}
	FT_ParallelExchCompGridArrayBuffer(source,front,NULL);

	for (j = 0; j <= ctop_gmax[1]; j++)
        for (i = 0; i <= ctop_gmax[0]; i++)
	{
	    index  = d_index2d(i,j,ctop_gmax);
	    if (!ifluid_comp(ctop_comp[index]))
		continue;
	    div_U[index] = source[index];
	    source[index] /= accum_dt;
	    array[index] = d_phi[index];
	}

	if(debugging("step_size"))
	{
	    for (j = cjmin; j <= cjmax; j++)
	    for (i = cimin; i <= cimax; i++)
	    {
		index = d_index2d(i,j,ctop_gmax);
	        if (!ifluid_comp(ctop_comp[index]))
		    continue;
	        value = fabs(div_U[index]);
		sum_div = sum_div + div_U[index];
		if (max_value < value)
		    max_value = value;
	    }
	    pp_global_sum(&sum_div,1);
	    (void) printf("\nThe summation of divergence of U is %.16g\n",
					sum_div);
	    pp_global_max(&max_value,1);
	    (void) printf("\nThe max value of divergence of U is %.16g\n",
					max_value);
	    max_value = 0.0;
	}
        dual_elliptic_solver.D = diff_coeff;
        dual_elliptic_solver.source = source;
        dual_elliptic_solver.ilower = cilower;
        dual_elliptic_solver.iupper = ciupper;
        dual_elliptic_solver.soln = array;
        dual_elliptic_solver.obst_comp = SOLID_COMP;
        dual_elliptic_solver.ij_to_I = cij_to_I;
	dual_elliptic_solver.set_solver_domain();
	dual_elliptic_solver.getStateVar = getStatePhi;
	dual_elliptic_solver.findStateAtCrossing = 
		ifluid_find_state_at_dual_crossing;
	dual_elliptic_solver.solve(array);

	for (j = 0; j <= ctop_gmax[1]; j++)
	for (i = 0; i <= ctop_gmax[0]; i++)
	{
	    index  = d_index2d(i,j,ctop_gmax);
	    d_phi[index] = array[index];
	}
}	/* end computeProjectionDual */

void Incompress_Solver_Smooth_2D_Cartesian::
        computeDiffusionImplicit(void)
{
        COMPONENT comp;
        int index,index_nb[4],size;
        int I,I_nb[4];
        int i,j,k,l,nb,icoords[MAXD];
        double coords[MAXD], crx_coords[MAXD];
        double coeff[4],mu[4],mu0,rho,rhs,U_nb[4];
        double *x;
        GRID_DIRECTION dir[4] = {WEST,EAST,SOUTH,NORTH};
        POINTER intfc_state;
        HYPER_SURF *hs;
        PetscInt num_iter;
        double rel_residual;
        double aII;
        double source[MAXD];
        double **vel = field->vel;
        double **f_surf = field->f_surf;
        INTERFACE *grid_intfc = front->grid_intfc;

        if (debugging("trace"))
            (void) printf("Entering Incompress_Solver_Smooth_3D_Cartesian::"
                        "computeDiffusionImplicit()\n");

        setIndexMap();

        size = iupper - ilower;
        FT_VectorMemoryAlloc((POINTER*)&x,size,sizeof(double));

        for (l = 0; l < dim; ++l)
        {
            PETSc solver;
            solver.Create(ilower, iupper-1, 5, 5);
            solver.Reset_A();
            solver.Reset_b();
            solver.Reset_x();

            for (j = jmin; j <= jmax; j++)
            for (i = imin; i <= imax; i++)
            {
                I  = ij_to_I[i][j];
                if (I == -1) continue;

                index  = d_index2d(i,j,top_gmax);
                index_nb[0] = d_index2d(i-1,j,top_gmax);
                index_nb[1] = d_index2d(i+1,j,top_gmax);
                index_nb[2] = d_index2d(i,j-1,top_gmax);
                index_nb[3] = d_index2d(i,j+1,top_gmax);

                icoords[0] = i;
                icoords[1] = j;
                comp = top_comp[index];

                I_nb[0] = ij_to_I[i-1][j]; //west
                I_nb[1] = ij_to_I[i+1][j]; //east
                I_nb[2] = ij_to_I[i][j-1]; //south
                I_nb[3] = ij_to_I[i][j+1]; //north


                mu0   = field->mu[index];
                rho   = field->rho[index];

                for (nb = 0; nb < 4; nb++)
                {
                    if (FT_StateStructAtGridCrossing(front,grid_intfc,icoords,
                                dir[nb],comp,&intfc_state,&hs,crx_coords) &&
                                wave_type(hs) != FIRST_PHYSICS_WAVE_TYPE)
                    {
                        if (wave_type(hs) == DIRICHLET_BOUNDARY &&
                            boundary_state_function(hs) &&
                            strcmp(boundary_state_function_name(hs),
                            "flowThroughBoundaryState") == 0)
                        {
                            U_nb[nb] = vel[l][index];
                        }
                        else
                        {
                            U_nb[nb] = getStateVel[l](intfc_state);
                        }
                        if (wave_type(hs) == DIRICHLET_BOUNDARY ||
                            neumann_type_bdry(wave_type(hs)))
                            mu[nb] = mu0;
                        else
                            mu[nb] = 1.0/2*(mu0 + field->mu[index_nb[nb]]);
                    }
                    else
                    {
                        U_nb[nb] = vel[l][index_nb[nb]];
                        mu[nb] = 1.0/2*(mu0 + field->mu[index_nb[nb]]);
                    }
                }

                coeff[0] = m_dt/rho*mu[0]/(top_h[0]*top_h[0]);
                coeff[1] = m_dt/rho*mu[1]/(top_h[0]*top_h[0]);
                coeff[2] = m_dt/rho*mu[2]/(top_h[1]*top_h[1]);
                coeff[3] = m_dt/rho*mu[3]/(top_h[1]*top_h[1]);

                getRectangleCenter(index, coords);
                computeSourceTerm(coords, source);

                aII = 1+coeff[0]+coeff[1]+coeff[2]+coeff[3];
                rhs = vel[l][index];

                for(nb = 0; nb < 4; nb++)
                {
                    if (!(*findStateAtCrossing)(front,icoords,dir[nb],comp,
                                &intfc_state,&hs,crx_coords))
                        solver.Set_A(I,I_nb[nb],-coeff[nb]);
                    else
                    {
                        if (wave_type(hs) == DIRICHLET_BOUNDARY &&
                            boundary_state_function(hs) &&
                            strcmp(boundary_state_function_name(hs),
                            "flowThroughBoundaryState") == 0)
                            aII -= coeff[nb];
                        else
                            rhs += coeff[nb]*U_nb[nb];
                    }
                }
                rhs += m_dt*source[l];
                rhs += m_dt*f_surf[l][index];
                //rhs -= m_dt*grad_q[l][index]/rho;
                solver.Set_A(I,I,aII);
                solver.Set_b(I, rhs);
            }

            solver.SetMaxIter(40000);
            solver.SetTol(1e-14);

            start_clock("Befor Petsc solve");
            //solver.Solve_GMRES();
            solver.Solve();
            solver.GetNumIterations(&num_iter);
            solver.GetFinalRelativeResidualNorm(&rel_residual);

            stop_clock("After Petsc solve");
            // get back the solution
            solver.Get_x(x);

            if (debugging("PETSc"))
                (void) printf("L_CARTESIAN::"
                        "computeDiffusionImplicit: "
                        "num_iter = %d, rel_residual = %g. \n",
                        num_iter,rel_residual);

            for (j = jmin; j <= jmax; j++)
            for (i = imin; i <= imax; i++)
            {
                I = ij_to_I[i][j];
                index = d_index2d(i,j,top_gmax);
                if (I >= 0)
                    vel[l][index] = x[I-ilower];
                else
                    vel[l][index] = 0.0;
            }
	}
	FT_ParallelExchGridVectorArrayBuffer(vel,front);
        FT_FreeThese(1,x);

        if (debugging("trace"))
            (void) printf("Leaving Incompress_Solver_Smooth_3D_Cartesian::"
                        "computeDiffusionImplicit()\n");
}       /* end computeDiffusionImplicit */

static int parab_find_state_at_crossing(
        Front *front,
        int *icoords,
        GRID_DIRECTION dir,
        int comp,
        POINTER *state,
        HYPER_SURF **hs,
        double *crx_coords)
{
        boolean status;
        INTERFACE *grid_intfc = front->grid_intfc;

        status = FT_StateStructAtGridCrossing(front,grid_intfc,icoords,dir,
                                comp,state,hs,crx_coords);
        if (status == NO) return NO_PDE_BOUNDARY;

        switch (wave_type(*hs))
        {
        case FIRST_PHYSICS_WAVE_TYPE:
            return NO_PDE_BOUNDARY;
        case DIRICHLET_BOUNDARY:
        case GROWING_BODY_BOUNDARY:
        case MOVABLE_BODY_BOUNDARY:
        case NEUMANN_BOUNDARY:
        case ICE_PARTICLE_BOUNDARY:
            return DIRICHLET_PDE_BOUNDARY;
        default:
            (void) printf("In parab_find_state_at_crossing()\n");
            (void) printf("Unknown wave type %s\n",
                                f_wave_type_as_string(wave_type(*hs)));
            clean_up(ERROR);
        }
}       /* end parab_find_state_at_crossing */

void Incompress_Solver_Smooth_2D_Cartesian::
        computeDiffusionParab(void)
{
        static PARABOLIC_SOLVER parab_solver(*front);

        COMPONENT comp;
        int index;
        int i,j,l,icoords[MAXD];
        double source[MAXD];
        double **vel = field->vel;
        double **f_surf = field->f_surf;
        static double *nu;
	static boolean first = YES;

        if (debugging("trace"))
            (void) printf("Entering Incompress_Solver_Smooth_3D_Cartesian::"
                        "computeDiffusionParab()\n");

        if (nu == NULL)
        {
            int size = 1;
            for (i = 0; i < dim; ++i)  size *= (top_gmax[i] + 1);
            FT_VectorMemoryAlloc((POINTER*)&nu,size,sizeof(double));
        }
        setIndexMap();
        for (j = jmin; j <= jmax; j++)
        for (i = imin; i <= imax; i++)
        {
            index  = d_index2d(i,j,top_gmax);
            nu[index] = field->mu[index]/field->rho[index];
        }
        FT_ParallelExchGridArrayBuffer(nu,front,NULL);

        parab_solver.soln_comp = LIQUID_COMP2;
        parab_solver.obst_comp = SOLID_COMP;
        parab_solver.ilower = ilower;
        parab_solver.iupper = iupper;
        parab_solver.dt = m_dt;
        parab_solver.order = 2;
        parab_solver.a = NULL;
        parab_solver.findStateAtCrossing = parab_find_state_at_crossing;
        parab_solver.first = first;
        parab_solver.set_solver_domain();
	first = NO;
        switch(dim)
        {
        case 2:
            parab_solver.ij_to_I = ij_to_I;
            break;
        case 3:
            parab_solver.ijk_to_I = ijk_to_I;
            break;
        }
        for (l = 0; l < dim; ++l)
        {
            parab_solver.var = vel[l];
            parab_solver.soln = vel[l];
            parab_solver.getStateVarFunc = getStateVel[l];
            parab_solver.solveIM();
        }

        if (debugging("trace"))
            (void) printf("Leaving Incompress_Solver_Smooth_3D_Cartesian::"
                        "computeDiffusionParab()\n");
}       /* end computeDiffusionParab */

void Incompress_Solver_Smooth_2D_Cartesian::
	computeDiffusionExplicit(void)
{
        COMPONENT comp;
        int index,index_nb[4],size;
        int I,I_nb[4];
	int i,j,l,nb,icoords[MAXD];
	double coords[MAXD], crx_coords[MAXD];
	double coeff[4],mu[4],mu0,rho,rhs,U_nb[4];
        double *x;
	GRID_DIRECTION dir[4] = {WEST,EAST,SOUTH,NORTH};
	POINTER intfc_state;
	HYPER_SURF *hs;
	PetscInt num_iter;
	double rel_residual;
	double aII;
	double source[MAXD];
	double **vel = field->vel;
	double **f_surf = field->f_surf;
	INTERFACE *grid_intfc = front->grid_intfc;

	if (debugging("trace"))
	    (void) printf("Entering Incompress_Solver_Smooth_2D_Cartesian::"
			"computeDiffusionExplicit()\n");

        setIndexMap();

        size = iupper - ilower;
        FT_VectorMemoryAlloc((POINTER*)&x,size,sizeof(double));

	for (l = 0; l < dim; ++l)
	{
            for (j = jmin; j <= jmax; j++)
            for (i = imin; i <= imax; i++)
            {
            	I  = ij_to_I[i][j];
            	if (I == -1) continue;

            	index  = d_index2d(i,j,top_gmax);
            	index_nb[0] = d_index2d(i-1,j,top_gmax);
            	index_nb[1] = d_index2d(i+1,j,top_gmax);
            	index_nb[2] = d_index2d(i,j-1,top_gmax);
            	index_nb[3] = d_index2d(i,j+1,top_gmax);

		icoords[0] = i;
		icoords[1] = j;
		comp = top_comp[index];

            	mu0   = field->mu[index];
            	rho   = field->rho[index];

            	for (nb = 0; nb < 4; nb++)
            	{
                    if (FT_StateStructAtGridCrossing(front,grid_intfc,icoords,
				dir[nb],comp,&intfc_state,&hs,crx_coords) &&
                                wave_type(hs) != FIRST_PHYSICS_WAVE_TYPE)
		    {
			if (wave_type(hs) == DIRICHLET_BOUNDARY &&
                    	    boundary_state_function(hs) &&
                    	    strcmp(boundary_state_function_name(hs),
                    	    "flowThroughBoundaryState") == 0)
                    	    U_nb[nb] = vel[l][index];
			else
			    U_nb[nb] = getStateVel[l](intfc_state);
			if (wave_type(hs) == DIRICHLET_BOUNDARY || 
			    neumann_type_bdry(wave_type(hs)))
			    mu[nb] = mu0;
			else
			    mu[nb] = 1.0/2*(mu0 + field->mu[index_nb[nb]]);
		    }
                    else
		    {
                    	U_nb[nb] = vel[l][index_nb[nb]];
			mu[nb] = 1.0/2*(mu0 + field->mu[index_nb[nb]]);
		    }
            	}

            	coeff[0] = m_dt/rho*mu[0]/(top_h[0]*top_h[0]);
            	coeff[1] = m_dt/rho*mu[1]/(top_h[0]*top_h[0]);
            	coeff[2] = m_dt/rho*mu[2]/(top_h[1]*top_h[1]);
            	coeff[3] = m_dt/rho*mu[3]/(top_h[1]*top_h[1]);

            	getRectangleCenter(index, coords);
            	computeSourceTerm(coords, source);

		rhs = (-coeff[0]-coeff[1]-coeff[2]-coeff[3])*vel[l][index];

		int num_nb = 0;
		for(nb = 0; nb < 4; nb++)
		{
		    rhs += coeff[nb]*U_nb[nb];
		}
		rhs += m_dt*source[l];
		rhs += m_dt*f_surf[l][index];
                x[I-ilower] = vel[l][index] + rhs;
            }

            for (j = jmin; j <= jmax; j++)
            for (i = imin; i <= imax; i++)
            {
                I = ij_to_I[i][j];
                index = d_index2d(i,j,top_gmax);
                if (I >= 0)
                    vel[l][index] = x[I-ilower];
                else
                    vel[l][index] = 0.0;
            }
	}
	FT_ParallelExchGridVectorArrayBuffer(vel,front);
        FT_FreeThese(1,x);

	if (debugging("trace"))
	    (void) printf("Leaving Incompress_Solver_Smooth_2D_Cartesian::"
			"computeDiffusionExplicit()\n");
}       /* end computeDiffusionExplicit */

void Incompress_Solver_Smooth_2D_Cartesian::vtk_plot_scalar(
        char *outname, const char* varname)
{
}	/* end vtk_plot_scalar */

void Incompress_Solver_Smooth_2D_Cartesian::extractFlowThroughVelocity()
{
	int index,index_nb,index_op;
	int i,j,k,l,nb,icoords[MAXD],icoords_nb[MAXD],icoords_op[MAXD];
	GRID_DIRECTION dir[3][2] = {{WEST,EAST},{SOUTH,NORTH},{LOWER,UPPER}};
	POINTER intfc_state;
        HYPER_SURF *hs;
	double crx_coords[MAXD];
	double vel_save,**vel = field->vel;
	double div;
	int status;

	/* Extract velocity for zero interior divergence */
        for (j = 0; j <= top_gmax[1]; ++j)
        for (i = 0; i <= top_gmax[0]; ++i)
	{
	    if (i >= imin && i <= imax && j >= jmin && j <= jmax)
		continue;

	    icoords[0] = i;
	    icoords[1] = j;
            index = d_index(icoords,top_gmax,dim);
            if (ifluid_comp(top_comp[index])) continue;
	    for (l = 0; l < dim; ++l)
		icoords_nb[l] = icoords_op[l] = icoords[l];
	    for (l = 0; l < dim; ++l)
	    for (nb = 0; nb < 2; ++nb)
	    {
		status = (*findStateAtCrossing)(front,icoords,dir[l][nb],
                                top_comp[index],&intfc_state,&hs,crx_coords);
		if (status == NO_PDE_BOUNDARY) continue;
		icoords_nb[l] = (nb == 0) ? icoords[l] - 1 : icoords[l] + 1;
		icoords_op[l] = (nb == 0) ? icoords[l] - 2 : icoords[l] + 2;
            	index_nb = d_index(icoords_nb,top_gmax,dim);
            	index_op = d_index(icoords_op,top_gmax,dim);
		if (!ifluid_comp(top_comp[index_nb]))
		    continue;
		vel[l][index] = 0.0;
                div = computeFieldPointDiv(icoords_nb,vel);
                vel[l][index] = (nb == 0) ? vel[l][index] - 2.0*top_h[l]*div
                                : vel[l][index] + 2.0*top_h[l]*div;
	    }
	}
}	/* end extractFlowThroughVelocity */


void Incompress_Solver_Smooth_2D_Cartesian::updateComponent(void)
{
	int i,j,l,icoords[MAXD];
	int index;

	/*update the component of pressure on dual grid*/
	for (j = jmin; j <= jmax; ++j)
        for (i = imin; i <= imax; ++i)
        {
            icoords[0] = i;
            icoords[1] = j;
            index = d_index(icoords,top_gmax,dim);
            if (ifluid_comp(top_comp[index]))
            {
                int cl[MAXD], cu[MAXD];
                for (l = 0; l < dim; ++l)
                    cl[l] = icoords[l] - offset[l];
                for (int m = 0; m < 2; ++m)
                for (int n = 0; n < 2; ++n)
                {
                    cu[0] = cl[0] + m;
                    cu[1] = cl[1] + n;
		    if (cu[0]<cimin||cu[0]>cimax||cu[1]<cjmin||cu[1]>cjmax)
			continue;
                    int index_tmp = d_index(cu,ctop_gmax,dim);
		    if (!ifluid_comp(ctop_comp[index_tmp]))
                    	ctop_comp[index_tmp] = top_comp[index];
                }
            }
        }

	/*Set rho for boundary layer on computational grid*/
	if (field->rho == NULL)
	    return;
	for (j = jmin; j <= jmax; ++j)
	for (i = imin; i <= imax; ++i)
	{
	    icoords[0] = i;
            icoords[1] = j;
            index = d_index(icoords,top_gmax,dim);
            if (!ifluid_comp(top_comp[index]))
	    {
		int cl[MAXD],cu[MAXD],indexl,indexu;
		boolean VelSet = NO;
		for (l = 0; l < dim && !VelSet; ++l)
		{
		    cl[l] = icoords[l]-1;
		    cu[l] = icoords[l]+1;
		    for (int n = -1; n < 2 && !VelSet; ++n)
		    {
			cl[(l+1)%dim] = cu[(l+1)%dim] = icoords[(l+1)%dim]+n;
			indexl = d_index(cl,top_gmax,dim);
                        if (ifluid_comp(top_comp[indexl]))
                        {
                            field->rho[index] = field->rho[indexl];
                            VelSet = YES;
                            continue;
                        }
                        indexu = d_index(cu,top_gmax,dim);
                        if (ifluid_comp(top_comp[indexu]))
                        {
                            field->rho[index] = field->rho[indexu];
                            VelSet = YES;
                            continue;
                        }
		    }
		}
	    }
	}
}	/* end updateComponent */

void Incompress_Solver_Smooth_2D_Cartesian::computeVarIncrement(
	double *var_old,
	double *var_new,
	boolean use_dual_grid)
{
	int i,j,k,index,size;
	double mag;
	double Lnorm[3];

	if (use_dual_grid)
	{
	    ;	// To add dual grid computation.
	}
	else
	{
	    Lnorm[0] = Lnorm[1] = Lnorm[2] = 0.0;
	    size = (imax - imin + 1)*(jmax - jmin + 1);
	    for (j = jmin; j <= jmax; j++)
            for (i = imin; i <= imax; i++)
            {
                index = d_index2d(i,j,top_gmax);
		mag = fabs(var_new[index] - var_old[index]);
		Lnorm[0] += mag;
		Lnorm[1] += sqr(mag);
		if (Lnorm[2] < mag)
		{
		    Lnorm[2] = mag;
		}
            }
	    Lnorm[0] /= (imax - imin + 1)*(jmax - jmin + 1);
	    Lnorm[1] = sqrt(Lnorm[1]/size);
	    (void) printf("L-1 norm = %20.14f  L-1/dt = %20.14f\n",
					Lnorm[0],Lnorm[0]/m_dt);
	    (void) printf("L-2 norm = %20.14f  L-2/dt = %20.14f\n",
					Lnorm[1],Lnorm[1]/m_dt);
	    (void) printf("L-I norm = %20.14f  L-I/dt = %20.14f\n",
					Lnorm[2],Lnorm[2]/m_dt);
	    Lnorm[0] = 0.0;
	    for (j = jmin; j <= jmax; j++)
            for (i = imin; i <= imax; i++)
	    {
		mag = fabs(var_new[index]);
		Lnorm[0] += mag;
	    }
	    Lnorm[0] /= (imax - imin + 1)*(jmax - jmin + 1);
	    (void) printf("L-1 norm of old variable = %20.14f\n",Lnorm[0]);
	}
}	/* end computeVarIncrement */	

void Incompress_Solver_Smooth_2D_Cartesian::computeVelDivergence()
{
	double *div_U = field->div_U;
	double **vel = field->vel;
	int i,j,index,icoords[MAXD];

	/* Compute velocity divergence */
	for (j = jmin; j <= jmax; j++)
        for (i = imin; i <= imax; i++)
	{
	    icoords[0] = i;
	    icoords[1] = j;
	    index  = d_index(icoords,top_gmax,dim);
	    if (!ifluid_comp(top_comp[index]))
		div_U[index] = 0.0;
	    div_U[index] = computeFieldPointDiv(icoords,vel);
	}
}	/* end computeVelDivergence */
