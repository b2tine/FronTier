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

#include "MonteCarlo.h"

static void initXgraphPlots(char*,int*,double**);
static void initMovieFrames(char*,int*,double*,double*);
static void makeXgraphPlot(char*,int,double*,double*,int,double);
static void makeSolnMovie(char*,double*,double**,int,int,double,double,char*);
static void readMonteCarloParams(PARAMS*,char*);
static void plotDistFunction(PARAMS,char*);
static void goMonteCarlo(PARAMS,char*);
#if defined(__GSL__)
static void makeGIGdensityPlot(PARAMS*,char*,int,double*,double*,int,double);
static void makeGHdensityPlot(PARAMS*,char*,int,double*,double*,int,double);
#endif /* if defined(__GSL__) */

double (*random_func)(POINTER,unsigned short int*);

char *in_name,*out_name;

int main(int argc, char **argv)
{
	static Front front;
	static F_BASIC_DATA f_basic;
	PARAMS params;

	FT_Init(argc,argv,&f_basic);
	f_basic.dim = 1;
	in_name                 = f_basic.in_name;
	out_name                = f_basic.out_name;

	FT_ReadSpaceDomain(in_name,&f_basic);
	params.XL = f_basic.L[0];
	params.XU = f_basic.U[0];
	params.nx = f_basic.gmax[0];
	FT_StartUp(&front,&f_basic);
	readMonteCarloParams(&params,in_name);
	plotDistFunction(params,out_name);
	goMonteCarlo(params,out_name);

	clean_up(0);
}

static void plotDistFunction(
	PARAMS params,
	char *outname)
{
        double mu,sigma;
        int i,ii,M,N;
	double XL,XU,dx;
        double *x,*f;
        double xrd;
        unsigned short int xsubi[3] = {271,6253,176};
	POINTER pdf_params = params.pdf_params;

	M = params.nx;
	XL = params.XL;
	XU = params.XU;
	FT_VectorMemoryAlloc((POINTER*)&x,M,sizeof(double));
	FT_VectorMemoryAlloc((POINTER*)&f,M,sizeof(double));

	dx = (XU-XL)/M;

        for (i = 0; i < M; ++i)
        {
            x[i] = XL + (i + 0.5)*dx;
            f[i] = 0.0;
        }

	/* Increasing N will make the profile convergent */
        N = params.num_samples; 

        for (i = 0; i < N; ++i)
        {
            xrd = random_func(pdf_params,xsubi);
            ii = (int)floor((xrd - XL + 0.5*dx)/dx);
            if (ii >= 0 && ii < M)
                f[ii] += 1.0/dx/N;
        }
	makeXgraphPlot(outname,M,x,f,0,0.0);

#if defined(__GSL__)
	if (params.rand_type == GIG)
	    makeGIGdensityPlot(&params,outname,M,x,f,0,0.0);
	else if (params.rand_type == GH)
	    makeGHdensityPlot(&params,outname,M,x,f,0,0.0);
#endif /* if defined(__GSL__) */

	FT_FreeThese(2,x,f);

}	/* end plotDistFunction */

static void makeXgraphPlot(
	char *out_name,
	int N,
	double *x,
	double *u,
	int l,
	double time)
{
	char xname[256];
	FILE *xfile;
	int i;

	sprintf(xname,"%s/soln-%d.xg",out_name,l);
	xfile = fopen(xname,"w");
	fprintf(xfile,"color=%s\n","red");
	fprintf(xfile,"thickness = 2.0\n");
	for (i = 0; i < N; ++i)
	    	fprintf(xfile,"%f %f\n",x[i],u[i]);
	fclose(xfile);
}	/* end makeXgraphPlot */

static void makeSolnMovie(
	char *out_name,
	double *x,
	double **u,
	int N,
	int M,
	double Ts,
	double Te,
	char *caption)
{
	int l,i;
	double umin = HUGE;
	double umax = -HUGE;
	double xmin,xmax;
	double height;
	static boolean first = YES;
	char time_label[100];
	char gd_name[256];
	double dt,time;

	sprintf(gd_name,"%s/function.gif",out_name);
	for (l = 0; l < M; ++l)
	for (i = 0; i <= N; ++i)
	{
	    if (umin > u[l][i]) umin = u[l][i];
	    if (umax < u[l][i]) umax = u[l][i];
	}
	height = umax - umin;
	umin -= 0.15*height;	umax += 0.15*height;
	xmin = x[0];		xmax = x[N];

	sprintf(time_label,"t = %6.3f",Ts);
	gd_initplot(gd_name,caption,xmin,xmax,umin,umax,3);
	gd_plotdata(N+1,x,u[0]);
	gd_plotframe(time_label);
	if (M == 1) return;
	dt = (Te - Ts)/(double)(M - 1);
	for (l = 1; l < M; ++l)
	{
	    time = Ts + l*dt;
	    sprintf(time_label,"t = %6.3f",time);
	    gd_appendplot(gd_name,caption,xmin,xmax,umin,umax,3);
	    gd_plotdata(N+1,x,u[l]);
	    gd_plotframe(time_label);
	}
}	/* end makeSolnMovie */

static void initXgraphPlots(
	char *inname,
	int *N,
	double **time)
{
	char string[100];
	FILE *infile = fopen(inname,"r");
	static double *tstore;
	int i,n;
	double T0,T1,dT;

	CursorAfterString(infile,"Enter number of time frames:");
	fscanf(infile,"%d",&n);
        (void) printf(" %d\n",n);
	FT_VectorMemoryAlloc((POINTER*)&tstore,n,sizeof(double));
	CursorAfterString(infile,"Enter yes for non-uniform time sequence:");
	fscanf(infile,"%s",string);
        (void) printf(" %s\n",string);
	if (string[0] == 'y' || string[0] == 'Y')
	{
	    for (i = 0; i < n; ++i)
	    {
		sprintf(string,"Enter time for %d-th frame:",i);
		CursorAfterString(infile,string);
		fscanf(infile,"%lf",&tstore[i]);
		printf(" %f\n",tstore[i]);
	    }
	}
	else
	{
	    CursorAfterString(infile,"Enter start time:");
	    fscanf(infile,"%lf",&T0);
	    printf("%f\n",T0);
	    CursorAfterString(infile,"Enter end time:");
	    fscanf(infile,"%lf",&T1);
	    printf("%f\n",T1);
	    dT = (T1 - T0)/(n-1);
	    for (i = 0; i < n; ++i)
		tstore[i] = T0 + n*dT;
	}
	fclose(infile);

	*N = n;
	*time = tstore;
}	/* initXgraphPlots */

static void initMovieFrames(
	char *inname,
	int *N,
	double *Ts,
	double *Te)
{
	char string[100];
	FILE *infile = fopen(inname,"r");

	CursorAfterString(infile,"Enter number of movie frames:");
	fscanf(infile,"%d",N);
        (void) printf(" %d\n",*N);
	CursorAfterString(infile,"Enter start time:");
	fscanf(infile,"%lf",Ts);
	printf("%f\n",*Ts);
	CursorAfterString(infile,"Enter end time:");
	fscanf(infile,"%lf",Te);
	printf("%f\n",*Te);
	fclose(infile);
}	/* initXgraphPlots */

static void readMonteCarloParams(
	PARAMS *params,
	char *inname)
{
	static GAUSS_PARAMS *gauss_params;
	static EXP_PARAMS *exp_params;
	static POWER_PARAMS *power_params;
	static UNIFORM_PARAMS *uniform_params;
	static STABLE_PARAMS *stable_params;
	static GIG_PARAMS *gig_params;
	static GH_PARAMS *gh_params;
	char string[100];
        FILE *infile = fopen(inname,"r");

	CursorAfterString(infile,"Enter name of distribution function:");
        fscanf(infile,"%s",string);
        (void) printf(" %s\n",string);
	switch (string[0])
	{
	case 'g':
	case 'G':
	    switch (string[1])
	    {
	    case 'a':
	    case 'A':
	    	FT_ScalarMemoryAlloc((POINTER*)&gauss_params,
			sizeof(GAUSS_PARAMS));
	    	switch (string[6])
	    	{
	    	case 'n':
	    	case 'N':
	    	    params->rand_type = GAUSS_NEWTON;
		    break;
	    	case 'b':
	    	case 'B':
	    	    params->rand_type = GAUSS_BM;
		    break;
	    	case 'c':
	    	case 'C':
	    	    params->rand_type = GAUSS_CL;
		    break;
	    	}
	    	CursorAfterString(infile,"Enter mathematical expectation:");
            	fscanf(infile,"%lf",&gauss_params->mu);
            	(void) printf(" %f\n",gauss_params->mu);
	    	CursorAfterString(infile,"Enter standard deviation:");
            	fscanf(infile,"%lf",&gauss_params->sigma);
            	(void) printf(" %f\n",gauss_params->sigma);
	    	params->pdf_params = (POINTER)gauss_params;
		break;
	    case 'i':
	    case 'I':
		FT_ScalarMemoryAlloc((POINTER*)&gig_params,sizeof(GIG_PARAMS));
		params->rand_type = GIG;
		CursorAfterString(infile,"Enter parameter lambda:");
                fscanf(infile,"%lf",&gig_params->lambda);
                (void) printf(" %f\n",gig_params->lambda);
                CursorAfterString(infile,"Enter parameter psi:");
                fscanf(infile,"%lf",&gig_params->psi);
                (void) printf(" %f\n",gig_params->psi);
		CursorAfterString(infile,"Enter parameter chi:");
                fscanf(infile,"%lf",&gig_params->chi);
                (void) printf(" %f\n",gig_params->chi);
                params->pdf_params = (POINTER)gig_params;
		break;
	    case 'h':
            case 'H':
                FT_ScalarMemoryAlloc((POINTER*)&gh_params,sizeof(GH_PARAMS));
		params->rand_type = GH;
                CursorAfterString(infile,"Enter parameter lambda:");
                fscanf(infile,"%lf",&gh_params->lambda);
                (void) printf(" %f\n",gh_params->lambda);
                CursorAfterString(infile,"Enter parameter alpha:");
                fscanf(infile,"%lf",&gh_params->alpha);
                (void) printf(" %f\n",gh_params->alpha);
                CursorAfterString(infile,"Enter parameter beta:");
                fscanf(infile,"%lf",&gh_params->beta);
                (void) printf(" %f\n",gh_params->beta);
		CursorAfterString(infile,"Enter parameter mu:");
                fscanf(infile,"%lf",&gh_params->mu);
                (void) printf(" %f\n",gh_params->mu);
		CursorAfterString(infile,"Enter parameter delta:");
                fscanf(infile,"%lf",&gh_params->delta);
                (void) printf(" %f\n",gh_params->delta);
                params->pdf_params = (POINTER)gh_params;
                break;
	    }
	    break;
	case 'e':
	case 'E':
	    FT_ScalarMemoryAlloc((POINTER*)&exp_params,sizeof(EXP_PARAMS));
	    params->rand_type = EXPONENTIAL;
	    CursorAfterString(infile,"Enter parameter lambda:");
            fscanf(infile,"%lf",&exp_params->lambda);
            (void) printf(" %f\n",exp_params->lambda);
	    params->pdf_params = (POINTER)exp_params;
	    break;
	case 'c':
	case 'C':
	    params->rand_type = CAUCHY;
	    FT_ScalarMemoryAlloc((POINTER*)&gauss_params,sizeof(GAUSS_PARAMS));
	    CursorAfterString(infile,"Enter mathematical expectation:");
            fscanf(infile,"%lf",&gauss_params->mu);
            (void) printf(" %f\n",gauss_params->mu);
	    CursorAfterString(infile,"Enter standard deviation:");
            fscanf(infile,"%lf",&gauss_params->sigma);
            (void) printf(" %f\n",gauss_params->sigma);
	    params->pdf_params = (POINTER)gauss_params;
	    break;
	case 'p':
	case 'P':
	    FT_ScalarMemoryAlloc((POINTER*)&power_params,sizeof(POWER_PARAMS));
	    params->rand_type = POWER;
	    CursorAfterString(infile,"Enter the power:");
            fscanf(infile,"%d",&power_params->power);
            (void) printf(" %d\n",power_params->power);
	    params->pdf_params = (POINTER)power_params;
	    break;
	case 'u':
	case 'U':
	    FT_ScalarMemoryAlloc((POINTER*)&uniform_params,
				sizeof(UNIFORM_PARAMS));
	    params->rand_type = UNIFORM;
	    CursorAfterString(infile,"Enter the lower and upper bounds:");
            fscanf(infile,"%lf %lf",&uniform_params->a,&uniform_params->b);
            (void) printf(" %f %f\n",uniform_params->a,uniform_params->b);
	    params->pdf_params = (POINTER)uniform_params;
	    break;
	case 's':
	case 'S':
	    FT_ScalarMemoryAlloc((POINTER*)&stable_params,
				sizeof(STABLE_PARAMS));
	    params->rand_type = STABLE;
	    CursorAfterString(infile,"Enter alpha:");
            fscanf(infile,"%lf",&stable_params->alpha);
            (void) printf(" %f\n",stable_params->alpha);
	    CursorAfterString(infile,"Enter beta:");
            fscanf(infile,"%lf",&stable_params->beta);
            (void) printf(" %f\n",stable_params->beta);
	    CursorAfterString(infile,"Enter sigma:");
            fscanf(infile,"%lf",&stable_params->sigma);
            (void) printf(" %f\n",stable_params->sigma);
	    CursorAfterString(infile,"Enter mu:");
            fscanf(infile,"%lf",&stable_params->mu);
            (void) printf(" %f\n",stable_params->mu);
	    params->pdf_params = (POINTER)stable_params;
	    break;
	case 'm':
	case 'M':
	    params->rand_type = MIDDLE;
	    break;
	default:
	    (void) printf("Unknown random variable type!\n");
	    clean_up(ERROR);
	}
	CursorAfterString(infile,"Enter number of samples:");
        fscanf(infile,"%d",&params->num_samples);
        (void) printf(" %d\n",params->num_samples);

	switch (params->rand_type)
	{
	case GAUSS_NEWTON:
	    random_func = gauss_newton;
	    break;
	case GAUSS_BM:
	    random_func = gauss_box_muller;
	    break;
	case GAUSS_CL:
	    random_func = gauss_center_limit;
	    break;
	case EXPONENTIAL:
	    random_func = dist_exponential;
	    break;
	case CAUCHY:
	    random_func = dist_cauchy;
	    break;
	case POWER:
	    random_func = dist_power;
	    break;
	case MIDDLE:
	    random_func = dist_middle;
	    break;
	case UNIFORM:
	    random_func = dist_uniform;
	    break;
	case STABLE:
	    random_func = dist_stable;
	    break;
	case GIG:
	    random_func = dist_gig;
	    break;
	case GH:
	    random_func = dist_gh;
	    break;
	default:
	    (void) printf("Unknown random variable\n");
	    clean_up(ERROR);
	}
	CursorAfterString(infile,"Enter type of random seed:");
        fscanf(infile,"%s",string);
        (void) printf(" %s\n",string);
	switch (string[0])
	{
	case 'f':
	case 'F':
	    params->seed_type = FIXED_SEED;
	    break;
	case 'r':
	case 'R':
	    params->seed_type = RANDOM_SEED;
	    break;
	case 'i':
	case 'I':
	    params->seed_type = INPUT_SEED;
	    break;
	}
	if (params->seed_type == INPUT_SEED)
	{
	    CursorAfterString(infile,"Enter three random seeds:");
            fscanf(infile,"%hu %hu %hu",&params->seeds[0],&params->seeds[1],
					&params->seeds[2]);
            (void) printf(" %d %d %d\n",params->seeds[0],params->seeds[1],
					params->seeds[2]);
	}
	/* The following is for stock simulation */
	params->do_monte_carlo = NO;
	params->do_option_price = NO;
	if (CursorAfterStringOpt(infile,"Enter yes for stock simulation:"))
	{
            fscanf(infile,"%s",string);
            (void) printf(" %s\n",string);
	    if (string[0] == 'y' || string[0] == 'Y')
	    {
		params->do_monte_carlo = YES;
	    	CursorAfterString(infile,"Enter market growth rate:");
            	fscanf(infile,"%lf",&params->mu);
		(void) printf(" %f\n",params->mu);
	    	CursorAfterString(infile,"Enter market volatility:");
            	fscanf(infile,"%lf",&params->sigma);
		(void) printf(" %f\n",params->sigma);
	    	CursorAfterString(infile,"Enter stock start price:");
            	fscanf(infile,"%lf",&params->S0);
		(void) printf(" %f\n",params->S0);
	    	CursorAfterString(infile,"Enter ending time:");
            	fscanf(infile,"%lf",&params->T);
		(void) printf(" %f\n",params->T);
	    	CursorAfterString(infile,"Enter number of time steps:");
            	fscanf(infile,"%d",&params->num_steps);
		(void) printf(" %d\n",params->num_steps);
	    	CursorAfterString(infile,"Enter number of simulations:");
            	fscanf(infile,"%d",&params->num_sims);
		(void) printf(" %d\n",params->num_sims);
		params->print_detail = NO;
	    	if (CursorAfterStringOpt(infile,"Enter yes to print details:"))
		{
            	    fscanf(infile,"%s",string);
		    (void) printf(" %s\n",string);
		    if (string[0] == 'y' || string[0] == 'Y')
		    	params->print_detail = YES;
	    	    CursorAfterString(infile,"Enter number of print cases:");
            	    fscanf(infile,"%d",&params->num_print_sims);
		    (void) printf(" %d\n",params->num_print_sims);
		}
		if (CursorAfterStringOpt(infile,
			"Enter yes to calculate option proce:"))
		{
            	    fscanf(infile,"%s",string);
            	    (void) printf(" %s\n",string);
	    	    if (string[0] == 'y' || string[0] == 'Y')
	    	    {
			params->do_option_price = YES;
	    		CursorAfterString(infile,"Enter strike price:");
            		fscanf(infile,"%lf",&params->E);
			(void) printf(" %f\n",params->E);
		    }
		}
	    }
	}

	fclose(infile);
}	/* end readMonteCarloParams */

static void assignRandomSeeds(
	PARAMS *params)
{
	time_t seconds;

	switch (params->seed_type)
	{
	case INPUT_SEED:
	    return;
	case FIXED_SEED:
	    params->seeds[0] = 271;
	    params->seeds[1] = 6253;
	    params->seeds[2] = 176;
	    return;
	case RANDOM_SEED:
	    seconds = time(NULL);
	    params->seeds[0] = seconds%9829;
	    params->seeds[1] = seconds%9883;
	    params->seeds[2] = seconds%9743;
	    return;
	}
}	/* end assignRandomSeeds */

static void goMonteCarlo(
	PARAMS params,
	char *outname)
{
	int n,N = params.num_steps;
	int i,num_sims = params.num_sims;
	double T = params.T;
	double *time,dX,dt = T/N;
	char xname[200];
	FILE *xfile;
	double dS,S0,S,*S2,*S_ave,*S_bank;
	double GS,*GS2,*GS_ave;
	double sigma = params.sigma;
	double *sample_S,*sample_GS;
	double mu = params.mu;
	POINTER pdf_params = params.pdf_params;
	boolean print_detail = params.print_detail;
	int num_print_sims = params.num_print_sims;
	double *C,*C2,*P,*P2,*sample_C,*sample_P;	// Option price
	double E = params.E;
	double *emut;			// exp(-mu*time)
	double randv;

	if (params.do_monte_carlo == NO) return;

	FT_VectorMemoryAlloc((POINTER*)&S_ave,N+1,sizeof(double));
	FT_VectorMemoryAlloc((POINTER*)&S2,N+1,sizeof(double));
	FT_VectorMemoryAlloc((POINTER*)&GS_ave,N+1,sizeof(double));
	FT_VectorMemoryAlloc((POINTER*)&GS2,N+1,sizeof(double));
	FT_VectorMemoryAlloc((POINTER*)&S_bank,N+1,sizeof(double));
	FT_VectorMemoryAlloc((POINTER*)&sample_S,N+1,sizeof(double));
	FT_VectorMemoryAlloc((POINTER*)&sample_GS,N+1,sizeof(double));
	FT_VectorMemoryAlloc((POINTER*)&time,N+1,sizeof(double));
	if (params.do_option_price == YES)
	{
	    FT_VectorMemoryAlloc((POINTER*)&C,N+1,sizeof(double));
	    FT_VectorMemoryAlloc((POINTER*)&C2,N+1,sizeof(double));
	    FT_VectorMemoryAlloc((POINTER*)&sample_C,N+1,sizeof(double));
	    FT_VectorMemoryAlloc((POINTER*)&P,N+1,sizeof(double));
	    FT_VectorMemoryAlloc((POINTER*)&P2,N+1,sizeof(double));
	    FT_VectorMemoryAlloc((POINTER*)&sample_P,N+1,sizeof(double));
	    FT_VectorMemoryAlloc((POINTER*)&emut,N+1,sizeof(double));
	}
	assignRandomSeeds(&params);
	for (n = 0; n <= N; ++n)
	{
	    time[n] = n*dt;
	    S_ave[n] = S2[n] = 0.0;
	    if (params.do_option_price == YES)
	    {
		C[n] = C2[n] = P[n] = P2[n] = 0.0;
		emut[n] = exp(-time[n]*mu);
	    }
	}
	for (i = 0; i < num_sims; ++i)
	{
	    GS = S = S0 = S_bank[0] = params.S0;
	    S_ave[0] += S0;
	    S2[0] += sqr(S0);
	    GS_ave[0] += S0;
	    GS2[0] += sqr(S0);
	    if (params.do_option_price == YES)
	    {
		if (S-E > 0.0)
		{
		    C[0] += emut[0]*(S - E);
		    C2[0] += sqr(emut[0]*(S-E));
		}
		else
		{
		    P[0] += emut[0]*(E-S);
		    P2[0] += sqr(emut[0]*(S-E));
		}
	    }
	    if (print_detail && i < num_print_sims)
	    {
	    	sprintf(xname,"%s/stock-%d.xg",out_name,i);
	    	xfile = fopen(xname,"w");
	    	fprintf(xfile,"color = %s\n",xcolor[(i+1)%max_color]);
	    	fprintf(xfile,"thickness = 1.5\n");
	    	fprintf(xfile,"%f %f\n",time[0],S);
	    }
	    for (n = 1; n <= N; ++n)
	    {
		randv = random_func(pdf_params,params.seeds);
		dX = sqrt(dt)*randv;
		dS = S*(sigma*dX + mu*dt);
		S += dS;
		S_ave[n] += S;
		S2[n] += sqr(S);

		GS = S0*exp((mu - 0.5*sqr(sigma))*time[n] + 
				sigma*randv*sqrt(time[n]));
		GS_ave[n] += GS;
		GS2[n] += sqr(GS);

		S_bank[n] = S_bank[n-1]*(1.0 + mu*dt);

	    	if (params.do_option_price == YES)
		{
		    if (S-E > 0.0)
		    {
			C[n] += emut[n]*(S - E);
			C2[n] += sqr(emut[n]*(S-E));
		    }
		    else
		    {
			P[n] += emut[n]*(E-S);
			P2[n] += sqr(emut[n]*(S-E));
		    }
		}
	    	if (print_detail && i < num_print_sims)
		    fprintf(xfile,"%f %f\n",time[n],S);
	    }
	    if (print_detail && i < num_print_sims)
	    {
	    	fprintf(xfile,"\n");
	    	fclose(xfile);
	    }
	}

	for (n = 0; n <= N; ++n)
	{
	    S_ave[n] /= num_sims;
	    GS_ave[n] /= num_sims;
	    sample_S[n] = S2[n] - num_sims*sqr(S_ave[n]);
	    sample_S[n] /= num_sims*(num_sims - 1);
	    sample_S[n] = sqrt(sample_S[n]);
	    sample_GS[n] = GS2[n] - num_sims*sqr(GS_ave[n]);
	    sample_GS[n] /= num_sims*(num_sims - 1);
	    sample_GS[n] = sqrt(sample_GS[n]);
	    fprintf(xfile,"%f %f\n",time[n],S_ave[n]);
	    if (params.do_option_price == YES)
	    {
		C[n] /= num_sims;
		P[n] /= num_sims;
		sample_C[n] = fabs(C2[n] - num_sims*sqr(C[n]));
		sample_C[n] /= num_sims*(num_sims - 1);
            	sample_C[n] = sqrt(sample_C[n]);
		sample_P[n] = fabs(P2[n] - num_sims*sqr(P[n]));
		sample_P[n] /= num_sims*(num_sims - 1);
            	sample_P[n] = sqrt(sample_P[n]);
	    }
	}
	sprintf(xname,"%s/stock-ave.xg",out_name);
	xfile = fopen(xname,"w");
	fprintf(xfile,"color = %s\n",xcolor[0]);
	fprintf(xfile,"thickness = 3.0\n");
	for (n = 0; n <= N; ++n)
	    fprintf(xfile,"%f %f\n",time[n],S_ave[n]);
	fprintf(xfile,"Next\n");
	fprintf(xfile,"color = %s\n",xcolor[0]);
	fprintf(xfile,"thickness = 2.0\n");
	for (n = 0; n <= N; ++n)
	    fprintf(xfile,"%f %f\n",time[n],S_ave[n]-sample_S[n]);
	fprintf(xfile,"Next\n");
	fprintf(xfile,"color = %s\n",xcolor[0]);
	fprintf(xfile,"thickness = 2.0\n");
	for (n = 0; n <= N; ++n)
	    fprintf(xfile,"%f %f\n",time[n],S_ave[n]+sample_S[n]);
	fprintf(xfile,"Next\n");
	fprintf(xfile,"color = %s\n",xcolor[0]);
	fprintf(xfile,"thickness = 2.0\n");
	for (n = 0; n <= N; ++n)
	    fprintf(xfile,"%f %f\n",time[n],S_ave[n]-2.0*sample_S[n]);
	fprintf(xfile,"Next\n");
	fprintf(xfile,"color = %s\n",xcolor[0]);
	fprintf(xfile,"thickness = 2.0\n");
	for (n = 0; n <= N; ++n)
	    fprintf(xfile,"%f %f\n",time[n],S_ave[n]+2.0*sample_S[n]);
	fclose(xfile);

	sprintf(xname,"%s/gb-stock-ave.xg",out_name);
	xfile = fopen(xname,"w");
	fprintf(xfile,"color = %s\n",xcolor[max_color-1]);
	fprintf(xfile,"thickness = 3.0\n");
	for (n = 0; n <= N; ++n)
	    fprintf(xfile,"%f %f\n",time[n],GS_ave[n]);
	fprintf(xfile,"Next\n");
	fprintf(xfile,"color = %s\n",xcolor[max_color-1]);
	fprintf(xfile,"thickness = 2.0\n");
	for (n = 0; n <= N; ++n)
	    fprintf(xfile,"%f %f\n",time[n],GS_ave[n]-sample_GS[n]);
	fprintf(xfile,"Next\n");
	fprintf(xfile,"color = %s\n",xcolor[max_color-1]);
	fprintf(xfile,"thickness = 2.0\n");
	for (n = 0; n <= N; ++n)
	    fprintf(xfile,"%f %f\n",time[n],GS_ave[n]+sample_GS[n]);
	fprintf(xfile,"Next\n");
	fprintf(xfile,"color = %s\n",xcolor[max_color-1]);
	fprintf(xfile,"thickness = 2.0\n");
	for (n = 0; n <= N; ++n)
	    fprintf(xfile,"%f %f\n",time[n],GS_ave[n]-2.0*sample_GS[n]);
	fprintf(xfile,"Next\n");
	fprintf(xfile,"color = %s\n",xcolor[max_color-1]);
	fprintf(xfile,"thickness = 2.0\n");
	for (n = 0; n <= N; ++n)
	    fprintf(xfile,"%f %f\n",time[n],GS_ave[n]+2.0*sample_GS[n]);
	fclose(xfile);

	if (params.do_option_price == YES)
	{
	    sprintf(xname,"%s/call-ave.xg",out_name);
	    xfile = fopen(xname,"w");
	    fprintf(xfile,"Next\n");
	    fprintf(xfile,"color = navy\n");
	    fprintf(xfile,"thickness = 3.0\n");
	    for (n = 0; n <= N; ++n)
	    	fprintf(xfile,"%f %f\n",time[n],C[n]);
	    fprintf(xfile,"Next\n");
	    fprintf(xfile,"color = navy\n");
	    fprintf(xfile,"thickness = 2.0\n");
	    for (n = 0; n <= N; ++n)
	    	fprintf(xfile,"%f %f\n",time[n],C[n]-sample_C[n]);
	    fprintf(xfile,"Next\n");
	    fprintf(xfile,"color = navy\n");
	    fprintf(xfile,"thickness = 2.0\n");
	    for (n = 0; n <= N; ++n)
	    	fprintf(xfile,"%f %f\n",time[n],C[n]+sample_C[n]);
	    fprintf(xfile,"Next\n");
	    fprintf(xfile,"color = navy\n");
	    fprintf(xfile,"thickness = 2.0\n");
	    for (n = 0; n <= N; ++n)
	    	fprintf(xfile,"%f %f\n",time[n],C[n]-2.0*sample_C[n]);
	    fprintf(xfile,"Next\n");
	    fprintf(xfile,"color = navy\n");
	    fprintf(xfile,"thickness = 2.0\n");
	    for (n = 0; n <= N; ++n)
	    	fprintf(xfile,"%f %f\n",time[n],C[n]+2.0*sample_C[n]);
	    fclose(xfile);

	    sprintf(xname,"%s/put-ave.xg",out_name);
	    xfile = fopen(xname,"w");
	    fprintf(xfile,"color = orange\n");
	    fprintf(xfile,"thickness = 3.0\n");
	    for (n = 0; n <= N; ++n)
	    	fprintf(xfile,"%f %f\n",time[n],P[n]);
	    fprintf(xfile,"Next\n");
	    fprintf(xfile,"color = orange\n");
	    fprintf(xfile,"thickness = 2.0\n");
	    for (n = 0; n <= N; ++n)
	    	fprintf(xfile,"%f %f\n",time[n],P[n]-sample_P[n]);
	    fprintf(xfile,"Next\n");
	    fprintf(xfile,"color = orange\n");
	    fprintf(xfile,"thickness = 2.0\n");
	    for (n = 0; n <= N; ++n)
	    	fprintf(xfile,"%f %f\n",time[n],P[n]+sample_P[n]);
	    fprintf(xfile,"Next\n");
	    fprintf(xfile,"color = orange\n");
	    fprintf(xfile,"thickness = 2.0\n");
	    for (n = 0; n <= N; ++n)
	    	fprintf(xfile,"%f %f\n",time[n],P[n]-2.0*sample_P[n]);
	    fprintf(xfile,"Next\n");
	    fprintf(xfile,"color = orange\n");
	    fprintf(xfile,"thickness = 2.0\n");
	    for (n = 0; n <= N; ++n)
	    	fprintf(xfile,"%f %f\n",time[n],P[n]+2.0*sample_P[n]);
	    fclose(xfile);
	}

	sprintf(xname,"%s/savings.xg",out_name);
	xfile = fopen(xname,"w");
	fprintf(xfile,"color = %s\n",xcolor[max_color/2]);
	fprintf(xfile,"thickness = 3.0\n");
	for (n = 0; n <= N; ++n)
	    fprintf(xfile,"%f %f\n",time[n],S_bank[n]);
	fclose(xfile);
	FT_FreeThese(2,S_ave,S_bank);
}	/* end goMonteCarlo */

#if defined(__GSL__)
static void makeGIGdensityPlot(
	PARAMS *params,
	char *out_name,
	int N,
	double *x,
	double *u,
	int l,
	double time)
{
	char xname[256];
	FILE *xfile;
	int i;

	GIG_PARAMS *gig_params = (GIG_PARAMS*)params->pdf_params;	
	double lambda = gig_params->lambda;
	double psi = gig_params->psi;
	double chi = gig_params->chi;
	double coeff = 0.5*pow(psi/chi,lambda/2.0)
			/gsl_sf_bessel_Knu(fabs(lambda),sqrt(chi*psi));

	sprintf(xname,"%s/exact-%d.xg",out_name,l);
	xfile = fopen(xname,"w");
	fprintf(xfile,"color=%s\n","blue");
	fprintf(xfile,"\"time = %6.3f\"\n",time);
	for (i = 0; i < N; ++i)
	{
	    u[i] = coeff*pow(x[i],lambda-1.0)*exp(-0.5*(chi/x[i]+psi*x[i]));
	    fprintf(xfile,"%f %f\n",x[i],u[i]);
	}
	fclose(xfile);
}	/* end makeGIGdensityPlot */

static void makeGHdensityPlot(
	PARAMS *params,
	char *out_name,
	int N,
	double *x,
	double *u,
	int l,
	double time)
{
	char xname[256];
	FILE *xfile;
	int i;

	GH_PARAMS *gh_params = (GH_PARAMS*)params->pdf_params;	
	double lambda = gh_params->lambda;
	double alpha = gh_params->alpha;
	double beta = gh_params->beta;
	double delta = gh_params->delta;
	double mu = gh_params->mu;
	double gamma = sqrt(alpha*alpha-beta*beta);
	double coeff = pow(gamma/delta,lambda)/sqrt(2*PI)
			/gsl_sf_bessel_Knu(fabs(lambda),gamma*delta);

	sprintf(xname,"%s/exact-%d.xg",out_name,l);
	xfile = fopen(xname,"w");
	fprintf(xfile,"color=%s\n","blue");
	fprintf(xfile,"\"time = %6.3f\"\n",time);
	for (i = 0; i < N; ++i)
	{
	    double temp = sqrt(delta*delta+(x[i]-mu)*(x[i]-mu));
	    u[i] = coeff*exp(beta*(x[i]-mu))
		/pow(temp/alpha,0.5-lambda)
		*gsl_sf_bessel_Knu(fabs(lambda-0.5),alpha*temp);
	    fprintf(xfile,"%f %f\n",x[i],u[i]);
	}
	fclose(xfile);
}	/* end makeGIGdensityPlot */
#endif /* if defined(__GSL__) */
