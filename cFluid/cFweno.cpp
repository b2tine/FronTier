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

#include <cFluid.h>
#include <ctime>
#include <time.h>

#if defined(__GPU__)
#include "cFweno_gpu.cuh"
#endif

static double weno5_scal(double *f);
static void matmvec(double *b, double L[5][5], double *x);
static void f2is(double *f, double *s);
static void u2f(double *u, double* f);
static void weno5_get_flux(POINTER,int,int,double**, std::vector<std::vector<double> >);
static void arti_compression(POINTER,double*,double*,double,double,double*,int,double &c);

extern void WENO_flux(
        POINTER params,
        SWEEP *vst,
        FSWEEP *vflux,
        int n,
        int sizevst)
{
    double* u_old[6];

	u_old[0] = vst->dens;
	u_old[1] = vst->momn[0];
	u_old[2] = vst->momn[1];
	u_old[3] = vst->momn[2];
	u_old[4] = vst->engy;
	u_old[5] = vst->pres;

    std::vector<std::vector<double> > flux(5);

    flux[0].assign(vflux->dens_flux, vflux->dens_flux + sizevst);
    flux[1].assign(vflux->momn_flux[0], vflux->momn_flux[0] + sizevst);
    flux[2].assign(vflux->momn_flux[1], vflux->momn_flux[1] + sizevst);
    flux[3].assign(vflux->momn_flux[2], vflux->momn_flux[2] + sizevst);
    flux[4].assign(vflux->engy_flux, vflux->engy_flux + sizevst);

    int ghost_size = 3;
	int extend_size = n + 2*ghost_size;

	SCHEME_PARAMS *scheme_params = (SCHEME_PARAMS*)params;
	double lambda = scheme_params->lambda;

#if defined(__GPU__)
	//startClock("Total_time_flux_gpu");
	weno5_get_flux_gpu(scheme_params->gamma,scheme_params->lambda,
				extend_size,ghost_size,u_old,flux);
	//stopClock("Total_time_flux_gpu");
#else
	weno5_get_flux(params,extend_size,ghost_size,u_old,flux);
#endif

	for (int i = ghost_size; i < n+ghost_size; ++i)
	{
	    vflux->dens_flux[i] = -lambda*(flux[0][i+1] - flux[0][i]);
	    vflux->momn_flux[0][i] = -lambda*(flux[1][i+1] - flux[1][i]);
	    vflux->momn_flux[1][i] = -lambda*(flux[2][i+1] - flux[2][i]);
	    vflux->momn_flux[2][i] = -lambda*(flux[3][i+1] - flux[3][i]);
	    vflux->engy_flux[i] = -lambda*(flux[4][i+1] - flux[4][i]);
	
        //TODO: put in debug function
        if (isnan(vflux->dens_flux[i]))
	    {
            for (int j = 0; j < extend_size; ++j)
                printf("u = %f %f %f %f %f\n",
                        u_old[0][j],u_old[1][j],u_old[2][j],
                        u_old[3][j],u_old[4][j]);
    
            clean_up(ERROR);
	    }
	}

}	/* end weno5_flux */

static void weno5_get_flux(
	POINTER params,
	int extend_size, 
	int ghost_size, 
	double **u_old,
    std::vector<std::vector<double> > flux)
{

	SCHEME_PARAMS *scheme_params = (SCHEME_PARAMS*)params;
	double gamma = scheme_params->gamma;
    double gm = gamma - 1.0;

    double** f;
    FT_MatrixMemoryAlloc((POINTER*)&f,extend_size,5,sizeof(double));

    double maxeig[5] = {0 ,0 ,0, 0, 0};
    for(int i = 0; i < extend_size; ++i)
    {
        double a = sqrt(gamma * u_old[5][i]/u_old[0][i]);
        double v = u_old[1][i]/u_old[0][i];

        maxeig[0] = std::max(maxeig[0], fabs(v - a));
        maxeig[1] = std::max(maxeig[1], fabs(v));
        maxeig[4] = std::max(maxeig[4], fabs(v + a));
        
	    double u[6];
        for (int j = 0; j < 6; ++j)
            u[j] = u_old[j][i];

        u2f(u,f[i]);
    }

    maxeig[2] = maxeig[1];
    maxeig[3] = maxeig[1];

    //#pragma omp parallel for num_threads(4)
    for(int i = ghost_size; i < extend_size - ghost_size + 1; ++i)
    {
        /*** Get u_1/2 ***/

        double u_mid[11];	// rho,mx,my,mz,e,u,v,w,p,u^2+v^2+w^2,a
        for(int j = 0; j < 5; ++j)
        {
            u_mid[j] = 0.5*(u_old[j][i-1] + u_old[j][i]);
        }

        u_mid[5] = u_mid[1]/u_mid[0];
        u_mid[6] = u_mid[2]/u_mid[0];
        u_mid[7] = u_mid[3]/u_mid[0];
        u_mid[9] = sqr(u_mid[5]) + sqr(u_mid[6]) + sqr(u_mid[7]);
        u_mid[8] = (u_mid[4] - 0.5*u_mid[0]*u_mid[9])*(gamma - 1.0);
        u_mid[10] = sqrt(gamma*u_mid[8]/u_mid[0]);

        /*** R(u_1/2) & R^-1(u_1/2) ***/

        //h enthalpy
        double h = 0.5*u_mid[9] + gamma*u_mid[8]/(gamma - 1.0)/u_mid[0];

        double R[5][5];
 
        R[0][0] = 1.0;
        R[0][1] = 1.0;
        R[0][2] = 0.0;
        R[0][3] = 0.0;
        R[0][4] = 1.0;

        R[1][0] = u_mid[5] - u_mid[10];
        R[1][1] = u_mid[5];
        R[1][2] = 0.0;
        R[1][3] = 0.0;
        R[1][4] = u_mid[5] + u_mid[10];

        R[2][0] = u_mid[6];
        R[2][1] = u_mid[6];
        R[2][2] = 1.0;
        R[2][3] = 0.0;
        R[2][4] = u_mid[6];

        R[3][0] = u_mid[7];
        R[3][1] = u_mid[7];
        R[3][2] = 0.0;
        R[3][3] = 1.0;
        R[3][4] = u_mid[7];

        R[4][0] = h - u_mid[5]*u_mid[10];
        R[4][1] = 0.5*u_mid[9];
        R[4][2] = u_mid[6];
        R[4][3] = u_mid[7];
        R[4][4] = h + u_mid[5]*u_mid[10];;

        double  L[5][5];

        L[0][0] = h + u_mid[10]*(u_mid[5] - u_mid[10]) / gm;
        L[0][1] = -1.0*(u_mid[5] + u_mid[10] / gm);
        L[0][2] = -1.0*u_mid[6];
        L[0][3] = -1.0*u_mid[7];
        L[0][4] = 1.0;

        L[1][0] = -2.0*h + 4.0*sqr(u_mid[10])/gm;
        L[1][1] = 2.0*u_mid[5];
        L[1][2] = 2.0*u_mid[6];
        L[1][3] = 2.0*u_mid[7];
        L[1][4] = -2.0;

        L[2][0] = -2.0*u_mid[6]*sqr(u_mid[10])/gm;
        L[2][1] = 0;
        L[2][2] = 2.0*sqr(u_mid[10])/gm;
        L[2][3] = 0;
        L[2][4] = 0;

        L[3][0] = -2.0*u_mid[7]*sqr(u_mid[10])/gm;
        L[3][1] = 0;
        L[3][2] = 0;
        L[3][3] = 2.0*sqr(u_mid[10])/gm;
        L[3][4] = 0;

        L[4][0] = h - u_mid[10]*(u_mid[5] + u_mid[10])/gm;
        L[4][1] = -1.0*u_mid[5] + u_mid[10]/gm;
        L[4][2] = -1.0*u_mid[6];
        L[4][3] = -1.0*u_mid[7];
        L[4][4] = 1.0;

        for(int j = 0; j < 5; ++j)
        {
            for(int k = 0; k < 5; ++k)
                L[j][k] *= gm/(2.0*sqr(u_mid[10])); 
        }

        /*** Get R^-1 * u and R^-1 * F ***/	    

        double sten_u[6][5];
        double sten_f[6][5]; //f_-2, f_-1, f_0, f_1, f_2, f_3
     
        for(int j = 0; j < 6; ++j)
        {
            double uu[5];
            for (int k = 0; k < 5; ++k)
                uu[k] = u_old[k][i - ghost_size + j];

            matmvec(sten_u[j], L, uu);
            matmvec(sten_f[j], L, f[i-ghost_size+j]);
        }

        double gfluxp[5][5];
        double gfluxm[5][5]; 
 
        for(int j = 0; j < 5; ++j)
        {
            for(int k = 0; k < 5; ++k)
            {
                gfluxp[j][k] = 0.5*(sten_f[k][j] + maxeig[j]*sten_u[k][j]);
                gfluxm[j][k] = 0.5*(sten_f[5-k][j] - maxeig[j]*sten_u[5-k][j]);
            }
        }

	    double f_prevp[5];
        double f_nowp[5];
        double f_prevm[5];
        double f_nowm[5];
        double f_tmp[5];
	
	    double gflux_tmp[5];
        double vecp[5][4];
        double vecm[5][4];
        
        for(int j = 0; j < 5; ++j)
        {
            f_prevp[j] = f_nowp[j];
            f_nowp[j] = weno5_scal(gfluxp[j]);
            f_tmp[j] = f_nowp[j];
      
            /* artificial compression */
            for(int k = 0; k < 5; ++k)
                gflux_tmp[k] = 0.5*(sten_f[5-k][j] + maxeig[j]*sten_u[5-k][j]);
                
            double c;
            arti_compression(params,gfluxp[j],gflux_tmp,f_prevp[j],f_nowp[j],vecp[j],1,c);
            f_tmp[j] += c;
            /* end of artificial compression */

            f_prevm[j] = f_nowm[j];
            f_nowm[j] = weno5_scal(gfluxm[j]);
            f_tmp[j] += f_nowm[j];
     
            /* artificial compression */
            for(int k = 0; k < 5; ++k)
                gflux_tmp[k] = 0.5*(sten_f[k][j] - maxeig[j]*sten_u[k][j]);
    
            arti_compression(params,gfluxm[j],gflux_tmp,f_prevm[j],f_nowm[j],vecm[j],-1,c);
            f_tmp[j] += c;
            /* end of artificial compression */
        }

        double ff[5];
        matmvec(ff, R, f_tmp);
        
        for(int j = 0; j < 5; ++j)
        {
            //TODO: put in debugging function
            if ( isnan(ff[j]) )
            {
                (void) printf("In weno5_flux(): flux[%d][%d] = %f\n",
                        j, i, ff[j]);

                for (int k = 0; k < extend_size; ++k)
                {
                    printf("u[%d] = %f %f %f %f %f\n",
                            k,u_old[0][k],u_old[1][k],
                            u_old[2][k],u_old[3][k],u_old[4][k]);
                }

                clean_up(ERROR);
            }

            flux[j][i] = ff[j];
        }
       
    }

    FT_FreeThese(1,f);
}

static double weno5_scal(double *f)
{
    const double eps = 1.e-40;
    const int p = 2;

    double d[3] = {0.1, 0.6, 0.3}; //*** Optimal weights C_k 

    const double a[3][5] = { {1.0/3, -7.0/6, 11.0/6, 0, 0},
                        {0, -1.0/6, 5.0/6, 1.0/3, 0}, 
                        {0, 0, 1.0/3, 5.0/6, -1.0/6} }; 
    
    //*** coefficients for 2nd-order ENO interpolation stencil
    double w[5]; //weight for every point

    
    double is[3]; //*** a smoothness measurement of the flux function
    f2is(f,is);


    double alpha[3];
    double sum = 0.0;
    for(int i = 0; i < 3; ++i)
    {
        alpha[i] = d[i]/pow(eps + is[i],p);
        sum += alpha[i];
    }
    
    double omega[3]; // stencil weights
    for(int i = 0; i < 3; ++i)
    {
        omega[i] = alpha[i] / sum;
    }

    double flux = 0.0;
    for(int i = 0; i < 5; ++i)
    {
        w[i] = 0.0;
        for(int j = 0; j < 3; ++j)
        {
            w[i] += omega[j] * a[j][i];
        }
        
        flux += w[i] * f[i];
    }

	return flux;
}

static void f2is(
	double *f, 
	double *s)
{
	s[0] = 13.0/12*sqr((f[0] - 2.0*f[1] + f[2])) +
                0.25*sqr((f[0] - 4.0*f[1] + 3.0*f[2]));
        s[1] = 13.0/12*sqr((f[1] - 2.0*f[2] + f[3])) +
                0.25*sqr((f[1] - f[3]));
        s[2] = 13.0/12*sqr((f[2] - 2.0*f[3] + f[4])) +
                0.25*sqr((3.0*f[2] - 4.0*f[3] + f[4]));
}

static void matmvec(
	double *b, 
	double L[5][5], 
	double *x)
{
    	int i, j;

    	for(i = 0; i < 5; ++i)
    	{
            b[i] = 0.0;
            for(j = 0; j < 5; ++j)
	        {
                b[i] += L[i][j] * x[j]; 
	        }
        }
}

static void u2f(
	double* u,
    double* f)
{
	double v = u[1]/u[0];

    	f[0] = u[1];
    	f[1] = v*u[1] + u[5];
    	f[2] = v*u[2];
    	f[3] = v*u[3];
    	f[4] = v*(u[4] + u[5]);
}


//writes to vec and c
static void arti_compression(
	POINTER params,
	double *sten,
	double *sten_tmp,
	double f_prev,
	double f_now,
	double *vec,
	int sign,
	double &c)
{
	double tmp,f_ac_tmp;
	double ac_alpha,ac_alpha_nume,ac_alpha_deno;
	SCHEME_PARAMS *scheme_params = (SCHEME_PARAMS*)params; 

	c = 0.0;
	if (scheme_params->artificial_compression == NO)
        return;
	
    ac_alpha_nume = fabs(sten[3] - 2.0*sten[2] + sten[1]);
	ac_alpha_deno = fabs(sten[3] - sten[2]) + fabs(sten[2] - sten[1]);
	f_ac_tmp = weno5_scal(sten_tmp);

    //data race
	if (sign == 1)
	{
	    vec[1] = vec[0];
	    vec[0] = f_ac_tmp - f_now;
	    vec[2] = sten[3] - f_now;
	    vec[3] = vec[1] + f_prev - sten[1];
	}
	else if (sign == -1)
	{
	    vec[0] = vec[1];
	    vec[1] = f_ac_tmp - f_now;
	    vec[2] = sten[1] - f_prev;
        vec[3] = vec[1] + f_now - sten[3];
	}

	if (ac_alpha_deno != 0.0)
	{
	    ac_alpha = 33.0 * pow(ac_alpha_nume/ac_alpha_deno,2);
	    
        if (ac_alpha > 1.0)
        {
            if (vec[0] > 0.0 && vec[1] > 0.0)
                tmp = (vec[0] < vec[1]) ? vec[0] : vec[1];
            else if (vec[0] < 0.0 && vec[1] < 0.0)
                tmp = (vec[0] > vec[1]) ? vec[0] : vec[1];
            else
                tmp = 0.0;

            tmp *= ac_alpha / 2.0;
            if (tmp > 0.0 && vec[2] > 0.0 && vec[3] > 0.0)
            {
                c = (tmp < vec[2]) ? tmp : vec[2];
                c = (c < vec[3]) ? c : vec[3];
            }
            else if (tmp < 0.0 && vec[2] < 0.0 && vec[3] < 0.0)
            {
                c = (tmp > vec[2]) ? tmp : vec[2];
                c = (c > vec[3]) ? c : vec[3];
            }
            else
                c = 0.0;
        }
	}
}
