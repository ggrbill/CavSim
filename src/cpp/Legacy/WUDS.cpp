#include <math.h>

#include "WUDS.hpp"

void Calc_WUDS_Coef_X(
    double rho,
    double mi,
    int nv,
    double dx,
    double dy,
    double** u,
    double** v,
    CVBoundaries** alpha_x,
    CVBoundaries** beta_x)
{

    CVBoundaries Re;
	// Left-bottom corner
	Re.e = (rho*((u[0][0]+u[1][0])/2.)*dx)/mi;
	Re.n = (rho*((v[0][0]+v[1][0])/2.)*dy)/mi;
	Re.w = (rho*(u[0][0]/2.)*dx)/mi;
	Re.s = 0.;
	alpha_x[0][0].e = pow(Re.e,2.)/(10.+2.*pow(Re.e,2.));
	alpha_x[0][0].w = pow(Re.w,2.)/(10.+2.*pow(Re.w,2.));
	alpha_x[0][0].n = pow(Re.n,2.)/(10.+2.*pow(Re.n,2.));
	alpha_x[0][0].s = pow(Re.s,2.)/(10.+2.*pow(Re.s,2.));
	beta_x[0][0].e = (1.+(0.005*pow(Re.e,2.)))/(1.+(0.05*pow(Re.e,2.)));
	beta_x[0][0].w = (1.+(0.005*pow(Re.w,2.)))/(1.+(0.05*pow(Re.w,2.)));
	beta_x[0][0].n = (1.+(0.005*pow(Re.n,2.)))/(1.+(0.05*pow(Re.n,2.)));
	beta_x[0][0].s = (1.+(0.005*pow(Re.s,2.)))/(1.+(0.05*pow(Re.s,2.)));

	if(Re.e < 0) { alpha_x[0][0].e = -alpha_x[0][0].e; }
	if(Re.w < 0) { alpha_x[0][0].w = -alpha_x[0][0].w; }
	if(Re.n < 0) { alpha_x[0][0].n = -alpha_x[0][0].n; }
	if(Re.s < 0) { alpha_x[0][0].s = -alpha_x[0][0].s; }

	// Right-bottom corner
	Re.e = (rho*((u[nv-2][0])/2.)*dx)/mi;
	Re.n = (rho*((v[nv-2][0]+v[nv-1][0])/2.)*dx)/mi;
	Re.w = (rho*((u[nv-2][0]+u[nv-3][0])/2.)*dx)/mi; 
	Re.s = 0.;
	alpha_x[nv-2][0].e = pow(Re.e,2.)/(10.+2.*pow(Re.e,2.));
	alpha_x[nv-2][0].w = pow(Re.w,2.)/(10.+2.*pow(Re.w,2.));
	alpha_x[nv-2][0].n = pow(Re.n,2.)/(10.+2.*pow(Re.n,2.));
	alpha_x[nv-2][0].s = pow(Re.s,2.)/(10.+2.*pow(Re.s,2.));
	beta_x[nv-2][0].e = (1.+(0.005*pow(Re.e,2.)))/(1.+(0.05*pow(Re.e,2.)));
	beta_x[nv-2][0].w = (1.+(0.005*pow(Re.w,2.)))/(1.+(0.05*pow(Re.w,2.)));
	beta_x[nv-2][0].n = (1.+(0.005*pow(Re.n,2.)))/(1.+(0.05*pow(Re.n,2.)));
	beta_x[nv-2][0].s = (1.+(0.005*pow(Re.s,2.)))/(1.+(0.05*pow(Re.s,2.)));

	if(Re.e < 0) { alpha_x[nv-2][0].e = -alpha_x[nv-2][0].e; }
	if(Re.w < 0) { alpha_x[nv-2][0].w = -alpha_x[nv-2][0].w; }
	if(Re.n < 0) { alpha_x[nv-2][0].n = -alpha_x[nv-2][0].n; }
	if(Re.s < 0) { alpha_x[nv-2][0].s = -alpha_x[nv-2][0].s; }

	// Left-up corner
	Re.e = (rho*((u[0][nv-1]+u[1][nv-1])/2.)*dx)/mi;
	Re.n = 0.; //(rho*(U)*dx)/mi;
	Re.w = (rho*(u[0][nv-1]/2.)*dx)/mi;
	Re.s = (rho*((v[0][nv-2]+v[1][nv-2])/2.)*dx)/mi;
	alpha_x[0][nv-1].e = pow(Re.e,2.)/(10.+2.*pow(Re.e,2.));
	alpha_x[0][nv-1].w = pow(Re.w,2.)/(10.+2.*pow(Re.w,2.));
	alpha_x[0][nv-1].n = pow(Re.n,2.)/(10.+2.*pow(Re.n,2.));
	alpha_x[0][nv-1].s = pow(Re.s,2.)/(10.+2.*pow(Re.s,2.));
	beta_x[0][nv-1].e = (1.+(0.005*pow(Re.e,2.)))/(1.+(0.05*pow(Re.e,2.)));
	beta_x[0][nv-1].w = (1.+(0.005*pow(Re.w,2.)))/(1.+(0.05*pow(Re.w,2.)));
	beta_x[0][nv-1].n = (1.+(0.005*pow(Re.n,2.)))/(1.+(0.05*pow(Re.n,2.)));
	beta_x[0][nv-1].s = (1.+(0.005*pow(Re.s,2.)))/(1.+(0.05*pow(Re.s,2.)));

	if(Re.e < 0) { alpha_x[0][nv-1].e = -alpha_x[0][nv-1].e; }
	if(Re.w < 0) { alpha_x[0][nv-1].w = -alpha_x[0][nv-1].w; }
	if(Re.n < 0) { alpha_x[0][nv-1].n = -alpha_x[0][nv-1].n; }
	if(Re.s < 0) { alpha_x[0][nv-1].s = -alpha_x[0][nv-1].s; }

	// Right-up corner
	Re.e = (rho*((u[nv-2][nv-1])/2.)*dx)/mi;
	Re.n = 0.; //(rho*(U)*dx)/mi;
	Re.w = (rho*((u[nv-2][nv-1]+u[nv-3][nv-1])/2.)*dx)/mi;
	Re.s = (rho*((v[nv-2][nv-2]+v[nv-1][nv-2])/2.)*dx)/mi;
	alpha_x[nv-2][nv-1].e = pow(Re.e,2.)/(10.+2.*pow(Re.e,2.));
	alpha_x[nv-2][nv-1].w = pow(Re.w,2.)/(10.+2.*pow(Re.w,2.));
	alpha_x[nv-2][nv-1].n = pow(Re.n,2.)/(10.+2.*pow(Re.n,2.));
	alpha_x[nv-2][nv-1].s = pow(Re.s,2.)/(10.+2.*pow(Re.s,2.));
	beta_x[nv-2][nv-1].e = (1.+(0.005*pow(Re.e,2.)))/(1.+(0.05*pow(Re.e,2.)));
	beta_x[nv-2][nv-1].w = (1.+(0.005*pow(Re.w,2.)))/(1.+(0.05*pow(Re.w,2.)));
	beta_x[nv-2][nv-1].n = (1.+(0.005*pow(Re.n,2.)))/(1.+(0.05*pow(Re.n,2.)));
	beta_x[nv-2][nv-1].s = (1.+(0.005*pow(Re.s,2.)))/(1.+(0.05*pow(Re.s,2.)));

	if(Re.e < 0) { alpha_x[nv-2][nv-1].e = -alpha_x[nv-2][nv-1].e; }
	if(Re.w < 0) { alpha_x[nv-2][nv-1].w = -alpha_x[nv-2][nv-1].w; }
	if(Re.n < 0) { alpha_x[nv-2][nv-1].n = -alpha_x[nv-2][nv-1].n; }
	if(Re.s < 0) { alpha_x[nv-2][nv-1].s = -alpha_x[nv-2][nv-1].s; }

	// North boundary volumes
	for(int i=1;i<(nv-2);i++)
	{
		Re.e = (rho*((u[i][nv-1]+u[i+1][nv-1])/2.)*dx)/mi;
		Re.s = (rho*((v[i][nv-2]+v[i+1][nv-2])/2.)*dx)/mi;
		Re.w = (rho*((u[i][nv-1]+u[i-1][nv-1])/2.)*dx)/mi;
		Re.n = 0.; 
        alpha_x[i][nv-1].e = pow(Re.e,2.)/(10.+2.*pow(Re.e,2.));
		alpha_x[i][nv-1].w = pow(Re.w,2.)/(10.+2.*pow(Re.w,2.));
		alpha_x[i][nv-1].n = pow(Re.n,2.)/(10.+2.*pow(Re.n,2.));
		alpha_x[i][nv-1].s = pow(Re.s,2.)/(10.+2.*pow(Re.s,2.));
		beta_x[i][nv-1].e = (1.+(0.005*pow(Re.e,2.)))/(1.+(0.05*pow(Re.e,2.)));
		beta_x[i][nv-1].w = (1.+(0.005*pow(Re.w,2.)))/(1.+(0.05*pow(Re.w,2.)));
		beta_x[i][nv-1].n = (1.+(0.005*pow(Re.n,2.)))/(1.+(0.05*pow(Re.n,2.)));
		beta_x[i][nv-1].s = (1.+(0.005*pow(Re.s,2.)))/(1.+(0.05*pow(Re.s,2.)));

		if(Re.e < 0) { alpha_x[i][nv-1].e = -alpha_x[i][nv-1].e; }
		if(Re.w < 0) { alpha_x[i][nv-1].w = -alpha_x[i][nv-1].w; }
		if(Re.n < 0) { alpha_x[i][nv-1].n = -alpha_x[i][nv-1].n; }
		if(Re.s < 0) { alpha_x[i][nv-1].s = -alpha_x[i][nv-1].s; }

	}

	// South boundary volumes
	for(int i=1;i<(nv-2);i++)
	{
		Re.e = (rho*((u[i][0]+u[i+1][0])/2.)*dx)/mi;
		Re.n = (rho*((v[i][0]+v[i+1][0])/2.)*dx)/mi;
		Re.w = (rho*((u[i][0]+u[i-1][0])/2.)*dx)/mi;
		Re.s = 0.;
		alpha_x[i][0].e = pow(Re.e,2.)/(10.+2.*pow(Re.e,2.));
		alpha_x[i][0].w = pow(Re.w,2.)/(10.+2.*pow(Re.w,2.));
		alpha_x[i][0].n = pow(Re.n,2.)/(10.+2.*pow(Re.n,2.));
		alpha_x[i][0].s = pow(Re.s,2.)/(10.+2.*pow(Re.s,2.));
		beta_x[i][0].e = (1.+(0.005*pow(Re.e,2.)))/(1.+(0.05*pow(Re.e,2.)));
		beta_x[i][0].w = (1.+(0.005*pow(Re.w,2.)))/(1.+(0.05*pow(Re.w,2.)));
		beta_x[i][0].n = (1.+(0.005*pow(Re.n,2.)))/(1.+(0.05*pow(Re.n,2.)));
		beta_x[i][0].s = (1.+(0.005*pow(Re.s,2.)))/(1.+(0.05*pow(Re.s,2.)));

		if(Re.e < 0) { alpha_x[i][0].e = -alpha_x[i][0].e; }
		if(Re.w < 0) { alpha_x[i][0].w = -alpha_x[i][0].w; }
		if(Re.n < 0) { alpha_x[i][0].n = -alpha_x[i][0].n; }
		if(Re.s < 0) { alpha_x[i][0].s = -alpha_x[i][0].s; }
		
	}
	
	// East boundary volumes
	for(int j=1;j<(nv-1);j++)
	{
		Re.e = (rho*((u[nv-2][j])/2.)*dx)/mi;
		Re.n = (rho*((v[nv-2][j]+v[nv-1][j])/2.)*dx)/mi;
		Re.w = (rho*((u[nv-2][j]+u[nv-3][j])/2.)*dx)/mi;
		Re.s = (rho*((v[nv-2][j-1]+v[nv-1][j-1])/2.)*dx)/mi;
		alpha_x[nv-2][j].e = pow(Re.e,2.)/(10.+2.*pow(Re.e,2.));
		alpha_x[nv-2][j].w = pow(Re.w,2.)/(10.+2.*pow(Re.w,2.));
		alpha_x[nv-2][j].n = pow(Re.n,2.)/(10.+2.*pow(Re.n,2.));
		alpha_x[nv-2][j].s = pow(Re.s,2.)/(10.+2.*pow(Re.s,2.));
		beta_x[nv-2][j].e = (1.+(0.005*pow(Re.e,2.)))/(1.+(0.05*pow(Re.e,2.)));
		beta_x[nv-2][j].w = (1.+(0.005*pow(Re.w,2.)))/(1.+(0.05*pow(Re.w,2.)));
		beta_x[nv-2][j].n = (1.+(0.005*pow(Re.n,2.)))/(1.+(0.05*pow(Re.n,2.)));
		beta_x[nv-2][j].s = (1.+(0.005*pow(Re.s,2.)))/(1.+(0.05*pow(Re.s,2.)));

		if(Re.e < 0) { alpha_x[nv-2][j].e = -alpha_x[nv-2][j].e; }
		if(Re.w < 0) { alpha_x[nv-2][j].w = -alpha_x[nv-2][j].w; }
		if(Re.n < 0) { alpha_x[nv-2][j].n = -alpha_x[nv-2][j].n; }
		if(Re.s < 0) { alpha_x[nv-2][j].s = -alpha_x[nv-2][j].s; }
	}
	
	// West boundary volumes
	for(int j=1;j<(nv-1);j++)
	{
		Re.e = (rho*((u[0][j]+u[1][j])/2.)*dx)/mi;
		Re.n = (rho*((v[0][j]+v[1][j])/2.)*dx)/mi;
		Re.w = (rho*((u[0][j])/2.)*dx)/mi;
		Re.s = (rho*((v[0][j-1]+v[1][j-1])/2.)*dx)/mi;
		alpha_x[0][j].e = pow(Re.e,2.)/(10.+2.*pow(Re.e,2.));
		alpha_x[0][j].w = pow(Re.w,2.)/(10.+2.*pow(Re.w,2.));
		alpha_x[0][j].n = pow(Re.n,2.)/(10.+2.*pow(Re.n,2.));
		alpha_x[0][j].s = pow(Re.s,2.)/(10.+2.*pow(Re.s,2.));
		beta_x[0][j].e = (1.+(0.005*pow(Re.e,2.)))/(1.+(0.05*pow(Re.e,2.)));
		beta_x[0][j].w = (1.+(0.005*pow(Re.w,2.)))/(1.+(0.05*pow(Re.w,2.)));
		beta_x[0][j].n = (1.+(0.005*pow(Re.n,2.)))/(1.+(0.05*pow(Re.n,2.)));
		beta_x[0][j].s = (1.+(0.005*pow(Re.s,2.)))/(1.+(0.05*pow(Re.s,2.)));

		if(Re.e < 0) { alpha_x[0][j].e = -alpha_x[0][j].e; }
		if(Re.w < 0) { alpha_x[0][j].w = -alpha_x[0][j].w; }
		if(Re.n < 0) { alpha_x[0][j].n = -alpha_x[0][j].n; }
		if(Re.s < 0) { alpha_x[0][j].s = -alpha_x[0][j].s; }
	}
	// Core volumes
	for(int i=1;i<(nv-2);i++)
	{
		for(int j=1;j<(nv-1);j++)
		{
			Re.e = (rho*((u[i][j]+u[i+1][j])/2.)*dx)/mi;
			Re.n = (rho*((v[i+1][j]+v[i][j])/2.)*dy)/mi;
			Re.w = (rho*((u[i][j]+u[i-1][j])/2.)*dx)/mi;
			Re.s = (rho*((v[i][j-1]+v[i+1][j-1])/2.)*dy)/mi;
			alpha_x[i][j].e = pow(Re.e,2.)/(10.+2.*pow(Re.e,2.));
			alpha_x[i][j].w = pow(Re.w,2.)/(10.+2.*pow(Re.w,2.));
			alpha_x[i][j].n = pow(Re.n,2.)/(10.+2.*pow(Re.n,2.));
			alpha_x[i][j].s = pow(Re.s,2.)/(10.+2.*pow(Re.s,2.));
			beta_x[i][j].e = (1.+(0.005*pow(Re.e,2.)))/(1.+(0.05*pow(Re.e,2.)));
			beta_x[i][j].w = (1.+(0.005*pow(Re.w,2.)))/(1.+(0.05*pow(Re.w,2.)));
			beta_x[i][j].n = (1.+(0.005*pow(Re.n,2.)))/(1.+(0.05*pow(Re.n,2.)));
			beta_x[i][j].s = (1.+(0.005*pow(Re.s,2.)))/(1.+(0.05*pow(Re.s,2.)));
			
			if(Re.e < 0) { alpha_x[i][j].e = -alpha_x[i][j].e; }
			if(Re.w < 0) { alpha_x[i][j].w = -alpha_x[i][j].w; }
			if(Re.n < 0) { alpha_x[i][j].n = -alpha_x[i][j].n; }
			if(Re.s < 0) { alpha_x[i][j].s = -alpha_x[i][j].s; }
		}
	}
}

void Calc_WUDS_Coef_Y(
    double rho,
    double mi,
    int nv,
    double dx,
    double dy,
    double** u,
    double** v,
    CVBoundaries** alpha_y,
    CVBoundaries** beta_y)
{

    CVBoundaries Re;
	// Left-bottom corner
	Re.e = (rho*((u[0][0]+u[0][1])/2.)*dy)/mi;
	Re.n = (rho*((v[0][1]+v[0][0])/2.)*dy)/mi;
	Re.w = 0.;
	Re.s = (rho*((v[0][0])/2.)*dy)/mi;
	alpha_y[0][0].e = pow(Re.e,2.)/(10.+2.*pow(Re.e,2.));
	alpha_y[0][0].w = pow(Re.w,2.)/(10.+2.*pow(Re.w,2.));
	alpha_y[0][0].n = pow(Re.n,2.)/(10.+2.*pow(Re.n,2.));
	alpha_y[0][0].s = pow(Re.s,2.)/(10.+2.*pow(Re.s,2.));
	beta_y[0][0].e = (1.+(0.005*pow(Re.e,2.)))/(1.+(0.05*pow(Re.e,2.)));
	beta_y[0][0].w = (1.+(0.005*pow(Re.w,2.)))/(1.+(0.05*pow(Re.w,2.)));
	beta_y[0][0].n = (1.+(0.005*pow(Re.n,2.)))/(1.+(0.05*pow(Re.n,2.)));
	beta_y[0][0].s = (1.+(0.005*pow(Re.s,2.)))/(1.+(0.05*pow(Re.s,2.)));
		
	if(Re.e < 0) { alpha_y[0][0].e = -alpha_y[0][0].e; }
	if(Re.w < 0) { alpha_y[0][0].w = -alpha_y[0][0].w; }
	if(Re.n < 0) { alpha_y[0][0].n = -alpha_y[0][0].n; }
	if(Re.s < 0) { alpha_y[0][0].s = -alpha_y[0][0].s; }
	
	// Right-bottom corner
	Re.e = 0.;
	Re.n = (rho*((v[nv-1][1]+v[nv-1][0])/2.)*dy)/mi;
	Re.w = (rho*((u[nv-2][0]+u[nv-2][1])/2.)*dy)/mi;
	Re.s = (rho*((v[nv-1][0])/2.)*dy)/mi;
	alpha_y[nv-1][0].e = pow(Re.e,2.)/(10.+2.*pow(Re.e,2.));
	alpha_y[nv-1][0].w = pow(Re.w,2.)/(10.+2.*pow(Re.w,2.));
	alpha_y[nv-1][0].n = pow(Re.n,2.)/(10.+2.*pow(Re.n,2.));
	alpha_y[nv-1][0].s = pow(Re.s,2.)/(10.+2.*pow(Re.s,2.));
	beta_y[nv-1][0].e = (1.+(0.005*pow(Re.e,2.)))/(1.+(0.05*pow(Re.e,2.)));
	beta_y[nv-1][0].w = (1.+(0.005*pow(Re.w,2.)))/(1.+(0.05*pow(Re.w,2.)));
	beta_y[nv-1][0].n = (1.+(0.005*pow(Re.n,2.)))/(1.+(0.05*pow(Re.n,2.)));
	beta_y[nv-1][0].s = (1.+(0.005*pow(Re.s,2.)))/(1.+(0.05*pow(Re.s,2.)));
			
	if(Re.e < 0) { alpha_y[nv-1][0].e = -alpha_y[nv-1][0].e; }
	if(Re.w < 0) { alpha_y[nv-1][0].w = -alpha_y[nv-1][0].w; }
	if(Re.n < 0) { alpha_y[nv-1][0].n = -alpha_y[nv-1][0].n; }
	if(Re.s < 0) { alpha_y[nv-1][0].s = -alpha_y[nv-1][0].s; }
			
	// Left-up corner
	Re.e = (rho*((u[0][nv-2]+u[0][nv-1])/2.)*dy)/mi;
	Re.n = (rho*((v[0][nv-2])/2.)*dy)/mi;
	Re.w = 0.;
	Re.s = (rho*((v[0][nv-3]+v[0][nv-2])/2.)*dy)/mi;
	alpha_y[0][nv-2].e = pow(Re.e,2.)/(10.+2.*pow(Re.e,2.));
	alpha_y[0][nv-2].w = pow(Re.w,2.)/(10.+2.*pow(Re.w,2.));
	alpha_y[0][nv-2].n = pow(Re.n,2.)/(10.+2.*pow(Re.n,2.));
	alpha_y[0][nv-2].s = pow(Re.s,2.)/(10.+2.*pow(Re.s,2.));
	beta_y[0][nv-2].e = (1.+(0.005*pow(Re.e,2.)))/(1.+(0.05*pow(Re.e,2.)));
	beta_y[0][nv-2].w = (1.+(0.005*pow(Re.w,2.)))/(1.+(0.05*pow(Re.w,2.)));
	beta_y[0][nv-2].n = (1.+(0.005*pow(Re.n,2.)))/(1.+(0.05*pow(Re.n,2.)));
	beta_y[0][nv-2].s = (1.+(0.005*pow(Re.s,2.)))/(1.+(0.05*pow(Re.s,2.)));
			
	if(Re.e < 0) { alpha_y[0][nv-2].e = -alpha_y[0][nv-2].e; }
	if(Re.w < 0) { alpha_y[0][nv-2].w = -alpha_y[0][nv-2].w; }
	if(Re.n < 0) { alpha_y[0][nv-2].n = -alpha_y[0][nv-2].n; }
	if(Re.s < 0) { alpha_y[0][nv-2].s = -alpha_y[0][nv-2].s; }

	// Right-up corner
	Re.e = 0.;
	Re.n = (rho*((v[nv-1][nv-2])/2.)*dy)/mi;
	Re.w = (rho*((u[nv-2][nv-2]+u[nv-2][nv-1])/2.)*dy)/mi;
	Re.s = (rho*((v[nv-1][nv-3]+v[nv-1][nv-2])/2.)*dy)/mi;
	alpha_y[nv-1][nv-2].e = pow(Re.e,2.)/(10.+2.*pow(Re.e,2.));
	alpha_y[nv-1][nv-2].w = pow(Re.w,2.)/(10.+2.*pow(Re.w,2.));
	alpha_y[nv-1][nv-2].n = pow(Re.n,2.)/(10.+2.*pow(Re.n,2.));
	alpha_y[nv-1][nv-2].s = pow(Re.s,2.)/(10.+2.*pow(Re.s,2.));
	beta_y[nv-1][nv-2].e = (1.+(0.005*pow(Re.e,2.)))/(1.+(0.05*pow(Re.e,2.)));
	beta_y[nv-1][nv-2].w = (1.+(0.005*pow(Re.w,2.)))/(1.+(0.05*pow(Re.w,2.)));
	beta_y[nv-1][nv-2].n = (1.+(0.005*pow(Re.n,2.)))/(1.+(0.05*pow(Re.n,2.)));
	beta_y[nv-1][nv-2].s = (1.+(0.005*pow(Re.s,2.)))/(1.+(0.05*pow(Re.s,2.)));
			
	if(Re.e < 0) { alpha_y[nv-1][nv-2].e = -alpha_y[nv-1][nv-2].e; }
	if(Re.w < 0) { alpha_y[nv-1][nv-2].w = -alpha_y[nv-1][nv-2].w; }
	if(Re.n < 0) { alpha_y[nv-1][nv-2].n = -alpha_y[nv-1][nv-2].n; }
	if(Re.s < 0) { alpha_y[nv-1][nv-2].s = -alpha_y[nv-1][nv-2].s; }
			
	// North boundary volumes
	for(int i=1;i<(nv-1);i++)
	{
		Re.e = (rho*((u[i][nv-2]+u[i][nv-1])/2.)*dy)/mi;
		Re.n = (rho*((v[i][nv-2])/2.)*dy)/mi;
		Re.w = (rho*((u[i-1][nv-1]+u[i-1][nv-2])/2.)*dy)/mi;
		Re.s = (rho*((v[i][nv-3]+v[i][nv-2])/2.)*dy)/mi;
		alpha_y[i][nv-2].e = pow(Re.e,2.)/(10.+2.*pow(Re.e,2.));
		alpha_y[i][nv-2].w = pow(Re.w,2.)/(10.+2.*pow(Re.w,2.));
		alpha_y[i][nv-2].n = pow(Re.n,2.)/(10.+2.*pow(Re.n,2.));
		alpha_y[i][nv-2].s = pow(Re.s,2.)/(10.+2.*pow(Re.s,2.));
		beta_y[i][nv-2].e = (1.+(0.005*pow(Re.e,2.)))/(1.+(0.05*pow(Re.e,2.)));
		beta_y[i][nv-2].w = (1.+(0.005*pow(Re.w,2.)))/(1.+(0.05*pow(Re.w,2.)));
		beta_y[i][nv-2].n = (1.+(0.005*pow(Re.n,2.)))/(1.+(0.05*pow(Re.n,2.)));
		beta_y[i][nv-2].s = (1.+(0.005*pow(Re.s,2.)))/(1.+(0.05*pow(Re.s,2.)));
			
		if(Re.e < 0) { alpha_y[i][nv-2].e = -alpha_y[i][nv-2].e; }
		if(Re.w < 0) { alpha_y[i][nv-2].w = -alpha_y[i][nv-2].w; }
		if(Re.n < 0) { alpha_y[i][nv-2].n = -alpha_y[i][nv-2].n; }
		if(Re.s < 0) { alpha_y[i][nv-2].s = -alpha_y[i][nv-2].s; }
	}
	// South boundary volumes
	for(int i=1;i<(nv-1);i++)
	{
		Re.e = (rho*((u[i][0]+u[i][1])/2.)*dy)/mi;
		Re.n = (rho*((v[i][0]+v[i][1])/2.)*dy)/mi;
		Re.w = (rho*((u[i-1][0]+u[i-1][1])/2.)*dy)/mi;
		Re.s = (rho*((v[i][0])/2.)*dy)/mi;
		alpha_y[i][0].e = pow(Re.e,2.)/(10.+2.*pow(Re.e,2.));
		alpha_y[i][0].w = pow(Re.w,2.)/(10.+2.*pow(Re.w,2.));
		alpha_y[i][0].n = pow(Re.n,2.)/(10.+2.*pow(Re.n,2.));
		alpha_y[i][0].s = pow(Re.s,2.)/(10.+2.*pow(Re.s,2.));
		beta_y[i][0].e = (1.+(0.005*pow(Re.e,2.)))/(1.+(0.05*pow(Re.e,2.)));
		beta_y[i][0].w = (1.+(0.005*pow(Re.w,2.)))/(1.+(0.05*pow(Re.w,2.)));
		beta_y[i][0].n = (1.+(0.005*pow(Re.n,2.)))/(1.+(0.05*pow(Re.n,2.)));
		beta_y[i][0].s = (1.+(0.005*pow(Re.s,2.)))/(1.+(0.05*pow(Re.s,2.)));
		
		if(Re.e < 0) { alpha_y[i][0].e = -alpha_y[i][0].e; }
		if(Re.w < 0) { alpha_y[i][0].w = -alpha_y[i][0].w; }
		if(Re.n < 0) { alpha_y[i][0].n = -alpha_y[i][0].n; }
		if(Re.s < 0) { alpha_y[i][0].s = -alpha_y[i][0].s; }
	}
	// East boundary volumes
	for(int j=1;j<(nv-2);j++)
	{
		Re.e = 0.;
		Re.n = (rho*((v[nv-1][j+1]+v[nv-1][j])/2.)*dy)/mi;
		Re.w = (rho*((u[nv-2][j+1]+u[nv-2][j])/2.)*dy)/mi;
		Re.s = (rho*((v[nv-1][j-1]+v[nv-1][j])/2.)*dy)/mi;
		alpha_y[nv-1][j].e = pow(Re.e,2.)/(10.+2.*pow(Re.e,2.));
		alpha_y[nv-1][j].w = pow(Re.w,2.)/(10.+2.*pow(Re.w,2.));
		alpha_y[nv-1][j].n = pow(Re.n,2.)/(10.+2.*pow(Re.n,2.));
		alpha_y[nv-1][j].s = pow(Re.s,2.)/(10.+2.*pow(Re.s,2.));
		beta_y[nv-1][j].e = (1.+(0.005*pow(Re.e,2.)))/(1.+(0.05*pow(Re.e,2.)));
		beta_y[nv-1][j].w = (1.+(0.005*pow(Re.w,2.)))/(1.+(0.05*pow(Re.w,2.)));
		beta_y[nv-1][j].n = (1.+(0.005*pow(Re.n,2.)))/(1.+(0.05*pow(Re.n,2.)));
		beta_y[nv-1][j].s = (1.+(0.005*pow(Re.s,2.)))/(1.+(0.05*pow(Re.s,2.)));
			
		if(Re.e < 0) { alpha_y[nv-1][j].e = -alpha_y[nv-1][j].e; }
		if(Re.w < 0) { alpha_y[nv-1][j].w = -alpha_y[nv-1][j].w; }
		if(Re.n < 0) { alpha_y[nv-1][j].n = -alpha_y[nv-1][j].n; }
		if(Re.s < 0) { alpha_y[nv-1][j].s = -alpha_y[nv-1][j].s; }
	}
	// West boundary volumes
	for(int j=1;j<(nv-2);j++)
	{
		Re.e = (rho*((u[0][j]+u[0][j+1])/2.)*dy)/mi;
		Re.n = (rho*((v[0][j]+v[0][j+1])/2.)*dy)/mi;
		Re.w = 0.;
		Re.s = (rho*((v[0][j]+v[0][j-1])/2.)*dy)/mi;
		alpha_y[0][j].e = pow(Re.e,2.)/(10.+2.*pow(Re.e,2.));
		alpha_y[0][j].w = pow(Re.w,2.)/(10.+2.*pow(Re.w,2.));
		alpha_y[0][j].n = pow(Re.n,2.)/(10.+2.*pow(Re.n,2.));
		alpha_y[0][j].s = pow(Re.s,2.)/(10.+2.*pow(Re.s,2.));
		beta_y[0][j].e = (1.+(0.005*pow(Re.e,2.)))/(1.+(0.05*pow(Re.e,2.)));
		beta_y[0][j].w = (1.+(0.005*pow(Re.w,2.)))/(1.+(0.05*pow(Re.w,2.)));
		beta_y[0][j].n = (1.+(0.005*pow(Re.n,2.)))/(1.+(0.05*pow(Re.n,2.)));
		beta_y[0][j].s = (1.+(0.005*pow(Re.s,2.)))/(1.+(0.05*pow(Re.s,2.)));
		
		if(Re.e < 0) { alpha_y[0][j].e = -alpha_y[0][j].e; }
		if(Re.w < 0) { alpha_y[0][j].w = -alpha_y[0][j].w; }
		if(Re.n < 0) { alpha_y[0][j].n = -alpha_y[0][j].n; }
		if(Re.s < 0) { alpha_y[0][j].s = -alpha_y[0][j].s; }		
	}
	// Core volummes
	for(int i=1;i<(nv-1);i++)
	{
		for(int j=1;j<(nv-2);j++)
		{
			Re.e = (rho*((u[i][j]+u[i][j+1])/2.)*dy)/mi;
			Re.n = (rho*((v[i][j+1]+v[i][j])/2.)*dy)/mi;
			Re.w = (rho*((u[i-1][j]+u[i-1][j+1])/2.)*dy)/mi;
			Re.s = (rho*((v[i][j-1]+v[i][j])/2.)*dy)/mi;
			alpha_y[i][j].e = pow(Re.e,2.)/(10.+2.*pow(Re.e,2.));
			alpha_y[i][j].w = pow(Re.w,2.)/(10.+2.*pow(Re.w,2.));
			alpha_y[i][j].n = pow(Re.n,2.)/(10.+2.*pow(Re.n,2.));
			alpha_y[i][j].s = pow(Re.s,2.)/(10.+2.*pow(Re.s,2.));
			beta_y[i][j].e = (1.+(0.005*pow(Re.e,2.)))/(1.+(0.05*pow(Re.e,2.)));
			beta_y[i][j].w = (1.+(0.005*pow(Re.w,2.)))/(1.+(0.05*pow(Re.w,2.)));
			beta_y[i][j].n = (1.+(0.005*pow(Re.n,2.)))/(1.+(0.05*pow(Re.n,2.)));
			beta_y[i][j].s = (1.+(0.005*pow(Re.s,2.)))/(1.+(0.05*pow(Re.s,2.)));
			
			if(Re.e < 0) { alpha_y[i][j].e = -alpha_y[i][j].e; }
			if(Re.w < 0) { alpha_y[i][j].w = -alpha_y[i][j].w; }
			if(Re.n < 0) { alpha_y[i][j].n = -alpha_y[i][j].n; }
			if(Re.s < 0) { alpha_y[i][j].s = -alpha_y[i][j].s; }
		}
	}
	
}