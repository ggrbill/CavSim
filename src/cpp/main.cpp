#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <math.h>

#include "Legacy/solver.hpp"
#include "Legacy/IO.hpp"
#include "Legacy/numeric.hpp"
#include "Legacy/Structures.hpp"
#include "Legacy/WUDS.hpp"

using namespace std;

// Resulta matrices
double **u;
double **v;
double **uOLD;
double **vOLD;
double **P;
double **Pn;
double **u_hat;
double **v_hat;

// Coefficients Matrices
// x-velocity - u
double **Ap_u;
double **Aw_u;
double **Ae_u;
double **An_u;
double **As_u;
double **B_u;
// y-velocity - v
double **Ap_v;
double **Aw_v;
double **Ae_v;
double **An_v;
double **As_v;
double **B_v;
// Pressure
double **Ap_p;
double **Aw_p;
double **Ae_p;
double **An_p;
double **As_p;
double **B_p;

CVBoundaries **alpha_x;
CVBoundaries **beta_x;
CVBoundaries **alpha_y;
CVBoundaries **beta_y;

// Input data
double L;	// Length 
double H;	// Height
int    nv;	// Number of rows and columns
double rho;	// density
double U;	// velocity at north boundary
double mi;	// viscosity

// Delta X e Delta Y
double dx = 0.;
double dy = 0.;

void Calc_Coef_NS_X()  //Calculate Coefficients for u - NS_x
{
	// Left-bottom corner
	Ae_u[0][0] = (((mi*beta_x[0][0].e*dy)/dx)-((rho*((u[0][0]+u[1][0])/2.)*dy)*(0.5-alpha_x[0][0].e)));
	An_u[0][0] = (((mi*beta_x[0][0].n*dx)/dy)-((rho*((v[1][0]+v[0][0])/2.)*dx)*(0.5-alpha_x[0][0].n)));
	Aw_u[0][0] = (((mi*beta_x[0][0].w*dy)/dx)+((rho*((u[0][0])/2.)*dy)*(0.5+alpha_x[0][0].w)));
	As_u[0][0] = ((2.*mi*beta_x[0][0].s*dx)/dy);
	Ap_u[0][0] = Ae_u[0][0] + Aw_u[0][0] + An_u[0][0] + As_u[0][0];
	B_u[0][0] = 0.;

	// Right-bottom corner
	Aw_u[nv-2][0] = (((mi*beta_x[nv-2][0].w*dy)/dx)+((rho*((u[nv-2][0]+u[nv-3][0])/2.)*dy)*(0.5+alpha_x[nv-2][0].w)));
	An_u[nv-2][0] = (((mi*beta_x[nv-2][0].n*dx)/dy)-((rho*((v[nv-2][0]+v[nv-1][0])/2.)*dx)*(0.5-alpha_x[nv-2][0].n)));
	Ae_u[nv-2][0] = (((mi*beta_x[nv-2][0].e*dy)/dx)-((rho*((u[nv-2][0])/2.)*dy)*(0.5-alpha_x[nv-2][0].e)));
	As_u[nv-2][0] = ((2.*mi*beta_x[nv-2][0].s*dx)/dy);
	Ap_u[nv-2][0] = Ae_u[nv-2][0] + Aw_u[nv-2][0] + An_u[nv-2][0] + As_u[nv-2][0];
	B_u[nv-2][0] = 0.;
	
	// Left-up corner
	Ae_u[0][nv-1] = (((mi*beta_x[0][nv-1].e*dy)/dx)-((rho*((u[0][nv-1]+u[1][nv-1])/2.)*dy)*(0.5-alpha_x[0][nv-1].e)));
	As_u[0][nv-1] = (((mi*beta_x[0][nv-1].s*dx)/dy)+((rho*((v[1][nv-2]+v[0][nv-2])/2.)*dx)*(0.5+alpha_x[0][nv-1].s)));
	Aw_u[0][nv-1] = (((mi*beta_x[0][nv-1].w*dy)/dx)+((rho*((u[0][nv-1])/2.)*dy)*(0.5+alpha_x[0][nv-1].w)));
	An_u[0][nv-1] = ((2.*mi*beta_x[0][nv-1].n*dx)/dy);
	Ap_u[0][nv-1] = Ae_u[0][nv-1] + Aw_u[0][nv-1] + An_u[0][nv-1] + As_u[0][nv-1];
	B_u[0][nv-1] = ((2.*mi*beta_x[0][nv-1].n*U*dx)/dy);
	
	// Right-up corner
	Aw_u[nv-2][nv-1] = (((mi*beta_x[nv-2][nv-1].w*dy)/dx)+((rho*((u[nv-2][nv-1]+u[nv-3][nv-1])/2.)*dy)*(0.5+alpha_x[nv-2][nv-1].w)));
	As_u[nv-2][nv-1] = (((mi*beta_x[nv-2][nv-1].s*dx)/dy)+((rho*((v[nv-2][nv-2]+v[nv-1][nv-2])/2.)*dx)*(0.5+alpha_x[nv-2][nv-1].s)));
	Ae_u[nv-2][nv-1] = (((mi*beta_x[nv-2][nv-1].e*dy)/dx)-((rho*((u[nv-2][nv-1])/2.)*dy)*(0.5-alpha_x[nv-2][nv-1].e)));
	An_u[nv-2][nv-1] = ((2.*mi*beta_x[nv-2][nv-1].n*dx)/dy);
	Ap_u[nv-2][nv-1] = Ae_u[nv-2][nv-1] + Aw_u[nv-2][nv-1] + An_u[nv-2][nv-1] + As_u[nv-2][nv-1];
	B_u[nv-2][nv-1] = ((2.*mi*beta_x[nv-2][nv-1].n*U*dx)/dy);

	// North boundary volumes
	for(int i=1;i<(nv-2);i++)
	{
		Ae_u[i][nv-1] = (((mi*beta_x[i][nv-1].e*dy)/dx)-((rho*((u[i][nv-1]+u[i+1][nv-1])/2.)*dy)*(0.5-alpha_x[i][nv-1].e)));
		As_u[i][nv-1] = (((mi*beta_x[i][nv-1].s*dx)/dy)+((rho*((v[i+1][nv-2]+v[i][nv-2])/2.)*dx)*(0.5+alpha_x[i][nv-1].s)));
		Aw_u[i][nv-1] = (((mi*beta_x[i][nv-1].w*dy)/dx)+((rho*((u[i-1][nv-1]+u[i][nv-1])/2.)*dy)*(0.5+alpha_x[i][nv-1].w)));
		An_u[i][nv-1] = ((2.*mi*beta_x[i][nv-1].n*dx)/dy);
		Ap_u[i][nv-1] = Ae_u[i][nv-1] + Aw_u[i][nv-1] + An_u[i][nv-1] + As_u[i][nv-1];
		B_u[i][nv-1] = ((2.*mi*beta_x[i][nv-1].n*U*dx)/dy);
	}

	// South boundary volumes
	for(int i=1;i<(nv-2);i++)
	{
		Ae_u[i][0] = (((mi*beta_x[i][0].e*dy)/dx)-((rho*((u[i][0]+u[i+1][0])/2.)*dy)*(0.5-alpha_x[i][0].e)));
		An_u[i][0] = (((mi*beta_x[i][0].n*dx)/dy)-((rho*((v[i+1][0]+v[i][0])/2.)*dx)*(0.5-alpha_x[i][0].n)));
		Aw_u[i][0] = (((mi*beta_x[i][0].w*dy)/dx)+((rho*((u[i-1][0]+u[i][0])/2.)*dy)*(0.5+alpha_x[i][0].w)));
		As_u[i][0] = ((2.*mi*beta_x[i][0].s*dx)/dy);
		Ap_u[i][0] = Ae_u[i][0] + Aw_u[i][0] + An_u[i][0] + As_u[i][0];
		B_u[i][0] = 0.;
	}
	
	// East boundary volumes
	for(int j=1;j<(nv-1);j++)
	{
		As_u[nv-2][j] = (((mi*beta_x[nv-2][j].s*dx)/dy)+((rho*((v[nv-2][j-1]+v[nv-1][j-1])/2.)*dx)*(0.5+alpha_x[nv-2][j].s)));
		An_u[nv-2][j] = (((mi*beta_x[nv-2][j].n*dx)/dy)-((rho*((v[nv-2][j]+v[nv-1][j])/2.)*dx)*(0.5-alpha_x[nv-2][j].n)));
		Aw_u[nv-2][j] = (((mi*beta_x[nv-2][j].w*dy)/dx)+((rho*((u[nv-2][j]+u[nv-3][j])/2.)*dy)*(0.5+alpha_x[nv-2][j].w)));
		Ae_u[nv-2][j] = ((mi*beta_x[nv-2][j].e*dx)/dy)-((rho*((u[nv-2][j]/2.))*dy)*(0.5-alpha_x[nv-2][j].e));
		Ap_u[nv-2][j] = Ae_u[nv-2][j] + Aw_u[nv-2][j] + An_u[nv-2][j] + As_u[nv-2][j];
		B_u[nv-2][j] = 0.;
	}
	
	// West boundary volumes
	for(int j=1;j<(nv-1);j++)
	{
		As_u[0][j] = (((mi*beta_x[0][j].s*dx)/dy)+((rho*((v[0][j-1]+v[1][j-1])/2.)*dx)*(0.5+alpha_x[0][j].s)));
		An_u[0][j] = (((mi*beta_x[0][j].n*dx)/dy)-((rho*((v[0][j]+v[1][j])/2.)*dx)*(0.5-alpha_x[0][j].n)));
		Ae_u[0][j] = (((mi*beta_x[0][j].e*dy)/dx)-((rho*((u[0][j]+u[1][j])/2.)*dy)*(0.5-alpha_x[0][j].e)));
		Aw_u[0][j] = (((mi*beta_x[0][j].w*dy)/dx)+((rho*((u[0][j]/2.))*dy)*(0.5+alpha_x[0][j].w)));
		Ap_u[0][j] = Ae_u[0][j] + Aw_u[0][j] + An_u[0][j] + As_u[0][j];
		B_u[0][j] = 0.;
	}
	// Core volumes
	for(int i=1;i<(nv-2);i++)
	{
		for(int j=1;j<(nv-1);j++)
		{
			As_u[i][j] = (((mi*beta_x[i][j].s*dx)/dy)+((rho*((v[i][j-1]+v[i+1][j-1])/2.)*dx)*(0.5+alpha_x[i][j].s)));
			An_u[i][j] = (((mi*beta_x[i][j].n*dx)/dy)-((rho*((v[i][j]+v[i+1][j])/2.)*dx)*(0.5-alpha_x[i][j].n)));
			Ae_u[i][j] = (((mi*beta_x[i][j].e*dy)/dx)-((rho*((u[i][j]+u[i+1][j])/2.)*dy)*(0.5-alpha_x[i][j].e)));
			Aw_u[i][j] = (((mi*beta_x[i][j].w*dy)/dx)+((rho*((u[i-1][j]+u[i][j])/2.)*dy)*(0.5+alpha_x[i][j].w)));
			Ap_u[i][j] = Ae_u[i][j] + Aw_u[i][j] + An_u[i][j] + As_u[i][j];
			B_u[i][j] = 0.;
		}
	}

}

void Calc_Coef_NS_Y()
{
	// Left-bottom corner
	As_v[0][0] = (((mi*beta_y[0][0].s*dx)/dy)+((rho*((v[0][0])/2.)*dx)*(0.5+alpha_y[0][0].s)));
	An_v[0][0] = (((mi*beta_y[0][0].n*dx)/dy)-((rho*((v[0][1]+v[0][0])/2.)*dx)*(0.5-alpha_y[0][0].n)));
	Ae_v[0][0] = (((mi*beta_y[0][0].e*dy)/dx)-((rho*((u[0][1]+u[0][0])/2.)*dy)*(0.5-alpha_y[0][0].e)));
	Aw_v[0][0] = ((2.*mi*beta_y[0][0].w*dy)/dx);
	Ap_v[0][0] = Ae_v[0][0] + Aw_v[0][0] + An_v[0][0] + As_v[0][0];
	B_v[0][0] = 0.;

	// Right-bottom corner
	As_v[nv-1][0] = (((mi*beta_y[nv-1][0].s*dx)/dy)+((rho*((v[nv-1][0])/2.)*dx)*(0.5+alpha_y[nv-1][0].s)));
	An_v[nv-1][0] = (((mi*beta_y[nv-1][0].n*dx)/dy)-((rho*((v[nv-1][0]+v[nv-1][1])/2.)*dx)*(0.5-alpha_y[nv-1][0].n)));
	Ae_v[nv-1][0] = ((2.*mi*beta_y[nv-1][0].e*dy)/dx);
	Aw_v[nv-1][0] = (((mi*beta_y[nv-1][0].w*dy)/dx)+((rho*((u[nv-2][0]+u[nv-2][1])/2.)*dy)*(0.5+alpha_y[nv-1][0].w)));
	Ap_v[nv-1][0] = Ae_v[nv-1][0] + Aw_v[nv-1][0] + An_v[nv-1][0] + As_v[nv-1][0];
	B_v[nv-1][0] = 0.;

	// Left-up corner
	As_v[0][nv-2] = (((mi*beta_y[0][nv-2].s*dx)/dy)+((rho*((v[0][nv-3]+v[0][nv-2])/2.)*dx)*(0.5+alpha_y[0][nv-2].s)));
	An_v[0][nv-2] = (((mi*beta_y[0][nv-2].n*dx)/dy)-((rho*((v[0][nv-2])/2.)*dx)*(0.5-alpha_y[0][nv-2].n)));
	Ae_v[0][nv-2] = (((mi*beta_y[0][nv-2].e*dy)/dx)-((rho*((u[0][nv-1]+u[0][nv-2])/2.)*dy)*(0.5-alpha_y[0][nv-2].e)));
	Aw_v[0][nv-2] = ((2.*mi*beta_y[0][nv-2].w*dy)/dx);
	Ap_v[0][nv-2] = Ae_v[0][nv-2] + Aw_v[0][nv-2] + An_v[0][nv-2] + As_v[0][nv-2];
	B_v[0][nv-2] = 0.;

	// Right-up corner
	As_v[nv-1][nv-2] = (((mi*beta_y[nv-1][nv-2].s*dx)/dy)+((rho*((v[nv-1][nv-3]+v[nv-1][nv-2])/2.)*dx)*(0.5+alpha_y[nv-1][nv-2].s)));
	An_v[nv-1][nv-2] = (((mi*beta_y[nv-1][nv-2].n*dx)/dy)-((rho*((v[nv-1][nv-2])/2.)*dx)*(0.5-alpha_y[nv-1][nv-2].n)));
	Ae_v[nv-1][nv-2] = ((2.*mi*beta_y[nv-1][nv-2].e*dy)/dx);
	Aw_v[nv-1][nv-2] = (((mi*beta_y[nv-1][nv-2].w*dy)/dx)+((rho*((u[nv-2][nv-1]+u[nv-2][nv-2])/2)*dy)*(0.5+alpha_y[nv-1][nv-2].w)));
	Ap_v[nv-1][nv-2] = Ae_v[nv-1][nv-2] + Aw_v[nv-1][nv-2] + An_v[nv-1][nv-2] + As_v[nv-1][nv-2];
	B_v[nv-1][nv-2] = 0.;

	// North boundary volumes
	for(int i=1;i<(nv-1);i++)
	{
		As_v[i][nv-2] = (((mi*beta_y[i][nv-2].s*dx)/dy)+((rho*((v[i][nv-3]+v[i][nv-2])/2.)*dx)*(0.5+alpha_y[i][nv-2].s)));
		An_v[i][nv-2] = (((mi*beta_y[i][nv-2].n*dx)/dy)-((rho*((v[i][nv-2])/2.)*dx)*(0.5-alpha_y[i][nv-2].n)));
		Ae_v[i][nv-2] = (((mi*beta_y[i][nv-2].e*dy)/dx)-((rho*((u[i][nv-1]+u[i][nv-2])/2.)*dy)*(0.5-alpha_y[i][nv-2].e)));
		Aw_v[i][nv-2] = (((mi*beta_y[i][nv-2].w*dy)/dx)+((rho*((u[i-1][nv-1]+u[i-1][nv-2])/2.)*dy)*(0.5+alpha_y[i][nv-2].w)));
		Ap_v[i][nv-2] = Ae_v[i][nv-2] + Aw_v[i][nv-2] + An_v[i][nv-2] + As_v[i][nv-2];
		B_v[i][nv-2] = 0.;
	}
	// South boundary volumes
	for(int i=1;i<(nv-1);i++)
	{
		As_v[i][0] = (((mi*beta_y[i][0].s*dx)/dy)+((rho*((v[i][0])/2.)*dx)*(0.5+alpha_y[i][0].s)));
		An_v[i][0] = (((mi*beta_y[i][0].n*dx)/dy)-((rho*((v[i][1]+v[i][0])/2.)*dx)*(0.5-alpha_y[i][0].n)));
		Ae_v[i][0] = (((mi*beta_y[i][0].e*dy)/dx)-((rho*((u[i][1]+u[i][0])/2.)*dy)*(0.5-alpha_y[i][0].e)));
		Aw_v[i][0] = (((mi*beta_y[i][0].w*dy)/dx)+((rho*((u[i-1][1]+u[i-1][0])/2.)*dy)*(0.5+alpha_y[i][0].w)));
		Ap_v[i][0] = Ae_v[i][0] + Aw_v[i][0] + An_v[i][0] + As_v[i][0];
		B_v[i][0] = 0.;
	}
	// East boundary volumes
	for(int j=1;j<(nv-2);j++)
	{
		As_v[nv-1][j] = (((mi*beta_y[nv-1][j].s*dx)/dy)+((rho*((v[nv-1][j-1]+v[nv-1][j])/2.)*dx)*(0.5+alpha_y[nv-1][j].s)));
		An_v[nv-1][j] = (((mi*beta_y[nv-1][j].n*dx)/dy)-((rho*((v[nv-1][j+1]+v[nv-1][j])/2.)*dx)*(0.5-alpha_y[nv-1][j].n)));
		Ae_v[nv-1][j] = ((2.*mi*beta_y[nv-1][j].e*dy)/dx);
		Aw_v[nv-1][j] = (((mi*beta_y[nv-1][j].w*dy)/dx)+((rho*((u[nv-2][j+1]+u[nv-2][j])/2.)*dy)*(0.5+alpha_y[nv-1][j].w)));
		Ap_v[nv-1][j] = Ae_v[nv-1][j] + Aw_v[nv-1][j] + An_v[nv-1][j] + As_v[nv-1][j];
		B_v[nv-1][j] = 0.;
	}
	// West boundary volumes
	for(int j=1;j<(nv-2);j++)
	{
		As_v[0][j] = (((mi*beta_y[0][j].s*dx)/dy)+((rho*((v[0][j-1]+v[0][j])/2.)*dx)*(0.5+alpha_y[0][j].s)));
		An_v[0][j] = (((mi*beta_y[0][j].n*dx)/dy)-((rho*((v[0][j+1]+v[0][j])/2.)*dx)*(0.5-alpha_y[0][j].n)));
		Ae_v[0][j] = (((mi*beta_y[0][j].e*dy)/dx)-((rho*((u[0][j+1]+u[0][j])/2.)*dy)*(0.5-alpha_y[0][j].e)));
		Aw_v[0][j] = ((2.*mi*beta_y[0][j].w*dy)/dx);
		Ap_v[0][j] = Ae_v[0][j] + Aw_v[0][j] + An_v[0][j] + As_v[0][j];
		B_v[0][j] = 0.;
	}
	// Core volummes
	for(int i=1;i<(nv-1);i++)
	{
		for(int j=1;j<(nv-2);j++)
		{
			As_v[i][j] = (((mi*beta_y[i][j].s*dx)/dy)+((rho*((v[i][j-1]+v[i][j])/2.)*dx)*(0.5+alpha_y[i][j].s)));
			An_v[i][j] = (((mi*beta_y[i][j].n*dx)/dy)-((rho*((v[i][j+1]+v[i][j])/2.)*dx)*(0.5-alpha_y[i][j].n)));
			Ae_v[i][j] = (((mi*beta_y[i][j].e*dy)/dx)-((rho*((u[i][j+1]+u[i][j])/2.)*dy)*(0.5-alpha_y[i][j].e)));
			Aw_v[i][j] = (((mi*beta_y[i][j].w*dy)/dx)+((rho*((u[i-1][j+1]+u[i-1][j])/2.)*dy)*(0.5+alpha_y[i][j].w)));
			Ap_v[i][j] = Ae_v[i][j] + Aw_v[i][j] + An_v[i][j] + As_v[i][j];
			B_v[i][j] = 0.;
		}
	}
}

void Calc_u_hat()
{
	// Left-up Corner
	u_hat[0][nv-1] = ((Ae_u[0][nv-1]*u[1][nv-1])+(As_u[0][nv-1]*u[0][nv-2])+(B_u[0][nv-1]))/Ap_u[0][nv-1];
	// Right-up corner
	u_hat[nv-2][nv-1] = ((Aw_u[nv-2][nv-1]*u[nv-3][nv-1])+(As_u[nv-2][nv-1]*u[nv-2][nv-2])+(B_u[nv-2][nv-1]))/Ap_u[nv-2][nv-1];
	// Left-bottom corner
	u_hat[0][0] = ((Ae_u[0][0]*u[1][0])+(An_u[0][0]*u[0][1])+(B_u[0][0]))/Ap_u[0][0];
	// Right-bottom corner
	u_hat[nv-2][0] = ((Aw_u[nv-2][0]*u[nv-3][0])+(An_u[nv-2][0]*u[nv-2][1])+(B_u[nv-2][0]))/Ap_u[nv-2][0];
	// East boundary
	for(int j=1;j<(nv-1);j++)
	{
		u_hat[nv-2][j] = ((Aw_u[nv-2][j]*u[nv-3][j])+(An_u[nv-2][j]*u[nv-2][j+1])+(As_u[nv-2][j]*u[nv-2][j-1])+(B_u[nv-2][j]))/Ap_u[nv-2][j];
	}
	// West boundary
	for(int j=1;j<(nv-1);j++)
	{
		u_hat[0][j] = ((Ae_u[0][j]*u[1][j])+(An_u[0][j]*u[0][j+1])+(As_u[0][j]*u[0][j-1])+(B_u[0][j]))/Ap_u[0][j];
	}
	// North boundary
	for(int i=1;i<(nv-2);i++)
	{
		u_hat[i][nv-1] = ((Ae_u[i][nv-1]*u[i+1][nv-1])+(Aw_u[i][nv-1]*u[i-1][nv-1])+(As_u[i][nv-1]*u[i][nv-2])+(B_u[i][nv-1]))/Ap_u[i][nv-1];
	}
	// South boundary
	for(int i=1;i<(nv-2);i++)
	{
		u_hat[i][0] = ((Ae_u[i][0]*u[i+1][0])+(Aw_u[i][0]*u[i-1][0])+(An_u[i][0]*u[i][1])+(B_u[i][0]))/Ap_u[i][0];
	}
	// Core volumes
	for(int i=1;i<(nv-2);i++)
	{
		for(int j=1;j<(nv-1);j++)
		{
			u_hat[i][j] = ((Ae_u[i][j]*u[i+1][j])+(Aw_u[i][j]*u[i-1][j])+(An_u[i][j]*u[i][j+1])+(As_u[i][j]*u[i][j-1])+(B_u[i][j]))/Ap_u[i][j];
		}
	}
}

void Calc_v_hat()
{
	// Left-up Corner
	v_hat[0][nv-2] = ((Ae_v[0][nv-2]*v[1][nv-2])+(As_v[0][nv-2]*v[0][nv-3])+(B_v[0][nv-2]))/Ap_v[0][nv-2];
	// Right-up corner
	v_hat[nv-1][nv-2] = ((Aw_v[nv-1][nv-2]*v[nv-2][nv-2])+(As_v[nv-1][nv-2]*v[nv-1][nv-3])+(B_v[nv-1][nv-2]))/Ap_v[nv-1][nv-2];
	// Left-bottom corner
	v_hat[0][0] = ((Ae_v[0][0]*v[1][0])+(An_v[0][0]*v[0][1])+(B_v[0][0]))/Ap_v[0][0];
	// Right-bottom corner
	v_hat[nv-1][0] = ((Aw_v[nv-1][0]*v[nv-2][0])+(An_v[nv-1][0]*v[nv-1][1])+(B_v[nv-1][0]))/Ap_v[nv-1][0];
	// East boundary
	for(int j=1;j<(nv-2);j++)
	{
		v_hat[nv-1][j] = ((Aw_v[nv-1][j]*v[nv-2][j])+(An_v[nv-1][j]*v[nv-1][j+1])+(As_v[nv-1][j]*v[nv-1][j-1])+(B_v[nv-1][j]))/Ap_v[nv-1][j];
	}
	// West boundary
	for(int j=1;j<(nv-2);j++)
	{
		v_hat[0][j] = ((Ae_v[0][j]*v[1][j])+(An_v[0][j]*v[0][j+1])+(As_v[0][j]*v[0][j-1])+(B_v[0][j]))/Ap_v[0][j];
	}
	// North boundary
	for(int i=1;i<(nv-1);i++)
	{
		v_hat[i][nv-2] = ((Ae_v[i][nv-2]*v[i+1][nv-2])+(Aw_v[i][nv-2]*v[i-1][nv-2])+(As_v[i][nv-2]*v[i][nv-3])+(B_v[i][nv-2]))/Ap_v[i][nv-2];
	}
	// South boundary
	for(int i=1;i<(nv-1);i++)
	{
		v_hat[i][0] = ((Ae_v[i][0]*v[i+1][0])+(Aw_v[i][0]*v[i-1][0])+(An_v[i][0]*v[i][1])+(B_v[i][0]))/Ap_v[i][0];
	}
	// Core volumes
	for(int i=1;i<(nv-1);i++)
	{
		for(int j=1;j<(nv-2);j++)
		{
			v_hat[i][j] = ((Ae_v[i][j]*v[i+1][j])+(Aw_v[i][j]*v[i-1][j])+(An_v[i][j]*v[i][j+1])+(As_v[i][j]*v[i][j-1])+(B_v[i][j]))/Ap_v[i][j];
		}
	}
}

void Calc_Coef_Pressao()
{
	// Left-up corner
	Ae_p[0][nv-1] = (rho*dy*dx*dy)/(dx*Ap_u[0][nv-1]);
	An_p[0][nv-1] = 0.;
	Aw_p[0][nv-1] = 0.;
	As_p[0][nv-1] = (rho*dx*dx*dy)/(dy*Ap_v[0][nv-2]);
	
	Ap_p[0][nv-1] = Ae_p[0][nv-1]+As_p[0][nv-1];
	B_p[0][nv-1] = -(rho*u_hat[0][nv-1]*dy)+(rho*v_hat[0][nv-2]*dx);
	// Right-up corner
	Ae_p[nv-1][nv-1] = 0.;
	An_p[nv-1][nv-1] = 0.;
	Aw_p[nv-1][nv-1] = (rho*dx*dx*dy)/(dy*Ap_u[nv-2][nv-1]);
	As_p[nv-1][nv-1] = (rho*dx*dx*dy)/(dy*Ap_v[nv-1][nv-2]);
	
	Ap_p[nv-1][nv-1] = Aw_p[nv-1][nv-1]+As_p[nv-1][nv-1];
	B_p[nv-1][nv-1] = (rho*u_hat[nv-2][nv-1]*dy)+(rho*v_hat[nv-1][nv-2]*dx);
	// Left-bottom corner
	Ae_p[0][0] = (rho*dy*dx*dy)/(dx*Ap_u[0][0]);
	An_p[0][0] = (rho*dx*dx*dy)/(dy*Ap_v[0][0]);
	Aw_p[0][0] = 0.;
	As_p[0][0] = 0.;
			
	Ap_p[0][0] = Ae_p[0][0]+An_p[0][0];
	B_p[0][0] = -(rho*u_hat[0][0]*dy)-(rho*v_hat[0][0]*dy);
	// Right-bottom corner
	Ae_p[nv-1][0] = 0.;
	An_p[nv-1][0] = (rho*dx*dx*dy)/(dy*Ap_v[nv-1][0]);
	Aw_p[nv-1][0] = (rho*dx*dx*dy)/(dy*Ap_u[nv-2][0]);
	As_p[nv-1][0] = 0.;
		
	Ap_p[nv-1][0] = Aw_p[nv-1][0]+An_p[nv-1][0];
	B_p[nv-1][0] = (rho*u_hat[nv-2][0]*dy)-(rho*v_hat[nv-1][0]*dy);
	// East boundary
	for(int j=1;j<(nv-1);j++)
	{
		Ae_p[nv-1][j] = 0.;
		An_p[nv-1][j] = (rho*dx*dx*dy)/(dy*Ap_v[nv-1][j]);
		Aw_p[nv-1][j] = (rho*dx*dx*dy)/(dy*Ap_u[nv-2][j]);
		As_p[nv-1][j] = (rho*dx*dx*dy)/(dy*Ap_v[nv-1][j-1]);
			
		Ap_p[nv-1][j] = Aw_p[nv-1][j]+An_p[nv-1][j]+As_p[nv-1][j];
		B_p[nv-1][j] = (rho*u_hat[nv-2][j]*dy)+(rho*v_hat[nv-1][j-1]*dx)-(rho*v_hat[nv-1][j]*dy);
	}
	// West boundary
	for(int j=1;j<(nv-1);j++)
	{
		Ae_p[0][j] = (rho*dy*dx*dy)/(dx*Ap_u[0][j]);
		An_p[0][j] = (rho*dx*dx*dy)/(dy*Ap_v[0][j]);
		Aw_p[0][j] = 0.;
		As_p[0][j] = (rho*dx*dx*dy)/(dy*Ap_v[0][j-1]);
		
		Ap_p[0][j] = Ae_p[0][j]+An_p[0][j]+As_p[0][j];
		B_p[0][j] = -(rho*u_hat[0][j]*dy)+(rho*v_hat[0][j-1]*dx)-(rho*v_hat[0][j]*dx);
	}
	// North boundary
	for(int i=1;i<(nv-1);i++)
	{
		Ae_p[i][nv-1] = (rho*dy*dx*dy)/(dx*Ap_u[i][nv-1]);
		An_p[i][nv-1] = 0.;
		Aw_p[i][nv-1] = (rho*dx*dx*dy)/(dy*Ap_u[i-1][nv-1]);
		As_p[i][nv-1] = (rho*dx*dx*dy)/(dy*Ap_v[i][nv-2]);
		
		Ap_p[i][nv-1] = Ae_p[i][nv-1]+Aw_p[i][nv-1]+As_p[i][nv-1];
		B_p[i][nv-1] = (rho*u_hat[i-1][nv-1]*dy)-(rho*u_hat[i][nv-1]*dy)+(rho*v_hat[i][nv-2]*dx);
	}
	// South boundary
	for(int i=1;i<(nv-1);i++)
	{
		Ae_p[i][0] = (rho*dy*dx*dy)/(dx*Ap_u[i][0]);
		An_p[i][0] = (rho*dx*dx*dy)/(dy*Ap_v[i][0]);
		Aw_p[i][0] = (rho*dx*dx*dy)/(dy*Ap_u[i-1][0]);
		As_p[i][0] = 0.;
		
		Ap_p[i][0] = Ae_p[i][0]+Aw_p[i][0]+An_p[i][0];
		B_p[i][0] = (rho*u_hat[i-1][0]*dy)-(rho*u_hat[i][0]*dy)-(rho*v_hat[i][0]*dy);
	}
	// Core volumes
	for(int i=1;i<(nv-1);i++)
	{
		for(int j=1;j<(nv-1);j++)
		{
			Ae_p[i][j] = (rho*dy*dx*dy)/(dx*Ap_u[i][j]);
			An_p[i][j] = (rho*dx*dx*dy)/(dy*Ap_v[i][j]);
			Aw_p[i][j] = (rho*dx*dx*dy)/(dy*Ap_u[i-1][j]);
			As_p[i][j] = (rho*dx*dx*dy)/(dy*Ap_v[i][j-1]);
			
			Ap_p[i][j] = Ae_p[i][j]+Aw_p[i][j]+An_p[i][j]+As_p[i][j];
			B_p[i][j] = (rho*u_hat[i-1][j]*dy)-(rho*u_hat[i][j]*dy)+(rho*v_hat[i][j-1]*dx)-(rho*v_hat[i][j]*dy);
		}
	}
}

void correct_u_v()
{
	// Correct x-velocity u
	for(int i=0;i<(nv-1);i++)
	{
		for(int j=0;j<nv;j++)
		{
			uOLD[i][j] = u[i][j];
			u[i][j]    = u_hat[i][j] - ((Pn[i+1][j]-Pn[i][j])*dx*dy)/(Ap_u[i][j]*dx);
		}
	}
	// Correct y-velocity v
	for(int i=0;i<nv;i++)
	{
		for(int j=0;j<(nv-1);j++)
		{
			vOLD[i][j] = v[i][j];
			v[i][j]    = v_hat[i][j] - ((Pn[i][j+1]-Pn[i][j])*dy*dx)/(Ap_v[i][j]*dy);
		}
	}
}

int main()
{
	std:: string filename_input = "./inCav.txt";
	std:: string filename_results = "./outCav.txt";
	
	std::tie(L, H, nv, rho, U, mi) = read_input_data(filename_input);

	// Calculate Delta X e Delta Y
	dx = (double)L/nv; 
	dy = (double)H/nv;

	int n_x = nv;
	int n_y = nv;

	allocate_vector_2d(u    , n_x-1, n_y);
	allocate_vector_2d(uOLD , n_x-1, n_y);
	allocate_vector_2d(u_hat, n_x-1, n_y);
	allocate_vector_2d(Ap_u , n_x-1, n_y);
	allocate_vector_2d(Aw_u , n_x-1, n_y);
	allocate_vector_2d(Ae_u , n_x-1, n_y);
	allocate_vector_2d(An_u , n_x-1, n_y);
	allocate_vector_2d(As_u , n_x-1, n_y);
	allocate_vector_2d(B_u  , n_x-1, n_y);

	allocate_cvbound_2d(alpha_x, n_x-1, n_y);
	allocate_cvbound_2d(beta_x , n_x-1, n_y);

	allocate_vector_2d(v    , n_x, n_y-1);
	allocate_vector_2d(vOLD , n_x, n_y-1);
	allocate_vector_2d(v_hat, n_x, n_y-1);
	allocate_vector_2d(Ap_v , n_x, n_y-1);
	allocate_vector_2d(Aw_v , n_x, n_y-1);
	allocate_vector_2d(Ae_v , n_x, n_y-1);
	allocate_vector_2d(An_v , n_x, n_y-1);
	allocate_vector_2d(As_v , n_x, n_y-1);
	allocate_vector_2d(B_v  , n_x, n_y-1);

	allocate_cvbound_2d(alpha_y, n_x, n_y-1);
	allocate_cvbound_2d(beta_y , n_x, n_y-1);

	allocate_vector_2d(P   , n_x, n_y);
	allocate_vector_2d(Pn  , n_x, n_y);
	allocate_vector_2d(Ap_p, n_x, n_y);
	allocate_vector_2d(Aw_p, n_x, n_y);
	allocate_vector_2d(Ae_p, n_x, n_y);
	allocate_vector_2d(An_p, n_x, n_y);
	allocate_vector_2d(As_p, n_x, n_y);
	allocate_vector_2d(B_p , n_x, n_y);

	int IT = 0;
	int TESTE =0;
	while ( TESTE == 0 )
	{
		IT++;
		cout << IT << " ";
		Calc_WUDS_Coef_X(rho, mi, nv, dx, dy, u, v, alpha_x, beta_x);
		Calc_Coef_NS_X(); 
		Calc_WUDS_Coef_Y(rho, mi, nv, dx, dy, u, v, alpha_y, beta_y);
		Calc_Coef_NS_Y();
		Calc_u_hat();
		Calc_v_hat();
		Calc_Coef_Pressao();
		SOR_structured(
			Ap_p, Aw_p, Ae_p, An_p, As_p,
			P, Pn, B_p, 
			nv, 50, 1.6
		);	
		correct_u_v();

		double error_u = calculate_vec_diff_L2_norm(u, uOLD, n_x-1, n_y);
		double error_v = calculate_vec_diff_L2_norm(v, vOLD, n_x, n_y-1);
		cout <<"error -u:" << setw(7) << setprecision(5) << error_u
			 << " -v:" << setw(7) << setprecision(5) << error_v << endl;
		if((IT % 500) == 0)
		{
			cout << endl << "......Saving Partial Solution....." << endl;
			save_results(filename_results, u, v, Pn, nv, dx, dy, U);
		}
		if( ( (error_u < (0.0001)) and (error_v < (0.0001) )) or (IT == 100000) )
		{
			TESTE = 1;
		}
	}
	IT++;
	cout << IT << " ";
	Calc_WUDS_Coef_X(rho, mi, nv, dx, dy, u, v, alpha_x, beta_x);
	Calc_Coef_NS_X(); 
	Calc_WUDS_Coef_Y(rho, mi, nv, dx, dy, u, v, alpha_y, beta_y);
	Calc_Coef_NS_Y();
	Calc_u_hat();
	Calc_v_hat();
	Calc_Coef_Pressao();
	SOR_structured(
		Ap_p, Aw_p, Ae_p, An_p, As_p,
		P, Pn, B_p,
		nv, 1000, 1.6
	);
	correct_u_v();

	double error_u = calculate_vec_diff_L2_norm(u, uOLD, n_x-1, n_y);
	double error_v = calculate_vec_diff_L2_norm(v, vOLD, n_x, n_y-1);
	cout << "error -u:" << setw(7) << setprecision(5) << error_u
		 << " -v:" << setw(7) << setprecision(5) << error_v << endl;
	save_results(filename_results, u, v, Pn, nv, dx, dy, U);

	deallocate_vector_2d(u    , n_x-1, n_y);
	deallocate_vector_2d(uOLD , n_x-1, n_y);
	deallocate_vector_2d(u_hat, n_x-1, n_y);
	deallocate_vector_2d(Ap_u , n_x-1, n_y);
	deallocate_vector_2d(Aw_u , n_x-1, n_y);
	deallocate_vector_2d(Ae_u , n_x-1, n_y);
	deallocate_vector_2d(An_u , n_x-1, n_y);
	deallocate_vector_2d(As_u , n_x-1, n_y);
	deallocate_vector_2d(B_u  , n_x-1, n_y);

	deallocate_cvbound_2d(alpha_x, n_x-1, n_y);
	deallocate_cvbound_2d(beta_x , n_x-1, n_y);

	deallocate_vector_2d(v    , n_x, n_y-1);
	deallocate_vector_2d(vOLD , n_x, n_y-1);
	deallocate_vector_2d(v_hat, n_x, n_y-1);
	deallocate_vector_2d(Ap_v , n_x, n_y-1);
	deallocate_vector_2d(Aw_v , n_x, n_y-1);
	deallocate_vector_2d(Ae_v , n_x, n_y-1);
	deallocate_vector_2d(An_v , n_x, n_y-1);
	deallocate_vector_2d(As_v , n_x, n_y-1);
	deallocate_vector_2d(B_v  , n_x, n_y-1);

	deallocate_cvbound_2d(alpha_y, n_x, n_y-1);
	deallocate_cvbound_2d(beta_y , n_x, n_y-1);

	deallocate_vector_2d(P   , n_x, n_y);
	deallocate_vector_2d(Pn  , n_x, n_y);
	deallocate_vector_2d(Ap_p, n_x, n_y);
	deallocate_vector_2d(Aw_p, n_x, n_y);
	deallocate_vector_2d(Ae_p, n_x, n_y);
	deallocate_vector_2d(An_p, n_x, n_y);
	deallocate_vector_2d(As_p, n_x, n_y);
	deallocate_vector_2d(B_p , n_x, n_y);
}
