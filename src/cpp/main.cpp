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
#include "Legacy/PRIMEcorrection.hpp"
#include "Legacy/MassEquation.hpp"

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

	double tol = 1.e-4;
	int MAX_IT = 100000;
	int IT = 1;
	while (true)
	{
		cout << IT << " ";
		calculate_WUDS_coefficients_X(rho, mi, nv, dx, dy, u, v, alpha_x, beta_x);
		Calc_Coef_NS_X(); 
		calculate_WUDS_coefficients_Y(rho, mi, nv, dx, dy, u, v, alpha_y, beta_y);
		Calc_Coef_NS_Y();
		calculate_u_hat(nv, Ap_u, Ae_u, Aw_u, As_u, An_u, B_u, u, u_hat);
		calculate_v_hat(nv, Ap_v, Ae_v, Aw_v, As_v, An_v, B_v, v, v_hat);
		calculate_pressure_coefficients(nv, dx, dy, rho, u_hat, v_hat, Ap_u, Ap_v,
										Ap_p, Ae_p, Aw_p, As_p, An_p, B_p);
		SOR_structured(
			Ap_p, Aw_p, Ae_p, An_p, As_p,
			P, Pn, B_p, 
			nv, 50, 1.6
		);	
		correct_u_v(nv, dx, dy, Pn, Ap_u, uOLD,	u_hat, u, Ap_v, vOLD, v_hat, v);

		double error_u = calculate_vec_diff_L2_norm(u, uOLD, n_x-1, n_y);
		double error_v = calculate_vec_diff_L2_norm(v, vOLD, n_x, n_y-1);
		cout <<"error -u:" << setw(7) << setprecision(5) << error_u
			 << " -v:" << setw(7) << setprecision(5) << error_v << endl;
		if((IT % 500) == 0)
		{
			cout << endl << "......Saving Partial Solution....." << endl;
			save_results(filename_results, u, v, Pn, nv, dx, dy, U);
		}
		IT++;
		if( ( (error_u <= (tol)) and (error_v <= (tol) )) or (IT >= MAX_IT + 1) )
		{
			break;
		}
	}
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
