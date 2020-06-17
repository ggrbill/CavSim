#include "MomentumEquation.hpp"

void calculate_velocity_coeficients_X(
    double U,
    double rho,
    double mi,
    int nv,
    double dx,
    double dy,
    double **u,
    double **v,
    CVBoundaries **alpha_x,
    CVBoundaries **beta_x,
    double **Ap_u,
    double **Ae_u,
    double **Aw_u,
    double **As_u,
    double **An_u,
    double **B_u
) 
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


void calculate_velocity_coeficients_Y(
    double rho,
    double mi,
    int nv,
    double dx,
    double dy,
    double **u,
    double **v,
    CVBoundaries **alpha_y,
    CVBoundaries **beta_y,
    double **Ap_v,
    double **Ae_v,
    double **Aw_v,
    double **As_v,
    double **An_v,
    double **B_v
)
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