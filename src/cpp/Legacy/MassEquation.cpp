#include "MassEquation.hpp"

void calculate_pressure_coefficients(
    int nv,
    double dx,
    double dy,
    double rho,
    double **u_hat,
    double **v_hat,
    double **Ap_u,
    double **Ap_v,
    double **Ap_p,
    double **Ae_p,
    double **Aw_p,
    double **As_p,
    double **An_p,
    double **B_p
)
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
