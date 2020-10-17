#include "PRIMEcorrection.hpp"

void calculate_u_hat(
    int nv,
    double** Ap_u,
    double** Ae_u,
    double** Aw_u,
    double** As_u,
    double** An_u,
    double** B_u,
    double** u,
    double** u_hat
)
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


void calculate_v_hat(
	int nv,
    double** Ap_v,
    double** Ae_v,
    double** Aw_v,
    double** As_v,
    double** An_v,
    double** B_v,
    double** v,
    double** v_hat
)
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


void correct_u_v(
	int nv,
	double dx,
	double dy,
	double **Pn,
	double **Ap_u,
	double **uOLD,
	double **u_hat,
	double **u,
	double **Ap_v,
	double **vOLD,
	double **v_hat,
	double **v
)
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
