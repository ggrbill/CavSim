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
#include "Legacy/MomentumEquation.hpp"

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

	int saving_interval = 500;
	double tol = 1.e-4;
	int MAX_IT = 100000;
	int IT = 1;
	while (true)
	{
		cout << IT << " ";
		calculate_WUDS_coefficients_X(rho, mi, nv, dx, dy, u, v, alpha_x, beta_x);
		calculate_WUDS_coefficients_Y(rho, mi, nv, dx, dy, u, v, alpha_y, beta_y);
		
		calculate_velocity_coeficients_X(U, rho, mi, nv, dx, dy, u, v, alpha_x, beta_x,
										 Ap_u, Ae_u, Aw_u, As_u, An_u, B_u); 
		calculate_velocity_coeficients_Y(rho, mi, nv, dx, dy, u, v, alpha_y, beta_y,
										 Ap_v, Ae_v, Aw_v, As_v, An_v, B_v); 
		
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
		
		if((IT % saving_interval) == 0) {
			cout << endl << "......Saving Partial Solution....." << endl;
			save_results(filename_results, u, v, Pn, nv, dx, dy, U);
		}
		
		IT++;
		bool velocity_error_condition = (error_u <= (tol)) and (error_v <= (tol));
		bool max_iteration_condition = IT >= MAX_IT + 1;
		bool stop_condition =  velocity_error_condition or max_iteration_condition;
		if (stop_condition) {
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
