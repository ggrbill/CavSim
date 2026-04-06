#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <math.h>

#include "classic/solver.hpp"
#include "classic/IO.hpp"
#include "classic/numeric.hpp"
#include "classic/Structures.hpp"
#include "classic/WUDS.hpp"
#include "classic/PRIMEcorrection.hpp"
#include "classic/MassEquation.hpp"
#include "classic/MomentumEquation.hpp"

using namespace std;

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
    std:: string filename_results_csv = "./outCav.csv";
	
	std::tie(L, H, nv, rho, U, mi) = read_input_data(filename_input);

	// Calculate Delta X e Delta Y
	dx = (double)L/nv; 
	dy = (double)H/nv;

	int n_x = nv;
	int n_y = nv;

	PrimeCoefficients p_coeffs(n_x, n_y);
	CavSimResult r(n_x, n_y);

	int saving_interval = 500;
	double tol = 1.e-4;
	int MAX_IT = 100000;
	int IT = 1;
	while (true)
	{
		cout << IT << " ";
		calculate_WUDS_coefficients_X(rho, mi, nv, dx, dy, r.u, r.v, p_coeffs.alpha_x, p_coeffs.beta_x);
		calculate_WUDS_coefficients_Y(rho, mi, nv, dx, dy, r.u, r.v, p_coeffs.alpha_y, p_coeffs.beta_y);
		
		calculate_velocity_coeficients_X(U, rho, mi, nv, dx, dy, r.u, r.v, p_coeffs.alpha_x, p_coeffs.beta_x,
										 p_coeffs.Ap_u, p_coeffs.Ae_u, p_coeffs.Aw_u, p_coeffs.As_u, p_coeffs.An_u, p_coeffs.B_u); 
		calculate_velocity_coeficients_Y(rho, mi, nv, dx, dy, r.u, r.v, p_coeffs.alpha_y, p_coeffs.beta_y,
										 p_coeffs.Ap_v, p_coeffs.Ae_v, p_coeffs.Aw_v, p_coeffs.As_v, p_coeffs.An_v, p_coeffs.B_v); 
		
		calculate_u_hat(nv, p_coeffs.Ap_u, p_coeffs.Ae_u, p_coeffs.Aw_u, p_coeffs.As_u, p_coeffs.An_u, p_coeffs.B_u, r.u, r.u_hat);
		calculate_v_hat(nv, p_coeffs.Ap_v, p_coeffs.Ae_v, p_coeffs.Aw_v, p_coeffs.As_v, p_coeffs.An_v, p_coeffs.B_v, r.v, r.v_hat);
		
		calculate_pressure_coefficients(nv, dx, dy, rho, r.u_hat, r.v_hat, p_coeffs.Ap_u, p_coeffs.Ap_v,
										p_coeffs.Ap_p, p_coeffs.Ae_p, p_coeffs.Aw_p, p_coeffs.As_p, p_coeffs.An_p, p_coeffs.B_p);
		SOR_structured(
			p_coeffs.Ap_p, p_coeffs.Aw_p, p_coeffs.Ae_p, p_coeffs.An_p, p_coeffs.As_p,
			r.P, r.Pn, p_coeffs.B_p, 
			nv, 50, 1.6
		);	
		correct_u_v(nv, dx, dy, r.Pn, p_coeffs.Ap_u, r.u_old, r.u_hat, r.u, p_coeffs.Ap_v, r.v_old, r.v_hat, r.v);

		double error_u = calculate_vec_diff_L2_norm(r.u, r.u_old, n_x-1, n_y);
		double error_v = calculate_vec_diff_L2_norm(r.v, r.v_old, n_x, n_y-1);
		cout <<"error -u:" << setw(7) << setprecision(5) << error_u
			 << " -v:" << setw(7) << setprecision(5) << error_v << endl;
		
		if((IT % saving_interval) == 0) {
			cout << endl << "......Saving Partial Solution....." << endl;
			save_results_tecplot(filename_results, r.u, r.v, r.Pn, nv, dx, dy, U);
		}
		
		IT++;
		bool velocity_error_condition = (error_u <= (tol)) and (error_v <= (tol));
		bool max_iteration_condition = IT >= MAX_IT + 1;
		bool stop_condition =  velocity_error_condition or max_iteration_condition;
		if (stop_condition) {
			break;
		}
	}
	save_results_tecplot(filename_results, r.u, r.v, r.Pn, nv, dx, dy, U);
	save_results_csv(filename_results_csv, r.u, r.v, r.Pn, nv, dx, dy, U);
}
