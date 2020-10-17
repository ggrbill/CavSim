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
	
	std::tie(L, H, nv, rho, U, mi) = read_input_data(filename_input);

	// Calculate Delta X e Delta Y
	dx = (double)L/nv; 
	dy = (double)H/nv;

	int n_x = nv;
	int n_y = nv;

	CavSimAux a(n_x, n_y);
	CavSimData d(n_x, n_y);
	CavSimResult r(n_x, n_y);

	int saving_interval = 500;
	double tol = 1.e-4;
	int MAX_IT = 100000;
	int IT = 1;
	while (true)
	{
		cout << IT << " ";
		calculate_WUDS_coefficients_X(rho, mi, nv, dx, dy, r.u, r.v, a.alpha_x, a.beta_x);
		calculate_WUDS_coefficients_Y(rho, mi, nv, dx, dy, r.u, r.v, a.alpha_y, a.beta_y);
		
		calculate_velocity_coeficients_X(U, rho, mi, nv, dx, dy, r.u, r.v, a.alpha_x, a.beta_x,
										 d.Ap_u, d.Ae_u, d.Aw_u, d.As_u, d.An_u, d.B_u); 
		calculate_velocity_coeficients_Y(rho, mi, nv, dx, dy, r.u, r.v, a.alpha_y, a.beta_y,
										 d.Ap_v, d.Ae_v, d.Aw_v, d.As_v, d.An_v, d.B_v); 
		
		calculate_u_hat(nv, d.Ap_u, d.Ae_u, d.Aw_u, d.As_u, d.An_u, d.B_u, r.u, r.u_hat);
		calculate_v_hat(nv, d.Ap_v, d.Ae_v, d.Aw_v, d.As_v, d.An_v, d.B_v, r.v, r.v_hat);
		
		calculate_pressure_coefficients(nv, dx, dy, rho, r.u_hat, r.v_hat, d.Ap_u, d.Ap_v,
										d.Ap_p, d.Ae_p, d.Aw_p, d.As_p, d.An_p, d.B_p);
		SOR_structured(
			d.Ap_p, d.Aw_p, d.Ae_p, d.An_p, d.As_p,
			r.P, r.Pn, d.B_p, 
			nv, 50, 1.6
		);	
		correct_u_v(nv, dx, dy, r.Pn, d.Ap_u, r.uOLD, r.u_hat, r.u, d.Ap_v, r.vOLD, r.v_hat, r.v);

		double error_u = calculate_vec_diff_L2_norm(r.u, r.uOLD, n_x-1, n_y);
		double error_v = calculate_vec_diff_L2_norm(r.v, r.vOLD, n_x, n_y-1);
		cout <<"error -u:" << setw(7) << setprecision(5) << error_u
			 << " -v:" << setw(7) << setprecision(5) << error_v << endl;
		
		if((IT % saving_interval) == 0) {
			cout << endl << "......Saving Partial Solution....." << endl;
			save_results(filename_results, r.u, r.v, r.Pn, nv, dx, dy, U);
		}
		
		IT++;
		bool velocity_error_condition = (error_u <= (tol)) and (error_v <= (tol));
		bool max_iteration_condition = IT >= MAX_IT + 1;
		bool stop_condition =  velocity_error_condition or max_iteration_condition;
		if (stop_condition) {
			break;
		}
	}
	save_results(filename_results, r.u, r.v, r.Pn, nv, dx, dy, U);
}
