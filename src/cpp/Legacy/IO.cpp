#include <fstream>

#include "IO.hpp"

std::tuple<double, double, int, double, double, double> read_input_data(std::string filename)
{
	double Length; 
	double Height;
	int    n_rc; // Number of rows and columns
	double density;
	double vel_top_bound;
	double viscosity;

	std::fstream fin(filename);
	fin >>  Length;
	fin >>  Height;
	fin >>  n_rc;
	fin >>  density;
	fin >>  vel_top_bound; 
	fin >>  viscosity;
	fin.close();

	return std::make_tuple(
		Length,
		Height,
		n_rc,
		density,
		vel_top_bound,
		viscosity);
}

void save_results(
	std::string filename, 
	std::shared_ptr<CavitySetup> cav_setup, 
	double** u,
	double** v,
	double** Pn)
{	
	// For while it is assuming that n_x is equal to n_y
	int nv = cav_setup->n_x; 
	std::ofstream fout;
	fout << std::scientific;
	fout.open(filename);
	fout << "TITLE = \"OutCav\" " << std::endl;
	fout << "VARIABLES = \"X\", \"Y\", \"U\", \"V\", \"P\" " << std::endl;
	fout <<"ZONE T=\"" << "OUTCAV" << "\", N=" << (nv+1)*(nv+1) << ", E=" << nv*nv  <<  
  	",  DATAPACKING=BLOCK, ZONETYPE=FEQUADRILATERAL, VARLOCATION=([3-5]=CELLCENTERED), SOLUTIONTIME=" << "1" << std::endl;

	// x coordinates
	int k = 0;
	for (int j = 1; j < (nv + 2); j++) {
		fout << 0. << "  ";
		k++;
		if (k == 9) {
			fout << std::endl;
			k = 0;
		}
		for (int i = 1; i < (nv + 1); i++) {
			fout << i * cav_setup->dx << "  ";
			k++;
			if (k == 9) {
				fout << std::endl;
				k = 0;
			}
		}
	}
	fout << std::endl;
	fout << std::endl;

	// y coordinates
	k = 0;
	for (int j = 1; j < (nv + 2); j++) {
		fout << 0. << "  ";
		k++;
		if (k == 9) {
			fout << std::endl;
			k = 0;
		}
	}
	for (int j = 1; j < (nv + 1); j++) {
		for (int i = 1; i < (nv + 2); i++) {
			fout << j * cav_setup->dy << "  ";
			k++;
			if (k == 9) {
				fout << std::endl;
				k = 0;
			}
		}
	}
	fout << std::endl;
	fout << std::endl;

	// x-velocity U
	k = 0;
	for (int j = 0; j < (nv); j++) {
		for (int i = 0; i < (nv - 1); i++) {
			fout << u[i][j] / cav_setup->U << "  ";
			k++;
			if (k == 9) {
				fout << std::endl;
				k = 0;
			}
		}
		fout << "0." << "  ";
		k++;
		if (k == 9) {
			fout << std::endl;
			k = 0;
		}
	}
	fout << std::endl;
	fout << std::endl;

	// y-velocity V
	k = 0;
	for (int j = 0; j < (nv - 1); j++) {
		for (int i = 0; i < (nv); i++) {
			fout << v[i][j] / cav_setup->U << "  ";
			k++;
			if (k == 9) {
				fout << std::endl;
				k = 0;
			}
		}
	}
	for (int i = 0; i < (nv); i++) {
		fout << "0." << "  ";
		k++;
		if (k == 9) {
			fout << std::endl;
			k = 0;
		}
	}
	fout << std::endl;
	fout << std::endl;

	// Pressure
	k = 0;
	for (int j = 0; j < (nv); j++) {
		for (int i = 0; i < (nv); i++) {
			fout << Pn[i][j] << "  ";
			k++;
			if (k == 9) {
				fout << std::endl;
				k = 0;
			}
		}
	}
	fout << std::endl;
	fout << std::endl;

	// Cells (Control Volumes)
	for (int i = 1; i < (nv + 1); i++) {
		fout << i + (nv + 1) << "\t" << i + (nv + 2) << "\t" << i + 1 << "\t" << i << std::endl;
	}
	for (int j = 1; j < (nv); j++) {
		for (int i = 1; i < (nv + 1); i++) {
			fout << i + (j * (nv + 1)) + (nv + 1) << "\t" << i + (j * (nv + 1)) + (nv + 2) << "\t" << i + (j * (nv + 1)) + 1 << "\t" << i + (j * (nv + 1)) << std::endl;
		}
	}

	fout << std::endl;
	fout.close();
}