#include <fstream>

#include "IO.hpp"

std::shared_ptr<CavitySetup> read_input_data(std::string filename) // Read input data
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

	return std::make_shared<CavitySetup>(
		Length, 
		Height, 
		n_rc, 
		n_rc, 
		density, 
		viscosity,
		vel_top_bound
	);
}