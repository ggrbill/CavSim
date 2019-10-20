#ifndef IO_HPP
#define IO_HPP

#include <memory>

#include "Cavity.hpp"

/*!
    Input data reader.
*/
std::shared_ptr<CavitySetup> read_input_data(std::string filename);

/*!
    Tecplot format Saver of results.
*/
void save_results(
    std::string filename, 
	std::shared_ptr<CavitySetup> cav_setup, 
	double** u,
	double** v,
	double** Pn);

#endif