#ifndef IO_HPP
#define IO_HPP

#include <memory>
#include <tuple>

/*!
    Input data reader.
*/
std::tuple<double, double, int, double, double, double> read_input_data(std::string filename);

/*!
    Tecplot format Saver of results.
*/
void save_results(
    std::string filename, 
	double** u,
	double** v,
	double** Pn,
	double nv,
	double dx,
	double dy,
	double U);

#endif