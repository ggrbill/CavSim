#ifndef IO_HPP
#define IO_HPP

#include <memory>
#include <tuple>
#include "Structures.hpp"

/*!
    Input data reader.
*/
std::tuple<double, double, int, double, double, double> read_input_data(std::string filename);

/*!
    Tecplot format Saver of results.
*/
void save_results_tecplot(
    std::string filename, 
	DoubleArray2D& u,
	DoubleArray2D& v,
	DoubleArray2D& Pn,
	double nv,
	double dx,
	double dy,
	double U);

/*!
    CSV format Saver of results.
*/
void save_results_csv(
    std::string filename, 
	DoubleArray2D& u,
	DoubleArray2D& v,
	DoubleArray2D& Pn,
	int nv,
	double dx,
	double dy,
	double U);

#endif
