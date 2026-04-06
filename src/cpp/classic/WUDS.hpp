#ifndef WUDS_HPP
#define WUDS_HPP

#include "Structures.hpp"

/*!
    Calculates Alpha and Beta Coefficients of WUDS scheme.
    
    Those coefficients are used to calculate the Advective and 
    Diffusive terms of Discrete Navier-Stokes Equation in respect 
    to the X-axis.  
*/
void calculate_WUDS_coefficients_X(
    double rho,
    double mi,
    int nv,
    double dx,
    double dy,
    DoubleArray2D& u,
    DoubleArray2D& v,
    FacesArray2D& alpha_x,
    FacesArray2D& beta_x);

/*!
    Calculates Alpha and Beta Coefficients of WUDS scheme.
    
    Those coefficients are used to calculate the Advective and 
    Diffusive terms of Discrete Navier-Stokes Equation in respect 
    to the Y-axis.  
*/
void calculate_WUDS_coefficients_Y(
    double rho,
    double mi,
    int nv,
    double dx,
    double dy,
    DoubleArray2D& u,
    DoubleArray2D& v,
    FacesArray2D& alpha_y,
    FacesArray2D& beta_y);


#endif
