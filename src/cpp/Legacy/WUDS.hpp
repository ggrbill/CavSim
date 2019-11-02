#ifndef WUDS_HPP
#define WUDS_HPP

#include "Structures.hpp"

/*!
    Calculates Alpha and Beta Coefficients of WUDS scheme.
    
    Those coefficients are used to calculate the Advective and 
    Diffusive terms of Discrete Navier-Stokes Equation in respect 
    to the X-axis.  
*/
void Calc_WUDS_Coef_X(
    double rho,
    double mi,
    int nv,
    double dx,
    double dy,
    double** u,
    double** v,
    CVBoundaries** alpha_x,
    CVBoundaries** beta_x);

/*!
    Calculates Alpha and Beta Coefficients of WUDS scheme.
    
    Those coefficients are used to calculate the Advective and 
    Diffusive terms of Discrete Navier-Stokes Equation in respect 
    to the Y-axis.  
*/
void Calc_WUDS_Coef_Y(
    double rho,
    double mi,
    int nv,
    double dx,
    double dy,
    double** u,
    double** v,
    CVBoundaries** alpha_y,
    CVBoundaries** beta_y);


#endif