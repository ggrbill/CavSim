#ifndef MASS_EQUATION_HPP
#define MASS_EQUATION_HPP

#include "Structures.hpp"

/*!
    Calculates the pressure coefficients.
    
    Those coefficients are related to the mass conservation equation and 
    generates a linear system which the pressure is the main variable.
*/
void calculate_pressure_coefficients(
    int nv,
    double dx,
    double dy,
    double rho,
    DoubleArray2D&u_hat,
    DoubleArray2D&v_hat,
    DoubleArray2D&Ap_u,
    DoubleArray2D&Ap_v,
    DoubleArray2D&Ap_p,
    DoubleArray2D&Ae_p,
    DoubleArray2D&Aw_p,
    DoubleArray2D&As_p,
    DoubleArray2D&An_p,
    DoubleArray2D&B_p
);

#endif
