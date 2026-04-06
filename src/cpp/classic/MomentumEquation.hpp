#ifndef MOMENTUM_EQUATION_HPP
#define MOMENTUM_EQUATION_HPP

#include "Structures.hpp"

/*!
    Calculates the velocities coefficients for direction X.
    
    Those coefficients are related to the momentum conservation equation
    (Navier-Stokes Equation) and it is used to compute the U-velocity 
    (component at X-direction).
*/
void calculate_velocity_coeficients_X(
    double U,
    double rho,
    double mi,
    int nv,
    double dx,
    double dy,
    DoubleArray2D&u,
    DoubleArray2D&v,
    FacesArray2D&alpha_x,
    FacesArray2D&beta_x,
    DoubleArray2D&Ap_u,
    DoubleArray2D&Ae_u,
    DoubleArray2D&Aw_u,
    DoubleArray2D&As_u,
    DoubleArray2D&An_u,
    DoubleArray2D&B_u
);

/*!
    Calculates the velocities coefficients for direction Y.
    
    Those coefficients are related to the momentum conservation equation
    (Navier-Stokes Equation) and it is used to compute the V-velocity 
    (component at Y-direction).
*/
void calculate_velocity_coeficients_Y(
    double rho,
    double mi,
    int nv,
    double dx,
    double dy,
    DoubleArray2D&u,
    DoubleArray2D&v,
    FacesArray2D&alpha_x,
    FacesArray2D&beta_x,
    DoubleArray2D&Ap_u,
    DoubleArray2D&Ae_u,
    DoubleArray2D&Aw_u,
    DoubleArray2D&As_u,
    DoubleArray2D&An_u,
    DoubleArray2D&B_u
);

#endif