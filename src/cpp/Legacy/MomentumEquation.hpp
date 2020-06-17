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
    double **u,
    double **v,
    CVBoundaries **alpha_x,
    CVBoundaries **beta_x,
    double **Ap_u,
    double **Ae_u,
    double **Aw_u,
    double **As_u,
    double **An_u,
    double **B_u
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
    double **u,
    double **v,
    CVBoundaries **alpha_x,
    CVBoundaries **beta_x,
    double **Ap_u,
    double **Ae_u,
    double **Aw_u,
    double **As_u,
    double **An_u,
    double **B_u
);

#endif