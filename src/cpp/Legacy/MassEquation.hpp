#ifndef MASS_EQUATION_HPP
#define MASS_EQUATION_HPP

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
    double **u_hat,
    double **v_hat,
    double **Ap_u,
    double **Ap_v,
    double **Ap_p,
    double **Ae_p,
    double **Aw_p,
    double **As_p,
    double **An_p,
    double **B_p
);

#endif
