#ifndef PRIME_CORRECTION_HPP
#define PRIME_CORRECTION_HPP

/*!
    Calculates correction velocity in x direction (u_hat).
    
    Those velocities will be used to correct real velocities
    in x direction (u).  
*/
void calculate_u_hat(
    int nv,
    double** Ap_u,
    double** Ae_u,
    double** Aw_u,
    double** As_u,
    double** An_u,
    double** B_u,
    double** u,
    double** u_hat
);

/*!
    Calculates correction velocity in y direction (v_hat).
    
    Those velocities will be used to correct real velocities
    in y direction (v).  
*/
void calculate_v_hat(
    int nv,
    double** Ap_v,
    double** Ae_v,
    double** Aw_v,
    double** As_v,
    double** An_v,
    double** B_v,
    double** v,
    double** v_hat
);

/*!
    Corrects the velocities u and v.
    
    The correction is preformed to obtain the solution 
    of u and v.
*/
void correct_u_v(
	int nv,
	double dx,
	double dy,
	double **Pn,
	double **Ap_u,
	double **uOLD,
	double **u_hat,
	double **u,
	double **Ap_v,
	double **vOLD,
	double **v_hat,
	double **v
);

#endif
