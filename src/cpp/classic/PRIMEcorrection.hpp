#ifndef PRIME_CORRECTION_HPP
#define PRIME_CORRECTION_HPP

#include "Structures.hpp"

/*!
    Calculates correction velocity in x direction (u_hat).
    
    Those velocities will be used to correct real velocities
    in x direction (u).  
*/
void calculate_u_hat(
    int nv,
    DoubleArray2D& Ap_u,
    DoubleArray2D& Ae_u,
    DoubleArray2D& Aw_u,
    DoubleArray2D& As_u,
    DoubleArray2D& An_u,
    DoubleArray2D& B_u,
    DoubleArray2D& u,
    DoubleArray2D& u_hat
);

/*!
    Calculates correction velocity in y direction (v_hat).
    
    Those velocities will be used to correct real velocities
    in y direction (v).  
*/
void calculate_v_hat(
    int nv,
    DoubleArray2D& Ap_v,
    DoubleArray2D& Ae_v,
    DoubleArray2D& Aw_v,
    DoubleArray2D& As_v,
    DoubleArray2D& An_v,
    DoubleArray2D& B_v,
    DoubleArray2D& v,
    DoubleArray2D& v_hat
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
	DoubleArray2D& Pn,
	DoubleArray2D& Ap_u,
	DoubleArray2D& u_old,
	DoubleArray2D& u_hat,
	DoubleArray2D& u,
	DoubleArray2D& Ap_v,
	DoubleArray2D& v_old,
	DoubleArray2D& v_hat,
	DoubleArray2D& v
);

#endif
