#ifndef NUMERIC_H
#define NUMERIC_H

#include "Structures.hpp"

/*!
    Calculates L2 norm of a difference between two vectors.
    
    This difference is used to check the convergence of velocities
    during the iteration process of PRIME (PRessure Implicit Momentum 
    Explicit) method.  

    result = |v1-v2|_l2
*/
double calculate_vec_diff_L2_norm(DoubleArray2D& v1, DoubleArray2D& v2, int nx, int ny);

#endif
