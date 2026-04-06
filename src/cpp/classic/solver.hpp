#ifndef CAVITY_SOLVER_H
#define CAVITY_SOLVER_H

#include "Structures.hpp"

/*!
    SOR(Successive Over Relaxation) Solver for Structured grid.
*/
void SOR_structured(
    DoubleArray2D& Ap, DoubleArray2D& Aw, DoubleArray2D& Ae, DoubleArray2D& An, DoubleArray2D& As, 
    DoubleArray2D& x, DoubleArray2D& xn, DoubleArray2D& b, 
    int const size,
    int const MAX_IT, float const w, // relaxation factor
    double const tol = 0.0001// convergence tolerance
);

#endif
