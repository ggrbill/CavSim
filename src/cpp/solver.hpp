

/*!
    SOR(Successive Over Relaxation) Solver for Structured grid.
*/
void SOR_structured(
    double** Ap, double** Aw, double** Ae, double** An, double** As, 
    double** x, double** xn, double** b, 
    int const size,
    int const MAX_IT, float const w, // relaxation factor
    double const tol = 0.0001// convergence tolerance
);