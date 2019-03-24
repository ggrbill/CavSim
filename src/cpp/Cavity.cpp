#include "Cavity.hpp"

Cavity::Cavity(
    double length, 
    double height, 
    int n_x, 
    int n_y,  
    double rho,   
    double mu,
    double U_lid 
)
: L(length)
, H(height)
, n_x(n_x)
, n_y(n_y)
, rho(rho)
, mu(mu)
, U(U_lid) 
, dx((double)L/n_x)
, dy((double)L/n_y)
{}