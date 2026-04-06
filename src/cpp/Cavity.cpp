#include "Cavity.hpp"
#include <stdexcept>

CavitySetup::CavitySetup(
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
, dy((double)H/n_y)
{}

Cavity::Cavity(std::shared_ptr<CavitySetup> setup)
: setup(setup)
, p_coeffs(std::make_shared<PrimeCoefficients>(setup->n_x, setup->n_y))
, results(std::make_shared<CavSimResult>(setup->n_x, setup->n_y))
{}

void Cavity::run_simulation()
{
    throw std::logic_error("Cavity::run_simulation is not implemented yet.");
}
