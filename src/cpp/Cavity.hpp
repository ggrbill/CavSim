#ifndef CAVITY_HPP
#define CAVITY_HPP

#include<memory>

class CavitySetup
{
public:
    CavitySetup(
        double length, 
        double height, 
        int n_x, 
        int n_y,  
        double rho,   // Density 
        double mu,    // Viscosity
        double U_lid  // Upper lid velocity
    );
    ~CavitySetup() {}

    double L, H;
    int n_x, n_y; 
    double rho, mu, U;
    double dx, dy;
}; // class Cavity Setup


class Cavity
{
public:
    Cavity(std::shared_ptr<CavitySetup> setup);
    ~Cavity() {}

protected:
    std::shared_ptr<CavitySetup> setup;

private:

}; // class Cavity


#endif