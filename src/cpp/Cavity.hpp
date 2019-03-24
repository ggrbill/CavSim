#ifndef CAVITY_H
#define CAVITY_H

class Cavity
{
public:
    Cavity(
        double length, 
        double height, 
        int n_x, 
        int n_y,  
        double rho,   // Density 
        double mi,    // Viscosity
        double U_lid  // Upper lid velocity
    );
    ~Cavity() {}

protected:
    double L, H;
    int n_x, n_y; 
    double rho, mu, U;
    double dx, dy;

private:

}; // class Cavity


#endif