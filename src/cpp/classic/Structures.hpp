#ifndef STRUCTURES_HPP
#define STRUCTURES_HPP

#include <vector>

using DoubleArray2D = std::vector<std::vector<double>>;

// Control Volume Boundaries interpolation
struct Faces{
	double e = 0.0;
	double w = 0.0;
	double n = 0.0;
	double s = 0.0;
};

using FacesArray2D = std::vector<std::vector<Faces>>;

struct CavSimAux{
	CavSimAux(int n_x, int n_y);
	~CavSimAux() = default;

	FacesArray2D alpha_x;
	FacesArray2D beta_x;
	FacesArray2D alpha_y;
	FacesArray2D beta_y;
};

struct CavSimData{
	CavSimData(int n_x, int n_y);
	~CavSimData() = default;

	// Coefficients Matrices
	// x-velocity - u
	DoubleArray2D Ap_u;
	DoubleArray2D Aw_u;
	DoubleArray2D Ae_u;
	DoubleArray2D An_u;
	DoubleArray2D As_u;
	DoubleArray2D B_u;
	// y-velocity - v
	DoubleArray2D Ap_v;
	DoubleArray2D Aw_v;
	DoubleArray2D Ae_v;
	DoubleArray2D An_v;
	DoubleArray2D As_v;
	DoubleArray2D B_v;
	// Pressure
	DoubleArray2D Ap_p;
	DoubleArray2D Aw_p;
	DoubleArray2D Ae_p;
	DoubleArray2D An_p;
	DoubleArray2D As_p;
	DoubleArray2D B_p;
};

struct CavSimResult{
	CavSimResult(int n_x, int n_y);
	~CavSimResult() = default;

	// Velocities
	DoubleArray2D u;
	DoubleArray2D v;
	DoubleArray2D u_old;
	DoubleArray2D v_old;
	DoubleArray2D u_hat;
	DoubleArray2D v_hat;
	// Pressure
	DoubleArray2D P;
	DoubleArray2D Pn;
};

#endif
