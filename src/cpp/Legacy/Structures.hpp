#ifndef STRUCTURES_HPP
#define STRUCTURES_HPP

// Control Volume Boundaries interpolation
struct CVBoundaries{
	double e;
	double w;
	double n;
	double s;
};

struct CavSimAux{
	CavSimAux(int n_x, int n_y);
	~CavSimAux();

	// Number of divisions
	int n_x;
	int n_y;

	CVBoundaries **alpha_x;
	CVBoundaries **beta_x;
	CVBoundaries **alpha_y;
	CVBoundaries **beta_y;
};

struct CavSimData{
	CavSimData(int n_x, int n_y);
	~CavSimData();

	// Number of divisions
	int n_x;
	int n_y;

	// Coefficients Matrices
	// x-velocity - u
	double **Ap_u;
	double **Aw_u;
	double **Ae_u;
	double **An_u;
	double **As_u;
	double **B_u;
	// y-velocity - v
	double **Ap_v;
	double **Aw_v;
	double **Ae_v;
	double **An_v;
	double **As_v;
	double **B_v;
	// Pressure
	double **Ap_p;
	double **Aw_p;
	double **Ae_p;
	double **An_p;
	double **As_p;
	double **B_p;
};

struct CavSimResult{
	CavSimResult(int n_x, int n_y);
	~CavSimResult();
	
	// Number of divisions
	int n_x;
	int n_y;

	// Velocities
	double **u;
	double **v;
	double **uOLD;
	double **vOLD;
	double **u_hat;
	double **v_hat;
	// Pressure
	double **P;
	double **Pn;
};

// Allocate Vectors 2d 
void allocate_vector_2d(double** &vec, int nx, int ny);
// Deallocate Vectors 2d 
void deallocate_vector_2d(double** &vec, int nx, int ny);

// Allocate CVBoundaries 2d 
void allocate_cvbound_2d(CVBoundaries** &vec, int nx, int ny);
// Deallocate CVBoundaries 2d 
void deallocate_cvbound_2d(CVBoundaries** &vec, int nx, int ny);

#endif
