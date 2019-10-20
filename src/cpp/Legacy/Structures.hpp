#ifndef STRUCTURES_HPP
#define STRUCTURES_HPP

// Control Volume Boundaries interpolation
struct CVBoundaries{
	double e;
	double w;
	double n;
	double s;
};

// Allocate Vectors 2d 
void allocate_vector_2d(double** &vec, int nx, int ny);
// deallocate Vectors 2d 
void deallocate_vector_2d(double** &vec, int nx, int ny);

// Allocate CVBoundaries 2d 
void allocate_cvbound_2d(CVBoundaries** &vec, int nx, int ny);
// deallocate CVBoundaries 2d 
void deallocate_cvbound_2d(CVBoundaries** &vec, int nx, int ny);

#endif
