#include "Structures.hpp"

void allocate_vector_2d(double** &vec, int nx, int ny){
    // Allocating
    vec = new double *[nx];
    for (int i = 0; i < nx; ++i){
        vec[i] = new double [ny];
    }
    // Initializing
    for (int i = 0; i < nx; ++i) {
        for (int j = 0; j < ny; ++j) {
            vec[i][j] = 0.0;
        }
    }
}

void deallocate_vector_2d(double** &vec, int nx, int ny){
    // Allocating
    for (int i = 0; i < nx; ++i){
        delete [] vec[i];
    }
    delete [] vec;
}

void allocate_cvbound_2d(CVBoundaries** &vec, int nx, int ny){
    // Allocating
    vec = new CVBoundaries *[nx];
    for (int i = 0; i < nx; ++i){
        vec[i] = new CVBoundaries [ny];
    }
    // Initializing
    for (int i = 0; i < nx; ++i) {
        for (int j = 0; j < ny; ++j) {
            vec[i][j].e = 0.0;
            vec[i][j].w = 0.0;
            vec[i][j].n = 0.0;
            vec[i][j].s = 0.0;
        }
    }
}

void deallocate_cvbound_2d(CVBoundaries** &vec, int nx, int ny){
    // Allocating
    for (int i = 0; i < nx; ++i){
        delete [] vec[i];
    }
    delete [] vec;
}