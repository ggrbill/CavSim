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


CavSimData::CavSimData(int n_x, int n_y)
    : n_x(n_x), n_y(n_y)
{
    allocate_vector_2d(this->Ap_u , this->n_x-1, this->n_y);
	allocate_vector_2d(this->Aw_u , this->n_x-1, this->n_y);
	allocate_vector_2d(this->Ae_u , this->n_x-1, this->n_y);
	allocate_vector_2d(this->An_u , this->n_x-1, this->n_y);
	allocate_vector_2d(this->As_u , this->n_x-1, this->n_y);
	allocate_vector_2d(this->B_u  , this->n_x-1, this->n_y);

    allocate_vector_2d(this->Ap_v , this->n_x, this->n_y-1);
	allocate_vector_2d(this->Aw_v , this->n_x, this->n_y-1);
	allocate_vector_2d(this->Ae_v , this->n_x, this->n_y-1);
	allocate_vector_2d(this->An_v , this->n_x, this->n_y-1);
	allocate_vector_2d(this->As_v , this->n_x, this->n_y-1);
	allocate_vector_2d(this->B_v  , this->n_x, this->n_y-1);

    allocate_vector_2d(this->Ap_p, this->n_x, this->n_y);
	allocate_vector_2d(this->Aw_p, this->n_x, this->n_y);
	allocate_vector_2d(this->Ae_p, this->n_x, this->n_y);
	allocate_vector_2d(this->An_p, this->n_x, this->n_y);
	allocate_vector_2d(this->As_p, this->n_x, this->n_y);
	allocate_vector_2d(this->B_p , this->n_x, this->n_y);
}

CavSimData::~CavSimData(){
    deallocate_vector_2d(Ap_u , n_x-1, n_y);
	deallocate_vector_2d(Aw_u , n_x-1, n_y);
	deallocate_vector_2d(Ae_u , n_x-1, n_y);
	deallocate_vector_2d(An_u , n_x-1, n_y);
	deallocate_vector_2d(As_u , n_x-1, n_y);
	deallocate_vector_2d(B_u  , n_x-1, n_y);

    deallocate_vector_2d(Ap_v , n_x, n_y-1);
	deallocate_vector_2d(Aw_v , n_x, n_y-1);
	deallocate_vector_2d(Ae_v , n_x, n_y-1);
	deallocate_vector_2d(An_v , n_x, n_y-1);
	deallocate_vector_2d(As_v , n_x, n_y-1);
	deallocate_vector_2d(B_v  , n_x, n_y-1);

    deallocate_vector_2d(Ap_p, n_x, n_y);
	deallocate_vector_2d(Aw_p, n_x, n_y);
	deallocate_vector_2d(Ae_p, n_x, n_y);
	deallocate_vector_2d(An_p, n_x, n_y);
	deallocate_vector_2d(As_p, n_x, n_y);
	deallocate_vector_2d(B_p , n_x, n_y);
}

CavSimResult::CavSimResult(int n_x, int n_y)
    : n_x(n_x), n_y(n_y)
{
    allocate_vector_2d(this->u    , this->n_x-1, this->n_y);
	allocate_vector_2d(this->uOLD , this->n_x-1, this->n_y);
	allocate_vector_2d(this->u_hat, this->n_x-1, this->n_y);

    allocate_vector_2d(this->v    , this->n_x, this->n_y-1);
	allocate_vector_2d(this->vOLD , this->n_x, this->n_y-1);
	allocate_vector_2d(this->v_hat, this->n_x, this->n_y-1);

    allocate_vector_2d(this->P   , this->n_x, this->n_y);
	allocate_vector_2d(this->Pn  , this->n_x, this->n_y);
}

CavSimResult::~CavSimResult(){
    deallocate_vector_2d(this->u    , this->n_x-1, this->n_y);
	deallocate_vector_2d(this->uOLD , this->n_x-1, this->n_y);
	deallocate_vector_2d(this->u_hat, this->n_x-1, this->n_y);
    
    deallocate_vector_2d(this->v    , this->n_x, this->n_y-1);
	deallocate_vector_2d(this->vOLD , this->n_x, this->n_y-1);
	deallocate_vector_2d(this->v_hat, this->n_x, this->n_y-1);

    deallocate_vector_2d(this->P   , this->n_x, this->n_y);
	deallocate_vector_2d(this->Pn  , this->n_x, this->n_y);
}

CavSimAux::CavSimAux(int n_x, int n_y)
    : n_x(n_x), n_y(n_y)
{
    allocate_cvbound_2d(alpha_x, n_x-1, n_y);
	allocate_cvbound_2d(beta_x , n_x-1, n_y);

	allocate_cvbound_2d(alpha_y, n_x, n_y-1);
	allocate_cvbound_2d(beta_y , n_x, n_y-1);
}

CavSimAux::~CavSimAux(){
    deallocate_cvbound_2d(alpha_x, n_x-1, n_y);
	deallocate_cvbound_2d(beta_x , n_x-1, n_y);
	
	deallocate_cvbound_2d(alpha_y, n_x, n_y-1);
	deallocate_cvbound_2d(beta_y , n_x, n_y-1);
}