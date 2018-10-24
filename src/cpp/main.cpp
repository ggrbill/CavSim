//Este Programa Simula o Problema da Cavidade Isolada a Oeste, Leste e Sul com uma Velocidade Tangencial U em Norte

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <math.h>

#include "solver.hpp"

using namespace std;

//MATRIZES DE RESULTADOS
double **u;
double **v;
double **uOLD;
double **vOLD;
double **P;
double **Pn;
double **u_hat;
double **v_hat;


// MATRIZES DE COEFICIENTES

// Velocidade em x - u
double **Ap_u;
double **Aw_u;
double **Ae_u;
double **An_u;
double **As_u;
double **B_u;
// Velocidade em y - v
double **Ap_v;
double **Aw_v;
double **Ae_v;
double **An_v;
double **As_v;
double **B_v;
// Press�o
double **Ap_p;
double **Aw_p;
double **Ae_p;
double **An_p;
double **As_p;
double **B_p;

// VARI�VEIS DE ENTRADA
double L;	// Largura 
double H;	// Altura
int    nv;	// Numero de Linhas e Colunas
double rho;	// Massa Especifica
double U;	// Velocidade na Fronteira Norte
double mi;	// Viscosidade

// VARI�VEIS DE INTERPOLA��O NAS FRONTEIRAS
struct fronteiras{
	double e;
	double w;
	double n;
	double s;
};
fronteiras **Re_x;
fronteiras **alpha_x;
fronteiras **beta_x;
fronteiras **Re_y;
fronteiras **alpha_y;
fronteiras **beta_y;

// Delta X e Delta Y
double dx = 0.;
double dy = 0.;

//Variaveis auxiliares
int verifica = 0;
int N_ITE = 0;
double T = 0;

void read_in() // L� os Dados de Entrada
{
	
	fstream fin("./inCav.txt");
	fin >>  L;
	fin >>  H;
	fin >>  nv;
	fin >>  rho;
	fin >>  U; 
	fin >>  mi;
	fin.close();
}

void Alocate() // Aloca e Inicializa Vari�veis
{
	
	//Alocando Mem�ria

	// u   Velocidade-x
	u     = new double *[nv-1];
	uOLD  = new double *[nv-1];
	u_hat = new double *[nv-1];
	Ap_u  = new double *[nv-1];
	Aw_u  = new double *[nv-1];
	Ae_u  = new double *[nv-1];
	An_u  = new double *[nv-1];
	As_u  = new double *[nv-1];
	B_u   = new double *[nv-1];
	for(int i=0;i<(nv-1);i++)
	{
		u[i]    = new double [nv];
		uOLD[i] = new double [nv];
		u_hat[i]= new double [nv];
		Ap_u[i] = new double [nv];
		Aw_u[i] = new double [nv];
		Ae_u[i] = new double [nv];
		An_u[i] = new double [nv];
		As_u[i] = new double [nv];
		B_u[i]  = new double [nv];
	}
	for(int i=0;i<(nv-1);i++)
	{
		for(int j=0;j<nv;j++)
		{
			u[i][j]    = 0.;
			uOLD[i][j] = 0.;
			u_hat[i][j]= 0.;
			Ap_u[i][j] = 0.;
			Aw_u[i][j] = 0.;
			Ae_u[i][j] = 0.;
			An_u[i][j] = 0.;
			As_u[i][j] = 0.;
			B_u[i][j]  = 0.;
		}
	}
	// v   Velocidade-y
	v     = new double *[nv];
	vOLD  = new double *[nv];
	v_hat = new double *[nv];
	Ap_v  = new double *[nv];
	Aw_v  = new double *[nv];
	Ae_v  = new double *[nv];
	An_v  = new double *[nv];
	As_v  = new double *[nv];
	B_v   = new double *[nv];
	for(int i=0;i<nv;i++)
	{
		v[i]    = new double [nv-1];
		vOLD[i] = new double [nv-1];
		v_hat[i]= new double [nv-1];
		Ap_v[i] = new double [nv-1];
		Aw_v[i] = new double [nv-1];
		Ae_v[i] = new double [nv-1];
		An_v[i] = new double [nv-1];
		As_v[i] = new double [nv-1];
		B_v[i]  = new double [nv-1];
	}
	for(int i=0;i<nv;i++)
	{
		for(int j=0;j<(nv-1);j++)
		{
			v[i][j]    = 0.;
			vOLD[i][j] = 0.;
			v_hat[i][j]= 0.;
			Ap_v[i][j] = 0.;
			Aw_v[i][j] = 0.;
			Ae_v[i][j] = 0.;
			An_v[i][j] = 0.;
			As_v[i][j] = 0.;
			B_v[i][j]  = 0.;
		}
	}
	// Press�o
	P     = new double *[nv];
	Pn    = new double *[nv];
	Ap_p  = new double *[nv];
	Aw_p  = new double *[nv];
	Ae_p  = new double *[nv];
	An_p  = new double *[nv];
	As_p  = new double *[nv];
	B_p   = new double *[nv];
	for(int i=0;i<nv;i++)
	{
		P[i]    = new double [nv];
		Pn[i]   = new double [nv];
		Ap_p[i] = new double [nv];
		Aw_p[i] = new double [nv];
		Ae_p[i] = new double [nv];
		An_p[i] = new double [nv];
		As_p[i] = new double [nv];
		B_p[i]  = new double [nv];
	}
	for(int i=0;i<nv;i++)
	{
		for(int j=0;j<nv;j++)
		{
			P[i][j]    = 0.;
			Pn[i][j]   = 0.;
			Ap_p[i][j] = 0.;
			Aw_p[i][j] = 0.;
			Ae_p[i][j] = 0.;
			An_p[i][j] = 0.;
			As_p[i][j] = 0.;
			B_p[i][j] = 0.;
		}
	}
	// Dados nas Fronteiras
	
	// Dire��o x - Malha de u
	Re_x     = new fronteiras *[nv-1];
	alpha_x  = new fronteiras *[nv-1];
	beta_x   = new fronteiras *[nv-1];
	for(int i=0;i<(nv-1);i++)
	{
		Re_x[i]     = new fronteiras [nv];
		alpha_x[i]  = new fronteiras [nv];
		beta_x[i]   = new fronteiras [nv];
	}
	for(int i=0;i<(nv-1);i++)
	{
		for(int j=0;j<nv;j++)
		{
			Re_x[i][j].e     = 0.;
			Re_x[i][j].w     = 0.;
			Re_x[i][j].n     = 0.;
			Re_x[i][j].s     = 0.;
			alpha_x[i][j].e  = 0.;
			alpha_x[i][j].w  = 0.;
			alpha_x[i][j].n  = 0.;
			alpha_x[i][j].s  = 0.;
			beta_x[i][j].e   = 0.;
			beta_x[i][j].w   = 0.;
			beta_x[i][j].n   = 0.;
			beta_x[i][j].s   = 0.;
		}
	}
	// Dire��o y - Malha de v
	Re_y     = new fronteiras *[nv];
	alpha_y  = new fronteiras *[nv];
	beta_y   = new fronteiras *[nv];
	for(int i=0;i<nv;i++)
	{
		Re_y[i]     = new fronteiras [nv-1];
		alpha_y[i]  = new fronteiras [nv-1];
		beta_y[i]   = new fronteiras [nv-1];
	}
	for(int i=0;i<nv;i++)
	{
		for(int j=0;j<(nv-1);j++)
		{
			Re_y[i][j].e     = 0.;
			Re_y[i][j].w     = 0.;
			Re_y[i][j].n     = 0.;
			Re_y[i][j].s     = 0.;
			alpha_y[i][j].e  = 0.;
			alpha_y[i][j].w  = 0.;
			alpha_y[i][j].n  = 0.;
			alpha_y[i][j].s  = 0.;
			beta_y[i][j].e   = 0.;
			beta_y[i][j].w   = 0.;
			beta_y[i][j].n   = 0.;
			beta_y[i][j].s   = 0.;
		}
	}
	
}

void Calc_Coef_NS_X()  //Calculo dos Coeficientes de u - NS_x
{
	//Volume do Canto Inferior Esquerdo
	Re_x[0][0].e = (rho*((u[0][0]+u[1][0])/2.)*dx)/mi;
	Re_x[0][0].n = (rho*((v[0][0]+v[1][0])/2.)*dy)/mi;
	Re_x[0][0].w = (rho*(u[0][0]/2.)*dx)/mi;
	Re_x[0][0].s = 0.;
	alpha_x[0][0].e = pow(Re_x[0][0].e,2.)/(10.+2.*pow(Re_x[0][0].e,2.));
	alpha_x[0][0].w = pow(Re_x[0][0].w,2.)/(10.+2.*pow(Re_x[0][0].w,2.));
	alpha_x[0][0].n = pow(Re_x[0][0].n,2.)/(10.+2.*pow(Re_x[0][0].n,2.));
	alpha_x[0][0].s = pow(Re_x[0][0].s,2.)/(10.+2.*pow(Re_x[0][0].s,2.));
	beta_x[0][0].e = (1.+(0.005*pow(Re_x[0][0].e,2.)))/(1.+(0.05*pow(Re_x[0][0].e,2.)));
	beta_x[0][0].w = (1.+(0.005*pow(Re_x[0][0].w,2.)))/(1.+(0.05*pow(Re_x[0][0].w,2.)));
	beta_x[0][0].n = (1.+(0.005*pow(Re_x[0][0].n,2.)))/(1.+(0.05*pow(Re_x[0][0].n,2.)));
	beta_x[0][0].s = (1.+(0.005*pow(Re_x[0][0].s,2.)))/(1.+(0.05*pow(Re_x[0][0].s,2.)));

	if(Re_x[0][0].e < 0) { alpha_x[0][0].e = -alpha_x[0][0].e; }
	if(Re_x[0][0].w < 0) { alpha_x[0][0].w = -alpha_x[0][0].w; }
	if(Re_x[0][0].n < 0) { alpha_x[0][0].n = -alpha_x[0][0].n; }
	if(Re_x[0][0].s < 0) { alpha_x[0][0].s = -alpha_x[0][0].s; }

	Ae_u[0][0] = (((mi*beta_x[0][0].e*dy)/dx)-((rho*((u[0][0]+u[1][0])/2.)*dy)*(0.5-alpha_x[0][0].e)));
	An_u[0][0] = (((mi*beta_x[0][0].n*dx)/dy)-((rho*((v[1][0]+v[0][0])/2.)*dx)*(0.5-alpha_x[0][0].n)));
	Aw_u[0][0] = (((mi*beta_x[0][0].w*dy)/dx)+((rho*((u[0][0])/2.)*dy)*(0.5+alpha_x[0][0].w)));
	As_u[0][0] = ((2.*mi*beta_x[0][0].s*dx)/dy);
	Ap_u[0][0] = Ae_u[0][0] + Aw_u[0][0] + An_u[0][0] + As_u[0][0];
	B_u[0][0] = 0.;

	//Volume do Canto Inferior Direito
	Re_x[nv-2][0].e = (rho*((u[nv-2][0])/2.)*dx)/mi;
	Re_x[nv-2][0].n = (rho*((v[nv-2][0]+v[nv-1][0])/2.)*dx)/mi;
	Re_x[nv-2][0].w = (rho*((u[nv-2][0]+u[nv-3][0])/2.)*dx)/mi; 
	Re_x[nv-2][0].s = 0.;
	alpha_x[nv-2][0].e = pow(Re_x[nv-2][0].e,2.)/(10.+2.*pow(Re_x[nv-2][0].e,2.));
	alpha_x[nv-2][0].w = pow(Re_x[nv-2][0].w,2.)/(10.+2.*pow(Re_x[nv-2][0].w,2.));
	alpha_x[nv-2][0].n = pow(Re_x[nv-2][0].n,2.)/(10.+2.*pow(Re_x[nv-2][0].n,2.));
	alpha_x[nv-2][0].s = pow(Re_x[nv-2][0].s,2.)/(10.+2.*pow(Re_x[nv-2][0].s,2.));
	beta_x[nv-2][0].e = (1.+(0.005*pow(Re_x[nv-2][0].e,2.)))/(1.+(0.05*pow(Re_x[nv-2][0].e,2.)));
	beta_x[nv-2][0].w = (1.+(0.005*pow(Re_x[nv-2][0].w,2.)))/(1.+(0.05*pow(Re_x[nv-2][0].w,2.)));
	beta_x[nv-2][0].n = (1.+(0.005*pow(Re_x[nv-2][0].n,2.)))/(1.+(0.05*pow(Re_x[nv-2][0].n,2.)));
	beta_x[nv-2][0].s = (1.+(0.005*pow(Re_x[nv-2][0].s,2.)))/(1.+(0.05*pow(Re_x[nv-2][0].s,2.)));

	if(Re_x[nv-2][0].e < 0) { alpha_x[nv-2][0].e = -alpha_x[nv-2][0].e; }
	if(Re_x[nv-2][0].w < 0) { alpha_x[nv-2][0].w = -alpha_x[nv-2][0].w; }
	if(Re_x[nv-2][0].n < 0) { alpha_x[nv-2][0].n = -alpha_x[nv-2][0].n; }
	if(Re_x[nv-2][0].s < 0) { alpha_x[nv-2][0].s = -alpha_x[nv-2][0].s; }

	Aw_u[nv-2][0] = (((mi*beta_x[nv-2][0].w*dy)/dx)+((rho*((u[nv-2][0]+u[nv-3][0])/2.)*dy)*(0.5+alpha_x[nv-2][0].w)));
	An_u[nv-2][0] = (((mi*beta_x[nv-2][0].n*dx)/dy)-((rho*((v[nv-2][0]+v[nv-1][0])/2.)*dx)*(0.5-alpha_x[nv-2][0].n)));
	Ae_u[nv-2][0] = (((mi*beta_x[nv-2][0].e*dy)/dx)-((rho*((u[nv-2][0])/2.)*dy)*(0.5-alpha_x[nv-2][0].e)));
	As_u[nv-2][0] = ((2.*mi*beta_x[nv-2][0].s*dx)/dy);
	Ap_u[nv-2][0] = Ae_u[nv-2][0] + Aw_u[nv-2][0] + An_u[nv-2][0] + As_u[nv-2][0];
	B_u[nv-2][0] = 0.;
	
	//Volume do Canto Superior Esquerdo
	Re_x[0][nv-1].e = (rho*((u[0][nv-1]+u[1][nv-1])/2.)*dx)/mi;
	Re_x[0][nv-1].n = 0.; //(rho*(U)*dx)/mi;
	Re_x[0][nv-1].w = (rho*(u[0][nv-1]/2.)*dx)/mi;
	Re_x[0][nv-1].s = (rho*((v[0][nv-2]+v[1][nv-2])/2.)*dx)/mi;
	alpha_x[0][nv-1].e = pow(Re_x[0][nv-1].e,2.)/(10.+2.*pow(Re_x[0][nv-1].e,2.));
	alpha_x[0][nv-1].w = pow(Re_x[0][nv-1].w,2.)/(10.+2.*pow(Re_x[0][nv-1].w,2.));
	alpha_x[0][nv-1].n = pow(Re_x[0][nv-1].n,2.)/(10.+2.*pow(Re_x[0][nv-1].n,2.));
	alpha_x[0][nv-1].s = pow(Re_x[0][nv-1].s,2.)/(10.+2.*pow(Re_x[0][nv-1].s,2.));
	beta_x[0][nv-1].e = (1.+(0.005*pow(Re_x[0][nv-1].e,2.)))/(1.+(0.05*pow(Re_x[0][nv-1].e,2.)));
	beta_x[0][nv-1].w = (1.+(0.005*pow(Re_x[0][nv-1].w,2.)))/(1.+(0.05*pow(Re_x[0][nv-1].w,2.)));
	beta_x[0][nv-1].n = (1.+(0.005*pow(Re_x[0][nv-1].n,2.)))/(1.+(0.05*pow(Re_x[0][nv-1].n,2.)));
	beta_x[0][nv-1].s = (1.+(0.005*pow(Re_x[0][nv-1].s,2.)))/(1.+(0.05*pow(Re_x[0][nv-1].s,2.)));

	if(Re_x[0][nv-1].e < 0) { alpha_x[0][nv-1].e = -alpha_x[0][nv-1].e; }
	if(Re_x[0][nv-1].w < 0) { alpha_x[0][nv-1].w = -alpha_x[0][nv-1].w; }
	if(Re_x[0][nv-1].n < 0) { alpha_x[0][nv-1].n = -alpha_x[0][nv-1].n; }
	if(Re_x[0][nv-1].s < 0) { alpha_x[0][nv-1].s = -alpha_x[0][nv-1].s; }

	Ae_u[0][nv-1] = (((mi*beta_x[0][nv-1].e*dy)/dx)-((rho*((u[0][nv-1]+u[1][nv-1])/2.)*dy)*(0.5-alpha_x[0][nv-1].e)));
	As_u[0][nv-1] = (((mi*beta_x[0][nv-1].s*dx)/dy)+((rho*((v[1][nv-2]+v[0][nv-2])/2.)*dx)*(0.5+alpha_x[0][nv-1].s)));
	Aw_u[0][nv-1] = (((mi*beta_x[0][nv-1].w*dy)/dx)+((rho*((u[0][nv-1])/2.)*dy)*(0.5+alpha_x[0][nv-1].w)));
	An_u[0][nv-1] = ((2.*mi*beta_x[0][nv-1].n*dx)/dy);
	Ap_u[0][nv-1] = Ae_u[0][nv-1] + Aw_u[0][nv-1] + An_u[0][nv-1] + As_u[0][nv-1];
	B_u[0][nv-1] = ((2.*mi*beta_x[0][nv-1].n*U*dx)/dy);
	
	//Volume do Canto Superior Direito
	Re_x[nv-2][nv-1].e = (rho*((u[nv-2][nv-1])/2.)*dx)/mi;
	Re_x[nv-2][nv-1].n = 0.; //(rho*(U)*dx)/mi;
	Re_x[nv-2][nv-1].w = (rho*((u[nv-2][nv-1]+u[nv-3][nv-1])/2.)*dx)/mi;
	Re_x[nv-2][nv-1].s = (rho*((v[nv-2][nv-2]+v[nv-1][nv-2])/2.)*dx)/mi;
	alpha_x[nv-2][nv-1].e = pow(Re_x[nv-2][nv-1].e,2.)/(10.+2.*pow(Re_x[nv-2][nv-1].e,2.));
	alpha_x[nv-2][nv-1].w = pow(Re_x[nv-2][nv-1].w,2.)/(10.+2.*pow(Re_x[nv-2][nv-1].w,2.));
	alpha_x[nv-2][nv-1].n = pow(Re_x[nv-2][nv-1].n,2.)/(10.+2.*pow(Re_x[nv-2][nv-1].n,2.));
	alpha_x[nv-2][nv-1].s = pow(Re_x[nv-2][nv-1].s,2.)/(10.+2.*pow(Re_x[nv-2][nv-1].s,2.));
	beta_x[nv-2][nv-1].e = (1.+(0.005*pow(Re_x[nv-2][nv-1].e,2.)))/(1.+(0.05*pow(Re_x[nv-2][nv-1].e,2.)));
	beta_x[nv-2][nv-1].w = (1.+(0.005*pow(Re_x[nv-2][nv-1].w,2.)))/(1.+(0.05*pow(Re_x[nv-2][nv-1].w,2.)));
	beta_x[nv-2][nv-1].n = (1.+(0.005*pow(Re_x[nv-2][nv-1].n,2.)))/(1.+(0.05*pow(Re_x[nv-2][nv-1].n,2.)));
	beta_x[nv-2][nv-1].s = (1.+(0.005*pow(Re_x[nv-2][nv-1].s,2.)))/(1.+(0.05*pow(Re_x[nv-2][nv-1].s,2.)));

	if(Re_x[nv-2][nv-1].e < 0) { alpha_x[nv-2][nv-1].e = -alpha_x[nv-2][nv-1].e; }
	if(Re_x[nv-2][nv-1].w < 0) { alpha_x[nv-2][nv-1].w = -alpha_x[nv-2][nv-1].w; }
	if(Re_x[nv-2][nv-1].n < 0) { alpha_x[nv-2][nv-1].n = -alpha_x[nv-2][nv-1].n; }
	if(Re_x[nv-2][nv-1].s < 0) { alpha_x[nv-2][nv-1].s = -alpha_x[nv-2][nv-1].s; }

	Aw_u[nv-2][nv-1] = (((mi*beta_x[nv-2][nv-1].w*dy)/dx)+((rho*((u[nv-2][nv-1]+u[nv-3][nv-1])/2.)*dy)*(0.5+alpha_x[nv-2][nv-1].w)));
	As_u[nv-2][nv-1] = (((mi*beta_x[nv-2][nv-1].s*dx)/dy)+((rho*((v[nv-2][nv-2]+v[nv-1][nv-2])/2.)*dx)*(0.5+alpha_x[nv-2][nv-1].s)));
	Ae_u[nv-2][nv-1] = (((mi*beta_x[nv-2][nv-1].e*dy)/dx)-((rho*((u[nv-2][nv-1])/2.)*dy)*(0.5-alpha_x[nv-2][nv-1].e)));
	An_u[nv-2][nv-1] = ((2.*mi*beta_x[nv-2][nv-1].n*dx)/dy);
	Ap_u[nv-2][nv-1] = Ae_u[nv-2][nv-1] + Aw_u[nv-2][nv-1] + An_u[nv-2][nv-1] + As_u[nv-2][nv-1];
	B_u[nv-2][nv-1] = ((2.*mi*beta_x[nv-2][nv-1].n*U*dx)/dy);

	//Volumes da Fronteira Norte
	for(int i=1;i<(nv-2);i++)
	{
		Re_x[i][nv-1].e = (rho*((u[i][nv-1]+u[i+1][nv-1])/2.)*dx)/mi;
		Re_x[i][nv-1].s = (rho*((v[i][nv-2]+v[i+1][nv-2])/2.)*dx)/mi;
		Re_x[i][nv-1].w = (rho*((u[i][nv-1]+u[i-1][nv-1])/2.)*dx)/mi;
		Re_x[i][nv-1].n = 0.; //(rho*(U)*dx)/mi;
		alpha_x[i][nv-1].e = pow(Re_x[i][nv-1].e,2.)/(10.+2.*pow(Re_x[i][nv-1].e,2.));
		alpha_x[i][nv-1].w = pow(Re_x[i][nv-1].w,2.)/(10.+2.*pow(Re_x[i][nv-1].w,2.));
		alpha_x[i][nv-1].n = pow(Re_x[i][nv-1].n,2.)/(10.+2.*pow(Re_x[i][nv-1].n,2.));
		alpha_x[i][nv-1].s = pow(Re_x[i][nv-1].s,2.)/(10.+2.*pow(Re_x[i][nv-1].s,2.));
		beta_x[i][nv-1].e = (1.+(0.005*pow(Re_x[i][nv-1].e,2.)))/(1.+(0.05*pow(Re_x[i][nv-1].e,2.)));
		beta_x[i][nv-1].w = (1.+(0.005*pow(Re_x[i][nv-1].w,2.)))/(1.+(0.05*pow(Re_x[i][nv-1].w,2.)));
		beta_x[i][nv-1].n = (1.+(0.005*pow(Re_x[i][nv-1].n,2.)))/(1.+(0.05*pow(Re_x[i][nv-1].n,2.)));
		beta_x[i][nv-1].s = (1.+(0.005*pow(Re_x[i][nv-1].s,2.)))/(1.+(0.05*pow(Re_x[i][nv-1].s,2.)));

		if(Re_x[i][nv-1].e < 0) { alpha_x[i][nv-1].e = -alpha_x[i][nv-1].e; }
		if(Re_x[i][nv-1].w < 0) { alpha_x[i][nv-1].w = -alpha_x[i][nv-1].w; }
		if(Re_x[i][nv-1].n < 0) { alpha_x[i][nv-1].n = -alpha_x[i][nv-1].n; }
		if(Re_x[i][nv-1].s < 0) { alpha_x[i][nv-1].s = -alpha_x[i][nv-1].s; }

		Ae_u[i][nv-1] = (((mi*beta_x[i][nv-1].e*dy)/dx)-((rho*((u[i][nv-1]+u[i+1][nv-1])/2.)*dy)*(0.5-alpha_x[i][nv-1].e)));
		As_u[i][nv-1] = (((mi*beta_x[i][nv-1].s*dx)/dy)+((rho*((v[i+1][nv-2]+v[i][nv-2])/2.)*dx)*(0.5+alpha_x[i][nv-1].s)));
		Aw_u[i][nv-1] = (((mi*beta_x[i][nv-1].w*dy)/dx)+((rho*((u[i-1][nv-1]+u[i][nv-1])/2.)*dy)*(0.5+alpha_x[i][nv-1].w)));
		An_u[i][nv-1] = ((2.*mi*beta_x[i][nv-1].n*dx)/dy);
		Ap_u[i][nv-1] = Ae_u[i][nv-1] + Aw_u[i][nv-1] + An_u[i][nv-1] + As_u[i][nv-1];
		B_u[i][nv-1] = ((2.*mi*beta_x[i][nv-1].n*U*dx)/dy);
	}

	//Volumes da Fronteira Sul
	for(int i=1;i<(nv-2);i++)
	{
		Re_x[i][0].e = (rho*((u[i][0]+u[i+1][0])/2.)*dx)/mi;
		Re_x[i][0].n = (rho*((v[i][0]+v[i+1][0])/2.)*dx)/mi;
		Re_x[i][0].w = (rho*((u[i][0]+u[i-1][0])/2.)*dx)/mi;
		Re_x[i][0].s = 0.;
		alpha_x[i][0].e = pow(Re_x[i][0].e,2.)/(10.+2.*pow(Re_x[i][0].e,2.));
		alpha_x[i][0].w = pow(Re_x[i][0].w,2.)/(10.+2.*pow(Re_x[i][0].w,2.));
		alpha_x[i][0].n = pow(Re_x[i][0].n,2.)/(10.+2.*pow(Re_x[i][0].n,2.));
		alpha_x[i][0].s = pow(Re_x[i][0].s,2.)/(10.+2.*pow(Re_x[i][0].s,2.));
		beta_x[i][0].e = (1.+(0.005*pow(Re_x[i][0].e,2.)))/(1.+(0.05*pow(Re_x[i][0].e,2.)));
		beta_x[i][0].w = (1.+(0.005*pow(Re_x[i][0].w,2.)))/(1.+(0.05*pow(Re_x[i][0].w,2.)));
		beta_x[i][0].n = (1.+(0.005*pow(Re_x[i][0].n,2.)))/(1.+(0.05*pow(Re_x[i][0].n,2.)));
		beta_x[i][0].s = (1.+(0.005*pow(Re_x[i][0].s,2.)))/(1.+(0.05*pow(Re_x[i][0].s,2.)));

		if(Re_x[i][0].e < 0) { alpha_x[i][0].e = -alpha_x[i][0].e; }
		if(Re_x[i][0].w < 0) { alpha_x[i][0].w = -alpha_x[i][0].w; }
		if(Re_x[i][0].n < 0) { alpha_x[i][0].n = -alpha_x[i][0].n; }
		if(Re_x[i][0].s < 0) { alpha_x[i][0].s = -alpha_x[i][0].s; }
		
		Ae_u[i][0] = (((mi*beta_x[i][0].e*dy)/dx)-((rho*((u[i][0]+u[i+1][0])/2.)*dy)*(0.5-alpha_x[i][0].e)));
		An_u[i][0] = (((mi*beta_x[i][0].n*dx)/dy)-((rho*((v[i+1][0]+v[i][0])/2.)*dx)*(0.5-alpha_x[i][0].n)));
		Aw_u[i][0] = (((mi*beta_x[i][0].w*dy)/dx)+((rho*((u[i-1][0]+u[i][0])/2.)*dy)*(0.5+alpha_x[i][0].w)));
		As_u[i][0] = ((2.*mi*beta_x[i][0].s*dx)/dy);
		Ap_u[i][0] = Ae_u[i][0] + Aw_u[i][0] + An_u[i][0] + As_u[i][0];
		B_u[i][0] = 0.;
	}
	
	//Volumes da Fronteira Leste
	for(int j=1;j<(nv-1);j++)
	{
		Re_x[nv-2][j].e = (rho*((u[nv-2][j])/2.)*dx)/mi;
		Re_x[nv-2][j].n = (rho*((v[nv-2][j]+v[nv-1][j])/2.)*dx)/mi;
		Re_x[nv-2][j].w = (rho*((u[nv-2][j]+u[nv-3][j])/2.)*dx)/mi;
		Re_x[nv-2][j].s = (rho*((v[nv-2][j-1]+v[nv-1][j-1])/2.)*dx)/mi;
		alpha_x[nv-2][j].e = pow(Re_x[nv-2][j].e,2.)/(10.+2.*pow(Re_x[nv-2][j].e,2.));
		alpha_x[nv-2][j].w = pow(Re_x[nv-2][j].w,2.)/(10.+2.*pow(Re_x[nv-2][j].w,2.));
		alpha_x[nv-2][j].n = pow(Re_x[nv-2][j].n,2.)/(10.+2.*pow(Re_x[nv-2][j].n,2.));
		alpha_x[nv-2][j].s = pow(Re_x[nv-2][j].s,2.)/(10.+2.*pow(Re_x[nv-2][j].s,2.));
		beta_x[nv-2][j].e = (1.+(0.005*pow(Re_x[nv-2][j].e,2.)))/(1.+(0.05*pow(Re_x[nv-2][j].e,2.)));
		beta_x[nv-2][j].w = (1.+(0.005*pow(Re_x[nv-2][j].w,2.)))/(1.+(0.05*pow(Re_x[nv-2][j].w,2.)));
		beta_x[nv-2][j].n = (1.+(0.005*pow(Re_x[nv-2][j].n,2.)))/(1.+(0.05*pow(Re_x[nv-2][j].n,2.)));
		beta_x[nv-2][j].s = (1.+(0.005*pow(Re_x[nv-2][j].s,2.)))/(1.+(0.05*pow(Re_x[nv-2][j].s,2.)));

		if(Re_x[nv-2][j].e < 0) { alpha_x[nv-2][j].e = -alpha_x[nv-2][j].e; }
		if(Re_x[nv-2][j].w < 0) { alpha_x[nv-2][j].w = -alpha_x[nv-2][j].w; }
		if(Re_x[nv-2][j].n < 0) { alpha_x[nv-2][j].n = -alpha_x[nv-2][j].n; }
		if(Re_x[nv-2][j].s < 0) { alpha_x[nv-2][j].s = -alpha_x[nv-2][j].s; }

		As_u[nv-2][j] = (((mi*beta_x[nv-2][j].s*dx)/dy)+((rho*((v[nv-2][j-1]+v[nv-1][j-1])/2.)*dx)*(0.5+alpha_x[nv-2][j].s)));
		An_u[nv-2][j] = (((mi*beta_x[nv-2][j].n*dx)/dy)-((rho*((v[nv-2][j]+v[nv-1][j])/2.)*dx)*(0.5-alpha_x[nv-2][j].n)));
		Aw_u[nv-2][j] = (((mi*beta_x[nv-2][j].w*dy)/dx)+((rho*((u[nv-2][j]+u[nv-3][j])/2.)*dy)*(0.5+alpha_x[nv-2][j].w)));
		Ae_u[nv-2][j] = ((mi*beta_x[nv-2][j].e*dx)/dy)-((rho*((u[nv-2][j]/2.))*dy)*(0.5-alpha_x[nv-2][j].e));
		Ap_u[nv-2][j] = Ae_u[nv-2][j] + Aw_u[nv-2][j] + An_u[nv-2][j] + As_u[nv-2][j];
		B_u[nv-2][j] = 0.;
	}
	
	//Volumes da Fronteira oeste
	for(int j=1;j<(nv-1);j++)
	{
		Re_x[0][j].e = (rho*((u[0][j]+u[1][j])/2.)*dx)/mi;
		Re_x[0][j].n = (rho*((v[0][j]+v[1][j])/2.)*dx)/mi;
		Re_x[0][j].w = (rho*((u[0][j])/2.)*dx)/mi;
		Re_x[0][j].s = (rho*((v[0][j-1]+v[1][j-1])/2.)*dx)/mi;
		alpha_x[0][j].e = pow(Re_x[0][j].e,2.)/(10.+2.*pow(Re_x[0][j].e,2.));
		alpha_x[0][j].w = pow(Re_x[0][j].w,2.)/(10.+2.*pow(Re_x[0][j].w,2.));
		alpha_x[0][j].n = pow(Re_x[0][j].n,2.)/(10.+2.*pow(Re_x[0][j].n,2.));
		alpha_x[0][j].s = pow(Re_x[0][j].s,2.)/(10.+2.*pow(Re_x[0][j].s,2.));
		beta_x[0][j].e = (1.+(0.005*pow(Re_x[0][j].e,2.)))/(1.+(0.05*pow(Re_x[0][j].e,2.)));
		beta_x[0][j].w = (1.+(0.005*pow(Re_x[0][j].w,2.)))/(1.+(0.05*pow(Re_x[0][j].w,2.)));
		beta_x[0][j].n = (1.+(0.005*pow(Re_x[0][j].n,2.)))/(1.+(0.05*pow(Re_x[0][j].n,2.)));
		beta_x[0][j].s = (1.+(0.005*pow(Re_x[0][j].s,2.)))/(1.+(0.05*pow(Re_x[0][j].s,2.)));

		if(Re_x[0][j].e < 0) { alpha_x[0][j].e = -alpha_x[0][j].e; }
		if(Re_x[0][j].w < 0) { alpha_x[0][j].w = -alpha_x[0][j].w; }
		if(Re_x[0][j].n < 0) { alpha_x[0][j].n = -alpha_x[0][j].n; }
		if(Re_x[0][j].s < 0) { alpha_x[0][j].s = -alpha_x[0][j].s; }
		
		As_u[0][j] = (((mi*beta_x[0][j].s*dx)/dy)+((rho*((v[0][j-1]+v[1][j-1])/2.)*dx)*(0.5+alpha_x[0][j].s)));
		An_u[0][j] = (((mi*beta_x[0][j].n*dx)/dy)-((rho*((v[0][j]+v[1][j])/2.)*dx)*(0.5-alpha_x[0][j].n)));
		Ae_u[0][j] = (((mi*beta_x[0][j].e*dy)/dx)-((rho*((u[0][j]+u[1][j])/2.)*dy)*(0.5-alpha_x[0][j].e)));
		Aw_u[0][j] = (((mi*beta_x[0][j].w*dy)/dx)+((rho*((u[0][j]/2.))*dy)*(0.5+alpha_x[0][j].w)));
		Ap_u[0][j] = Ae_u[0][j] + Aw_u[0][j] + An_u[0][j] + As_u[0][j];
		B_u[0][j] = 0.;
	}
	//Volumes Centrais
	for(int i=1;i<(nv-2);i++)
	{
		for(int j=1;j<(nv-1);j++)
		{
			Re_x[i][j].e = (rho*((u[i][j]+u[i+1][j])/2.)*dx)/mi;
			Re_x[i][j].n = (rho*((v[i+1][j]+v[i][j])/2.)*dy)/mi;
			Re_x[i][j].w = (rho*((u[i][j]+u[i-1][j])/2.)*dx)/mi;
			Re_x[i][j].s = (rho*((v[i][j-1]+v[i+1][j-1])/2.)*dy)/mi;
			alpha_x[i][j].e = pow(Re_x[i][j].e,2.)/(10.+2.*pow(Re_x[i][j].e,2.));
			alpha_x[i][j].w = pow(Re_x[i][j].w,2.)/(10.+2.*pow(Re_x[i][j].w,2.));
			alpha_x[i][j].n = pow(Re_x[i][j].n,2.)/(10.+2.*pow(Re_x[i][j].n,2.));
			alpha_x[i][j].s = pow(Re_x[i][j].s,2.)/(10.+2.*pow(Re_x[i][j].s,2.));
			beta_x[i][j].e = (1.+(0.005*pow(Re_x[i][j].e,2.)))/(1.+(0.05*pow(Re_x[i][j].e,2.)));
			beta_x[i][j].w = (1.+(0.005*pow(Re_x[i][j].w,2.)))/(1.+(0.05*pow(Re_x[i][j].w,2.)));
			beta_x[i][j].n = (1.+(0.005*pow(Re_x[i][j].n,2.)))/(1.+(0.05*pow(Re_x[i][j].n,2.)));
			beta_x[i][j].s = (1.+(0.005*pow(Re_x[i][j].s,2.)))/(1.+(0.05*pow(Re_x[i][j].s,2.)));
			
			if(Re_x[i][j].e < 0) { alpha_x[i][j].e = -alpha_x[i][j].e; }
			if(Re_x[i][j].w < 0) { alpha_x[i][j].w = -alpha_x[i][j].w; }
			if(Re_x[i][j].n < 0) { alpha_x[i][j].n = -alpha_x[i][j].n; }
			if(Re_x[i][j].s < 0) { alpha_x[i][j].s = -alpha_x[i][j].s; }
			
			As_u[i][j] = (((mi*beta_x[i][j].s*dx)/dy)+((rho*((v[i][j-1]+v[i+1][j-1])/2.)*dx)*(0.5+alpha_x[i][j].s)));
			An_u[i][j] = (((mi*beta_x[i][j].n*dx)/dy)-((rho*((v[i][j]+v[i+1][j])/2.)*dx)*(0.5-alpha_x[i][j].n)));
			Ae_u[i][j] = (((mi*beta_x[i][j].e*dy)/dx)-((rho*((u[i][j]+u[i+1][j])/2.)*dy)*(0.5-alpha_x[i][j].e)));
			Aw_u[i][j] = (((mi*beta_x[i][j].w*dy)/dx)+((rho*((u[i-1][j]+u[i][j])/2.)*dy)*(0.5+alpha_x[i][j].w)));
			Ap_u[i][j] = Ae_u[i][j] + Aw_u[i][j] + An_u[i][j] + As_u[i][j];
			B_u[i][j] = 0.;
		}
	}

}

void Calc_Coef_NS_Y()
{
	// Volume do Canto inferior Esquerdo
	Re_y[0][0].e = (rho*((u[0][0]+u[0][1])/2.)*dy)/mi;
	Re_y[0][0].n = (rho*((v[0][1]+v[0][0])/2.)*dy)/mi;
	Re_y[0][0].w = 0.;
	Re_y[0][0].s = (rho*((v[0][0])/2.)*dy)/mi;
	alpha_y[0][0].e = pow(Re_y[0][0].e,2.)/(10.+2.*pow(Re_y[0][0].e,2.));
	alpha_y[0][0].w = pow(Re_y[0][0].w,2.)/(10.+2.*pow(Re_y[0][0].w,2.));
	alpha_y[0][0].n = pow(Re_y[0][0].n,2.)/(10.+2.*pow(Re_y[0][0].n,2.));
	alpha_y[0][0].s = pow(Re_y[0][0].s,2.)/(10.+2.*pow(Re_y[0][0].s,2.));
	beta_y[0][0].e = (1.+(0.005*pow(Re_y[0][0].e,2.)))/(1.+(0.05*pow(Re_y[0][0].e,2.)));
	beta_y[0][0].w = (1.+(0.005*pow(Re_y[0][0].w,2.)))/(1.+(0.05*pow(Re_y[0][0].w,2.)));
	beta_y[0][0].n = (1.+(0.005*pow(Re_y[0][0].n,2.)))/(1.+(0.05*pow(Re_y[0][0].n,2.)));
	beta_y[0][0].s = (1.+(0.005*pow(Re_y[0][0].s,2.)))/(1.+(0.05*pow(Re_y[0][0].s,2.)));
		
	if(Re_y[0][0].e < 0) { alpha_y[0][0].e = -alpha_y[0][0].e; }
	if(Re_y[0][0].w < 0) { alpha_y[0][0].w = -alpha_y[0][0].w; }
	if(Re_y[0][0].n < 0) { alpha_y[0][0].n = -alpha_y[0][0].n; }
	if(Re_y[0][0].s < 0) { alpha_y[0][0].s = -alpha_y[0][0].s; }
	
	As_v[0][0] = (((mi*beta_y[0][0].s*dx)/dy)+((rho*((v[0][0])/2.)*dx)*(0.5+alpha_y[0][0].s)));
	An_v[0][0] = (((mi*beta_y[0][0].n*dx)/dy)-((rho*((v[0][1]+v[0][0])/2.)*dx)*(0.5-alpha_y[0][0].n)));
	Ae_v[0][0] = (((mi*beta_y[0][0].e*dy)/dx)-((rho*((u[0][1]+u[0][0])/2.)*dy)*(0.5-alpha_y[0][0].e)));
	Aw_v[0][0] = ((2.*mi*beta_y[0][0].w*dy)/dx);
	Ap_v[0][0] = Ae_v[0][0] + Aw_v[0][0] + An_v[0][0] + As_v[0][0];
	B_v[0][0] = 0.;

	// Volume do Canto Inferior Direito
	Re_y[nv-1][0].e = 0.;
	Re_y[nv-1][0].n = (rho*((v[nv-1][1]+v[nv-1][0])/2.)*dy)/mi;
	Re_y[nv-1][0].w = (rho*((u[nv-2][0]+u[nv-2][1])/2.)*dy)/mi;
	Re_y[nv-1][0].s = (rho*((v[nv-1][0])/2.)*dy)/mi;
	alpha_y[nv-1][0].e = pow(Re_y[nv-1][0].e,2.)/(10.+2.*pow(Re_y[nv-1][0].e,2.));
	alpha_y[nv-1][0].w = pow(Re_y[nv-1][0].w,2.)/(10.+2.*pow(Re_y[nv-1][0].w,2.));
	alpha_y[nv-1][0].n = pow(Re_y[nv-1][0].n,2.)/(10.+2.*pow(Re_y[nv-1][0].n,2.));
	alpha_y[nv-1][0].s = pow(Re_y[nv-1][0].s,2.)/(10.+2.*pow(Re_y[nv-1][0].s,2.));
	beta_y[nv-1][0].e = (1.+(0.005*pow(Re_y[nv-1][0].e,2.)))/(1.+(0.05*pow(Re_y[nv-1][0].e,2.)));
	beta_y[nv-1][0].w = (1.+(0.005*pow(Re_y[nv-1][0].w,2.)))/(1.+(0.05*pow(Re_y[nv-1][0].w,2.)));
	beta_y[nv-1][0].n = (1.+(0.005*pow(Re_y[nv-1][0].n,2.)))/(1.+(0.05*pow(Re_y[nv-1][0].n,2.)));
	beta_y[nv-1][0].s = (1.+(0.005*pow(Re_y[nv-1][0].s,2.)))/(1.+(0.05*pow(Re_y[nv-1][0].s,2.)));
			
	if(Re_y[nv-1][0].e < 0) { alpha_y[nv-1][0].e = -alpha_y[nv-1][0].e; }
	if(Re_y[nv-1][0].w < 0) { alpha_y[nv-1][0].w = -alpha_y[nv-1][0].w; }
	if(Re_y[nv-1][0].n < 0) { alpha_y[nv-1][0].n = -alpha_y[nv-1][0].n; }
	if(Re_y[nv-1][0].s < 0) { alpha_y[nv-1][0].s = -alpha_y[nv-1][0].s; }
			
	As_v[nv-1][0] = (((mi*beta_y[nv-1][0].s*dx)/dy)+((rho*((v[nv-1][0])/2.)*dx)*(0.5+alpha_y[nv-1][0].s)));
	An_v[nv-1][0] = (((mi*beta_y[nv-1][0].n*dx)/dy)-((rho*((v[nv-1][0]+v[nv-1][1])/2.)*dx)*(0.5-alpha_y[nv-1][0].n)));
	Ae_v[nv-1][0] = ((2.*mi*beta_y[nv-1][0].e*dy)/dx);
	Aw_v[nv-1][0] = (((mi*beta_y[nv-1][0].w*dy)/dx)+((rho*((u[nv-2][0]+u[nv-2][1])/2.)*dy)*(0.5+alpha_y[nv-1][0].w)));
	Ap_v[nv-1][0] = Ae_v[nv-1][0] + Aw_v[nv-1][0] + An_v[nv-1][0] + As_v[nv-1][0];
	B_v[nv-1][0] = 0.;

	//Volume do Canto Superior Esquerdo
	Re_y[0][nv-2].e = (rho*((u[0][nv-2]+u[0][nv-1])/2.)*dy)/mi;
	Re_y[0][nv-2].n = (rho*((v[0][nv-2])/2.)*dy)/mi;
	Re_y[0][nv-2].w = 0.;
	Re_y[0][nv-2].s = (rho*((v[0][nv-3]+v[0][nv-2])/2.)*dy)/mi;
	alpha_y[0][nv-2].e = pow(Re_y[0][nv-2].e,2.)/(10.+2.*pow(Re_y[0][nv-2].e,2.));
	alpha_y[0][nv-2].w = pow(Re_y[0][nv-2].w,2.)/(10.+2.*pow(Re_y[0][nv-2].w,2.));
	alpha_y[0][nv-2].n = pow(Re_y[0][nv-2].n,2.)/(10.+2.*pow(Re_y[0][nv-2].n,2.));
	alpha_y[0][nv-2].s = pow(Re_y[0][nv-2].s,2.)/(10.+2.*pow(Re_y[0][nv-2].s,2.));
	beta_y[0][nv-2].e = (1.+(0.005*pow(Re_y[0][nv-2].e,2.)))/(1.+(0.05*pow(Re_y[0][nv-2].e,2.)));
	beta_y[0][nv-2].w = (1.+(0.005*pow(Re_y[0][nv-2].w,2.)))/(1.+(0.05*pow(Re_y[0][nv-2].w,2.)));
	beta_y[0][nv-2].n = (1.+(0.005*pow(Re_y[0][nv-2].n,2.)))/(1.+(0.05*pow(Re_y[0][nv-2].n,2.)));
	beta_y[0][nv-2].s = (1.+(0.005*pow(Re_y[0][nv-2].s,2.)))/(1.+(0.05*pow(Re_y[0][nv-2].s,2.)));
			
	if(Re_y[0][nv-2].e < 0) { alpha_y[0][nv-2].e = -alpha_y[0][nv-2].e; }
	if(Re_y[0][nv-2].w < 0) { alpha_y[0][nv-2].w = -alpha_y[0][nv-2].w; }
	if(Re_y[0][nv-2].n < 0) { alpha_y[0][nv-2].n = -alpha_y[0][nv-2].n; }
	if(Re_y[0][nv-2].s < 0) { alpha_y[0][nv-2].s = -alpha_y[0][nv-2].s; }

	As_v[0][nv-2] = (((mi*beta_y[0][nv-2].s*dx)/dy)+((rho*((v[0][nv-3]+v[0][nv-2])/2.)*dx)*(0.5+alpha_y[0][nv-2].s)));
	An_v[0][nv-2] = (((mi*beta_y[0][nv-2].n*dx)/dy)-((rho*((v[0][nv-2])/2.)*dx)*(0.5-alpha_y[0][nv-2].n)));
	Ae_v[0][nv-2] = (((mi*beta_y[0][nv-2].e*dy)/dx)-((rho*((u[0][nv-1]+u[0][nv-2])/2.)*dy)*(0.5-alpha_y[0][nv-2].e)));
	Aw_v[0][nv-2] = ((2.*mi*beta_y[0][nv-2].w*dy)/dx);
	Ap_v[0][nv-2] = Ae_v[0][nv-2] + Aw_v[0][nv-2] + An_v[0][nv-2] + As_v[0][nv-2];
	B_v[0][nv-2] = 0.;

	//Volume do Canto Superior Direito
	Re_y[nv-1][nv-2].e = 0.;
	Re_y[nv-1][nv-2].n = (rho*((v[nv-1][nv-2])/2.)*dy)/mi;
	Re_y[nv-1][nv-2].w = (rho*((u[nv-2][nv-2]+u[nv-2][nv-1])/2.)*dy)/mi;
	Re_y[nv-1][nv-2].s = (rho*((v[nv-1][nv-3]+v[nv-1][nv-2])/2.)*dy)/mi;
	alpha_y[nv-1][nv-2].e = pow(Re_y[nv-1][nv-2].e,2.)/(10.+2.*pow(Re_y[nv-1][nv-2].e,2.));
	alpha_y[nv-1][nv-2].w = pow(Re_y[nv-1][nv-2].w,2.)/(10.+2.*pow(Re_y[nv-1][nv-2].w,2.));
	alpha_y[nv-1][nv-2].n = pow(Re_y[nv-1][nv-2].n,2.)/(10.+2.*pow(Re_y[nv-1][nv-2].n,2.));
	alpha_y[nv-1][nv-2].s = pow(Re_y[nv-1][nv-2].s,2.)/(10.+2.*pow(Re_y[nv-1][nv-2].s,2.));
	beta_y[nv-1][nv-2].e = (1.+(0.005*pow(Re_y[nv-1][nv-2].e,2.)))/(1.+(0.05*pow(Re_y[nv-1][nv-2].e,2.)));
	beta_y[nv-1][nv-2].w = (1.+(0.005*pow(Re_y[nv-1][nv-2].w,2.)))/(1.+(0.05*pow(Re_y[nv-1][nv-2].w,2.)));
	beta_y[nv-1][nv-2].n = (1.+(0.005*pow(Re_y[nv-1][nv-2].n,2.)))/(1.+(0.05*pow(Re_y[nv-1][nv-2].n,2.)));
	beta_y[nv-1][nv-2].s = (1.+(0.005*pow(Re_y[nv-1][nv-2].s,2.)))/(1.+(0.05*pow(Re_y[nv-1][nv-2].s,2.)));
			
	if(Re_y[nv-1][nv-2].e < 0) { alpha_y[nv-1][nv-2].e = -alpha_y[nv-1][nv-2].e; }
	if(Re_y[nv-1][nv-2].w < 0) { alpha_y[nv-1][nv-2].w = -alpha_y[nv-1][nv-2].w; }
	if(Re_y[nv-1][nv-2].n < 0) { alpha_y[nv-1][nv-2].n = -alpha_y[nv-1][nv-2].n; }
	if(Re_y[nv-1][nv-2].s < 0) { alpha_y[nv-1][nv-2].s = -alpha_y[nv-1][nv-2].s; }
			
	As_v[nv-1][nv-2] = (((mi*beta_y[nv-1][nv-2].s*dx)/dy)+((rho*((v[nv-1][nv-3]+v[nv-1][nv-2])/2.)*dx)*(0.5+alpha_y[nv-1][nv-2].s)));
	An_v[nv-1][nv-2] = (((mi*beta_y[nv-1][nv-2].n*dx)/dy)-((rho*((v[nv-1][nv-2])/2.)*dx)*(0.5-alpha_y[nv-1][nv-2].n)));
	Ae_v[nv-1][nv-2] = ((2.*mi*beta_y[nv-1][nv-2].e*dy)/dx);
	Aw_v[nv-1][nv-2] = (((mi*beta_y[nv-1][nv-2].w*dy)/dx)+((rho*((u[nv-2][nv-1]+u[nv-2][nv-2])/2)*dy)*(0.5+alpha_y[nv-1][nv-2].w)));
	Ap_v[nv-1][nv-2] = Ae_v[nv-1][nv-2] + Aw_v[nv-1][nv-2] + An_v[nv-1][nv-2] + As_v[nv-1][nv-2];
	B_v[nv-1][nv-2] = 0.;

	//Volumes da Fronteira Norte
	for(int i=1;i<(nv-1);i++)
	{
		Re_y[i][nv-2].e = (rho*((u[i][nv-2]+u[i][nv-1])/2.)*dy)/mi;
		Re_y[i][nv-2].n = (rho*((v[i][nv-2])/2.)*dy)/mi;
		Re_y[i][nv-2].w = (rho*((u[i-1][nv-1]+u[i-1][nv-2])/2.)*dy)/mi;
		Re_y[i][nv-2].s = (rho*((v[i][nv-3]+v[i][nv-2])/2.)*dy)/mi;
		alpha_y[i][nv-2].e = pow(Re_y[i][nv-2].e,2.)/(10.+2.*pow(Re_y[i][nv-2].e,2.));
		alpha_y[i][nv-2].w = pow(Re_y[i][nv-2].w,2.)/(10.+2.*pow(Re_y[i][nv-2].w,2.));
		alpha_y[i][nv-2].n = pow(Re_y[i][nv-2].n,2.)/(10.+2.*pow(Re_y[i][nv-2].n,2.));
		alpha_y[i][nv-2].s = pow(Re_y[i][nv-2].s,2.)/(10.+2.*pow(Re_y[i][nv-2].s,2.));
		beta_y[i][nv-2].e = (1.+(0.005*pow(Re_y[i][nv-2].e,2.)))/(1.+(0.05*pow(Re_y[i][nv-2].e,2.)));
		beta_y[i][nv-2].w = (1.+(0.005*pow(Re_y[i][nv-2].w,2.)))/(1.+(0.05*pow(Re_y[i][nv-2].w,2.)));
		beta_y[i][nv-2].n = (1.+(0.005*pow(Re_y[i][nv-2].n,2.)))/(1.+(0.05*pow(Re_y[i][nv-2].n,2.)));
		beta_y[i][nv-2].s = (1.+(0.005*pow(Re_y[i][nv-2].s,2.)))/(1.+(0.05*pow(Re_y[i][nv-2].s,2.)));
			
		if(Re_y[i][nv-2].e < 0) { alpha_y[i][nv-2].e = -alpha_y[i][nv-2].e; }
		if(Re_y[i][nv-2].w < 0) { alpha_y[i][nv-2].w = -alpha_y[i][nv-2].w; }
		if(Re_y[i][nv-2].n < 0) { alpha_y[i][nv-2].n = -alpha_y[i][nv-2].n; }
		if(Re_y[i][nv-2].s < 0) { alpha_y[i][nv-2].s = -alpha_y[i][nv-2].s; }
			
		As_v[i][nv-2] = (((mi*beta_y[i][nv-2].s*dx)/dy)+((rho*((v[i][nv-3]+v[i][nv-2])/2.)*dx)*(0.5+alpha_y[i][nv-2].s)));
		An_v[i][nv-2] = (((mi*beta_y[i][nv-2].n*dx)/dy)-((rho*((v[i][nv-2])/2.)*dx)*(0.5-alpha_y[i][nv-2].n)));
		Ae_v[i][nv-2] = (((mi*beta_y[i][nv-2].e*dy)/dx)-((rho*((u[i][nv-1]+u[i][nv-2])/2.)*dy)*(0.5-alpha_y[i][nv-2].e)));
		Aw_v[i][nv-2] = (((mi*beta_y[i][nv-2].w*dy)/dx)+((rho*((u[i-1][nv-1]+u[i-1][nv-2])/2.)*dy)*(0.5+alpha_y[i][nv-2].w)));
		Ap_v[i][nv-2] = Ae_v[i][nv-2] + Aw_v[i][nv-2] + An_v[i][nv-2] + As_v[i][nv-2];
		B_v[i][nv-2] = 0.;
	}
	//Volumes da Fronteira Sul
	for(int i=1;i<(nv-1);i++)
	{
		Re_y[i][0].e = (rho*((u[i][0]+u[i][1])/2.)*dy)/mi;
		Re_y[i][0].n = (rho*((v[i][0]+v[i][1])/2.)*dy)/mi;
		Re_y[i][0].w = (rho*((u[i-1][0]+u[i-1][1])/2.)*dy)/mi;
		Re_y[i][0].s = (rho*((v[i][0])/2.)*dy)/mi;
		alpha_y[i][0].e = pow(Re_y[i][0].e,2.)/(10.+2.*pow(Re_y[i][0].e,2.));
		alpha_y[i][0].w = pow(Re_y[i][0].w,2.)/(10.+2.*pow(Re_y[i][0].w,2.));
		alpha_y[i][0].n = pow(Re_y[i][0].n,2.)/(10.+2.*pow(Re_y[i][0].n,2.));
		alpha_y[i][0].s = pow(Re_y[i][0].s,2.)/(10.+2.*pow(Re_y[i][0].s,2.));
		beta_y[i][0].e = (1.+(0.005*pow(Re_y[i][0].e,2.)))/(1.+(0.05*pow(Re_y[i][0].e,2.)));
		beta_y[i][0].w = (1.+(0.005*pow(Re_y[i][0].w,2.)))/(1.+(0.05*pow(Re_y[i][0].w,2.)));
		beta_y[i][0].n = (1.+(0.005*pow(Re_y[i][0].n,2.)))/(1.+(0.05*pow(Re_y[i][0].n,2.)));
		beta_y[i][0].s = (1.+(0.005*pow(Re_y[i][0].s,2.)))/(1.+(0.05*pow(Re_y[i][0].s,2.)));
		
		if(Re_y[i][0].e < 0) { alpha_y[i][0].e = -alpha_y[i][0].e; }
		if(Re_y[i][0].w < 0) { alpha_y[i][0].w = -alpha_y[i][0].w; }
		if(Re_y[i][0].n < 0) { alpha_y[i][0].n = -alpha_y[i][0].n; }
		if(Re_y[i][0].s < 0) { alpha_y[i][0].s = -alpha_y[i][0].s; }
				
		As_v[i][0] = (((mi*beta_y[i][0].s*dx)/dy)+((rho*((v[i][0])/2.)*dx)*(0.5+alpha_y[i][0].s)));
		An_v[i][0] = (((mi*beta_y[i][0].n*dx)/dy)-((rho*((v[i][1]+v[i][0])/2.)*dx)*(0.5-alpha_y[i][0].n)));
		Ae_v[i][0] = (((mi*beta_y[i][0].e*dy)/dx)-((rho*((u[i][1]+u[i][0])/2.)*dy)*(0.5-alpha_y[i][0].e)));
		Aw_v[i][0] = (((mi*beta_y[i][0].w*dy)/dx)+((rho*((u[i-1][1]+u[i-1][0])/2.)*dy)*(0.5+alpha_y[i][0].w)));
		Ap_v[i][0] = Ae_v[i][0] + Aw_v[i][0] + An_v[i][0] + As_v[i][0];
		B_v[i][0] = 0.;
	}
	//Volumes da Fronteira Leste
	for(int j=1;j<(nv-2);j++)
	{
		Re_y[nv-1][j].e = 0.;
		Re_y[nv-1][j].n = (rho*((v[nv-1][j+1]+v[nv-1][j])/2.)*dy)/mi;
		Re_y[nv-1][j].w = (rho*((u[nv-2][j+1]+u[nv-2][j])/2.)*dy)/mi;
		Re_y[nv-1][j].s = (rho*((v[nv-1][j-1]+v[nv-1][j])/2.)*dy)/mi;
		alpha_y[nv-1][j].e = pow(Re_y[nv-1][j].e,2.)/(10.+2.*pow(Re_y[nv-1][j].e,2.));
		alpha_y[nv-1][j].w = pow(Re_y[nv-1][j].w,2.)/(10.+2.*pow(Re_y[nv-1][j].w,2.));
		alpha_y[nv-1][j].n = pow(Re_y[nv-1][j].n,2.)/(10.+2.*pow(Re_y[nv-1][j].n,2.));
		alpha_y[nv-1][j].s = pow(Re_y[nv-1][j].s,2.)/(10.+2.*pow(Re_y[nv-1][j].s,2.));
		beta_y[nv-1][j].e = (1.+(0.005*pow(Re_y[nv-1][j].e,2.)))/(1.+(0.05*pow(Re_y[nv-1][j].e,2.)));
		beta_y[nv-1][j].w = (1.+(0.005*pow(Re_y[nv-1][j].w,2.)))/(1.+(0.05*pow(Re_y[nv-1][j].w,2.)));
		beta_y[nv-1][j].n = (1.+(0.005*pow(Re_y[nv-1][j].n,2.)))/(1.+(0.05*pow(Re_y[nv-1][j].n,2.)));
		beta_y[nv-1][j].s = (1.+(0.005*pow(Re_y[nv-1][j].s,2.)))/(1.+(0.05*pow(Re_y[nv-1][j].s,2.)));
			
		if(Re_y[nv-1][j].e < 0) { alpha_y[nv-1][j].e = -alpha_y[nv-1][j].e; }
		if(Re_y[nv-1][j].w < 0) { alpha_y[nv-1][j].w = -alpha_y[nv-1][j].w; }
		if(Re_y[nv-1][j].n < 0) { alpha_y[nv-1][j].n = -alpha_y[nv-1][j].n; }
		if(Re_y[nv-1][j].s < 0) { alpha_y[nv-1][j].s = -alpha_y[nv-1][j].s; }
		
		As_v[nv-1][j] = (((mi*beta_y[nv-1][j].s*dx)/dy)+((rho*((v[nv-1][j-1]+v[nv-1][j])/2.)*dx)*(0.5+alpha_y[nv-1][j].s)));
		An_v[nv-1][j] = (((mi*beta_y[nv-1][j].n*dx)/dy)-((rho*((v[nv-1][j+1]+v[nv-1][j])/2.)*dx)*(0.5-alpha_y[nv-1][j].n)));
		Ae_v[nv-1][j] = ((2.*mi*beta_y[nv-1][j].e*dy)/dx);
		Aw_v[nv-1][j] = (((mi*beta_y[nv-1][j].w*dy)/dx)+((rho*((u[nv-2][j+1]+u[nv-2][j])/2.)*dy)*(0.5+alpha_y[nv-1][j].w)));
		Ap_v[nv-1][j] = Ae_v[nv-1][j] + Aw_v[nv-1][j] + An_v[nv-1][j] + As_v[nv-1][j];
		B_v[nv-1][j] = 0.;
	}
	//Volumes da Fronteira Oeste
	for(int j=1;j<(nv-2);j++)
	{
		Re_y[0][j].e = (rho*((u[0][j]+u[0][j+1])/2.)*dy)/mi;
		Re_y[0][j].n = (rho*((v[0][j]+v[0][j+1])/2.)*dy)/mi;
		Re_y[0][j].w = 0.;
		Re_y[0][j].s = (rho*((v[0][j]+v[0][j-1])/2.)*dy)/mi;
		alpha_y[0][j].e = pow(Re_y[0][j].e,2.)/(10.+2.*pow(Re_y[0][j].e,2.));
		alpha_y[0][j].w = pow(Re_y[0][j].w,2.)/(10.+2.*pow(Re_y[0][j].w,2.));
		alpha_y[0][j].n = pow(Re_y[0][j].n,2.)/(10.+2.*pow(Re_y[0][j].n,2.));
		alpha_y[0][j].s = pow(Re_y[0][j].s,2.)/(10.+2.*pow(Re_y[0][j].s,2.));
		beta_y[0][j].e = (1.+(0.005*pow(Re_y[0][j].e,2.)))/(1.+(0.05*pow(Re_y[0][j].e,2.)));
		beta_y[0][j].w = (1.+(0.005*pow(Re_y[0][j].w,2.)))/(1.+(0.05*pow(Re_y[0][j].w,2.)));
		beta_y[0][j].n = (1.+(0.005*pow(Re_y[0][j].n,2.)))/(1.+(0.05*pow(Re_y[0][j].n,2.)));
		beta_y[0][j].s = (1.+(0.005*pow(Re_y[0][j].s,2.)))/(1.+(0.05*pow(Re_y[0][j].s,2.)));
		
		if(Re_y[0][j].e < 0) { alpha_y[0][j].e = -alpha_y[0][j].e; }
		if(Re_y[0][j].w < 0) { alpha_y[0][j].w = -alpha_y[0][j].w; }
		if(Re_y[0][j].n < 0) { alpha_y[0][j].n = -alpha_y[0][j].n; }
		if(Re_y[0][j].s < 0) { alpha_y[0][j].s = -alpha_y[0][j].s; }		
	
		As_v[0][j] = (((mi*beta_y[0][j].s*dx)/dy)+((rho*((v[0][j-1]+v[0][j])/2.)*dx)*(0.5+alpha_y[0][j].s)));
		An_v[0][j] = (((mi*beta_y[0][j].n*dx)/dy)-((rho*((v[0][j+1]+v[0][j])/2.)*dx)*(0.5-alpha_y[0][j].n)));
		Ae_v[0][j] = (((mi*beta_y[0][j].e*dy)/dx)-((rho*((u[0][j+1]+u[0][j])/2.)*dy)*(0.5-alpha_y[0][j].e)));
		Aw_v[0][j] = ((2.*mi*beta_y[0][j].w*dy)/dx);
		Ap_v[0][j] = Ae_v[0][j] + Aw_v[0][j] + An_v[0][j] + As_v[0][j];
		B_v[0][j] = 0.;
	}
	//Volumes Centrais
	for(int i=1;i<(nv-1);i++)
	{
		for(int j=1;j<(nv-2);j++)
		{
			Re_y[i][j].e = (rho*((u[i][j]+u[i][j+1])/2.)*dy)/mi;
			Re_y[i][j].n = (rho*((v[i][j+1]+v[i][j])/2.)*dy)/mi;
			Re_y[i][j].w = (rho*((u[i-1][j]+u[i-1][j+1])/2.)*dy)/mi;
			Re_y[i][j].s = (rho*((v[i][j-1]+v[i][j])/2.)*dy)/mi;
			alpha_y[i][j].e = pow(Re_y[i][j].e,2.)/(10.+2.*pow(Re_y[i][j].e,2.));
			alpha_y[i][j].w = pow(Re_y[i][j].w,2.)/(10.+2.*pow(Re_y[i][j].w,2.));
			alpha_y[i][j].n = pow(Re_y[i][j].n,2.)/(10.+2.*pow(Re_y[i][j].n,2.));
			alpha_y[i][j].s = pow(Re_y[i][j].s,2.)/(10.+2.*pow(Re_y[i][j].s,2.));
			beta_y[i][j].e = (1.+(0.005*pow(Re_y[i][j].e,2.)))/(1.+(0.05*pow(Re_y[i][j].e,2.)));
			beta_y[i][j].w = (1.+(0.005*pow(Re_y[i][j].w,2.)))/(1.+(0.05*pow(Re_y[i][j].w,2.)));
			beta_y[i][j].n = (1.+(0.005*pow(Re_y[i][j].n,2.)))/(1.+(0.05*pow(Re_y[i][j].n,2.)));
			beta_y[i][j].s = (1.+(0.005*pow(Re_y[i][j].s,2.)))/(1.+(0.05*pow(Re_y[i][j].s,2.)));
			
			if(Re_y[i][j].e < 0) { alpha_y[i][j].e = -alpha_y[i][j].e; }
			if(Re_y[i][j].w < 0) { alpha_y[i][j].w = -alpha_y[i][j].w; }
			if(Re_y[i][j].n < 0) { alpha_y[i][j].n = -alpha_y[i][j].n; }
			if(Re_y[i][j].s < 0) { alpha_y[i][j].s = -alpha_y[i][j].s; }
			
			As_v[i][j] = (((mi*beta_y[i][j].s*dx)/dy)+((rho*((v[i][j-1]+v[i][j])/2.)*dx)*(0.5+alpha_y[i][j].s)));
			An_v[i][j] = (((mi*beta_y[i][j].n*dx)/dy)-((rho*((v[i][j+1]+v[i][j])/2.)*dx)*(0.5-alpha_y[i][j].n)));
			Ae_v[i][j] = (((mi*beta_y[i][j].e*dy)/dx)-((rho*((u[i][j+1]+u[i][j])/2.)*dy)*(0.5-alpha_y[i][j].e)));
			Aw_v[i][j] = (((mi*beta_y[i][j].w*dy)/dx)+((rho*((u[i-1][j+1]+u[i-1][j])/2.)*dy)*(0.5+alpha_y[i][j].w)));
			Ap_v[i][j] = Ae_v[i][j] + Aw_v[i][j] + An_v[i][j] + As_v[i][j];
			B_v[i][j] = 0.;
		}
	}
}

void Calc_u_hat()
{
	// Canto Superior Esquerdo
	u_hat[0][nv-1] = ((Ae_u[0][nv-1]*u[1][nv-1])+(As_u[0][nv-1]*u[0][nv-2])+(B_u[0][nv-1]))/Ap_u[0][nv-1];
	// Canto Superior Direito
	u_hat[nv-2][nv-1] = ((Aw_u[nv-2][nv-1]*u[nv-3][nv-1])+(As_u[nv-2][nv-1]*u[nv-2][nv-2])+(B_u[nv-2][nv-1]))/Ap_u[nv-2][nv-1];
	// Canto Inferior Esquerdo
	u_hat[0][0] = ((Ae_u[0][0]*u[1][0])+(An_u[0][0]*u[0][1])+(B_u[0][0]))/Ap_u[0][0];
	// Canto Inferior Direito
	u_hat[nv-2][0] = ((Aw_u[nv-2][0]*u[nv-3][0])+(An_u[nv-2][0]*u[nv-2][1])+(B_u[nv-2][0]))/Ap_u[nv-2][0];
	// Fronteira Leste
	for(int j=1;j<(nv-1);j++)
	{
		u_hat[nv-2][j] = ((Aw_u[nv-2][j]*u[nv-3][j])+(An_u[nv-2][j]*u[nv-2][j+1])+(As_u[nv-2][j]*u[nv-2][j-1])+(B_u[nv-2][j]))/Ap_u[nv-2][j];
	}
	// Fronteira Oeste
	for(int j=1;j<(nv-1);j++)
	{
		u_hat[0][j] = ((Ae_u[0][j]*u[1][j])+(An_u[0][j]*u[0][j+1])+(As_u[0][j]*u[0][j-1])+(B_u[0][j]))/Ap_u[0][j];
	}
	// Fronteira Norte
	for(int i=1;i<(nv-2);i++)
	{
		u_hat[i][nv-1] = ((Ae_u[i][nv-1]*u[i+1][nv-1])+(Aw_u[i][nv-1]*u[i-1][nv-1])+(As_u[i][nv-1]*u[i][nv-2])+(B_u[i][nv-1]))/Ap_u[i][nv-1];
	}
	// Fronteira Sul
	for(int i=1;i<(nv-2);i++)
	{
		u_hat[i][0] = ((Ae_u[i][0]*u[i+1][0])+(Aw_u[i][0]*u[i-1][0])+(An_u[i][0]*u[i][1])+(B_u[i][0]))/Ap_u[i][0];
	}
	// Volumes Centrais
	for(int i=1;i<(nv-2);i++)
	{
		for(int j=1;j<(nv-1);j++)
		{
			u_hat[i][j] = ((Ae_u[i][j]*u[i+1][j])+(Aw_u[i][j]*u[i-1][j])+(An_u[i][j]*u[i][j+1])+(As_u[i][j]*u[i][j-1])+(B_u[i][j]))/Ap_u[i][j];
		}
	}
}

void Calc_v_hat()
{
	// Canto Superior Esquerdo
	v_hat[0][nv-2] = ((Ae_v[0][nv-2]*v[1][nv-2])+(As_v[0][nv-2]*v[0][nv-3])+(B_v[0][nv-2]))/Ap_v[0][nv-2];
	// Canto Superior Direito
	v_hat[nv-1][nv-2] = ((Aw_v[nv-1][nv-2]*v[nv-2][nv-2])+(As_v[nv-1][nv-2]*v[nv-1][nv-3])+(B_v[nv-1][nv-2]))/Ap_v[nv-1][nv-2];
	// Canto Inferior Esquerdo
	v_hat[0][0] = ((Ae_v[0][0]*v[1][0])+(An_v[0][0]*v[0][1])+(B_v[0][0]))/Ap_v[0][0];
	// Canto Inferior Direito
	v_hat[nv-1][0] = ((Aw_v[nv-1][0]*v[nv-2][0])+(An_v[nv-1][0]*v[nv-1][1])+(B_v[nv-1][0]))/Ap_v[nv-1][0];
	// Fronteira Leste
	for(int j=1;j<(nv-2);j++)
	{
		v_hat[nv-1][j] = ((Aw_v[nv-1][j]*v[nv-2][j])+(An_v[nv-1][j]*v[nv-1][j+1])+(As_v[nv-1][j]*v[nv-1][j-1])+(B_v[nv-1][j]))/Ap_v[nv-1][j];
	}
	// Fronteira Oeste
	for(int j=1;j<(nv-2);j++)
	{
		v_hat[0][j] = ((Ae_v[0][j]*v[1][j])+(An_v[0][j]*v[0][j+1])+(As_v[0][j]*v[0][j-1])+(B_v[0][j]))/Ap_v[0][j];
	}
	// Fronteira Norte
	for(int i=1;i<(nv-1);i++)
	{
		v_hat[i][nv-2] = ((Ae_v[i][nv-2]*v[i+1][nv-2])+(Aw_v[i][nv-2]*v[i-1][nv-2])+(As_v[i][nv-2]*v[i][nv-3])+(B_v[i][nv-2]))/Ap_v[i][nv-2];
	}
	// Fronteira Sul
	for(int i=1;i<(nv-1);i++)
	{
		v_hat[i][0] = ((Ae_v[i][0]*v[i+1][0])+(Aw_v[i][0]*v[i-1][0])+(An_v[i][0]*v[i][1])+(B_v[i][0]))/Ap_v[i][0];
	}
	// Volumes Centrais
	for(int i=1;i<(nv-1);i++)
	{
		for(int j=1;j<(nv-2);j++)
		{
			v_hat[i][j] = ((Ae_v[i][j]*v[i+1][j])+(Aw_v[i][j]*v[i-1][j])+(An_v[i][j]*v[i][j+1])+(As_v[i][j]*v[i][j-1])+(B_v[i][j]))/Ap_v[i][j];
		}
	}
}

void Calc_Coef_Pressao()
{
	// Canto Superior Esquerdo
	Ae_p[0][nv-1] = (rho*dy*dx*dy)/(dx*Ap_u[0][nv-1]);
	An_p[0][nv-1] = 0.;
	Aw_p[0][nv-1] = 0.;
	As_p[0][nv-1] = (rho*dx*dx*dy)/(dy*Ap_v[0][nv-2]);
	
	Ap_p[0][nv-1] = Ae_p[0][nv-1]+As_p[0][nv-1];
	B_p[0][nv-1] = -(rho*u_hat[0][nv-1]*dy)+(rho*v_hat[0][nv-2]*dx);
	// Canto Superior Direito
	Ae_p[nv-1][nv-1] = 0.;
	An_p[nv-1][nv-1] = 0.;
	Aw_p[nv-1][nv-1] = (rho*dx*dx*dy)/(dy*Ap_u[nv-2][nv-1]);
	As_p[nv-1][nv-1] = (rho*dx*dx*dy)/(dy*Ap_v[nv-1][nv-2]);
	
	Ap_p[nv-1][nv-1] = Aw_p[nv-1][nv-1]+As_p[nv-1][nv-1];
	B_p[nv-1][nv-1] = (rho*u_hat[nv-2][nv-1]*dy)+(rho*v_hat[nv-1][nv-2]*dx);
	// Canto Inferior Esquerdo
	Ae_p[0][0] = (rho*dy*dx*dy)/(dx*Ap_u[0][0]);
	An_p[0][0] = (rho*dx*dx*dy)/(dy*Ap_v[0][0]);
	Aw_p[0][0] = 0.;
	As_p[0][0] = 0.;
			
	Ap_p[0][0] = Ae_p[0][0]+An_p[0][0];
	B_p[0][0] = -(rho*u_hat[0][0]*dy)-(rho*v_hat[0][0]*dy);
	// Canto Inferior Direito
	Ae_p[nv-1][0] = 0.;
	An_p[nv-1][0] = (rho*dx*dx*dy)/(dy*Ap_v[nv-1][0]);
	Aw_p[nv-1][0] = (rho*dx*dx*dy)/(dy*Ap_u[nv-2][0]);
	As_p[nv-1][0] = 0.;
		
	Ap_p[nv-1][0] = Aw_p[nv-1][0]+An_p[nv-1][0];
	B_p[nv-1][0] = (rho*u_hat[nv-2][0]*dy)-(rho*v_hat[nv-1][0]*dy);
	//Fronteira Leste
	for(int j=1;j<(nv-1);j++)
	{
		Ae_p[nv-1][j] = 0.;
		An_p[nv-1][j] = (rho*dx*dx*dy)/(dy*Ap_v[nv-1][j]);
		Aw_p[nv-1][j] = (rho*dx*dx*dy)/(dy*Ap_u[nv-2][j]);
		As_p[nv-1][j] = (rho*dx*dx*dy)/(dy*Ap_v[nv-1][j-1]);
			
		Ap_p[nv-1][j] = Aw_p[nv-1][j]+An_p[nv-1][j]+As_p[nv-1][j];
		B_p[nv-1][j] = (rho*u_hat[nv-2][j]*dy)+(rho*v_hat[nv-1][j-1]*dx)-(rho*v_hat[nv-1][j]*dy);
	}
	//Fronteira Oeste
	for(int j=1;j<(nv-1);j++)
	{
		Ae_p[0][j] = (rho*dy*dx*dy)/(dx*Ap_u[0][j]);
		An_p[0][j] = (rho*dx*dx*dy)/(dy*Ap_v[0][j]);
		Aw_p[0][j] = 0.;
		As_p[0][j] = (rho*dx*dx*dy)/(dy*Ap_v[0][j-1]);
		
		Ap_p[0][j] = Ae_p[0][j]+An_p[0][j]+As_p[0][j];
		B_p[0][j] = -(rho*u_hat[0][j]*dy)+(rho*v_hat[0][j-1]*dx)-(rho*v_hat[0][j]*dx);
	}
	//Fronteira Norte
	for(int i=1;i<(nv-1);i++)
	{
		Ae_p[i][nv-1] = (rho*dy*dx*dy)/(dx*Ap_u[i][nv-1]);
		An_p[i][nv-1] = 0.;
		Aw_p[i][nv-1] = (rho*dx*dx*dy)/(dy*Ap_u[i-1][nv-1]);
		As_p[i][nv-1] = (rho*dx*dx*dy)/(dy*Ap_v[i][nv-2]);
		
		Ap_p[i][nv-1] = Ae_p[i][nv-1]+Aw_p[i][nv-1]+As_p[i][nv-1];
		B_p[i][nv-1] = (rho*u_hat[i-1][nv-1]*dy)-(rho*u_hat[i][nv-1]*dy)+(rho*v_hat[i][nv-2]*dx);
	}
	//Fronteira Sul
	for(int i=1;i<(nv-1);i++)
	{
		Ae_p[i][0] = (rho*dy*dx*dy)/(dx*Ap_u[i][0]);
		An_p[i][0] = (rho*dx*dx*dy)/(dy*Ap_v[i][0]);
		Aw_p[i][0] = (rho*dx*dx*dy)/(dy*Ap_u[i-1][0]);
		As_p[i][0] = 0.;
		
		Ap_p[i][0] = Ae_p[i][0]+Aw_p[i][0]+An_p[i][0];
		B_p[i][0] = (rho*u_hat[i-1][0]*dy)-(rho*u_hat[i][0]*dy)-(rho*v_hat[i][0]*dy);
	}
	//volumes Centrais
	for(int i=1;i<(nv-1);i++)
	{
		for(int j=1;j<(nv-1);j++)
		{
			Ae_p[i][j] = (rho*dy*dx*dy)/(dx*Ap_u[i][j]);
			An_p[i][j] = (rho*dx*dx*dy)/(dy*Ap_v[i][j]);
			Aw_p[i][j] = (rho*dx*dx*dy)/(dy*Ap_u[i-1][j]);
			As_p[i][j] = (rho*dx*dx*dy)/(dy*Ap_v[i][j-1]);
			
			Ap_p[i][j] = Ae_p[i][j]+Aw_p[i][j]+An_p[i][j]+As_p[i][j];
			B_p[i][j] = (rho*u_hat[i-1][j]*dy)-(rho*u_hat[i][j]*dy)+(rho*v_hat[i][j-1]*dx)-(rho*v_hat[i][j]*dy);
		}
	}
}

void Correcao_u_v()
{
	//Corrigindo Velocidade u
	for(int i=0;i<(nv-1);i++)
	{
		for(int j=0;j<nv;j++)
		{
			uOLD[i][j] = u[i][j];
			u[i][j]    = u_hat[i][j] - ((Pn[i+1][j]-Pn[i][j])*dx*dy)/(Ap_u[i][j]*dx);
		}
	}
	//Corrigindo Velocidade v
	for(int i=0;i<nv;i++)
	{
		for(int j=0;j<(nv-1);j++)
		{
			vOLD[i][j] = v[i][j];
			v[i][j]    = v_hat[i][j] - ((Pn[i][j+1]-Pn[i][j])*dy*dx)/(Ap_v[i][j]*dy);
		}
	}
}

double Erro_u()
{
	double erro = 0.;
	for(int i=0;i<(nv-1);i++)
	{
		for(int j=0;j<nv;j++)
		{
			erro += pow((u[i][j] - uOLD[i][j]),2.);
		}
	}
	erro = sqrt(erro);
	return erro;
}
double Erro_v()
{	
	double erro = 0.;
	for(int i=0;i<nv;i++)
	{
		for(int j=0;j<(nv-1);j++)
		{
			erro += pow((v[i][j] - vOLD[i][j]),2.);
		}
	}
	erro = sqrt(erro);
	return erro;	
}

void Destroy()  //Limpando Mem�ria
{
	delete [] u;
	delete [] v;
	delete [] u_hat;
	delete [] v_hat;
	delete [] P;

	delete [] Ap_u;
   	delete [] Aw_u;
	delete [] Ae_u;
	delete [] An_u;
	delete [] As_u;
	delete [] B_u;

	delete [] Ap_v;
   	delete [] Aw_v;
	delete [] Ae_v;
	delete [] An_v;
	delete [] As_v;
	delete [] B_v;
	
	delete [] Ap_p;
   	delete [] Aw_p;
	delete [] Ae_p;
	delete [] An_p;
	delete [] As_p;
	delete [] B_p;

	delete [] Re_x;
	delete [] alpha_x;
	delete [] beta_x;

	delete [] Re_y;
	delete [] alpha_y;
	delete [] beta_y;
}

void GravaArquivo()
{
	
	ofstream fout;
	fout << std::scientific;
	fout.open("outCav.txt");
	fout << "TITLE = \"OutCav\" " << endl;
	fout << "VARIABLES = \"X\", \"Y\", \"U\", \"V\", \"P\" " << endl;
	fout <<"ZONE T=\"" << "OUTCAV" << "\", N=" << (nv+1)*(nv+1) << ", E=" << nv*nv  <<  
  	",  DATAPACKING=BLOCK, ZONETYPE=FEQUADRILATERAL, VARLOCATION=([3-5]=CELLCENTERED), SOLUTIONTIME=" << "1" << endl;
  
	  // x coordinates
	  int k=0;
	  for (int j=1; j<(nv+2); j++)
	  {
	  	fout << 0. << "  ";
	  	k++;
		if( k == 9 ) { fout << endl; k=0; }
	  	for(int i=1; i < (nv+1); i++ )
	  	{
			fout << i*dx << "  ";
			k++;
			if( k == 9 ) { fout << endl; k=0; }
		}
	  }
	  fout << endl;
	  fout << endl;
	  
	  // y coordinates
	  k = 0;	  
	  for(int j=1; j < (nv+2); j++ )
	  {
		fout << 0. << "  ";
		k++;	 
		if( k == 9 ) { fout << endl; k=0; }
	  }
	  for(int j=1; j < (nv+1); j++ )
	  {
	          for(int i=1; i < (nv+2); i++ )
		  {
			  fout << j*dy << "  ";
			  k++;
			  if( k == 9 ) { fout << endl; k=0; }
		  }
	  }
	  fout << endl;
	  fout << endl;
	  
	  // Velocidade U
	  k=0;
          for(int j=0; j < (nv); j++ )
	  {
	  	for(int i=0; i < (nv-1); i++ )
	  	{
	
			fout << u[i][j]/U << "  ";
			k++;
	  		if( k == 9 ) { fout << endl; k=0;}
		}
		fout << "0." << "  ";
		k++;
		if( k == 9 ) { fout << endl; k=0; }
	  }
	  fout << endl;
	  fout << endl;
	  
	  // Velocidade V
	  k=0;
          for(int j=0; j < (nv-1); j++ )
	  {
		for(int i=0; i < (nv); i++ )
		{
			fout << v[i][j]/U << "  ";
			k++;
	  		if( k == 9 ) { fout << endl; k=0; }
		}
	  }
	  for(int i=0; i < (nv); i++ )
	  {
	  	fout << "0." << "  ";
		k++;
	  	if( k == 9 ) { fout << endl; k=0; }
	  }
	  fout << endl;
	  fout << endl;

	  // PRESS�O
	  k=0;
          for(int j=0; j < (nv); j++ )
	  {
		for(int i=0; i < (nv); i++ )
		{
			fout << Pn[i][j] << "  ";
			k++;
	  		if( k == 9 ) { fout << endl; k=0; }
		}
	  }
	  fout << endl;
	  fout << endl;

	  //Elementos
	  for(int i=1; i < (nv+1); i++ )
	  {
		 fout << i+(nv+1) << "\t" << i+(nv+2) << "\t" << i+1 << "\t" << i << endl;
	  }
	  for(int j=1; j < (nv); j++ )
	  {
		  for(int i=1; i < (nv+1); i++ )
		  {
			fout << i+(j*(nv+1))+(nv+1) << "\t" << i+(j*(nv+1))+(nv+2) << "\t" << i+(j*(nv+1))+1 << "\t" << i+(j*(nv+1))   << endl;
		  }
	  }

	  fout << endl;
	  fout.close();
} 

int main()
{
	
	read_in(); // Lendo Dados de Entrada
	
	// Calculando Delta X e Delta Y
	dx = L/nv; 
	dy = H/nv;
	
	Alocate(); // Alocando Vari�veis na Mem�ria
	int IT = 0;
	int TESTE =0;
	double eu = 0.;
	double ev = 0.;
	while ( TESTE == 0 )
	{
		IT++;
		cout << IT << " ";
		Calc_Coef_NS_X(); 
		Calc_Coef_NS_Y();
		Calc_u_hat();
		Calc_v_hat();
		Calc_Coef_Pressao();
		SOR_structured(Ap_p, Aw_p, Ae_p, An_p, As_p, 
			P, Pn, B_p, 
			nv, 50, 1.6
		);	
		Correcao_u_v();
		cout <<"erro-u:" << setw(7) << setprecision(5) << Erro_u() << " -v:" << setw(7) << setprecision(5) << Erro_u() << endl;
		if((IT % 500) == 0)
		{
			cout << endl << "......Gravando Solu��o Parcial....." << endl;
			GravaArquivo();
		}
		if( ((Erro_u() < (0.0001)) and ( Erro_v() < (0.0001) )) or (IT == 100000) )
		{
			TESTE = 1;
		}
	}
	IT++;
	cout << IT << " ";
	Calc_Coef_NS_X(); 
	Calc_Coef_NS_Y();
	Calc_u_hat();
	Calc_v_hat();
	Calc_Coef_Pressao();
	SOR_structured(Ap_p, Aw_p, Ae_p, An_p, As_p, 
			P, Pn, B_p, 
			nv, 1000, 1.6
		);	
	Correcao_u_v();
	cout <<"erro-u:" << setw(7) << setprecision(5) << Erro_u() << " -v:" << setw(7) << setprecision(5) << Erro_u() << endl;
	GravaArquivo();

	Destroy(); // Deletando Vari�veis da Mem�ria
}

