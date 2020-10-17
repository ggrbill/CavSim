#include <iostream> 
#include <math.h>

#include "solver.hpp"

void SOR_structured(
    double** Ap,
    double** Aw,
    double** Ae,
    double** An,
    double** As, 
    double** x,
    double** xn,
    double** b, 
    int const size,
    int const MAX_IT, 
    float const w, // relaxation factor
    double const tol // convergence tolerance
)
{
    std::cout << "Solve Pressure ";
    int N_ITE = 0;
    double T = 0.0;
    while(true) {
        N_ITE++;   
        T = 0.0;
        for(int i=0; i<size; i++)
        {
            for(int j=0; j<size; j++)
            {
                x[i][j]  = xn[i][j];
            }
        }

        // Left-Up corner
        xn[0][size-1] = b[0][size-1]; 
        xn[0][size-1] += (As[0][size-1]*xn[0][size-2]+ Ae[0][size-1]*xn[1][size-1]);
        xn[0][size-1] /= Ap[0][size-1];
        xn[0][size-1] = w*xn[0][size-1] + (1-w)*x[0][size-1];

        // Up boundary
        for(int i=1; i<=(size-2); i++)
        {
                xn[i][size-1] = b[i][size-1]; 
                xn[i][size-1] += (Aw[i][size-1]*xn[i-1][size-1] + Ae[i][size-1]*xn[i+1][size-1] + As[i][size-1]*xn[i][size-2]);
                xn[i][size-1] /= Ap[i][size-1];
                xn[i][size-1] = w*xn[i][size-1] + (1-w)*x[i][size-1];
        }

        // Right-Up corner
        xn[size-1][size-1] = b[size-1][size-1]; 
        xn[size-1][size-1] += (As[size-1][size-1]*xn[size-1][size-2]+ Aw[size-1][size-1]*xn[size-2][size-1]);
        xn[size-1][size-1] /= Ap[size-1][size-1];
        xn[size-1][size-1] = w*xn[size-1][size-1] + (1-w)*x[size-1][size-1];

        for(int j=(size-2);j>=1;j--)
        {
            // West boundary
            xn[0][j] = b[0][j]; 
            xn[0][j] += (As[0][j]*xn[0][j-1] + Ae[0][j]*xn[1][j] + An[0][j]*xn[0][j+1]);
            xn[0][j] /= Ap[0][j];
            xn[0][j] = w*xn[0][j] + (1-w)*x[0][j];
            // Core control volumes
            for(int i=1; i<=(size-2); i++)
            { 
                xn[i][j] = b[i][j]; 
                xn[i][j] += (As[i][j]*xn[i][j-1] + Aw[i][j]*xn[i-1][j] + Ae[i][j]*xn[i+1][j] + An[i][j]*xn[i][j+1]);
                xn[i][j] /= Ap[i][j];
                xn[i][j] = w*xn[i][j] + (1-w)*x[i][j];
            }
            // East boundary
            xn[size-1][j] = b[size-1][j]; 
            xn[size-1][j] += (As[size-1][j]*xn[size-1][j-1] + Aw[size-1][j]*xn[size-2][j] + An[size-1][j]*xn[size-1][j+1]);
            xn[size-1][j] /= Ap[size-1][j];
            xn[size-1][j] = w*xn[size-1][j] + (1-w)*x[size-1][j];
        }

        // Left-Bottom corner	
        xn[0][0] = b[0][0]; 
        xn[0][0] += (An[0][0]*xn[0][1]+ Ae[0][0]*xn[1][0]);
        xn[0][0] /= Ap[0][0];
        xn[0][0] = w*xn[0][0] + (1-w)*x[0][0];	
        // South boundary 
        for(int i=1; i<=(size-2); i++)
        {
                xn[i][0] = b[i][0]; 
                xn[i][0] += (Aw[i][0]*xn[i-1][0] + Ae[i][0]*xn[i+1][0] + An[i][0]*xn[i][1]);
                xn[i][0] /= Ap[i][0];
                xn[i][0] = w*xn[i][0] + (1-w)*x[i][0];
        }
        // Right-Bottom Corner
        xn[size-1][0] = b[size-1][0]; 
        xn[size-1][0] += (An[size-1][0]*xn[size-1][1]+ Aw[size-1][0]*xn[size-2][0]);
        xn[size-1][0] /= Ap[size-1][0];
        xn[size-1][0] = w*xn[size-1][0] + (1-w)*x[size-1][0];

        for(int i=0; i<size; i++)
        {
            for(int j=0; j<size; j++)
            {
                T += pow((x[i][j] - xn[i][j]),2.0);
            }
        }
        T = sqrt(T);
        if ( (T < (tol)) or (N_ITE == MAX_IT) )
        {
            break;
        }	
    }
    std::cout << N_ITE << " Iterations ";
}