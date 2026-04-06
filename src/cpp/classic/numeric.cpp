#include "numeric.hpp"

#include <math.h>

double calculate_vec_diff_L2_norm(DoubleArray2D& v1, DoubleArray2D& v2, int nx, int ny)
{
    double norm = 0.;
	for(int i=0;i<nx;i++)
	{
		for(int j=0;j<ny;j++)
		{
			norm += pow((v1[i][j] - v2[i][j]),2.);
		}
	}
	norm = sqrt(norm);
	return norm;
}
