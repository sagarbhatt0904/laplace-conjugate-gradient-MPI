#include <iostream>
#include <vector>
#include <math.h>
#include "contour.h"

using namespace std;

int main(int argc, char const *argv[])
{
	int N=1000;
	double dx=1/((double)N);
	vector<double> x;
	vector<double> y;
	vector<vector<double> > u (N,vector<double>(N, 0));
	
	for (int i = 0; i < N; ++i)
	{
		x.push_back(i*dx);
		y.push_back(i*dx);
	}
	for (int i = 0; i < N; ++i)
	{
		for (int j = 0; j < N; ++j)
		{
			u[i][j]=sin(M_PI*x[i])*exp(-M_PI*y[j]);
		}
	}
	
	// Plotting
	char title[]="u_analytical.vtk";
	contour(u,N,title);
	
	return 0;
}
