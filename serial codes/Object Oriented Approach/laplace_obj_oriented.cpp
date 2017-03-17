#include <iostream>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <cblas.h>
#include "contour.h"
#include "CGKRYLOV.h"

using namespace std;

int main(int argc, char const *argv[])
{
	int N=100;
	int N1=N-2;
	int N2=pow(N1,2);
	int k, k_m;    // dummy counter
	double dx=1/((double)N-1), dy=1/((double)N-1);
	double bet=dx/dy;
	double alf=-2*(1+pow(bet,2));
	vector<double> x;    // Grid points
	vector<double> y;
	vector<vector<double> > u (N,vector<double>(N, 0)); // 2D u matrix
	vector<double> B(N2,0);								// RHS
	vector<double> U(N2,0);								// 1D u vector	for operations
	std::vector<std::vector<double> > A(N2,vector<double>(N2, 0));
	std::vector<std::vector<double> > A_m(N2,vector<double>(N2, 0));;


	// Storing A in a Sparse format
	for (int i = 0; i < N2; ++i)
	{	
		for (int j = 0; j < N2; ++j)
		{
			if (i==j)
			{
				A[i][j]=-4;
				A_m[i][j]=alf;
			}
			else if (i==j+1)
			{
				if (i!=0 && i%(N1)==0)
				{
					A[i][j]=0;
				}
				else
				{
					A[i][j]=1;
				}
			}
			else if (i+1==j)
			{
				if (j!=0 && j%(N1)==0)
				{
					A[i][j]=0;
				}
				else
				{
					A[i][j]=1;
				}
			}
			else if (i==j+N1)
			{
				A[i][j]=pow(bet,2);
			}
			else if (i+N1==j)
			{
				A[i][j]=pow(bet,2);
			}		
		}			
	}

	// Grid generation and Boundary conditions
	for (int i = 0; i < N; ++i)
	{
		x.push_back(i*dx);
		y.push_back(i*dx);
		u[i][0]=sin(M_PI*x[i]);
		u[i][N-1]=sin(M_PI*x[i])*exp(-M_PI);
		u[0][i]=0;
		u[N-1][i]=0;
	}
	
	// Initializing B, U 
	k=0;
	for (int i = 1; i < N-1; ++i)
	{
		for (int j = 1; j < N-1; ++j)
		{
			B[k]=-u[i+1][j]-u[i-1][j]-pow(bet,2)*u[i][j+1]-pow(bet,2)*u[i][j-1];			
			k+=1;

		}
	}
	/* Using Conjugate Gradient method to solve AU=B*/
	CGKRYLOV<double> solution(N2);
	U=solution.SolvePreCondCRS(A, A_m, B); 
  	// Updating values of u from computed results
    k=0;
    for (int i = 1; i < N-1; ++i)
    {
    	for (int j = 1; j < N-1; ++j)
	    {
	    	u[i][j]=U[k];
	    	k+=1; 
	    }
    }
		
	// Plotting
	char title[]="u_CG.vtk";
	char var[]="U";
	contour(u,N,title,var);

   return 0;
}