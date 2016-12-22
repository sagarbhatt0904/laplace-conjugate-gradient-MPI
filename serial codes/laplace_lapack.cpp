#include <iostream>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <lapacke.h>
#include "contour.h"


using namespace std;

int main(int argc, char const *argv[])
{
	int N=100, nrhs=1;
	int N1=N-2;
	int N2=pow(N1,2);
	int k;    // dummy counter
	int info; // info variable for dgesv
	double dx=1/((double)N-1), dy=1/((double)N-1);
	double conv=1,sum1=0,sum2=0, bet=dx/dy;
	double alf=-2*(1+pow(bet,2));
	vector<double> x;
	vector<double> y;

	vector<vector<double> > u (N,vector<double>(N, 0));
	vector<double> A(N2*N2);
	vector<double> B(N2);
	vector<int> U(N2);



	//Intializing  u and A
	for (int i = 0; i < N; ++i)
	{
		for (int j = 0; j < N; ++j)
		{ 
			u[i][j]=0;
		}
	}
	k=0;
	for (int i = 0; i < N2; ++i)
	{
		for (int j = 0; j < N2; ++j)
		{
			if (i==j)
			{
				A[k]=alf;
			}
			else if (i==j+1)
			{
				if (i!=0 && i%(N1)==0)
				{
					A[k]=0;
				}
				else
				{
					A[k]=1;
				}
			}
			else if (i+1==j)
			{
				if (j!=0 && j%(N1)==0)
				{
					A[k]=0;
				}
				else
				{
					A[k]=1;
				}
			}
			else if (i==j+N1)
			{
				A[k]=pow(bet,2);
			}
			else if (i+N1==j)
			{
				A[k]=pow(bet,2);
			}
			else
			{
				A[k]=0;
			}			
			k+=1;
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
	
	// Initializing B
	k=0;
	for (int i = 1; i < N-1; ++i)
	{
		for (int j = 1; j < N-1; ++j)
		{
			B[k]=-u[i+1][j]-u[i-1][j]-pow(bet,2)*u[i][j+1]-pow(bet,2)*u[i][j-1];
			k+=1;

		}
	}
	// Using lapack to solve AU=B
	
	dgesv_(&N2, &nrhs,A.data(), &N2, U.data(), B.data(), &N2, &info);

	if(info == 0)
	{
	  // Updating values of u from computed results
	    k=0;
	    for (int i = 1; i < N-1; ++i)
	    {
	    	for (int j = 1; j < N-1; ++j)
		    {
		    	u[i][j]=B[k];
		    	k+=1; 
		    }
	    }
			
		// Plotting
		char title[]="u_lapack.vtk";
		contour(u,N,title);
	}	
    else
   {
      // Write an error message
      std::cerr << "dgesv returned error " << info << "\n";
   }
   return 0;
}