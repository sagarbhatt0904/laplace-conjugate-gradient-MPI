#include <iostream>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <cblas.h>
#include "contour.h"

using namespace std;

vector<double> Sparse_CRS_MVMult(int N, vector<int> &row, vector<int> &col, vector<double> &val, vector<double> &v)
{
	vector<double> result;
	result.resize(N);	
	for (int i = 0; i < N; ++i)
	{	
		result[i]=0;
		for (int k = row[i]; k < row[i+1]; ++k)
		{
			result[i]=result[i]+val[k]*v[col[k]];			
		}
	}	
	return result;
}

int main(int argc, char const *argv[])
{
	int N=100;	
	int N1=N-2;
	int N2=pow(N1,2);
	int k, k_m;    // dummy counter
	double dx=1/((double)N-1), dy=1/((double)N-1);
	double bet=dx/dy;
	double alpha, rho, rho_new;
	double alf=-2*(1+pow(bet,2));
	vector<double> x;    // Grid points
	vector<double> y;
	vector<vector<double> > u (N,vector<double>(N, 0)); // 2D u matrix
	vector<double> B(N2,0);								// RHS
	vector<double> U(N2,0);								// 1D u vector	for operations
	vector<double> s(N2,0);				// A orthogonal vector
	vector<double> r(N2,0);				// residual
	vector<double> z(N2,0);				//preconditioned residual
	vector<double> p(N2,0);
	vector<double> a(N2,0);				//a=A*p

	vector<double> val;			// value vector for A matrix
	vector<double> val_m;		// value vector for preconditioner
	vector<int> row;			//row counter for A matrix
	vector<int> col;			// column indices for A matrix
	vector<int> row_m;			// row countrer for preconditioner
	vector<int> col_m;			// column indices for preconditioner


	k=0;
	k_m=0;

	// Storing A in a Sparse format
	for (int i = 0; i < N2; ++i)
	{	
		row.push_back(k);
		row_m.push_back(k_m);
		for (int j = 0; j < N2; ++j)
		{
			if (i==j)
			{
				val.push_back(-4);
				col.push_back(j);
				k++;
				val_m.push_back(alf);
				col_m.push_back(j);
				k_m++;
			}
			else if (i==j+1)
			{
				if (i!=0 && i%(N1)==0)
				{
					val.push_back(0);
					col.push_back(j);
					k++;
				}
				else
				{
					val.push_back(1);
					col.push_back(j);
					k++;
				}
			}
			else if (i+1==j)
			{
				if (j!=0 && j%(N1)==0)
				{
					val.push_back(0);
					col.push_back(j);
					k++;
				}
				else
				{
					val.push_back(1);
					col.push_back(j);
					k++;
				}
			}
			else if (i==j+N1)
			{
				val.push_back(pow(bet,2));
				col.push_back(j);
				k++;
			}
			else if (i+N1==j)
			{
				val.push_back(pow(bet,2));
				col.push_back(j);
				k++;
			}		
		}			
	}
	row.push_back(k);
	row_m.push_back(k_m);
	
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

	// Computing first residual
	s=Sparse_CRS_MVMult(N2, row, col, val, s); // computing using the first guess as 0
	for (int i = 0; i < N2; ++i)
	{
		r[i]=B[i]-s[i];
	}
	
	z=Sparse_CRS_MVMult(N2, row_m, col_m, val_m, r);
	cblas_dcopy(N2,z.data(),1,p.data(),1);
	rho=cblas_ddot(N2,r.data(),1,z.data(),1);

	int niter = 0;     // Init counter for number of iterations
	int flag = 0;      // Init break flag

	double tol=pow(10,-16);			// Convergence criteria
	int maxiter=10000;

	while (cblas_dnrm2(N2,r.data(),1) > tol)   // Test break condition
	{
		a=Sparse_CRS_MVMult(N2, row, col, val, p);
		alpha = rho/cblas_ddot(N2,a.data(),1,p.data(),1);
		cblas_daxpy(N2,alpha, p.data(),1,U.data(),1);
		cblas_daxpy(N2,-alpha, a.data(),1,r.data(),1);
		z=Sparse_CRS_MVMult(N2, row_m, col_m, val_m, r);
		rho_new=cblas_ddot(N2,r.data(),1,z.data(),1);
		
		for (int i = 0; i < N2; ++i)
		{
			p[i]=z[i]+(rho_new/rho)*p[i];
		}
		rho=rho_new;
		niter+=1;
		if (niter == maxiter)         // if max. number of iterations  is reached, break.
	    {
	    	flag = 1;                   
	        cout<<"ERROR!!!! \n Max Iterations Reached!!!!\n";
	        break;
	    }  
	}
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
	char title[]="u_CG_sparse.vtk";
	contour(u,N,title);

   return 0;
}