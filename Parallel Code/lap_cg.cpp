#include <iostream>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <mpi.h>
#include <mkl.h>
#include "contour.h"


using namespace std;

#define root 0

template <class T>

vector<T> Sparse_CRS_MVMult(int N, vector<int> &row, vector<int> &col, vector<T> &val, vector<T> &v)
{
	vector<T> result;
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
void proc_indices(int N, int np, int my_p, int count, int start, int end, int ind[3])
{
	count = N/np;
	int remainder= N-count*np;
	if (my_p<=remainder)
	{
		start=(my_p)*count;
    	end=(my_p+1)*count;	
	}
    else
    {
    	start=(my_p)*count + remainder;
    	end=(my_p+1)*count + remainder;	
    }

    
    ind[0]=start; ind[1]=end; ind[2]=count;
}


int main(int argc, char *argv[])
{
	int N=500;   // Grid points	
	int N1=N-2;
	int N2=pow(N1,2);
	double dx=1/((double)N-1), dy=1/((double)N-1);
	double bet=dx/dy;
	double alpha, alpha_g, rho, rho_new, rho_g, rho_new_g, ap, ap_g;
	double alf=-0.25;
	int k, k_m; // Dummy counter
	vector<double> x(N,0);    // Grid points
	vector<vector<double> > u (N,vector<double>(N, 0)); // 2D u matrix

	vector<double> val;			// value vector for A matrix
	vector<double> val_m;		// value vector for preconditioner
	vector<int> row;			//row counter for A matrix
	vector<int> col;			// column indices for A matrix
	vector<int> row_m;			// row countrer for preconditioner
	vector<int> col_m;			// column indices for preconditioner

	vector<double> r;				// residual
	vector<double> z;				//preconditioned residual
	vector<double> p;
	vector<double> a;				//a=A*p

	vector<double> B;
	vector<double> U;


	int my_p, np, ierr;
	ierr=MPI_Init(&argc, &argv);
	ierr=MPI_Comm_rank(MPI_COMM_WORLD, &my_p);
	ierr=MPI_Comm_size(MPI_COMM_WORLD, &np);

	/*MPI Success check*/
    if (ierr!= MPI_SUCCESS) 
    {
    	cout<<"Error starting MPI program.  Terminating. \n";
    	MPI_Abort(MPI_COMM_WORLD, ierr);
  	}

  
  	int start, end, count, ind[3];
  	proc_indices(N2, np, my_p, count, start, end, ind);
  	start=ind[0]; end=ind[1]; count=ind[2];
  

    int count1,start1,end1,ind1[3];
	proc_indices(N, np, my_p, count1, start1, end1,ind1);    
	start1=ind1[0]; end1=ind1[1]; count1=ind1[2];

    // Storing A in a Sparse format
    k=0;
    k_m=0;
	for (int i = start; i < end; ++i)
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



	for (int i = 0; i < N; ++i)
	{
		x[i]=i*dx;
		u[i][0]=sin(M_PI*x[i]);
		u[i][N-1]=sin(M_PI*x[i])*exp(-M_PI);
	}
	
	// Initializing B 
	
	for (int i = start1; i < end1; ++i)
	{
		for (int j = 1; j < N-1; ++j)
		{	
			if(i==N-1)
				break;
			else if (i==0)
				break;
			else
				B.push_back(-u[i+1][j]-u[i-1][j]-pow(bet,2)*u[i][j+1]-pow(bet,2)*u[i][j-1]);	

		}
	}
	
	
	// Computing first residual
	
	for (int i = 0; i < B.size(); ++i)
	{
		r.push_back(B[i]);
	}
	
	
	
	
	z=Sparse_CRS_MVMult(r.size(), start, end, row_m, col_m, val_m, r);
	
	p.resize(z.size());
	cblas_dcopy(z.size(),z.data(),1,p.data(),1);
	

	
	rho=cblas_ddot(r.size(),r.data(),1,z.data(),1);
	ierr=MPI_Allreduce(&rho, &rho_g, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	rho=rho_g;
	
	int niter = 0;     // Init counter for number of iterations

	double tol=pow(10,-16);			// Convergence criteria
	int maxiter=10000;
	double norm=1, norm_g;

	U.resize(p.size());
	while (norm>tol)   // Test break condition
	{		
		a=Sparse_CRS_MVMult(p.size(), start, end, row, col, val, p);
		ap=	cblas_ddot(a.size(),a.data(),1,p.data(),1);
		ierr=MPI_Allreduce(&ap, &ap_g, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
		ap=ap_g;
		alpha = rho/ap;
		ierr=MPI_Allreduce(&alpha, &alpha_g, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
		alpha=alpha_g;
		cblas_daxpy(p.size(),alpha, p.data(),1,U.data(),1); 
		cblas_daxpy(r.size(),-alpha, a.data(),1,r.data(),1); 
		z=Sparse_CRS_MVMult(r.size(), start, end, row_m, col_m, val_m, r);
		rho_new=cblas_ddot(r.size(),r.data(),1,z.data(),1);
		ierr=MPI_Allreduce(&rho_new, &rho_new_g, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
		rho_new=rho_new_g;
		for (int i = 0; i < p.size(); ++i)
		{
			p[i]=z[i]+(rho_new/rho)*p[i];
		}
		rho=rho_new;
		niter+=1;
		if (niter == maxiter)         // if max. number of iterations  is reached, break.
	    {
	        cout<<"ERROR!!!! \n Max Iterations Reached!!!!\n";
	        break;
	    }  
	    
	    norm=cblas_dnrm2(r.size(),r.data(),1);
	    ierr=MPI_Allreduce(&norm, &norm_g, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
		norm=norm_g;
		
	}
	
		
	int dm=a.size(); int dm_g;					
	ierr=MPI_Allreduce(&dm, &dm_g, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
	std::vector<double> U_g(dm_g);
	ierr=MPI_Gather(U.data(), U.size(), MPI_DOUBLE, U_g.data(), U.size(), MPI_DOUBLE, root, MPI_COMM_WORLD);
	

  	if (my_p==root)
  	{
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
		char title[]="u_CG_mpi.vtk";
		contour(u,N,title,"U");
  	}

  	ierr=MPI_Finalize();
  	return ierr;
 }