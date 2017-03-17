/* 

	Using Conjugate Gradient method to solve AU=B
	Created on: Feb 20, 2017
    Author: Sagar Bhatt

*/


#ifndef CGKRYLOV_h
#define CGKRYLOV_h

#include "CRS.h"
#include <cblas.h>
#include <vector>

using namespace std;

template <typename CG>
class CGKRYLOV
{
	int N;							// Grid Size (Excluding boundary nodes) 
	double alpha, rho, rho_new;
	vector<double> U;				// 1D u vector	for operations
	vector<double> s;				// A orthogonal vector
	vector<double> r;				// residual
	vector<double> z;				//preconditioned residual
	vector<double> p;
	vector<double> a;				//a=A*p

public:
	CGKRYLOV()
	{
		cout<<"\n\n\nCG Krylov : Problem Size not given during object intialization\n\n\n";
	}
	CGKRYLOV(const int& N)
	{	
		this->N=N;
		U.resize(N,0);
		s.resize(N,0);
		r.resize(N,0);
		z.resize(N,0);
		p.resize(N,0);
		a.resize(N,0);
		cout<<"\n\n\n \t\t\t Starting Conjugate Gradient Method to solve AX=B system\n\n\n";
	}
	void zeros(vector<double> &v)
	{
		for (int i = 0; i < v.size(); ++i)
		{
			v[i]=0;
		}
	}
	vector<CG> Solve(vector<vector<CG> > &, const  vector<CG> &);
	vector<CG> SolveCRS(vector<vector<CG> > &, const  vector<CG> &);
	vector<CG> SolvePreCond(vector<vector<CG> > &, vector<vector<CG> > &,const  vector<CG> &);
	vector<CG> SolvePreCondCRS(vector<vector<CG> > &, vector<vector<CG> > &,const  vector<CG> &);

};


/*

	For Sparse Matrices


*/

/*
	For pre conditioned system:

*/

template<typename CG>
vector<CG> CGKRYLOV<CG>::SolvePreCondCRS(vector<vector<CG> > &A, vector<vector<CG> > &A_m,const  vector<CG> &B)
{
	cout<<"Solving preconditioned system using Compressed Row Storage Format....\n";
	CRS<CG> A_crs(A);
	CRS<CG> A_m_crs(A_m);
	static_cast<vector<double>>(B);
	// Computing first residual
	
	s=A_crs.MVMult(s); // computing using the first guess as 0
	for (int i = 0; i < N; ++i)
	{
		r[i]=B[i]-s[i];
	}

	z=A_m_crs.MVMult(r);
	cblas_dcopy(N,z.data(),1,p.data(),1);
	rho=cblas_ddot(N,r.data(),1,z.data(),1);

	int niter = 0;     // Init counter for number of iterations
	int flag = 0;      // Init break flag

	double tol=pow(10,-16);			// Convergence criteria
	int maxiter=10000;

	while (cblas_dnrm2(N,r.data(),1) > tol)   // Test break condition
	{
		a=A_crs.MVMult(p);
		alpha = rho/cblas_ddot(N,a.data(),1,p.data(),1);
		cblas_daxpy(N,alpha, p.data(),1,U.data(),1);
		cblas_daxpy(N,-alpha, a.data(),1,r.data(),1);
		z=A_m_crs.MVMult(r);
		rho_new=cblas_ddot(N,r.data(),1,z.data(),1);
		
		for (int i = 0; i < N; ++i)
		{
			p[i]=z[i]+(rho_new/rho)*p[i];
		}
		rho=rho_new;
		niter+=1;
		if (niter == maxiter)         // if max. number of iterations  is reached, break.
	    {
	    	flag = 1;                   
	        cout<<"ERROR!!!! \n Max Iterations Reached!!!!\n"<<"Flag= "<<flag;
	        break;
	    }  
	}
	cout<<"Completed in "<<niter<<" iterations \n";
	return static_cast<vector<CG>>(U);		
}

/*
	For non pre conditioned system:

*/

template<typename CG>
vector<CG> CGKRYLOV<CG>::SolveCRS(vector<vector<CG> > &A, const  vector<CG> &B)
{
	cout<<"Solving using Compressed Row Storage Format....\n";
	CRS<CG> A_crs(A);
	
	static_cast<vector<double>>(B);
	// Computing first residual
	
	s=A_crs.MVMult(s); // computing using the first guess as 0
	for (int i = 0; i < N; ++i)
	{
		r[i]=B[i]-s[i];
	}

	cblas_dcopy(N,r.data(),1,p.data(),1);
	rho=cblas_ddot(N,r.data(),1,r.data(),1);

	int niter = 0;     // Init counter for number of iterations
	int flag = 0;      // Init break flag

	double tol=pow(10,-16);			// Convergence criteria
	int maxiter=10000;

	while (cblas_dnrm2(N,r.data(),1) > tol)   // Test break condition
	{
		a=A_crs.MVMult(p);
		alpha = rho/cblas_ddot(N,a.data(),1,p.data(),1);
		cblas_daxpy(N,alpha, p.data(),1,U.data(),1);
		cblas_daxpy(N,-alpha, a.data(),1,r.data(),1);
		rho_new=cblas_ddot(N,r.data(),1,r.data(),1);
		
		for (int i = 0; i < N; ++i)
		{
			p[i]=r[i]+(rho_new/rho)*p[i];
		}
		rho=rho_new;
		niter+=1;
		if (niter == maxiter)         // if max. number of iterations  is reached, break.
	    {
	    	flag = 1;                   
	        cout<<"ERROR!!!! \n Max Iterations Reached!!!!\n"<<"Flag= "<<flag;
	        break;
	    }  
	}
	cout<<"Completed in "<<niter<<" iterations \n";
	return static_cast<vector<CG>>(U);		
}


/*

	For Dense Matrices

*/


/*
	For pre conditioned system:

*/
template<typename CG>
vector<CG> CGKRYLOV<CG>::SolvePreCond(vector<vector<CG> > &A, vector<vector<CG> > &A_m,const  vector<CG> &B)
{
	cout<<"Solving preconditioned system....\n";
	vector<CG> A_oneD;
	vector<CG> A_m_oneD;
    for (int i = 0; i < N; ++i)
    {
    	for (int j = 0; j < N; ++j)
	    {
	    	A_oneD.push_back(A[i][j]);
	    	A_m_oneD.push_back(A_m[i][j]);
	    }
    }

	vector<double> dummy(N,0);
	static_cast<vector<double>>(B);
	// Computing first residual
	
	cblas_dcopy(N,U.data(),1,s.data(),1);

	// Computing first residual
	cblas_dgemv(CblasRowMajor,CblasNoTrans, N, N,1,A_oneD.data(), N, s.data(), 1, 1, dummy.data(), 1);
	cblas_dcopy(N,dummy.data(),1,s.data(),1);
	zeros(dummy);
	for (int i = 0; i < N; ++i)
	{
		r[i]=B[i]-s[i];
	}
	cblas_dcopy(N,r.data(),1,z.data(),1);
	cblas_dgemv(CblasRowMajor,CblasNoTrans, N, N,1,A_m_oneD.data(), N, z.data(), 1, 1, dummy.data(), 1);
	cblas_dcopy(N,dummy.data(),1,z.data(),1);
	zeros(dummy);
	cblas_dcopy(N,z.data(),1,p.data(),1);
	rho=cblas_ddot(N,r.data(),1,z.data(),1);

	int niter = 0;     // Init counter for number of iterations
	int flag = 0;      // Init break flag

	double tol=pow(10,-16);			// Convergence criteria
	int maxiter=10000;

	while (cblas_dnrm2(N,r.data(),1) > tol)   // Test break condition
	{
		cblas_dgemv(CblasRowMajor,CblasNoTrans, N, N,1,A_oneD.data(), N, p.data(), 1, 1, dummy.data(), 1);
		cblas_dcopy(N,dummy.data(),1,a.data(),1);
		zeros(dummy);		
		alpha = rho/cblas_ddot(N,a.data(),1,p.data(),1);
		cblas_daxpy(N,alpha, p.data(),1,U.data(),1);
		cblas_daxpy(N,-alpha, a.data(),1,r.data(),1);
		cblas_dcopy(N,r.data(),1,z.data(),1);
		cblas_dgemv(CblasRowMajor,CblasNoTrans, N, N,1,A_m_oneD.data(), N, z.data(), 1, 1, dummy.data(), 1);
		cblas_dcopy(N,dummy.data(),1,z.data(),1);
		zeros(dummy);
		rho_new=cblas_ddot(N,r.data(),1,z.data(),1);
		
		for (int i = 0; i < N; ++i)
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
	cout<<"Completed in "<<niter<<" iterations \n";
	return static_cast<vector<CG>>(U);		

}


/*
	For non pre conditioned system:

*/

template<typename CG>
vector<CG> CGKRYLOV<CG>::Solve(vector<vector<CG> > &A, const  vector<CG> &B)
{
	cout<<"Solving....\n";
	vector<CG> A_oneD;

    for (int i = 0; i < N; ++i)
    {
    	for (int j = 0; j < N; ++j)
	    {
	    	A_oneD.push_back(A[i][j]);
	    }
    }

	vector<double> dummy(N,0);
	static_cast<vector<double>>(B);
	// Computing first residual
	
	cblas_dcopy(N,U.data(),1,s.data(),1);

	// Computing first residual
	cblas_dgemv(CblasRowMajor,CblasNoTrans, N, N,1,A_oneD.data(), N, s.data(), 1, 1, dummy.data(), 1);
	cblas_dcopy(N,dummy.data(),1,s.data(),1);
	zeros(dummy);
	for (int i = 0; i < N; ++i)
	{
		r[i]=B[i]-s[i];
	}
	
	cblas_dcopy(N,r.data(),1,p.data(),1);
	rho=cblas_ddot(N,r.data(),1,r.data(),1);

	int niter = 0;     // Init counter for number of iterations
	int flag = 0;      // Init break flag

	double tol=pow(10,-16);			// Convergence criteria
	int maxiter=10000;

	while (cblas_dnrm2(N,r.data(),1) > tol)   // Test break condition
	{
		cblas_dgemv(CblasRowMajor,CblasNoTrans, N, N,1,A_oneD.data(), N, p.data(), 1, 1, dummy.data(), 1);
		cblas_dcopy(N,dummy.data(),1,a.data(),1);
		zeros(dummy);		
		alpha = rho/cblas_ddot(N,a.data(),1,p.data(),1);
		cblas_daxpy(N,alpha, p.data(),1,U.data(),1);
		cblas_daxpy(N,-alpha, a.data(),1,r.data(),1);
		rho_new=cblas_ddot(N,r.data(),1,r.data(),1);
		
		for (int i = 0; i < N; ++i)
		{
			p[i]=r[i]+(rho_new/rho)*p[i];
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
	cout<<"Completed in "<<niter<<" iterations \n";
	return static_cast<vector<CG>>(U);		

}

#endif