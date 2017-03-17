# ABOUT THE CODE:

* Folder:
	* Parallel Code: This folder contains parallel code with MPI. Use ```make``` command to build the code 
	* serial codes: This folder contains:
		* analytical.cpp : analytical solution to laplace equation
		* contour.h : header needed to generate VTK file
		* laplace.py : python solution to the problem using conjugate gradient method
		* laplace_krylov_cg_sparse.cpp : solution to the problem using conjugate gradient method
		* laplace_lapack.cpp : Direct solution to the problem. 
		* Object-Oriented Approach: This folder contains object oriented approach to the problem. Refer its README for further information

### NOTE: 
		Parallel codes will require mkl library to be installed.

	  	Serial codes will require CBLAS & LAPACK to be installed. Serial Codes can be
	  	compiled using:
			g++ -O2 <code_name>.cpp -lblas -o <code_name>

		For serial implementation with lapack:
			g++ -O2 <code_name>.cpp -lblas -llapack -o <code_name>

	  	Python implementation requires scipy to be installed	    
