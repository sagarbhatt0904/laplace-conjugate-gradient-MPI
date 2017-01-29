INSTRUCTIONS FOR BUILDING THE CODE:

1. Launch terminal in the 'Parallel Code' directory and type 'make' 

2. execute 'lap', it will generate u_CG_mpi.vtk file which can be used to plot the contour in Paraview.

NOTE: Parallel codes will require mkl library to be installed.

	  Serial codes will require CBLAS & LAPACK to be installed. Serial Codes can be compiled using:
	  g++ -O2 <code_name>.cpp -lblas -o <code_name>

	  For serial implementation with lapack:
	  g++ -O2 <code_name>.cpp -lblas -llapack -o <code_name>

	  Python implementation requires scipy to be installed	    