// Written by Sagar Bhatt
#include <fstream>

using namespace std;

// Function to plot contours
void contour(vector<vector<double> > &x, int N, char* str)
{
	ofstream fout(str);
	fout<<"# vtk DataFile Version 2.0\nVTK by Sagar Bhatt\n"<<"ASCII\n"<<"DATASET STRUCTURED_POINTS\n"<<"DIMENSIONS "<<N<<" "<<N<<" "<<1<<"\nSPACING 1 1 1 \n"<<"ORIGIN 0 0 0\n"<<"POINT_DATA "<<N*N<<"\nSCALARS U float 1\n"<<"LOOKUP_TABLE default\n";
	for (int i = 0; i < N; ++i)
	{
		for (int j = 0; j < N; ++j)
		{
			fout<<x[i][j]<<" ";
		}
	}
	
	fout.close();
}
