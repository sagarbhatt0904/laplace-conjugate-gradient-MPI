/*

Written by Sagar Bhatt


*/

#include <fstream>

using namespace std;

void contour(vector<vector<double> > &x, int N, char* str, char* str2)
{
	ofstream fout(str);
	fout<<"# vtk DataFile Version 2.0\nVTK from matlab\nASCII\nDATASET STRUCTURED_POINTS\nDIMENSIONS "<<N<<" "<<N<<" "<<1<<"\nSPACING 1 1 1 \nORIGIN 0 0 0\nPOINT_DATA "<<N*N<<"\nSCALARS "<<str2<<" float 1\n"<<"LOOKUP_TABLE default\n";
	for (int i = 0; i < N; ++i)
	{
		for (int j = 0; j < N; ++j)
		{
			fout<<x[i][j]<<" ";
		}
	}
	
	fout.close();
}
