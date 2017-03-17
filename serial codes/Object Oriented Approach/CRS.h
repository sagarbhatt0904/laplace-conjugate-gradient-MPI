/*

	Storing Sparse Matrix in Compressed Row Storage format
	Created on: Feb 20, 2017
    Author: Sagar Bhatt

*/

#ifndef CRS_h
#define CRS_h
#include <vector>

using namespace std;

template <typename T>
class CRS
{
	vector<T> val;
	vector<int> row;
	vector<int> col;
	int temp=0;

public:
	CRS(){}
	CRS(vector<vector<T> > &);
	vector<T> MVMult(vector<T> &);
	
};

template <typename T>
CRS<T>::CRS(vector<vector<T> > &matrix)
{
	for (int i = 0; i != matrix.size(); ++i)
	{
		row.push_back(temp);
		for (int j = 0; j != matrix[i].size(); ++j)
		{
			if (matrix[i][j]!=0)
			{
				val.push_back(matrix[i][j]);
				col.push_back(j);
				temp++;
			}
		}		
	}
	row.push_back(temp);
	for (int i = 0; i != matrix.size(); ++i)
	{
		matrix[i].clear();
	}
}

template <typename T> 
vector<T> CRS<T>::MVMult(vector<T> &v)
{
	vector<T> result;
	result.resize(v.size());	
	for (int i = 0; i < v.size(); ++i)
	{	
		result[i]=0;
		for (int k = row[i]; k < row[i+1]; ++k)
		{
			result[i]=result[i]+val[k]*v[col[k]];			
		}
	}	
	return result;
}

#endif