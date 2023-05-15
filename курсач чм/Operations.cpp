#include "Operations.h"
#include "Grid.h"
vector<vector<double>> Operations::MultMatrixonValue(vector<vector<double>> matrix, double value)
{
	for (int i = 0; i < matrix.size(); i++)
	{
		for (int j = 0; j < matrix[i].size(); j++)
		{
			matrix[i][j] *= value;
		}
	}
	return matrix;
}

vector<double> Operations::MultVectoronValue(vector<double> vector, double value)
{
	for (int i = 0; i < vector.size(); i++)
	{
		vector[i] *= value;
	}
	return vector;
}

vector<vector<double>> Operations::CumMatrixMatrix(vector<vector<double>> matrix1, vector<vector<double>> matrix2)
{
	vector<vector<double>> matrix3;
	if (matrix1.size() == matrix2.size())
		matrix3.resize(matrix1.size());
	for (int i = 0; i < matrix3.size(); i++)
	{
		if (matrix1[i].size() == matrix2[i].size())
			matrix3[i].resize(matrix1[i].size());
	}
	for (int i = 0; i < matrix3.size(); i++)
	{
		for (int j = 0; j < matrix3[i].size(); j++)
		{
			matrix3[i][j] = matrix1[i][j] + matrix2[i][j];
		}
	}
	return matrix3;
}

vector<double> Operations::CumVectorVector(vector<double>vector1, vector<double> vector2)
{
	vector<double> vector3;
	if (vector1.size() == vector2.size())
		vector3.resize(vector1.size());
	for (int i = 0; i < vector3.size(); i++)
	{
		vector3[i] = vector1[i] + vector2[i];
	}
	return vector3;
}

vector<double> Operations::MultMatrixonVector(vector<int> NodeIndex, vector<vector<double>> matrix, vector<double> _vector, int k)
{
	vector<double> result;
	result.resize(_vector.size());
	for (int i = 0; i < _vector.size(); i++)
	{
		result[i] = 0;
	}
	for (int i = 0; i < matrix.size(); i++)											/////////!!!!!!!!!!!!!!!!!!!!!!!!
	{
		for (int j = 0; j < matrix.size(); j++)
		{
			if (i < j)
			{
				result[NodeIndex[i]] += matrix[j][i] * _vector[NodeIndex[j]];
			}
			result[NodeIndex[i]] += matrix[i][j] * _vector[NodeIndex[j]];
		}
	}
	return result;
}

vector<double> Operations::Calc_q_j1(vector<int> NodeIndex, vector<vector<double>> matrix1, vector<vector<double>> matrix2, vector<double> vector1, double dt, double Nnode, int k)
{
	vector<double> vector2;
	vector2.resize(vector1.size());
	vector<vector<double>> matrix3;
	matrix3.resize(matrix1.size());
	for (int i = 0; i < matrix3.size(); i++)
	{
		matrix3[i].resize(matrix1[i].size());
	}
	matrix2 = Operations::MultMatrixonValue(matrix2, -0.5);
	matrix1 = Operations::MultMatrixonValue(matrix1, 1 / dt);
	matrix3 = Operations::CumMatrixMatrix(matrix1, matrix2);

	vector2 = Operations::MultMatrixonVector( NodeIndex , matrix3, vector1, k);
	return vector2;
}