#pragma once
#include <vector>
#include "Grid.h"
using namespace std;

class Operations
{
public:
	static vector<vector<double>> MultMatrixonValue(vector<vector<double>> matrix, double value);
	static vector<vector<double>> CumMatrixMatrix(vector<vector<double>> matrix1, vector<vector<double>> matrix2);
	static vector<double> MultVectoronValue(vector<double> vector, double value);
	static vector<double> CumVectorVector(vector<double> vector1, vector<double> vector2);
	static vector<double> MultMatrixonVector(vector<int> NodeIndex, vector<vector<double>> matrix, vector<double> vector, int k);
	static vector<double> Calc_q_j1(vector<int> NodeIndex, vector<vector<double>> matrix1, vector<vector<double>> matrix2, vector<double> vector1, double dt, double Nnode, int k);
};

