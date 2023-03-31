#pragma once
#include <vector>
#include <iostream>
#include "Grid.h"
#include "functionB.h"

class VectorB
{
public:
	VectorB(Grid& _Grid)
	{
		GIGAGRID = &_Grid;
		CalcB();
	}

	void CalcF(int i)
	{
		for (int j = 0; j < 3; j++)
		{
			F[j] = FunctionBclass::FunctionB(GIGAGRID->Nodes[GIGAGRID->Elems[i].NodeIndex[j]].r, GIGAGRID->Nodes[GIGAGRID->Elems[i].NodeIndex[j]].z);
		}
	}

	void CalcB()
	{
		B.resize(3);
		for (int k = 0; k < GIGAGRID->Elems.size(); k++)
		{
			for (int i = 0; i < 3; i++)
			{
				B[i] = 0;
			}
			CalcF(k);
			for (int i = 0; i < 3; i++)
			{
				for (int j = 0; j < 3 && i != j; j++)
				{
					B[i] += F[j] * GIGAGRID->Elems[k].MatrixMMMfinale[i][j]*GIGAGRID->Mats[GIGAGRID->Elems[k].MatIndex].gamma;
					B[j] += F[i] * GIGAGRID->Elems[k].MatrixMMMfinale[i][j] * GIGAGRID->Mats[GIGAGRID->Elems[k].MatIndex].gamma;
				}
				B[i] += F[i] * GIGAGRID->Elems[k].MatrixMMMfinale[i][i] * GIGAGRID->Mats[GIGAGRID->Elems[k].MatIndex].gamma;
			}
			GIGAGRID->Elems[k].VectorBBBfinale = B;
		}
	}
private:
	vector<double> B;
	Grid* GIGAGRID;
	double F[3];
};