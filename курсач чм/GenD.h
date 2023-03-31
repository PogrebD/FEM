#pragma once
#include <vector>
#include <iostream>
#include "Grid.h"

using namespace std;
typedef double type;

class GenD
{
public:
	void D(Grid& _Grid, int NElem)
	{

		D_1.resize(3);
		for (int i = 0; i < 3; i++)
		{
			D_1[i].resize(3);
		}
		for (int i = 0; i < NElem; i++)
		{
			type x21 = _Grid.Nodes[_Grid.Elems[i].NodeIndex[1]].r - _Grid.Nodes[_Grid.Elems[i].NodeIndex[0]].r;
			type y31 = _Grid.Nodes[_Grid.Elems[i].NodeIndex[2]].z - _Grid.Nodes[_Grid.Elems[i].NodeIndex[0]].z;
			type x31 = _Grid.Nodes[_Grid.Elems[i].NodeIndex[2]].r - _Grid.Nodes[_Grid.Elems[i].NodeIndex[0]].r;
			type y21 = _Grid.Nodes[_Grid.Elems[i].NodeIndex[1]].z - _Grid.Nodes[_Grid.Elems[i].NodeIndex[0]].z;

			detD = (x21 * y31) - (x31 * y21);

			D_1[0][0] = (_Grid.Nodes[_Grid.Elems[i].NodeIndex[1]].r * _Grid.Nodes[_Grid.Elems[i].NodeIndex[2]].z) -
				(_Grid.Nodes[_Grid.Elems[i].NodeIndex[2]].r * _Grid.Nodes[_Grid.Elems[i].NodeIndex[1]].z);
			D_1[1][0] = (_Grid.Nodes[_Grid.Elems[i].NodeIndex[2]].r * _Grid.Nodes[_Grid.Elems[i].NodeIndex[0]].z) -
				(_Grid.Nodes[_Grid.Elems[i].NodeIndex[0]].r * _Grid.Nodes[_Grid.Elems[i].NodeIndex[2]].z);
			D_1[2][0] = (_Grid.Nodes[_Grid.Elems[i].NodeIndex[0]].r * _Grid.Nodes[_Grid.Elems[i].NodeIndex[1]].z) -
				(_Grid.Nodes[_Grid.Elems[i].NodeIndex[1]].r * _Grid.Nodes[_Grid.Elems[i].NodeIndex[0]].z);
			D_1[0][1] = _Grid.Nodes[_Grid.Elems[i].NodeIndex[1]].z - _Grid.Nodes[_Grid.Elems[i].NodeIndex[2]].z;
			D_1[1][1] = y31;
			D_1[2][1] = y21 * (-1);
			D_1[0][2] = _Grid.Nodes[_Grid.Elems[i].NodeIndex[2]].r - _Grid.Nodes[_Grid.Elems[i].NodeIndex[1]].r;
			D_1[1][2] = x31 * (-1);
			D_1[2][2] = x21;
			for (int j = 0; j < 3; j++)
			{
				for (int k = 0; k < 3; k++)
				{
					D_1[j][k] /= detD;
				}
			}
			for (int j = 0; j < 3; j++)
			{
				for (int k = 0; k < 3; k++)
				{
					_Grid.Elems[i].BF[j].koef[k] = D_1[j][k];
				}
			}
		}
	}
private:
	vector<vector<type>> D_1;
	type detD;
};