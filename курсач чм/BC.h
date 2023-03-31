#pragma once
#include <vector>
#include "GlobalMatrixes.h"
#include "Grid.h"
using namespace std;

class BC2
{
public:
	int versh[2];
	double theta[2];
};

class BC3
{
public:
	int versh[2];
	double Beta, U[2];
};

class BC1
{
public:
	int versh[2];
	double U;
};


class BC
{
public:
	vector<BC2> BC2vector;
	vector<BC3> BC3vector;
	vector<BC1> BC1vector;
	void primeniaemKraevble(Grid _Grid, GlobalMatrix& _GlobalMatrix, int Nbc2, int Nbc3, int Nbc1, int Nnode)
	{
		double h;
		// 2
		for (int i = 0; i < Nbc2; i++)
		{
			if (BC2vector[i].versh[1] - BC2vector[i].versh[0] == 1)
			{
				h = _Grid.width;
				_GlobalMatrix.Globalvector[BC2vector[i].versh[0]] += ((h * _Grid.Nodes[BC2vector[i].versh[0]].r * (BC2vector[i].theta[0] * 2 + BC2vector[i].theta[1])) / 6) + ((h * h * (BC2vector[i].theta[0] + BC2vector[i].theta[1])) / 12);
				_GlobalMatrix.Globalvector[BC2vector[i].versh[1]] += ((h * _Grid.Nodes[BC2vector[i].versh[0]].r * (BC2vector[i].theta[1] * 2 + BC2vector[i].theta[0])) / 6) + ((h * h * (BC2vector[i].theta[0] + BC2vector[i].theta[1] * 3)) / 12);
			}
			else
			{
				h = _Grid.height;

				_GlobalMatrix.Globalvector[BC2vector[i].versh[0]] += (h * _Grid.Nodes[BC2vector[i].versh[0]].r * (BC2vector[i].theta[0] * 2 + BC2vector[i].theta[1])) / 6;
				_GlobalMatrix.Globalvector[BC2vector[i].versh[1]] += (h * _Grid.Nodes[BC2vector[i].versh[0]].r * (BC2vector[i].theta[1] * 2 + BC2vector[i].theta[0])) / 6;
			}
		}
		// 3 

		for (int i = 0; i < Nbc3; i++)
		{
			if (BC3vector[i].versh[1] - BC3vector[i].versh[0] == 1)
			{
				h = _Grid.width;
				_GlobalMatrix.Globalvector[BC3vector[i].versh[0]] += ((h * _Grid.Nodes[BC3vector[i].versh[0]].r * (BC3vector[i].U[0] * 2 + BC3vector[i].U[1]) * BC3vector[i].Beta) / 6) + ((h * h * (BC3vector[i].U[0] + BC3vector[i].U[1]) * BC3vector[i].Beta) / 12);
				_GlobalMatrix.Globalvector[BC3vector[i].versh[1]] += ((h * _Grid.Nodes[BC3vector[i].versh[0]].r * (BC3vector[i].U[1] * 2 + BC3vector[i].U[0]) * BC3vector[i].Beta) / 6) + ((h * h * (BC3vector[i].U[0] + BC3vector[i].U[1] * 3) * BC3vector[i].Beta) / 12);
				for (int j = _GlobalMatrix.ig[BC3vector[i].versh[1] - 1]; j < _GlobalMatrix.ig[BC3vector[i].versh[1]]; j++)
				{
					if (_GlobalMatrix.jg[j] == BC3vector[i].versh[0])
						_GlobalMatrix.GlobalGMtriangle[j] += ((h * _Grid.Nodes[BC3vector[i].versh[0]].r * BC3vector[i].Beta) / 6) + ((h * h * BC3vector[i].Beta) / 12);
				}

				_GlobalMatrix.GlobalGMdiag[BC3vector[i].versh[0]] += ((h * _Grid.Nodes[BC3vector[i].versh[0]].r * BC3vector[i].Beta) / 3) + ((h * h * BC3vector[i].Beta) / 12);
				_GlobalMatrix.GlobalGMdiag[BC3vector[i].versh[1]] += ((h * _Grid.Nodes[BC3vector[i].versh[0]].r * BC3vector[i].Beta) / 3) + ((h * h * BC3vector[i].Beta) / 4);

			}
			else
			{
				h = _Grid.height;
				_GlobalMatrix.Globalvector[BC3vector[i].versh[0]] += (h * _Grid.Nodes[BC3vector[i].versh[0]].r * (BC3vector[i].U[0] * 2 + BC3vector[i].U[1]) * BC3vector[i].Beta) / 6;
				_GlobalMatrix.Globalvector[BC3vector[i].versh[1]] += (h * _Grid.Nodes[BC3vector[i].versh[0]].r * (BC3vector[i].U[1] * 2 + BC3vector[i].U[0]) * BC3vector[i].Beta) / 6;
				for (int j = _GlobalMatrix.ig[BC3vector[i].versh[1] - 1]; j < _GlobalMatrix.ig[BC3vector[i].versh[1]]; j++)
				{
					if (_GlobalMatrix.jg[j] == BC3vector[i].versh[0])
						_GlobalMatrix.GlobalGMtriangle[j] += (h * _Grid.Nodes[BC3vector[i].versh[0]].r * BC3vector[i].Beta) / 6;
				}

				for (int j = 0; j < 2; j++)
				{
					_GlobalMatrix.GlobalGMdiag[BC3vector[i].versh[j]] += (h * _Grid.Nodes[BC3vector[i].versh[0]].r * BC3vector[i].Beta) / 3;
				}
			}



		}

		// 1
		for (int i = 0; i < Nbc1; i++)
		{
			for (int k = 0; k < 2; k++)
			{
				if (BC1vector[i].versh[k] != 0)
					for (int j = _GlobalMatrix.ig[BC1vector[i].versh[k] - 1]; j < _GlobalMatrix.ig[BC1vector[i].versh[k]]; j++)
					{
						_GlobalMatrix.Globalvector[_GlobalMatrix.jg[j]] -= BC1vector[i].U * _GlobalMatrix.GlobalGMtriangle[j];
						_GlobalMatrix.GlobalGMtriangle[j] = 0;
					}
				_GlobalMatrix.Globalvector[BC1vector[i].versh[k]] = BC1vector[i].U;
				_GlobalMatrix.GlobalGMdiag[BC1vector[i].versh[k]] = 1;
				vector<double> ystal;
				for (int j = BC1vector[i].versh[k] + 1; j < Nnode; j++)
				{
					for (int h = _GlobalMatrix.ig[j - 1]; h < _GlobalMatrix.ig[j]; h++)
					{
						if (_GlobalMatrix.jg[h] == BC1vector[i].versh[k])
						{
							_GlobalMatrix.Globalvector[j] -= BC1vector[i].U * _GlobalMatrix.GlobalGMtriangle[h];
							_GlobalMatrix.GlobalGMtriangle[h] = 0;
						}

					}
				}
			}
		}
	}
private:
	double b[2];
};
