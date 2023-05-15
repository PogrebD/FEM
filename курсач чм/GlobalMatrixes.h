#pragma once
#include"MatrixGM.h"

class GlobalMatrix
{
public:
	GlobalMatrix(Grid _Grid, int NElem, int Nnode)
	{
		portrait(_Grid, NElem, Nnode);
		GlobalGMtriangle.resize(ig[Nnode - 1]);
		GlobalGMdiag.resize(Nnode);
		Globalvector.resize(Nnode);
		MGlobal.resize(Nnode);
		GGlobal.resize(Nnode);
		for (int i = 0; i < Nnode; i++)
		{
			GGlobal[i].resize(Nnode);
			MGlobal[i].resize(Nnode);
		}

	}
	vector<double> GlobalGMtriangle, GlobalGMdiag, Globalvector;
	vector<vector<double>> MGlobal, GGlobal;
	vector<int> ig, jg;
	void portrait(Grid& _Grid, int NElem, int Nnode)
	{

		versh.resize(Nnode);
		ig.resize(Nnode);

		for (int i = 0; i < NElem; i++)
		{
			versh[_Grid.Elems[i].NodeIndex[1]].insert(_Grid.Elems[i].NodeIndex[0]);
			for (int j = 0; j < 2; j++)
			{
				versh[_Grid.Elems[i].NodeIndex[2]].insert(_Grid.Elems[i].NodeIndex[j]);
			}
		}
		ig[0] = 0;
		for (int i = 1; i < Nnode; i++)
		{
			ig[i] = ig[i - 1] + versh[i].size();

		}
		for (int i = 0; i < versh.size(); i++)
		{
			for (set<int>::iterator it = versh[i].begin(); it != versh[i].end(); ++it)
				jg.push_back(*it);
		}
	}

	//void vstavkaGM2(Grid _Grid, int NElem, int Nnode)
	//{
	//	for (int k = 0; k < NElem; k++)
	//	{
	//		for (int i = 0; i < 3; i++)
	//		{
	//			for (int j = 0; j < 3; j++)
	//			{
	//				MGlobal[i + k][j + k] += _Grid.Elems[k].MatrixMMMfinale[i][j];
	//			}
	//		}
	//	}
	//}

	void vstavkaGM(Grid _Grid, int NElem)
	{
		for (int i = 0; i < NElem; i++)
		{
			for (int j = 0; j < 3; j++)
			{
				GlobalGMdiag[_Grid.Elems[i].NodeIndex[j]] += _Grid.Elems[i].GIGAMATRIX[j][j];
			}
		}
		for (int i = 0; i < NElem; i++)
		{
			for (int k = ig[_Grid.Elems[i].NodeIndex[1] - 1]; k < ig[_Grid.Elems[i].NodeIndex[1]]; k++)
			{
				if (jg[k] == _Grid.Elems[i].NodeIndex[0])
					GlobalGMtriangle[k] += _Grid.Elems[i].GIGAMATRIX[1][0];
			}
			for (int j = 0; j < 2; j++)
			{
				for (int k = ig[_Grid.Elems[i].NodeIndex[2] - 1]; k < ig[_Grid.Elems[i].NodeIndex[2]]; k++)
				{
					if (jg[k] == _Grid.Elems[i].NodeIndex[j]) {
						GlobalGMtriangle[k] += _Grid.Elems[i].GIGAMATRIX[2][j];
						k += 3;
					}
				}
			}
		}

	}

	void vstavkaB(Grid _Grid, int NElem)
	{
		int n = Globalvector.size();
		Globalvector.clear();
		Globalvector.resize(n);

		for (int i = 0; i < NElem; i++)
		{
			for (int j = 0; j < 3; j++)
			{
				Globalvector[_Grid.Elems[i].NodeIndex[j]] += _Grid.Elems[i].VectorBBBfinale[j]; 
			}
		}
	}

private:
	vector<set<int>> versh;
};
