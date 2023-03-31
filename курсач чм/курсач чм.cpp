#include <fstream>
#include <vector>
#include <iostream>
#include <string>
#include <set>
#include <math.h>
#include "Config.h"
#include "ElemGenerator.h"
#include "MatGenerator.h"
#include "NodeGenerator.h"
#include "GenD.h"
#include "BasicFunction.h"
#include "Grid.h"
#include "MatrixGM.h"
#include "VectorB.h"
#include "GlobalMatrixes.h"
#include "BC.h"
#include "MSG.h"
#include "Solve.h"

using namespace std;
typedef double type;


int Nnode, NElem, NMat, Nbc2, Nbc3, Nbc1;
Grid _Grid;
vector<int> ig, jg;
GenD _GenD;
BC _BC;
Solve _Solve;



void Gen()
{
	GM BLIAT;
	BLIAT.GenMat();
	GE BLIAT1;
	BLIAT1.GenElem();
	GN BLIAT2;
	BLIAT2.GenNode();
}

void inputbc(BC& _BC, int& Nbc2, int& Nbc3, int& Nbc1)
{
	ifstream fin2("bc2.txt");
	ifstream fin3("bc3.txt");
	ifstream fin1("bc1.txt");
	fin2 >> Nbc2;
	_BC.BC2vector.resize(Nbc2);
	for (int i = 0; i < Nbc2; i++)
	{
		fin2 >> _BC.BC2vector[i].versh[0] >> _BC.BC2vector[i].versh[1] >> _BC.BC2vector[i].theta[0] >> _BC.BC2vector[i].theta[1];
	}
	fin3 >> Nbc3;
	_BC.BC3vector.resize(Nbc3);
	for (int i = 0; i < Nbc3; i++)
	{
		fin3 >> _BC.BC3vector[i].versh[0] >> _BC.BC3vector[i].versh[1] >> _BC.BC3vector[i].Beta >> _BC.BC3vector[i].U[0] >> _BC.BC3vector[i].U[1];
	}
	fin1 >> Nbc1;
	_BC.BC1vector.resize(Nbc1);
	for (int i = 0; i < Nbc1; i++)
	{
		fin1 >> _BC.BC1vector[i].versh[0] >> _BC.BC1vector[i].versh[1] >> _BC.BC1vector[i].U;
	}
}

void InputGrid(Grid& _Grid, int& Nnode, int& NElem, int& NMat)
{
	ifstream finN("Node.txt");
	ifstream finE("Elem.txt");
	ifstream finM("Mat.txt");
	finN >> Nnode;
	_Grid.Nodes.resize(Nnode);
	for (size_t i = 0; i < Nnode; i++)
	{
		finN >> _Grid.Nodes[i].r >> _Grid.Nodes[i].z;
	}
	finE >> NElem;
	_Grid.Elems.resize(NElem);
	for (size_t i = 0; i < NElem; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			finE >> _Grid.Elems[i].NodeIndex[j];
		}
		finE >> _Grid.Elems[i].MatIndex;
	}
	finM >> NMat;
	_Grid.Mats.resize(NMat);
	for (size_t i = 0; i < NMat; i++)
	{
		finM >> _Grid.Mats[i].gamma >> _Grid.Mats[i].L;
	}
	_Grid.height = _Grid.Nodes[_Grid.Elems[0].NodeIndex[2]].z - _Grid.Nodes[0].z;
	_Grid.width = _Grid.Nodes[1].r - _Grid.Nodes[0].r;
}

int main()
{

	Gen();
	InputGrid(_Grid, Nnode, NElem, NMat);
	_GenD.D(_Grid, NElem);
	MatrixMG _MatrixMG(_Grid);
	VectorB _VectorB(_Grid);
	GlobalMatrix _GlobalMatrix(_Grid, NElem, Nnode);
	inputbc(_BC, Nbc2, Nbc3, Nbc1);
	_BC.primeniaemKraevble(_Grid, _GlobalMatrix, Nbc2, Nbc3, Nbc1, Nnode);
	SLAU _SLAU(_GlobalMatrix, Nnode);
	_Solve.finale(_Grid, _SLAU);

}


