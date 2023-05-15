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
#include "Operations.h"


typedef double type;
vector<vector<double>> qT;
vector<double> q0, q_j1, b_j1, Mq_j1;
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
		fin1 >> _BC.BC1vector[i].versh[0] >> _BC.BC1vector[i].versh[1];
	}
}

void InputGrid(Grid& _Grid, int& Nnode, int& NElem, int& NMat)
{
	ifstream finN("Node.txt");
	ifstream finE("Elem.txt");
	ifstream finM("Mat.txt");
	ifstream finT("Time.txt");
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

	_Grid._Time.Tparse();
}

int main()
{
	cout.imbue(std::locale("German_germany"));
	vector<int> indx;
	indx.resize(3);
	Gen();
	InputGrid(_Grid, Nnode, NElem, NMat);
	inputbc(_BC, Nbc2, Nbc3, Nbc1);
	_GenD.D(_Grid, NElem);

	q0.resize(Nnode);
	q_j1.resize(Nnode);
	b_j1.resize(Nnode);
	Mq_j1.resize(Nnode);
	for (int i = 0; i < Nnode; i++)
	{
		q0[i] = Cfg::U(_Grid.Nodes[i].r, _Grid.Nodes[i].z, _Grid._Time.timeSloy[0]);
	}
	qT.push_back(q0);
	MatrixMG _MatrixMG(_Grid);

	VectorB _VectorB(_Grid);
	_VectorB.CalcB(_Grid._Time.timeSloy[0]);


	GlobalMatrix _GlobalMatrix(_Grid, NElem, Nnode);
	_GlobalMatrix.vstavkaB(_Grid, NElem);
	_GlobalMatrix.vstavkaGM(_Grid, NElem);
//	_BC.primeniaemKraevble1(_Grid, _GlobalMatrix, Nbc2, Nbc3, Nbc1, Nnode);

	b_j1 = _GlobalMatrix.Globalvector;

	vector<double> GlobalGMtriangleSave = _GlobalMatrix.GlobalGMtriangle;
	vector<double> GlobalGMdiagSave = _GlobalMatrix.GlobalGMdiag;

	_VectorB.CalcB(_Grid._Time.timeSloy[1]);
	_GlobalMatrix.vstavkaB(_Grid, NElem);
//	_BC.primeniaemKraevble1(_Grid, _GlobalMatrix, Nbc2, Nbc3, Nbc1, Nnode);

	vector<double> buf1;
	buf1.resize(Nnode);

	buf1 = Operations::CumVectorVector(b_j1, _GlobalMatrix.Globalvector);
	buf1 = Operations::MultVectoronValue(buf1, 0.5);
	b_j1 = _GlobalMatrix.Globalvector;

	for (int k = 0; k < NElem; k++)
	{
		vector<double> buf;
		buf.resize(Nnode);
		buf = Operations::Calc_q_j1({ _Grid.Elems[k].NodeIndex[0], _Grid.Elems[k].NodeIndex[1], _Grid.Elems[k].NodeIndex[2]}, _Grid.Elems[k].MatrixMMMfinale, _Grid.Elems[k].MatrixGGGfinale, q0, _Grid._Time.dt, Nnode, k);
		Mq_j1 = Operations::CumVectorVector(Mq_j1, buf);
	}
	_GlobalMatrix.Globalvector = Operations::CumVectorVector(Mq_j1, buf1);

	_BC.primeniaemKraevble1(_Grid, _GlobalMatrix, Nbc2, Nbc3, Nbc1, Nnode, 1);

	SLAU _SLAU(_GlobalMatrix, Nnode);

	_SLAU.MSG_no(Nnode);
	q_j1 = _SLAU.q;
	qT.push_back(q_j1);


	for (int i = 2; i < _Grid._Time.Ntime + 1; i++)
	{
		Mq_j1.clear();
		Mq_j1.resize(Nnode);
		for (int k = 0; k < NElem; k++)
		{
			vector<double> buf;
			buf.resize(Nnode);
			buf = Operations::Calc_q_j1({ _Grid.Elems[k].NodeIndex[0], _Grid.Elems[k].NodeIndex[1], _Grid.Elems[k].NodeIndex[2] }, _Grid.Elems[k].MatrixMMMfinale, _Grid.Elems[k].MatrixGGGfinale, q_j1, _Grid._Time.dt, Nnode, k);
			Mq_j1 = Operations::CumVectorVector(Mq_j1, buf);
		}

		_VectorB.CalcB(_Grid._Time.timeSloy[i]);
		_GlobalMatrix.vstavkaB(_Grid, NElem);

		buf1 = Operations::CumVectorVector(b_j1, _GlobalMatrix.Globalvector);
		buf1 = Operations::MultVectoronValue(buf1, 0.5);
		b_j1 = _GlobalMatrix.Globalvector;

		_GlobalMatrix.Globalvector = Operations::CumVectorVector(Mq_j1, buf1);

		_BC.primeniaemKraevble2(_Grid, _GlobalMatrix, GlobalGMtriangleSave, GlobalGMdiagSave, Nbc2, Nbc3, Nbc1, Nnode, i);

		_SLAU.MSG_no(Nnode);
		q_j1 = _SLAU.q;
		qT.push_back(q_j1);
		auto g = (double)_Grid._Time.Ntime / i;
		if ((((double)_Grid._Time.Ntime / i) == 1) || (((double)_Grid._Time.Ntime / i) == 1.5) || (((double)_Grid._Time.Ntime / i) == 3))
		{			
			cout << setprecision(15) << q_j1[4] << endl;
		}

	}
	cout << Cfg::U(453435, 2.5, _Grid._Time.timeSloy[10]) << endl << Cfg::U(453435, 2.5, _Grid._Time.timeSloy[20]) << endl << Cfg::U(453435, 2.5, _Grid._Time.timeSloy[30]) << endl;


	_Solve.finale(_Grid, _SLAU);

}


