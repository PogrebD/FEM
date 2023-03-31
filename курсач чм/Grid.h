#pragma once
#include <vector>
#include "BasicFunction.h"
#include "GenD.h"



class Mat
{
public:
	double gamma, L;
};

class Node
{
public:
	double r, z;
};

class Elem
{
public:
	int NodeIndex[3];
	double MatIndex;
	Bfunction BF[3];
	vector<vector<double>> MatrixGGGfinale;
	vector<vector<double>> MatrixMMMfinale;
	vector<double> VectorBBBfinale;
	vector<vector<double>> GIGAMATRIX;
};

class Grid
{
public:
	vector<Node> Nodes;
	vector<Elem> Elems;
	vector<Mat> Mats;
	double width = 0;
	double height = 0;


};