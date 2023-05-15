#pragma once
#include <vector>
#include <iostream>
#include <fstream>
#include "BasicFunction.h"
using namespace std;

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

class Time
{
public:
	double dt;
	vector<double> timeSloy;
	int Ntime;
	void Tparse()
	{
		ifstream fin("Time.txt");
		fin >> Ntime;
		timeSloy.resize(Ntime+1);
		fin >> timeSloy[0] >> timeSloy[Ntime];
		dt = (timeSloy[Ntime] - timeSloy[0]) / Ntime;
		for (int i = 0; i < Ntime; i++)
		{
			timeSloy[i + 1] = timeSloy[i] + dt;
		}
	}
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
	int N;
	Time _Time;

};