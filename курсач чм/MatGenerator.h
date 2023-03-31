#pragma once
#include <iostream>
#include <fstream>
#include "Config.h"

using namespace std;

class GM 
{
public:
	void GenMat()
	{
		ofstream fout("Mat.txt");
		fout << Cfg::Nmat << endl;
		for (int i = 0; i < Cfg::Nmat; i++)
		{
			fout << 1 << " " << 1 << endl;
		}
		fout.close();
	}
};