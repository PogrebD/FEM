#pragma once
#include <iostream>
#include <fstream>
#include "Config.h"

using namespace std;
typedef double type;
class GN
{
public:
	void GenNode() {

		ofstream fout("Node.txt");
		int Nnode = (Cfg::xElem + 1) * (Cfg::yElem + 1);
		fout << Nnode << endl;
		type stepx = (Cfg::x2 - Cfg::x1) / Cfg::xElem;
		type stepy = (Cfg::y2 - Cfg::y1) / Cfg::yElem;

		type ycord = Cfg::y1;
		for (int i = 0; i < Cfg::yElem + 1; i++)
		{
			type xcord = Cfg::x1;
			for (int j = 0; j < Cfg::xElem + 1; j++)
			{
				fout << xcord << " " << ycord << " " << endl;
				xcord += stepx;
			}
			ycord += stepy;
		}
		fout.close();
	}
};