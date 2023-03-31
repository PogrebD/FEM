#pragma once
#include <iostream>
#include <fstream>
#include "Config.h"

using namespace std;

class GE {
public:
	void GenElem()
	{
		ofstream fout("Elem.txt");
		int NElem = Cfg::xElem * Cfg::yElem * 2, v1 = 0, v2 = 0, v3 = 0;
		fout << NElem << endl;
		for (int i = 0; i < Cfg::yElem; i++)
		{
			for (int j = 0; j < Cfg::xElem; j++)
			{
				v1 = j + (Cfg::xElem + 1) * i;
				v2 = v1 + 1;
				v3 = j + (Cfg::xElem + 1) * (i + 1);
				fout << v1 << " " << v2 << " " << v3 << " " << 0 << " " << endl;
				fout << v2 << " " << v3 << " " << v3 + 1 << " " << 0 << " " << endl;
			}
		}
		fout.close();
	}
};