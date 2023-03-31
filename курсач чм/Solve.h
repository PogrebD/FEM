#pragma once
#include "Grid.h"
#include "Config.h"
#include <iomanip>


class Solve
{
public:
	double soooooolve;
	void finale(Grid _Grid, SLAU  _SLAU)
	{
		setlocale(LC_ALL, "Russian");
		cout << "Введите координату r" << endl;
		cin >> findR;
		cout << "Введите координату z" << endl;
		cin >> findZ;
		koefZ = findZ / _Grid.height;
		koefR = findR / _Grid.width;
		if (koefR > (Cfg::xElem - 1))
			koefR = koefR - 1;
		if (koefZ > (Cfg::yElem - 1))
			koefZ = koefZ - 1;
		koefdot = ((koefR + 1) * _Grid.width - (findR - Cfg::x1)) / ((findZ - Cfg::y1) - koefZ * _Grid.height);
		if (findZ - Cfg::y1 - koefZ * _Grid.height <= 0.0001)
			koefdot = 1;

		koefdiag = _Grid.width / _Grid.height;

		if (koefdot > koefdiag)
			IndexElem = ((koefZ * Cfg::xElem) + koefR) * 2;
		else
			IndexElem = (((koefZ * Cfg::xElem) + koefR) * 2) + 1;
		soooooolve = 0;
		for (int i = 0; i < 3; i++)
		{
			soooooolve += _SLAU.q[_Grid.Elems[IndexElem].NodeIndex[i]] * _Grid.Elems[IndexElem].BF[i].functionIn(findR, findZ);
		}
		cout << "ВООООТ: " << setprecision(15) << soooooolve;
		for (int i = 0; i < (Cfg::xElem + 1) * (Cfg::yElem + 1); i++)
		{
			cout << endl;
			cout << setprecision(15) <<  _SLAU.q[i] << " ";
		}

	}
private:
	double findR, findZ, koefdiag, koefdot;
	int koefZ, koefR, IndexElem;

};