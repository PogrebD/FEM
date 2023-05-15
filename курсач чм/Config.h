#pragma once
class Cfg
{
public:
	typedef double type;
	static const int Split = 2;
	static const int Nmat = 1;
	static const int xElem = 1 * Split;
	static const int yElem = 1 * Split;
	static constexpr type x1 = 1, x2 = 4, y1 = 1, y2 = 4;
	static double U(double r, double z, double t)
	{
		return z + t*t;
	}
	static double bc2(double r, double z, double t)
	{
		return 1;
	}

};