#pragma once
class Bfunction
{
public:
	double koef[3];
	double const delta = 0.0001;
	double functionIn(double r, double z)
	{
		return koef[0] + koef[1] * r + koef[2] * z;
	}
	double DerivativeZZZfunctionIn(double r, double z)
	{
		return (functionIn(r, z - delta) - functionIn(r , z + delta)) / (2 * delta);
	}
	double DerivativeRRRfunctionIn(double r, double z)
	{
		return (functionIn(r - delta, z) - functionIn(r + delta, z)) / (2 * delta);
	}
};