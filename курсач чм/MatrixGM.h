#pragma once
#include <vector>
#include <iostream>
#include "Grid.h"

class MatrixMG
{
public:

	MatrixMG(Grid& _Grid)
	{
		GIGAGRID = &_Grid;
		_heightStep = GIGAGRID->height / _numberOfSplits;
		_widthStep = GIGAGRID->width / _numberOfSplits;
		generatorGM();
		SumMatrixMG();
	}

	void generatorGM()
	{
		GenMatrixG();
		GenMatrixM();
		MdivGamma();
	}
private:
	vector<double>_rootsLegendrePolynomial = { -0.5773503, 0.5773503 };
	vector<double> _weights = { 1.0, 1.0 };

	int _gaussSplitsNumber = 2;
	int _numberOfSplits = 512;
	double _heightStep;
	double _widthStep;

	Grid* GIGAGRID;
	void GenMatrixG()
	{
		for (int evenElem = 0; evenElem < GIGAGRID->Elems.size(); evenElem += 2)
		{
			GIGAGRID->Elems[evenElem].MatrixGGGfinale = CalcMatrixGEVEN(evenElem);
		}
		for (int oddElem = 1; oddElem < GIGAGRID->Elems.size(); oddElem += 2)
		{
			GIGAGRID->Elems[oddElem].MatrixGGGfinale = CalcMatrixGOOD(oddElem);
		}
	}

	void GenMatrixM()
	{
		for (int evenElem = 0; evenElem < GIGAGRID->Elems.size(); evenElem += 2)
		{
			GIGAGRID->Elems[evenElem].MatrixMMMfinale = CalcMatrixMEVEN(evenElem);
		}
		for (int oddElem = 1; oddElem < GIGAGRID->Elems.size(); oddElem += 2)
		{
			GIGAGRID->Elems[oddElem].MatrixMMMfinale = CalcMatrixMODD(oddElem);
		}
	}

	vector<vector<double>> CalcMatrixGEVEN(int elemIndex)
	{
		vector<vector<double>> stiffnessMatrix(3, vector<double>(3));
		for (int q = 0; q < stiffnessMatrix.size(); q++)
		{
			for (int p = q; p < stiffnessMatrix.size(); p++)
			{
				double outerIntegralValue = 0.0;
				for (int i = 0; i < _gaussSplitsNumber; i++)
				{
					double sumOfOuterIntergal = 0.0;
					for (int r = 0; r < _numberOfSplits; r++)
					{
						double rI = (GIGAGRID->Nodes[GIGAGRID->Elems[elemIndex].NodeIndex[0]].r + r * _widthStep +
							GIGAGRID->Nodes[GIGAGRID->Elems[elemIndex].NodeIndex[0]].r +
							(r + 1) * _widthStep) / 2.0 + _rootsLegendrePolynomial[i] * _widthStep / 2;

						double innerIntergalValue = 0.0;

						for (int j = 0; j < _gaussSplitsNumber; j++)
						{
							double sumOfInnerIntegral = 0.0;
							for (int z = 0; z < _numberOfSplits - r; z++)
							{
								double zJ = (GIGAGRID->Nodes[GIGAGRID->Elems[elemIndex].NodeIndex[0]].z +
									z * _heightStep +
									GIGAGRID->Nodes[GIGAGRID->Elems[elemIndex].NodeIndex[0]].z +
									(z + 1) * _heightStep) / 2.0 +
									_rootsLegendrePolynomial[j] * _heightStep / 2;

								sumOfInnerIntegral +=
									rI * _heightStep *
									(GIGAGRID->Elems[elemIndex].BF[p].DerivativeRRRfunctionIn(rI, zJ) *
										GIGAGRID->Elems[elemIndex].BF[q].DerivativeRRRfunctionIn(rI, zJ) +
										(GIGAGRID->Elems[elemIndex].BF[p].DerivativeZZZfunctionIn(rI, zJ) *
											GIGAGRID->Elems[elemIndex].BF[q].DerivativeZZZfunctionIn(rI, zJ)));
							}
							innerIntergalValue += sumOfInnerIntegral * _weights[j] / 2.0;
						}
						sumOfOuterIntergal += _widthStep * innerIntergalValue;
					}
					outerIntegralValue += _weights[i] / 2.0 * sumOfOuterIntergal;
				}
				stiffnessMatrix[p][q] = outerIntegralValue;
			}
		}
		return stiffnessMatrix;
	}
	vector<vector<double>> CalcMatrixGOOD(int elemIndex)
	{
		vector<vector<double>> stiffnessMatrix(3, vector<double>(3));
		for (int q = 0; q < stiffnessMatrix.size(); q++)
		{
			for (int p = q; p < stiffnessMatrix.size(); p++)
			{
				double outerIntegralValue = 0.0;

				for (int i = 0; i < _gaussSplitsNumber; i++)
				{
					double sumOfOuterIntergal = 0.0;

					for (int r = 0; r < _numberOfSplits; r++)
					{
						double rI = (GIGAGRID->Nodes[GIGAGRID->Elems[elemIndex].NodeIndex[0]].r + r * _widthStep +
							GIGAGRID->Nodes[GIGAGRID->Elems[elemIndex].NodeIndex[0]].r +
							(r + 1) * _widthStep) / 2.0 + _rootsLegendrePolynomial[i] * _widthStep / 2;

						double innerIntergalValue = 0.0;

						for (int j = 0; j < _gaussSplitsNumber; j++)
						{
							double sumOfInnerIntegral = 0.0;
							for (int z = _numberOfSplits - r; z < _numberOfSplits; z++)
							{
								double zJ = (GIGAGRID->Nodes[GIGAGRID->Elems[elemIndex].NodeIndex[0]].z +
									z * _heightStep +
									GIGAGRID->Nodes[GIGAGRID->Elems[elemIndex].NodeIndex[0]].z +
									(z + 1) * _heightStep) / 2.0 +
									_rootsLegendrePolynomial[j] * _heightStep / 2;

								sumOfInnerIntegral +=
									rI * _heightStep *
									(GIGAGRID->Elems[elemIndex].BF[p].DerivativeRRRfunctionIn(rI, zJ) *
										GIGAGRID->Elems[elemIndex].BF[q].DerivativeRRRfunctionIn(rI, zJ) +
										(GIGAGRID->Elems[elemIndex].BF[p].DerivativeZZZfunctionIn(rI, zJ) *
											GIGAGRID->Elems[elemIndex].BF[q].DerivativeZZZfunctionIn(rI, zJ)));
							}

							innerIntergalValue += sumOfInnerIntegral * _weights[j] / 2.0;
						}

						sumOfOuterIntergal += _widthStep * innerIntergalValue;
					}

					outerIntegralValue += _weights[i] / 2.0 * sumOfOuterIntergal;
				}

				stiffnessMatrix[p][q] = outerIntegralValue;
			}
		}

		return stiffnessMatrix;
	}
	vector<vector<double>> CalcMatrixMEVEN(int elemIndex)
	{
		vector<vector<double>> massMatrix(3, vector<double>(3));
		for (int q = 0; q < massMatrix.size(); q++)
		{
			for (int p = q; p < massMatrix.size(); p++)
			{
				double outerIntegralValue = 0.0;

				for (int i = 0; i < _gaussSplitsNumber; i++)
				{
					double sumOfOuterIntergal = 0.0;

					for (int r = 0; r < _numberOfSplits; r++)
					{
						double rI = (GIGAGRID->Nodes[GIGAGRID->Elems[elemIndex].NodeIndex[0]].r + r * _widthStep +
							GIGAGRID->Nodes[GIGAGRID->Elems[elemIndex].NodeIndex[0]].r +
							(r + 1) * _widthStep) / 2.0 + _rootsLegendrePolynomial[i] * _widthStep / 2;

						double innerIntergalValue = 0.0;

						for (int j = 0; j < _gaussSplitsNumber; j++)
						{
							double sumOfInnerIntegral = 0.0;
							for (int z = 0; z < _numberOfSplits - r; z++)
							{
								double zJ = (GIGAGRID->Nodes[GIGAGRID->Elems[elemIndex].NodeIndex[0]].z +
									z * _heightStep +
									GIGAGRID->Nodes[GIGAGRID->Elems[elemIndex].NodeIndex[0]].z +
									(z + 1) * _heightStep) / 2.0 +
									_rootsLegendrePolynomial[j] * _heightStep / 2;

								sumOfInnerIntegral += rI * _heightStep *
									(GIGAGRID->Elems[elemIndex].BF[p].functionIn(rI, zJ) *
										GIGAGRID->Elems[elemIndex].BF[q].functionIn(rI, zJ));
							}

							innerIntergalValue += sumOfInnerIntegral * _weights[j] / 2.0;
						}

						sumOfOuterIntergal += _widthStep * innerIntergalValue;
					}

					outerIntegralValue += _weights[i] / 2.0 * sumOfOuterIntergal;
				}

				massMatrix[p][q] = outerIntegralValue;
			}
		}

		return massMatrix;
	}
	vector<vector<double>> CalcMatrixMODD(int elemIndex)
	{
		vector<vector<double>> massMatrix(3, vector<double>(3));
		for (int q = 0; q < massMatrix.size(); q++)
		{
			for (int p = q; p < massMatrix.size(); p++)
			{
				double outerIntegralValue = 0.0;

				for (int i = 0; i < _gaussSplitsNumber; i++)
				{
					double sumOfOuterIntergal = 0.0;

					for (int r = 0; r < _numberOfSplits; r++)
					{
						double rI = (GIGAGRID->Nodes[GIGAGRID->Elems[elemIndex].NodeIndex[0]].r + r * _widthStep +
							GIGAGRID->Nodes[GIGAGRID->Elems[elemIndex].NodeIndex[0]].r +
							(r + 1) * _widthStep) / 2.0 + _rootsLegendrePolynomial[i] * _widthStep / 2;

						double innerIntergalValue = 0.0;

						for (int j = 0; j < _gaussSplitsNumber; j++)
						{
							double sumOfInnerIntegral = 0.0;
							for (int z = _numberOfSplits - r; z < _numberOfSplits; z++)
							{
								double zJ = (GIGAGRID->Nodes[GIGAGRID->Elems[elemIndex].NodeIndex[0]].z +
									z * _heightStep +
									GIGAGRID->Nodes[GIGAGRID->Elems[elemIndex].NodeIndex[0]].z +
									(z + 1) * _heightStep) / 2.0 +
									_rootsLegendrePolynomial[j] * _heightStep / 2;

								sumOfInnerIntegral += rI * _heightStep *
									(GIGAGRID->Elems[elemIndex].BF[p].functionIn(rI, zJ) *
										GIGAGRID->Elems[elemIndex].BF[q].functionIn(rI, zJ));
							}

							innerIntergalValue += sumOfInnerIntegral * _weights[j] / 2.0;
						}

						sumOfOuterIntergal += _widthStep * innerIntergalValue;
					}

					outerIntegralValue += _weights[i] / 2.0 * sumOfOuterIntergal;
				}

				massMatrix[p][q] = outerIntegralValue;
			}
		}

		return massMatrix;
	}

	void MdivGamma()
	{
		for (int k = 0; k < GIGAGRID->Elems.size(); k++) {
			for (int i = 0; i < 3; i++) {
				for (int j = 0; j <= i; j++) {
					GIGAGRID->Elems[k].MatrixMMMfinale[i][j] = GIGAGRID->Elems[k].MatrixMMMfinale[i][j] * GIGAGRID->Mats[GIGAGRID->Elems[k].MatIndex].gamma;/////?
				}
			}
		}
	}

	void SumMatrixMG()
	{
		double sum;

		for (int k = 0; k < GIGAGRID->Elems.size(); k++) 
		{
			GIGAGRID->Elems[k].GIGAMATRIX.resize(3);
				for (int i = 0; i < 3; i++)
				{
					GIGAGRID->Elems[k].GIGAMATRIX[i].resize(3);
				}
			for (int i = 0; i < 3; i++) {
				for (int j = 0; j <= i; j++) {
					GIGAGRID->Elems[k].GIGAMATRIX[i][j] = GIGAGRID->Elems[k].MatrixGGGfinale[i][j] + GIGAGRID->Elems[k].MatrixMMMfinale[i][j];
				}
			}
		}
	}



};