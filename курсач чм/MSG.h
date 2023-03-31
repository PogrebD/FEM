#pragma once
#include "GlobalMatrixes.h"
class SLAU
{
public:

	SLAU(GlobalMatrix& _GlobalMatrix, int Nnode)
	{
		GIGAGLOBALMATRIX = &_GlobalMatrix;
		n = Nnode;
		MSG_no(Nnode);
	}

	vector<double> q;
	void MSG_no(int Nnode)
	{
		vector<type> z, r, y;
		z.resize(Nnode); r.resize(Nnode); y.resize(Nnode);
		int k_iter;
		int maxiter = 100000;
		double e = 0.00000000001;
		q.resize(Nnode);
		double alpha, betta;
		double rk_1_rk_1;
		multiplication_matrix_on_vector(q, y);
		summ(GIGAGLOBALMATRIX->Globalvector, -1, y, r);	
		for (int i = 0; i < Nnode; i++)	
			z[i] = r[i];

		for (int k = 1; k < maxiter; k++)
		{
			rk_1_rk_1 = vector_multiplication(r, r); 
			multiplication_matrix_on_vector(z, y); 
			alpha = rk_1_rk_1 / vector_multiplication(y, z); 
			summ(q, alpha, z, q); 
			multiplication_matrix_on_vector(z, y); 
			summ(r, -alpha, y, r); 
			betta = vector_multiplication(r, r) / rk_1_rk_1; 
			summ(r, betta, z, z);

			if (Norm(r) / Norm(GIGAGLOBALMATRIX->Globalvector) < e) 
			{
				k_iter = k;
				return;
			}
		}
		k_iter = maxiter;

	}
	vector<type> multiplication_matrix_on_vector(vector<type>& a, vector<type>& b)
	{
		for (int i = 0; i < n; i++)
			b[i] = GIGAGLOBALMATRIX->GlobalGMdiag[i] * a[i];

		for (int i = 1; i < n; i++)
		{
			int i0 = GIGAGLOBALMATRIX->ig[i - 1];
			int i1 = GIGAGLOBALMATRIX->ig[i];
			for (int j = 0; j < (i1 - i0); j++)
			{
				b[i] += GIGAGLOBALMATRIX->GlobalGMtriangle[i0 + j] * a[GIGAGLOBALMATRIX->jg[i0 + j]]; 
				b[GIGAGLOBALMATRIX->jg[i0 + j]] += GIGAGLOBALMATRIX->GlobalGMtriangle[i0 + j] * a[i]; 
			}
		}
		return b;
	}
	type Norm(vector<type>& y)
	{
		type norma = 0;
		for (int i = 0; i < n; i++)
			norma += y[i] * y[i];
		return sqrt(norma);
	}

	type vector_multiplication(vector<type>& a, vector<type>& b)
	{
		type s = 0;
		for (int i = 0; i < n; i++)
			s += a[i] * b[i];
		return s;
	}
	vector<type> summ(vector<type>& a, type b, vector<type>& c, vector<type>& d)
	{
		for (int i = 0; i < n; i++)
			d[i] = a[i] + b * c[i];
		return d;
	}

private:
	GlobalMatrix* GIGAGLOBALMATRIX;
	int n;
};
