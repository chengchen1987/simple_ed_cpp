#pragma once

#include <cstdlib>
#include <cstdio>
#include <iostream> 
#include <fstream>
#include <iomanip>
#include <cmath>
#include <cstring>
#include <vector>
#include <algorithm>

#include "basis.h" 

#define PI 3.14159265358979323846

class Ham_1dXXZ
{
private:
	Basis* basis;

public:
	Ham_1dXXZ(Basis* _basis, const double& _g, const int& _bc);
	~Ham_1dXXZ();

	int LatticeSize;
	l_int Dim;

	int bc;	//boundary condition: 0->pbc; 1->obc;
	double g;

	// Build Hamiltonian matrix for 
	// 1d XXZ model
    void DenseMat_Build();
	double* DMat;
	double* spec;

	void DenseMat_Eig();

	void DenseMat_SpecificHeat();

	// sparse matrix 
	l_int SparaseMat_Counts();
	void SparseMat_Build();
	double* SMat_vals;
	l_int* SMat_cols;
	l_int* SMat_PointerBE;

	void SparseMat_Eigs();

	void SparseMat_SpecificHeat();	// finite temperature Lanczos method 
// debug
	//inline void PrintHMatrix(int aindex);
};