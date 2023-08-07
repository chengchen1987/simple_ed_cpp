#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdio>
#include <ctime>
#include <cmath>
using namespace std;
#include "common.h"
#include "basis.h"
#include "hmatrix.h"

int main() {

	cout << "Exact diagonalization for 1D XXZ 1/2-spins." << endl;

	// test run without symmetry
	{
		// read lattice informationi from input files
		int LatticeSize;
		int bc;	//boundary condition: 0->pbc; 1->obc;
		double g;
		ifstream ifin("input.in", ios::in);
		ifin >> LatticeSize;
		ifin >> g;
		ifin >> bc;
		ifin.close();

		cout << "Number of sites: " << LatticeSize << endl;
		cout << "g = " << g << endl;
		cout << "bc = " << bc << endl;

		// dense matrix algorithm, no symmtry
		if (1)
		{
			//
			cout << "===================================================================" << endl;
			cout << "Dense matrix diagonalization, no symmtry" << endl;
			// generate basis
			int time0 = time(0);
			Basis basis(LatticeSize);
			int time1 = time(0);
			cout << "Time(s) to Build Basis: " << time1 - time0 << endl << endl;

			Ham_1dXXZ ham(&basis, g, bc);
			ham.DenseMat_Build();
			int time2 = time(0);
			cout << "Time(s) to Build DenseMat: " << time2 - time1 << endl << endl;

			ham.DenseMat_Eig();
			int time3 = time(0);
			cout << "Time(s) to  Diagonalize DenseMat: " << time3 - time2 << endl << endl;
			cout << "spec[0]: " << setprecision(14) << ham.spec[0] << endl << endl;

			// compute specific heat
			ham.DenseMat_SpecificHeat();
			Vec_fwrite_double("spec.bin", ham.spec, ham.Dim);
		}

		// dense matrix algorithm, conserved S^z_tot = 0
		if (1)
		{
			//
			cout << "===================================================================" << endl;
			cout << "Dense matrix diagonalization, conserved S^z_tot = 0" << endl;
			int nop = LatticeSize / 2;
			// generate basis
			int time0 = time(0);
			Basis basis(LatticeSize, nop);
			int time1 = time(0);
			cout << "Time(s) to Build Basis: " << time1 - time0 << endl << endl;

			Ham_1dXXZ ham(&basis, g, bc);
			ham.DenseMat_Build();
			int time2 = time(0);
			cout << "Time(s) to Build DenseMat: " << time2 - time1 << endl << endl;

			ham.DenseMat_Eig();
			int time3 = time(0);
			cout << "Time(s) to  Diagonalize DenseMat: " << time3 - time2 << endl << endl;
			cout << "spec[0]: " << setprecision(14) << ham.spec[0] << endl << endl;

			// compute specific heat
			//ham.DenseMat_SpecificHeat();
		}

		// sparse matrix algorithm, no symmtry
		if (1)
		{
			//
			cout << "===================================================================" << endl;
			cout << "Sparse matrix diagonalization, no symmtry" << endl;
			// generate basis
			int time0 = time(0);
			Basis basis(LatticeSize);
			int time1 = time(0);
			cout << "Time(s) to Build Basis: " << time1 - time0 << endl << endl;

			Ham_1dXXZ ham(&basis, g, bc);
			ham.SparseMat_Build();
			int time2 = time(0);
			cout << "Time(s) to Build SparseMat: " << time2 - time1 << endl << endl;

			ham.SparseMat_Eigs();
			int time3 = time(0);
			cout << "Time(s) to  Diagonalize SparseMat: " << time3 - time2 << endl << endl;

			ham.SparseMat_SpecificHeat();
			int time4 = time(0);
			cout << "Time(s) to do FTLM: " << time4 - time3 << endl << endl;
		}

		// sparse matrix algorithm, conserved S^z_tot = 0
		if (1)
		{
			//
			cout << "===================================================================" << endl;
			cout << "Sparse matrix diagonalization, conserved S^z_tot = 0" << endl;
			int nop = LatticeSize / 2;
			// generate basis
			int time0 = time(0);
			Basis basis(LatticeSize,nop);
			int time1 = time(0);
			cout << "Time(s) to Build Basis: " << time1 - time0 << endl << endl;

			Ham_1dXXZ ham(&basis, g, bc);
			ham.SparseMat_Build();
			int time2 = time(0);
			cout << "Time(s) to Build SparseMat: " << time2 - time1 << endl << endl;

			ham.SparseMat_Eigs();
			int time3 = time(0);
			cout << "Time(s) to  Diagonalize SparseMat: " << time3 - time2 << endl << endl;
		}

	}
	return 0;
}
