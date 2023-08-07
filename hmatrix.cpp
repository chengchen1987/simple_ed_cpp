using namespace std;
#include <omp.h>

#include "hmatrix.h"



Ham_1dXXZ::Ham_1dXXZ(Basis* _basis, const double& _g, const int& _bc) :
	basis(_basis),
	LatticeSize(basis->get_L()),
	Dim(basis->get_Dim()),
	g(_g),
	bc(_bc)
{
	cout << "Dim = " << Dim << endl;
}

Ham_1dXXZ::~Ham_1dXXZ()
{
}

// build dense matrix of H
void Ham_1dXXZ::DenseMat_Build()
{
	size_t matlen = (size_t)Dim * (size_t)Dim;
	cout << "Dim*Dim = " << matlen << endl;
	cout << "Estimated memory cost of the dense matrix: " << matlen * 8 / 1e9 << " GB" << endl;

	DMat = new double[matlen];
	spec = new double[Dim];

	for (size_t k = 0; k < matlen; k++) DMat[k] = 0;

	int L1 = LatticeSize - bc;
	//omp_set_num_threads(NUMTHREADS);
 //   #pragma omp parallel for schedule(dynamic)
	/*-------------------------------------------------------------------
	  H = sum_{i < j}J_{ij}[b_i^+ b_j^- + b_i^- b_j^+]
	  + sum_{i} (h_i + dh_i) (b_i^+b_i^- - 1/2)
	  --------------------------------------------------------------------*/
	for (l_int k = 0; k < Dim; k++) {
		l_int Num_k = basis->get_state(k);
		//l_int arrayindex = size_t(k) + size_t(k) * size_t(Dim);
		l_int arrayindex = (size_t)k + (size_t)k * (size_t)Dim;
		for (int i = 0; i < L1; i++)
		{
			int j = (i + 1) % LatticeSize;
			int k_i = ((Num_k >> i) & 1);
			int k_j = ((Num_k >> j) & 1);
			// diagonal term ----------------------------------------------	
			DMat[arrayindex] += g * (k_i - 0.5) * (k_j - 0.5);
			// off-diagonal term ------------------------------------------
			if (k_i != k_j)
			{
				// then flip the states on i and j 
				l_int Num_l = Num_k ^ (1 << i) ^ (1 << j);
				l_int state_l = basis->get_index(Num_l);
				//DMat[size_t(k) * size_t(Dim) + size_t(state_l)] += 0.5;
				DMat[(size_t)k * (size_t)Dim + (size_t)state_l] += 0.5;
			}
		}
	}
}

void Ham_1dXXZ::DenseMat_Eig()
{
	// 'V'/'N' with/without eigenvectors
	char job_vec = 'V';
	DenseMatrixEigenSolver(LAPACK_ROW_MAJOR, job_vec, 'U', Dim, DMat, Dim, spec);
}


// 
void Ham_1dXXZ::DenseMat_SpecificHeat()
{
	double dT = 0.05;
	int nT = 50;
	double one = 1;

	ofstream ofc("Dense_T_lnZ_H_H2_C.dat");

	for (int iT = 0; iT < nT; iT++)
	{
		double T = (iT + 1) * dT;
		double beta = 1 / T;
		double* aux = new double[Dim];
		cblas_dcopy(Dim, spec, 1, aux, 1);
		cblas_dscal(Dim, -beta, aux, 1);
		vdExp(Dim, aux, aux);
		double Z = cblas_ddot(Dim, aux, 1, &one, 0);
		double ob_H = cblas_ddot(Dim, aux, 1, spec, 1) / Z;
		vdMul(Dim, aux, spec, aux);
		double ob_H2 = cblas_ddot(Dim, aux, 1, spec, 1) / Z;
		double C = (ob_H2 - ob_H * ob_H) / T / T / LatticeSize;
		//
		ofc << setprecision(4) << T << "  ";
		ofc << setprecision(14) << log(Z) << "  ";
		ofc << setprecision(14) << ob_H << "  ";
		ofc << setprecision(14) << ob_H2 << "  ";
		ofc << setprecision(14) << C << endl;
	}

	ofc.close();
}