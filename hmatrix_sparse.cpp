using namespace std;
#include <omp.h>

#include "hmatrix.h"

l_int Ham_1dXXZ::SparaseMat_Counts()
{
	l_int counts = 0;
	int L1 = LatticeSize - bc;
	for (l_int k = 0; k < Dim; k++) {
		l_int Num_k = basis->get_state(k);
		// diagonal, always considered as non-zero
		counts++;
		// off-diagonal
		for (int i = 0; i < LatticeSize - 1; i++)
		{
			int j = i + 1;
			int k_i = ((Num_k >> i) & 1);
			int k_j = ((Num_k >> j) & 1);
			if (1 == k_i && 0 == k_j)
			{
				counts++;
			}
		}
		if (0 == bc)
		{
			int i = LatticeSize - 1, j = 0;
			int k_i = ((Num_k >> i) & 1);
			int k_j = ((Num_k >> j) & 1);
			if (0 == k_i && 1 == k_j)
			{
				counts++;
			}
		}
	}
	return counts;
}

// build sparse matrix of H 
void Ham_1dXXZ::SparseMat_Build()
{
	// get the number of nonzero elements
	l_int nnz = SparaseMat_Counts();

	SMat_vals = new double[nnz];
	SMat_cols = new l_int[nnz];
	SMat_PointerBE = new l_int[Dim + 1];

	l_int counts = 0;
	int L1 = LatticeSize - bc;
	for (l_int k = 0; k < Dim; k++) {
		l_int Num_k = basis->get_state(k);
		// diagonal 
		double aux_diag = 0;
		for (int i = 0; i < L1; i++)
		{
			int j = (i + 1) % LatticeSize;
			int k_i = ((Num_k >> i) & 1);
			int k_j = ((Num_k >> j) & 1);
			aux_diag += (k_i ^ k_j) ? -0.25 * g : 0.25 * g;
		}
		SMat_vals[counts] = aux_diag;
		SMat_cols[counts] = k;
		SMat_PointerBE[k] = counts;
		counts++;
		// off-diagonal
		vector <pair<l_int, double> > nonzero_pair;
		for (int i = 0; i < LatticeSize - 1; i++)
		{
			int j = i + 1;
			int k_i = ((Num_k >> i) & 1);
			int k_j = ((Num_k >> j) & 1);
			if (1 == k_i && 0 == k_j)
			{
				// then flip the states on i and j 
				l_int Num_l = Num_k ^ (1 << i) ^ (1 << j);
				l_int state_l = basis->get_index(Num_l);
				SMat_vals[counts] = 0.5;
				SMat_cols[counts] = state_l;
				counts++;
				//pair<l_int, double> nonzero_aux = make_pair(state_l, 0.5);
				//nonzero_pair.push_back(nonzero_aux);
			}
		}
		if (0 == bc)
		{
			int i = LatticeSize - 1, j = 0;
			int k_i = ((Num_k >> i) & 1);
			int k_j = ((Num_k >> j) & 1);
			if (0 == k_i && 1 == k_j)
			{
				// then flip the states on i and j
				l_int Num_l = Num_k ^ (1 << i) ^ (1 << j);
				l_int state_l = basis->get_index(Num_l);
				SMat_vals[counts] = 0.5;
				SMat_cols[counts] = state_l;
				counts++;
				//pair<l_int, double> nonzero_aux = make_pair(state_l, 0.5);
				//nonzero_pair.push_back(nonzero_aux);
			}
		}
		SMat_PointerBE[k + 1] = counts;
	}
}

void Ham_1dXXZ::SparseMat_Eigs()
{
	// 'V'/'I' with/without eigenvectors
	char job_vec = 'V';

	int M = 100;

	if ('N' == job_vec)
	{
		double* l_vecs = new double[Dim];
		rd_wf(Dim, l_vecs);
		//uni_wf(Dim, l_vecs);
		double* alpha = new double[M];
		Lanczos(Dim, SMat_vals, SMat_cols, SMat_PointerBE, M, job_vec, alpha, l_vecs, NULL);
		cout << "M = " << M << ", E[0] = " << setprecision(14) << alpha[0] << endl;
		delete[]alpha;
		delete[]l_vecs;
	}
	else if ('V' == job_vec)
	{
		double* l_vecs = new double[Dim * (M + 1)];		// Lanczos vectors
		double* l_eigvecs = new double[M * M];			// eigenvectors in Lanczos basis
		rd_wf(Dim, l_vecs);
		//uni_wf(Dim, l_vecs);
		double* alpha = new double[M];
		Lanczos(Dim, SMat_vals, SMat_cols, SMat_PointerBE, M, job_vec, alpha, l_vecs, l_eigvecs);
		cout << "M = " << M << ", E[0] = " << setprecision(14) << alpha[0] << endl;
		delete[]alpha;

		// groundstate wavefunction
		double* wf0 = new double[Dim];
		for (int i = 0; i < Dim; i++)
		{
			wf0[i] = cblas_ddot(M, &l_eigvecs[0], M, &l_vecs[i], Dim);
		}
		// check wf0
		//cout << "<wf0|U[:,0]> = " << cblas_ddot(Dim, wf0, 1, &DMat[0], Dim) << endl;
		delete[]wf0;

		delete[]l_vecs;
	}
	else
	{
		cout << "Error! job_vec for Lanczos method must be either \'N\' or \'V\'!" << endl;
		exit(666);
	}
}