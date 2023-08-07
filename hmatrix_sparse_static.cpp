using namespace std;

#include "hmatrix.h"

void Ham_1dXXZ::SparseMat_SpecificHeat()
{
	// fixed parameters
	int M = 20;
	double dT = 0.05;
	int nT = 50;
	double one = 1;

	int R = 2000;
	// sampling 
	// compute Z
	double* Z = new double[nT];
	double* ob_H = new double[nT];
	double* ob_H2 = new double[nT];
	double* C = new double[nT];
	for (int iT = 0; iT < nT; iT++)
	{
		Z[iT] = 0;
		ob_H[iT] = 0;
		ob_H2[iT] = 0;
		C[iT] = 0;
	}

	// save all <r,m|H|r,m> in file for test
	double* rspec = new double[R * M];
	// save obseracbles for each realization for test (get mean and standard deviation later)
	double* r_Z = new double[R * nT];
	double* r_H = new double[R * nT];
	double* r_H2 = new double[R * nT];

	for (int r = 0; r < R; r++)
	{
		

		// for a sample, random state, or a random product state?
		double* l_vecs = new double[Dim * (M + 1)];
		double* l_eigvecs = new double[M * M];			// eigenvectors in Lanczos basis
		double* l_eigvals = new double[M];
		
		rd_wf(Dim, l_vecs);
		// how about random product states?
		//rd_wf_product(Dim, l_vecs);

		Lanczos(Dim, SMat_vals, SMat_cols, SMat_PointerBE, M, 'V', l_eigvals, l_vecs, l_eigvecs);

		cblas_dcopy(M, l_eigvals, 1, &rspec[r * M], 1);

		for (int iT = 0; iT < nT; iT++)
		{
			double T = (iT + 1) * dT;
			double beta = 1 / T;
			double* aux = new double[M];
			cblas_dcopy(M, l_eigvals, 1, aux, 1);
			cblas_dscal(M, -beta, aux, 1);
			vdExp(M, aux, aux);

			double* psir2 = new double[M];
			cblas_dcopy(M, &l_eigvecs[0], 1, psir2, 1);
			vdSqr(M, psir2, psir2);
			double out0 = cblas_ddot(M, aux, 1, psir2, 1);
			Z[iT] += out0;
			r_Z[r * nT + iT] = out0;


			vdMul(M, psir2, l_eigvals, psir2);
			out0 = cblas_ddot(M, aux, 1, psir2, 1);
			ob_H[iT] += out0;
			r_H[iT] = out0;
			
			vdMul(M, psir2, l_eigvals, psir2);
			out0 = cblas_ddot(M, aux, 1, psir2, 1);
			ob_H2[iT] += out0;
			r_H2[iT] = out0;

			delete[]aux;
			delete[]psir2;
		}

		delete[]l_vecs;
		delete[]l_eigvecs;
		delete[]l_eigvals;
	}

	Vec_fwrite_double("r_spec.bin", rspec, R*M);
	Vec_fwrite_double("r_Z.bin", r_Z, R * nT);
	Vec_fwrite_double("r_H.bin", r_H, R * nT);
	Vec_fwrite_double("r_H2.bin", r_H2, R * nT);

	for (int iT = 0; iT < nT; iT++)
	{
		double T = (iT + 1) * dT;
		Z[iT] = Z[iT] * Dim / R;
		ob_H[iT] = ob_H[iT] * Dim / R / Z[iT];
		ob_H2[iT] = ob_H2[iT] * Dim / R / Z[iT];
		C[iT] = (ob_H2[iT] - ob_H[iT] * ob_H[iT]) / T / T / LatticeSize;
	}

	ofstream ofc("Sparse_T_lnZ_H_H2_C.dat");
	for (int iT = 0; iT < nT; iT++)
	{
		double T = (iT + 1) * dT;
		ofc << setw(8) << setprecision(4) << T;
		ofc << setw(24) << setprecision(14) << log(Z[iT]);
		ofc << setw(24) << setprecision(14) << ob_H[iT];
		ofc << setw(24) << setprecision(14) << ob_H2[iT];
		ofc << setw(24) << setprecision(14) << C[iT] << endl;
	}
	ofc.close();

	delete[]Z;
	delete[]C;
	delete[]ob_H;
	delete[]ob_H2;

	delete[]r_Z;
	delete[]r_H;
	delete[]r_H2;
	delete[]rspec;


// random state test
	sparse_matrix_t A;
	mkl_sparse_d_create_csr(&A, SPARSE_INDEX_BASE_ZERO, Dim, Dim, SMat_PointerBE, SMat_PointerBE + 1, SMat_cols, SMat_vals);
	struct matrix_descr descrA;
	descrA.type = SPARSE_MATRIX_TYPE_SYMMETRIC;
	descrA.mode = SPARSE_FILL_MODE_UPPER;
	descrA.diag = SPARSE_DIAG_NON_UNIT;
	mkl_sparse_optimize(A);

	int nrd = 10000;
	double* rd_spec = new double[nrd];
	
	double* wf = new double[Dim];
	double* wf1 = new double[Dim];
	for (int ir = 0; ir < nrd; ir++)
	{
		rd_wf(Dim, wf);
		mkl_sparse_d_mv(SPARSE_OPERATION_NON_TRANSPOSE, 1, A, descrA, wf, 0, wf1);
		rd_spec[ir] = cblas_ddot(Dim, wf, 1, wf1, 1);
	}
	Vec_fwrite_double("rd_spec.bin", rd_spec, nrd);

	// how about random product states?
	for (int ir = 0; ir < nrd; ir++)
	{
		rd_wf_product(Dim, wf);
		mkl_sparse_d_mv(SPARSE_OPERATION_NON_TRANSPOSE, 1, A, descrA, wf, 0, wf1);
		rd_spec[ir] = cblas_ddot(Dim, wf, 1, wf1, 1);
	}
	Vec_fwrite_double("rd_spec_product.bin", rd_spec, nrd);

	delete[]wf;
	delete[]wf1;
	//
	delete[]rd_spec;
}