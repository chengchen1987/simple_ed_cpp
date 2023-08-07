#include "common.h"

#include <iostream>
#include <iomanip>
#include <random>
std::random_device rd;  //Will be used to obtain a seed for the random number engine
std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
std::uniform_real_distribution<> uni_rd(0.0, 1.0); // call uni_rd to get a random number in [0,1)


//#include <vector>
using namespace std;

//
void rd_wf(const l_int& dim, double* wf)
{
	for (l_int i = 0; i < dim; i++)
	{
		wf[i] = (uni_rd(gen) - 0.5);
	}
	// normalize
	cblas_dscal(dim, 1.0 / cblas_dnrm2(dim, wf, 1), wf, 1);
}

void rd_wf_product(const l_int& dim, double* wf)
{
	for (l_int i = 0; i < dim; i++)
	{
		wf[i] = 0;
	}
	l_int idx = l_int(uni_rd(gen) * dim);
	wf[idx] = 1.0;
	//cout << "norm of random product states"  << cblas_dnrm2(dim, wf, 1) << endl;
}

void uni_wf(const l_int& dim, double* wf)
{
	for (l_int i = 0; i < dim; i++)
	{
		wf[i] = 1;
	}
	// normalize
	cblas_dscal(dim, 1.0 / cblas_dnrm2(dim, wf, 1), wf, 1);
}

// bitwise operations
int numOfBit1(const l_int& b)
{
	l_int a = b;
	int cnt = 0;
	while (a != 0)
	{
		++cnt;
		a &= (a - 1);
	}
	return cnt;
}

void print_binary(const l_int& a, const int& n)
{
	cout << "  ";
	for (int ix = 0; ix < n; ix++) cout << ((a >> (n - ix - 1)) & 1);
	cout << "  ";
}

l_int nchoosek(const int& n, const int& _k)
{
	if (_k > n) return 0;
	int k = _k < (n - _k) ? _k : (n - _k);
	if (0 == k) return 1;
	if (1 == k) return n;
	double aux = 1.0;
	for (int i = 0; i < k; i++)
	{
		aux *= double(n - i) / double(k - i);
	}
	return (l_int)(aux + 1e-2);
}

void Debug_Check_Conj_dense(double* mat, l_int dim)
{
	double aux = 0;
	for (l_int i = 0; i < dim; i++)
	{
		for (l_int j = i + 1; j < dim; j++)
		{
			aux += abs(mat[i * dim + j] - mat[j * dim + i]);
		}
	}
	cout << "abs(H - H^T) = " << aux << endl;
}

void Vec_fwrite_double(const char* fname, double* data, const int& dsize)
{
	FILE* f_out;
	f_out = fopen(fname, "wb");
	fwrite(data, sizeof(double), dsize, f_out);
	fclose(f_out);
}

void Vec_fread_double(const char* fname, double* data, const int& dsize)
{
	FILE* f_in;
	f_in = fopen(fname, "rb");
	fread(data, sizeof(double), dsize, f_in);
	fclose(f_in);
}


// ===============================================================================
// Matrix evd (use mkl lapacke function) 
void DenseMatrixEigenSolver(int matrix_layout, char jobz, char uplo, lapack_int n, double* a, lapack_int lda, double* w)
{
	int eigtype = 1;	// 0, dsyevd; 1, dsyev; 2 dsyevr;
	if (0 == eigtype)
	{
		LAPACKE_dsyev(matrix_layout, jobz, uplo, n, a, lda, w);
	}
	if (1 == eigtype)
	{
		LAPACKE_dsyevd(matrix_layout, jobz, uplo, n, a, lda, w);
	}
	if (2 == eigtype)
	{
		char range = 'A'; // all eigenvalues
		lapack_int* isuppz = new lapack_int[2 * n];
		double abstol = 0;
		LAPACKE_dsyevr(matrix_layout, jobz, range, uplo,
			n, a, lda, NULL, NULL, NULL, NULL,
			abstol, &n, w, a, n, isuppz);
	}
}

void DenseMatrixEigenSolver_FInterface(int matrix_layout, char jobz, char uplo, lapack_int n, double* a, lapack_int lda, double* w)
{
	int eigtype = 0;	// 0, dsyevd; 1, dsyev; 2 dsyevr;
	lapack_int info;
	// dsyev
	if (0 == eigtype)
	{
		lapack_int lwork;
		double wkopt;
		double* work;
		lwork = -1;
		dsyev(&jobz, &uplo, &n, a, &lda, w, &wkopt, &lwork, &info);
		lwork = (int)wkopt;
		cout << "lwork = " << lwork;
		work = (double*)malloc(lwork * sizeof(double));
		// Solve eigenproblem
		dsyev(&jobz, &uplo, &n, a, &lda, w, work, &lwork, &info);
	}
	// dsyevd
	if (1 == eigtype)
	{
		lapack_int lwork, liwork;
		lapack_int* iwork;
		double* work;
		/*
			int iwkopt;
			double wkopt;
			lwork = -1;
			liwork = -1;
			dsyevd(&jobz, &uplo, &n, a, &lda, w, &wkopt, &lwork, &iwkopt, &liwork, &info );
			lwork = (int)wkopt;
			liwork = iwkopt;
		*/
		lwork = 2 * n * n + 6 * n + 1;
		work = (double*)malloc(lwork * sizeof(double));
		liwork = 5 * n + 3;
		iwork = (lapack_int*)malloc(liwork * sizeof(int));
		// Solve eigenproblem
		dsyevd(&jobz, &uplo, &n, a, &lda, w, work, &lwork, iwork, &liwork, &info);
	}

	// dsyevr 
	if (2 == eigtype)
	{
		lapack_int il, iu, ldz, info, lwork, liwork;
		ldz = lda;
		double abstol, vl, vu;
		lapack_int iwkopt;
		lapack_int* iwork;
		double wkopt;
		double* work;
		lapack_int* isuppz = new lapack_int[n];

		abstol = -1;
		il = 1;
		iu = n;
		lwork = -1;
		liwork = -1;
		dsyevr("Vectors", "Indices", "Upper", &n, a, &lda, &vl, &vu, &il, &iu,
			&abstol, &n, w, a, &ldz, isuppz, &wkopt, &lwork, &iwkopt, &liwork,
			&info);
		lwork = (int)wkopt;
		work = (double*)malloc(lwork * sizeof(double));
		liwork = iwkopt;
		iwork = (lapack_int*)malloc(liwork * sizeof(int));
		cout << "lwork = " << lwork << ", liwork = " << liwork << endl;
		// Solve eigenproblem
		dsyevr("Vectors", "Indices", "Upper", &n, a, &lda, &vl, &vu, &il, &iu,
			&abstol, &n, w, a, &ldz, isuppz, work, &lwork, iwork, &liwork,
			&info);
	}

	/* Check for convergence */
	if (info > 0) {
		printf("The algorithm failed to compute eigenvalues.\n");
		exit(1);
	}
}

void MatrixSvd(int matrix_layout, char jobu, char jobvt, lapack_int m, lapack_int n, double* a, lapack_int lda, double* s, double* u, lapack_int ldu, double* vt, lapack_int ldvt) {
	double* superb = new double[m < n ? m : n];
	LAPACKE_dgesvd(matrix_layout, jobu, jobvt, m, n, a, lda, s, u, ldu, vt, ldvt, superb);
}

void Lanczos(const l_int& dim, double* vals, l_int* cols, l_int* PointerBE, const l_int& M, char job_vec, double* alpha, double* l_vecs, double* l_eigvecs)
{
	sparse_matrix_t A;
	mkl_sparse_d_create_csr(&A, SPARSE_INDEX_BASE_ZERO, dim, dim, PointerBE, PointerBE + 1, cols, vals);
	struct matrix_descr descrA;
	descrA.type = SPARSE_MATRIX_TYPE_SYMMETRIC;
	descrA.mode = SPARSE_FILL_MODE_UPPER;
	descrA.diag = SPARSE_DIAG_NON_UNIT;
	mkl_sparse_optimize(A);
	// 
	int steps = M;
	// Lanczos matrix elements
	double* beta = new double[steps];		// beta[0] is unused
	
	// if we do not need lanzcos vecters, only 3 vectors are used
	if ('N' == job_vec)
	{
		double* phi0 = new double[dim];
		double* phi1 = new double[dim];
		double* phi2 = new double[dim];
		// copy initial vector form 1st dim elements form l_vecs 
		cblas_dcopy(dim, l_vecs, 1, phi0, 1);
		
		// a[0],b[1]
		mkl_sparse_d_mv(SPARSE_OPERATION_NON_TRANSPOSE, 1, A, descrA, phi0, 0, phi1);
		alpha[0] = cblas_ddot(dim, phi0, 1, phi1, 1);
		cblas_daxpy(dim, -alpha[0], phi0, 1, phi1, 1);
		beta[1] = cblas_dnrm2(dim, phi1, 1);
		cblas_dscal(dim, 1.0 / beta[1], phi1, 1);

		//cout << setw(24) << setprecision(14) << alpha[0] << setw(24) << setprecision(14) << beta[1] << endl;

		//a[1 ~ M],b[1 ~ M-1]
		for (int m = 1; m < steps; m++)
		{
			mkl_sparse_d_mv(SPARSE_OPERATION_NON_TRANSPOSE, 1, A, descrA, phi1, 0, phi2);
			alpha[m] = cblas_ddot(dim, phi1, 1, phi2, 1);
			if (m < (steps - 1))
			{
				cblas_daxpy(dim, -alpha[m], phi1, 1, phi2, 1);
				cblas_daxpy(dim, -beta[m], phi0, 1, phi2, 1);
				beta[m + 1] = cblas_dnrm2(dim, phi2, 1);
				cblas_dscal(dim, 1.0 / beta[m + 1], phi2, 1);
				cblas_dcopy(dim, phi1, 1, phi0, 1);
				cblas_dcopy(dim, phi2, 1, phi1, 1);
				
				//cout << setw(24) << setprecision(14) << alpha[m] << setw(24) << setprecision(14) << beta[m+1] << endl;
			}
		}
		LAPACKE_dsteqr(LAPACK_ROW_MAJOR, 'N', steps, alpha, &beta[0]+1, NULL, steps);

		delete[]phi0;
		delete[]phi1;
		delete[]phi2;
	}

	// lanzcos vecters are required 
	if ('V' == job_vec)
	{
		// a[0],b[1]
		mkl_sparse_d_mv(SPARSE_OPERATION_NON_TRANSPOSE, 1, A, descrA, &l_vecs[0], 0, &l_vecs[dim]);
		alpha[0] = cblas_ddot(dim, &l_vecs[0], 1, &l_vecs[dim], 1);
		cblas_daxpy(dim, -alpha[0], &l_vecs[0], 1, &l_vecs[dim], 1);
		beta[1] = cblas_dnrm2(dim, &l_vecs[dim], 1);
		cblas_dscal(dim, 1.0 / beta[1], &l_vecs[dim], 1);

		//cout << setw(24) << setprecision(14) << alpha[0] << setw(24) << setprecision(14) << beta[1] << endl;

		//a[1 ~ M],b[1 ~ M-1]
		for (int m = 1; m < steps; m++)
		{
			mkl_sparse_d_mv(SPARSE_OPERATION_NON_TRANSPOSE, 1, A, descrA, &l_vecs[m*dim], 0, &l_vecs[(m+1) * dim]);
			alpha[m] = cblas_ddot(dim, &l_vecs[m * dim], 1, &l_vecs[(m+1) * dim], 1);
			if (m < (steps - 1))
			{
				cblas_daxpy(dim, -alpha[m], &l_vecs[m * dim], 1, &l_vecs[(m+1) * dim], 1);
				cblas_daxpy(dim, -beta[m], &l_vecs[(m-1) * dim], 1, &l_vecs[(m+1) * dim], 1);
				beta[m + 1] = cblas_dnrm2(dim, &l_vecs[(m+1) * dim], 1);
				cblas_dscal(dim, 1.0 / beta[m + 1], &l_vecs[(m + 1) * dim], 1);

				//cout << setw(24) << setprecision(14) << alpha[m] << setw(24) << setprecision(14) << beta[m+1] << endl;
			}
		}
		LAPACKE_dsteqr(LAPACK_ROW_MAJOR, 'I', steps, alpha, &beta[0] + 1, l_eigvecs, steps);
	}
	delete[]beta;
	mkl_sparse_destroy(A);
}