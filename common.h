#pragma once
#ifndef COMMON_H
#define COMMON_H

#define _CRT_SECURE_NO_DEPRECATE
#include <cstdint>

//#define l_int int64_t
//#define lapack_int int64_t
//#define lapack_uint uint64_t
//#define MKL_INT int64_t
//#define MKL_UINT uint64_t
#include <mkl.h>
#define l_int MKL_INT

void rd_wf(const l_int& dim, double* wf);
void uni_wf(const l_int& dim, double* wf);
void rd_wf_product(const l_int& dim, double* wf);

// some extra bitwise operations
int numOfBit1(const l_int& a);
void print_binary(const l_int& a, const int& n);

// 
l_int nchoosek(const int& k, const int& n);

void Vec_fwrite_double(const char* fname, double* data, const int& dsize);
void Vec_fread_double(const char* fname, double* data, const int& dsize);
//
void Debug_Check_Conj_dense(double* mat, l_int dim);

// 
void DenseMatrixEigenSolver(int matrix_layout, char jobz, char uplo, lapack_int n, double* a, lapack_int lda, double* w);
void DenseMatrixEigenSolver_FInterface(int matrix_layout, char jobz, char uplo, lapack_int n, double* a, lapack_int lda, double* w);
void MatrixSvd(int matrix_layout, char jobu, char jobvt, lapack_int m, lapack_int n, double* a, lapack_int lda, double* s, double* u, lapack_int ldu, double* vt, lapack_int ldvt);

void Lanczos(const l_int &Dim, double* vals, l_int* cols, l_int* PointerBE, const l_int& M, char job_vec, double* l_vals, double* l_vecs, double* l_eigvecs); // job_vec: ('N')'I' (no) Lanzcos vectors 
//void Lanczos_timestep(const l_int& Dim, double* vals, l_int* cols, l_int* PointerBE, const l_int& M, const double &dt, double *wft);

#endif // !COMMON_H