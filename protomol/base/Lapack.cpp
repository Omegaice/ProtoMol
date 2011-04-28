#include "protomol/base/Lapack.h"

#if defined(HAVE_LAPACK)
	#include <protomol/integrator/hessian/LapackProtomol.h>
#elif defined(HAVE_SIMTK_LAPACK)
	#include <SimTKlapack.h>
#elif defined(HAVE_MKL_LAPACK)
	#include <mkl_blas.h>
	#include <mkl_lapack.h>
#endif

bool Lapack::isEnabled(){
#if defined(HAVE_LAPACK)
	return true;
#elif defined(HAVE_SIMTK_LAPACK)
	return true;
#elif defined(HAVE_MKL_LAPACK)
	return true;
#endif
	return false;
}

void Lapack::dgemv(char *transA, int *m, int *n, double *alpha, double *A, int *lda, double *x, int *incx, double *beta, double *Y, int *incY){
#if defined(HAVE_LAPACK)
	dgemv_(transA, m, n, alpha, A, lda, x, incx, beta, Y, incY );
#elif defined(HAVE_SIMTK_LAPACK)
	dgemv_(*transA, *m, *n, *alpha, A, *lda, x, *incx, *beta, Y, *incY, 1);
#elif defined(HAVE_MKL_LAPACK)
	DGEMV(transA, m, n, alpha, A, lda, x, incx, beta, Y, incY );
#endif
}

void Lapack::dsyev(char *jobz, char *uplo, int *n, double *a, int *lda, double *w, double *work, int *lwork, int *info){
#if defined(HAVE_LAPACK)
	dsyev_(jobz, uplo, n, a, lda, w, work, lwork, info);
#elif defined(HAVE_SIMTK_LAPACK)
	dsyev_(*jobz, *uplo, *n, a, *lda, w, work, *lwork, *info);
#elif defined(HAVE_MKL_LAPACK)
	DSYEV(jobz, uplo, n, a, lda, w, work, lwork, info);
#endif
}

void Lapack::dsyevr(char *jobz, char *range, char *uplo, int *n, double *a, int *lda, double *vl, double *vu, int *il, int *iu, double *abstol, int *m, double *w, double *z,  int *ldz, int *isuppz, double *work, int *lwork, int *iwork, int *liwork, int *info){
#if defined(HAVE_LAPACK)
	dsyevr_(jobz, range, uplo, n, a, lda, vl, vu, il, iu, abstol, m, w, z, ldz, isuppz, work, lwork, iwork, liwork, info);
#elif defined(HAVE_SIMTK_LAPACK)
	dsyevr_( *jobz, *range, *uplo, *n, a, *lda, vl, vu, il, iu, abstol, *m, w, z, *ldz, isuppz, work, *lwork, iwork, liwork, *info, 1, 1, 1);
#elif defined(HAVE_MKL_LAPACK)
	DSYEVR(jobz, range, uplo, n, a, lda, vl, vu, il, iu, abstol, m, w, z, ldz, isuppz, work, lwork, iwork, liwork, info);
#endif
}

double Lapack::dlamch(char *cmach){
#if defined(HAVE_LAPACK)
	return dlamch_(cmach);
#elif defined(HAVE_SIMTK_LAPACK)
	return dlamch_( *cmach, 1 );
#elif defined(HAVE_MKL_LAPACK)
	return DLAMCH(cmach);
#endif
	return 0.0;
}

void Lapack::dgemm(char *transA, char *transB, int *m, int *n, int *k, double *alpha, double *A, int *lda, double *B, int *ldb, double *beta, double *C, int *l){
#if defined(HAVE_LAPACK)
	dgemm_( transA, transB, m, n, k, alpha, A, lda, B, ldb, beta, C, l );
#elif defined(HAVE_SIMTK_LAPACK)
	dgemm_( *transA, *transB, *m, *n, *k, *alpha, A, *lda, B, *ldb, *beta, C, *l );
#elif defined(HAVE_MKL_LAPACK)
	DGEMM( transA, transB, m, n, k, alpha, A, lda, B, ldb, beta, C, l );
#endif
}

double Lapack::ddot(int *n, double *x, int *incx, double *y, int *incy){
#if defined(HAVE_LAPACK)
	return ddot_( n, x, incx, y, incy );
#elif defined(HAVE_SIMTK_LAPACK)
	return ddot_( *n, x, *incx, y, *incy );
#elif defined(HAVE_MKL_LAPACK)
	return DDOT( n, x, incx, y, incy );
#endif
	return 0.0;
}

double Lapack::dnrm2(int *n, double *x, int *incx){
#if defined(HAVE_LAPACK)
	return dnrm2_( n, x, incx );
#elif defined(HAVE_SIMTK_LAPACK)
	return dnrm2_( *n, x, *incx );
#elif defined(HAVE_MKL_LAPACK)
	return DNRM2( n, x, incx );
#endif
	return 0.0;
}

void Lapack::dpotri(char *transA, int *n, double *A, int *lda, int *info){
#if defined(HAVE_LAPACK)
	dpotri_(transA, n, A, lda, info);
#elif defined(HAVE_SIMTK_LAPACK)
#elif defined(HAVE_MKL_LAPACK)
	DPOTRI(transA, n, A, lda, info);
#endif
}

void Lapack::dpotrf(char *transA, int *n, double *A, int *lda, int *info){
#if defined(HAVE_LAPACK)
	dpotrf_(transA, n, A, lda, info);
#elif defined(HAVE_SIMTK_LAPACK)
#elif defined(HAVE_MKL_LAPACK)
	DPOTRF(transA, n, A, lda, info);
#endif
}

void Lapack::dposv(char *transA, int *n, int *nrhs, double *a, int *lda, double *b, int *ldb, int *info){
#if defined(HAVE_LAPACK)
	dposv_(transA, n, nrhs, a, lda, b, ldb, info);
#elif defined(HAVE_SIMTK_LAPACK)
#elif defined(HAVE_MKL_LAPACK)
	DPOSV(transA, n, nrhs, a, lda, b, ldb, info);
#endif
}

void Lapack::dtrmm(char *sideA, char *ulA, char *transA, char *diagA, int *m, int *n, double *alpha, double *A, int *lda, double *B, int *ldb){
#if defined(HAVE_LAPACK)
	dtrmm_(sideA, ulA, transA, diagA, m, n, alpha, A, lda, B, ldb);
#elif defined(HAVE_SIMTK_LAPACK)
#elif defined(HAVE_MKL_LAPACK)
	DTRMM(sideA, ulA, transA, diagA, m, n, alpha, A, lda, B, ldb);
#endif
}

void Lapack::dtrsm(char *sideA, char *ulA, char *transA, char *diagA, int *m, int *n, double *alpha, double *A, int *lda, double *B, int *ldb){
#if defined(HAVE_LAPACK)
	dtrsm_(sideA, ulA, transA, diagA, m, n, alpha, A, lda, B, ldb);
#elif defined(HAVE_SIMTK_LAPACK)
#elif defined(HAVE_MKL_LAPACK)
	DTRSM(sideA, ulA, transA, diagA, m, n, alpha, A, lda, B, ldb);
#endif
}
