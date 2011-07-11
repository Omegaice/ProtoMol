#ifndef LAPACK_H_
#define LAPACK_H_

struct Lapack {
	static bool isEnabled();
	static void dgemv(char *transA, int *m, int *n, double *alpha, double *A,
                      int *lda, double *x, int *incx, double *beta, double *Y,
                      int *incY);
	static void dsyev(char *jobz, char *uplo, int *n, double *a, int *lda,
                      double *w, double *work, int *lwork, int *info);
	static void dsyevr(char *jobz, char *range, char *uplo, int *n, double *a,
                       int *lda, double *vl, double *vu, int *il, int *iu,
                       double *abstol, int *m, double *w, double *z,  int *ldz,
                       int *isuppz, double *work, int *lwork, int *iwork,
                       int *liwork, int *info);
	static double dlamch(char *cmach);
	static void dgemm(char *transA, char *transB, int *m, int *n, int *k,
                      double *alpha, double *A, int *lda, double *B, int *ldb,
                      double *beta, double *C, int *l);
	static double ddot(int *n, double *x, int *incx, double *y, int *incy);
	static double dnrm2(int *n, double *x, int *incx);
	static void dpotri(char *transA, int *n, double *A, int *lda, int *info);
	static void dpotrf(char *transA, int *n, double *A, int *lda, int *info);
	static void dposv(char *transA, int *n, int *nrhs, double *a, int *lda,
                      double *b, int *ldb,int *info);
	static void dtrmm(char *sideA, char *ulA, char *transA, char *diagA,
                      int *m, int *n, double *alpha, double *A, int *lda,
                      double *B, int *ldb);
	static void dtrsm(char *sideA, char *ulA, char *transA, char *diagA,
                      int *m, int *n, double *alpha, double *A, int *lda,
                      double *B, int *ldb);
};

#endif /* LAPACK_H_ */
