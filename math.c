#include "angular_power_spectrum.h"

/// Calculate the Legendre polynomial P_n(x) using the recurrence relation.
double Legendre(double x, int n)
{
  int i;
  double hold1, hold2, new;

  hold2 = 1.0;
  hold1 = x;
  if (n==0) return hold2;
  else if (n==1) return hold1;
  else {
    for (i=2; i<=n; i++){
      new = ((2.0*i-1.0)*x*hold1 - (i-1.0)*hold2)/(1.0*i);
      hold2 = hold1;
      hold1 = new;
    }
    return new;
  }

  printf("#Error in Legendre polynomial calculation.\n");
  return 0.0;
}


/** Calculate and return the trace of the matrix multiplication AB, where A is N x P and B is P x N.
  * NOTE: This is for matrices in column major order, which is the standard in Fortran not C. */
double trace_multiply(double *A, double *B, int N, int P)
{
  int i, j;
  double sum, trace;

  trace = 0.0;
  for (i=0; i<N; i++){
    sum = 0.0;
    for (j=0; j<P; j++){
      sum += A[i+j*N]*B[j+i*P];
    }
    trace += sum;
  }

  return trace;
}

/// In place inversion of a square size by size matrix
double invert_matrix(double *A, long size)
{
  int M, N, LDA, order, INFO, LWORK, *IPIV;
  double inversetime, *WORK;
  time_t t0, t1;

  order = M = N = LDA = LWORK = size;

  IPIV = (int *)malloc(N*sizeof(int));
  WORK = (double *)malloc(LWORK*LWORK*sizeof(double));

  time(&t0);
  dgetrf(&M, &N, A, &LDA, IPIV, &INFO);
  assert(INFO==0 && fflush(NULL)==0);
  dgetri(&order, A, &LDA, IPIV, WORK, &LWORK, &INFO);
  assert(INFO==0 && fflush(NULL)==0);
  time(&t1);
  inversetime = difftime(t1, t0);
  printf("#Done calculating inverse matrix. Time = %lf seconds.\n", inversetime);

  return inversetime;
}

/// Multiplication of two matrices of size (MxN) = (MxK)(KxN): C = AB
double multiply_matrices(double *A, double *B, double *C, long m, long k, long n)
{
  int LDA, LDB, LDC, M, N, K;
  double ALPHA, BETA, multiplicationtime;
  char TRANSA, TRANSB;
  time_t t0, t1;

  LDA = M = m;
  LDB = K = k;
  LDC = M = m;
  N = n;
  ALPHA = 1.0;
  BETA = 0.0;
  TRANSA = TRANSB = 'N';

  time(&t0);
  dgemm(&TRANSA, &TRANSB, &M, &N, &K, &ALPHA, A, &LDA, B, &LDB, &BETA, C, &LDC);
  time(&t1);
  multiplicationtime = difftime(t1, t0);
  printf("#Done calculating matrix multiplication. Time = %lf seconds.\n", multiplicationtime);

  return multiplicationtime;
}

/// Transpose m x n matrix A into n x m matrix B
void matrix_transpose(double *A, double *B, long m, long n)
{
  int i, j;

  for (i=0; i<m; i++) {  
    for (j=0; j<n; j++) {
      B[j+i*n] = A[i+j*m];
    }
  }

  return;
}

void matrix_square_root(double *A, double *B, long size)
     /* Calculate the square root of a diagonalizable, real, symmetric square size by size matrix A and store it in B */
{
  int INFO, LDA, N, LWORK;
  long i;
  double *W, *WORK, *inverse, *save, *C, *D;
  char UPLO, JOBZ;

  N = LDA = size;
  LWORK = 3*N - 1;
  JOBZ = 'V';  // 'N' for eigenvalues only, 'V' for eigenvalues and eigenvectors.
  UPLO = 'U';  // 'U' for upper triangle of B, 'L' for lower triangle of B.

  W = (double *)malloc(N*sizeof(double));
  WORK = (double *)malloc(LWORK*sizeof(double));
  inverse = (double *)malloc(size*size*sizeof(double));
  save = (double *)malloc(size*size*sizeof(double));
  C = (double *)malloc(size*size*sizeof(double));
  D = (double *)malloc(size*size*sizeof(double));

  for (i=0; i<size*size; i++) {
    save[i] = A[i];
    B[i] = 0.0;
  }
  
  dsyev(&JOBZ, &UPLO, &N, save, &LDA, W, WORK, &LWORK, &INFO);
  assert(INFO==0);

  for (i=0; i<size*size; i++) {
    inverse[i] = save[i];
  }

  invert_matrix(inverse, size);

  multiply_matrices(inverse, A, C, size, size, size);
  multiply_matrices(C, save, D, size, size, size);

  for (i=0; i<size; i++) {
    B[i+i*size] = sqrt(D[i+i*size]);
  }

  multiply_matrices(save, B, C, size, size, size);
  multiply_matrices(C, inverse, B, size, size, size);
  
  return;
}
