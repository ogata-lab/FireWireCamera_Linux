
//
// MTX::Matrix<> class template
//
// Jaeil Choi
// last modified in Aug, 2006
//
// matrix multiplication
//

#ifndef MATRIX_MULT_HPP
#define MATRIX_MULT_HPP

// ===================================================================
// 
// ===================================================================

template <class T>
void MatrixMultAB(T R[], int i, T A[], int k, T B[], int j)
{
  // calculate  R = A * B,  with A(i x k) and B(k x j)
  int r, c, t;
  T   *Arow = A, *Bcol, *Rp = R;
  for (r = 0; r < i; r++) {
    for (c = 0; c < j; c++) {
      Bcol = B + c;		// r-th row of A  .  c-th col of B
      for (*Rp = t = 0; t < k; t++) { *Rp += Arow[t] * Bcol[0];  Bcol += j; }
      Rp++;
    }
    Arow += k;
  }
}

template <class T>
void MatrixMultABt(T R[], int i, T A[], int k, T B[], int j)
{
  // calculate  R = A * B^T,  with A(i x k) and B(j x k)
  int r, c, t;
  T   *Arow = A, *Brow = B, *Rp = R;
  for (r = 0; r < i; r++) {
    for (c = 0; c < j; c++) {
      Brow = B + c * k;	// r-th row of A  .  c-th row of B
      for (*Rp = t = 0; t < k; t++) { *Rp += Arow[t] * Brow[t]; }
      Rp++;
    }
    Arow += k;
  }
}

template <class T>
void MatrixMultAtB(T R[], int i, T A[], int k, T B[], int j)
{
  // calculate  R = A^T * B,  with A(k x i) and B(k x j)
  int r, c, t;
  T   *Acol, *Bcol, *Rp = R;
  for (r = 0; r < i; r++) {
    for (c = 0; c < j; c++) {
      Acol = A + r;
      Bcol = B + c;		// r-th col of A  .  c-th col of B
      for (*Rp = t = 0; t < k; t++) { *Rp += Acol[0] * Bcol[0];  Acol += i; Bcol += j; }
      Rp++;
    }
  }
}

template <class T>
void MatrixMultAtBt(T R[], int i, T A[], int k, T B[], int j)
{
  // calculate  R = A^T * B^T,  with A(k x i) and B(j x k)
  int r, c, t;
  T   *Acol, *Brow = B, *Rp = R;
  for (r = 0; r < i; r++) {
    Brow = B;
    for (c = 0; c < j; c++) {
      Acol = A + r;		// r-th col of A  .  c-th row of B
      for (*Rp = t = 0; t < k; t++) { *Rp += Acol[0] * Brow[t];  Acol += i; }
      Brow += k;
      Rp++;
    }
  }
}

template <class T>
void MatrixMultABAt(T R[], int i, T A[], int j, T B[])
{
  // calculate  R = A * B * A^T,  with A(i x j) B(j x j)
  T *M = (T*) malloc( i * j * sizeof(T) );
  MatrixMultAB ( M, i, A, j, B, j );
  MatrixMultABt( R, i, M, j, A, i );
  free( M );
}

template <class T>
void MatrixMultABC(T R[], int i, T A[], int m, T B[], int n, T C[], int j)
{
  // calculate  R = A * B * C,  with A(i x m) B(m x n) and C(n x j)
  T *M = NULL;
  if ((m+j) > (n+i)) {		// Rij = ( Aim * Bmn ) * Cnj
    M = (T*) malloc( i * n * sizeof(T) );
    MatrixMultAB( M, i, A, m, B, n );
    MatrixMultAB( R, i, M, n, C, j );
  } else {			// Rij = Aim * ( Bmn * Cnj )
    M = (T*) malloc( m * j * sizeof(T) );
    MatrixMultAB( M, m, B, n, C, j );
    MatrixMultAB( R, i, A, m, M, j );
  }
  free( M );
}

template <class T>
void MatrixMultABCt(T R[], int i, T A[], int m, T B[], int n, T C[], int j)
{
  // calculate  R = A * B * C^T,  with A(i x m) B(m x n) and C(j x n)
  T *M = NULL;
  if ((m+j) > (n+i)) {		// Rij = ( Aim * Bmn ) * Cjn
    M = (T*) malloc( i * n * sizeof(T) );
    MatrixMultAB ( M, i, A, m, B, n );
    MatrixMultABt( R, i, M, n, C, j );
  } else {			// Rij = Aim * ( Bmn * Cjn )
    M = (T*) malloc( m * j * sizeof(T) );
    MatrixMultABt( M, m, B, n, C, j );
    MatrixMultAB ( R, i, A, m, M, j );
  }
  free( M );
}

template <class T>
void MatrixMultABtCt(T R[], int i, T A[], int m, T B[], int n, T C[], int j)
{
  // calculate  R = A * B^T * C^T,  with A(i x m) B(n x m) and C(j x n)
  T *M = NULL;
  if ((m+j) > (n+i)) {		// Rij = ( Aim * Bnm ) * Cjn
    M = (T*) malloc( i * n * sizeof(T) );
    MatrixMultABt( M, i, A, m, B, n );
    MatrixMultABt( R, i, M, n, C, j );
  } else {			// Rij = Aim * ( Bnm * Cjn )
    M = (T*) malloc( m * j * sizeof(T) );
    MatrixMultAtBt( M, m, B, n, C, j );
    MatrixMultAB  ( R, i, A, m, M, j );
  }
  free( M );
}

template <class T>
void MatrixMultABtC(T R[], int i, T A[], int m, T B[], int n, T C[], int j)
{
  // calculate  R = A * B^T * C,  with A(i x m) B(n x m) and C(n x j)
  T *M = NULL;
  if ((m+j) > (n+i)) {		// Rij = ( Aim * Bnm ) * Cjn
    M = (T*) malloc( i * n * sizeof(T) );
    MatrixMultABt( M, i, A, m, B, n );
    MatrixMultAB ( R, i, M, n, C, j );
  } else {			// Rij = Aim * ( Bnm * Cjn )
    M = (T*) malloc( m * j * sizeof(T) );
    MatrixMultAtB( M, m, B, n, C, j );
    MatrixMultAB ( R, i, A, m, M, j );
  }
  free( M );
}


#endif // MATRIX_MULT_HPP



#if 0 // =============================================================
#include <iostream>
#include "mtx_matrix.hpp"
#include "mth_matrix_mult.hpp"
using namespace std;
int main(void)
{
  MTX::Matrix<double> A, B, C, D, E, F, X, Y, Z;
  A.assign( 4, 3,
	    11.0, 12.0, 13.0,
	    21.0, 22.0, 23.0,
	    31.0, 32.0, 33.0,
	    41.0, 42.0, 43.0 );
  B.assign( 3, 5,
	    11.0, 12.0, 13.0, 14.0, 15.0,
	    21.0, 22.0, 23.0, 24.0, 25.0, 
	    31.0, 32.0, 33.0, 34.0, 35.0 );
  C.transpose( B );   D.transpose( A );
  X.mult(A, B);       Y.set( 4, 5 );
  MatrixMultAB  (Y.data, 4, A.data, 3, B.data, 5);  if (X!=Y) cout << "Wrong!" << endl;
  MatrixMultABt (Y.data, 4, A.data, 3, C.data, 5);  if (X!=Y) cout << "Wrong!" << endl;
  MatrixMultAtB (Y.data, 4, D.data, 3, B.data, 5);  if (X!=Y) cout << "Wrong!" << endl;
  MatrixMultAtBt(Y.data, 4, D.data, 3, C.data, 5);  if (X!=Y) cout << "Wrong!" << endl;
  for (int i = 0; i < C.nRow * C.nCol; i++) C.data[i] += 1.0;
  X.mult(A, B, C);
  Y.set( 4, 3 );
  MatrixMultABC( Y.data, 4, A.data, 3, B.data, 5, C.data, 3 );
  if (X != Y) cout << "Wrong!" << endl;
  D.mult( B, C );
  X.multABAt( A, D );
  Y.set( 4, 4 );
  MatrixMultABAt( Y.data,  4, A.data, 3, D.data );
  if (X != Y) cout << "Wrong!" << endl;
  return EXIT_SUCCESS;
}
#endif // ============================================================
