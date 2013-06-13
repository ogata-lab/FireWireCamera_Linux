
//
// Matrix<> class template
//
// Jaeil Choi
// last modified in Feb, 2008
//

// Usage:
//
// Matrix<float> m, m1, m2, m3, m4, m5;
// float values[] = {1.0, 2.0, 3.0,  2.0, 4.0, 7.0,  1.0, 3.0, 4.0};
// --- setting -------------------------------------------------------
// m.setMemory (array);		// let it use the given memory
// m.set (3, 4, array);		// let it use the given memory
// m.set (3, 4);		// set it as 3 x 4 matrix (reallocates memory, only if necessary)
// m.setIdentity(n);		//  "
// m.setZero(m,n);		//  "
// m1 = m2;			// copy matrix (reallocate memory, only if it's necessary)
// m(2,1) = 0;			// element at row 2, col 1 (index starts from 0)
// m(1) = 0;			// element on the diagonal, or of a vector (Nx1, 1xN)
// --- assigning values ----------------------------------------------
// m.assign(2, 3,  9.0, 8.0, 7.0,  6.0, 5.0, 4.0);     // not 9, 8, ..
// 				// Note that values must be given as double type.
// m.assignIdentity();		
// m.assignZero();		
// m.transpose();		// transpose itself
// m.transpose(m1);		// get transpose of m1
// m.copyValues(m2);		// copy values from another matrix
// m.copyValues(values);	// copy values from array
// m.clearNearZeroValues();
// --- test ----------------------------------------------------------
// if (m.isIdentity()) ;
// if (m.isZero()) ;
// if (m1 == m2) ;		if (m1 != m2) ;
// if (m1 == 0) ;		if (m1 != 0) ;
// --- algebraic calculation -----------------------------------------
// m.add(m1);			// m += m1
// m.add(A, B);			// m = A + B
// m.addWithScale(A, s, B)	// m = A + s * B
// m.add (m1, m2, m, m3, m4);	// m = m1 + m2 + m + m3 + m4
// m.sub(m2);			// m -= m2
// m.sub(A, B);			// m = A - B;
// m.mult(A, B);		// m = A * B;
// m.mult(m1, m2, m3, m4, m);	// m = m1 * m2 * m3 * m4 * m
// m.multVMV(v);		// vt * m * v
// m.multValue(2);		// m *= 2
// m.divValue(2);		// m /= 2
// --- miscellaneous -------------------------------------------------
// void appendMatrix(Matrix<T> &m);
// void changeColumnSize(int ncols);
// int  findMaxColInRow(int row, T *value, int except=-1, int except2=-1);
// int  findMaxRowInCol(int col, T *value, int except=-1, int except2=-1);
// T    getSubMatrixSum(int row, int col, int rsize, int csize);
// void sumOverCol(Matrix<T> &m);
// void sumOverRow(Matrix<T> &m);
// --- vector operations ---------------------------------------------
// value = v1.dot(v2);		// v1.v2   (v1 and v2 are Nx1)
// value = v.QuadraticMult(m);	// vt*m*v  (v is Nx1, m is NxN)
// --- file I/O ------------------------------------------------------
// m.writeFile("matrix_file.mtx");
// m.readFile("matrix_file.mtx");
// --- printing matrix content ---------------------------------------
// A.printInfo("size of matrix A");
// A.printMatrix("%12g ", "matrix A");
// A.printSubMatrix( 0,0, 3,4, "%12g ", "(3x4) submatrix starting at (0,0)" );
// A.printRow( 2, "%12g ", "2nd row    of matrix A");
// A.printCol( 3, "%12g ", "3rd column of matrix A");
// A.printAlongWith( x, b, "%12g ", "A * x = b" );
//

#ifndef MTX_MATRIX_HPP
#define MTX_MATRIX_HPP
#define USE_MTX_MATRIX

#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <cassert>
#include <string.h>
#include <stdlib.h>
#include <stdarg.h>
#define  SMALL_VALUE  (1e-6)
#define  ABSV(a)   ((a) < 0 ? (a) * -1 : (a))

namespace MTX {
  
using namespace std;


// ===================================================================
// Matrix
// ===================================================================

template <class T>
class Matrix {
public:
  int nRow, nCol;
  T   *data;
private:
  bool bMemFixed;

public:
  Matrix() : nRow(0), nCol(0), data(NULL), bMemFixed(false) {}
  Matrix(int row, int col) : nRow(0), nCol(0), data(NULL), bMemFixed(false) { set(row, col); }
  Matrix(int row, int col, bool init) : nRow(0), nCol(0), data(NULL), bMemFixed(false) { set(row, col, init); }
  Matrix(int row, int col, T *values) : nRow(0), nCol(0), data(NULL), bMemFixed(false) { set(row, col, values); }
  ~Matrix() { if (!bMemFixed && data) free(data); }

  // setting --------------------------------------------------------

  void clear(bool bFree=true) {
    if (bFree && !bMemFixed && data) free(data);
    nRow = nCol = 0;
    bMemFixed = false;
    data = NULL;
  }
  void setMemory(T* values) {
    clear(true);
    bMemFixed = true;
    data = values;
  }
  bool set(int row, int col, bool init=false) {
    if (!bMemFixed) {
      if (data && (nRow * nCol) != row * col) {	free(data); data = NULL; }
      if (data == NULL) data = (T*)calloc( row * col, sizeof(T) );
      if (data == NULL) return false;
    }
    nRow = row;    nCol = col;
    if (init && data) memset(data, 0, row*col*sizeof(T));
    return true;
  }
  bool set(int row, int col, T *values, bool init=false) {
    if (values) setMemory(values);
    return set(row, col, init);
  }
  bool setKeepingOldData(int row, int col, bool init=false) {
    int  oldsize = (data ? nRow*nCol : 0);
    if (!bMemFixed) {
      data = (T*)realloc( data, row * col * sizeof(T) );
      if (data == NULL) return false;
    }
    nRow = row;    nCol = col;
    if (init && data) memset(data+oldsize, 0, (row*col-oldsize)*sizeof(T));
    return true;
  }
  inline bool setIdentity(int d = 0) {
    if (d > 0 && set(d, d) == false) return false;
    assignIdentity();
    return true;
  }
  inline bool setZero(int row = 0, int col = 0) {
    if (row > 0 && set(row, col)==false) return false;
    assignZero();
    return true;
  }
  inline T& operator() (int idx) { return data[ idx ]; }
  inline T& operator() (int row, int col) { return data[ row * nCol + col ]; }
  
  // assigning values -----------------------------------------------
  
  bool assign(int row, int col, ... ) {
    // IMPORTANT!! - All values must be given as double types
    //   Incorrect! : M.assign(3, 3,  11,12,0, 21,22,0,  0,0,1);
    //   Correct    : M.assign(3, 3,  11.0, 12.0, 0.0,  21.0, 22.0, 0.0,  0.0, 0.0, 1.0); 
    if (set(row, col) == false) return false;
    va_list valist;
    va_start( valist, col );
    for (int i = 0; i < row*col; i++)
      data[i] = (T)va_arg(valist, double);
    va_end( valist );
    return true;
  }
  void assignZero(void) { memset(data, 0, nRow * nCol * sizeof(T)); }
  void assignOne(void) { for (int i = 0; i < nRow*nCol; i++) data[i] = 1; }
  
  void assignIdentity(void) {
    if (nRow != nCol) return;
    memset(data, 0, nRow * nCol * sizeof(T));
    for (int i = 0; i < nCol; i++) data [ i * nCol + i ] = 1;
  }
  Matrix<T> &transpose(void) {		// transpose itself
    Matrix<T> temp(nCol, nRow);
    for (int i = 0; i < nRow; i++)
      for (int j = 0; j < nCol; j++)
	temp.data[ j * nRow + i ] = data[ i * nCol + j ];
    this->moveDataFrom(temp);
    return *this;
  }
  Matrix<T> &transpose(Matrix<T> &m) {	// get transpose of m
    set(m.nCol, m.nRow);
    for (int i = 0; i < nRow; i++)
      for (int j = 0; j < nCol; j++)
	data[ i * nCol + j ] = m.data[ j * m.nCol + i ];
    return *this;
  }
  void flipElements(bool horizontally) {
    if (horizontally) {
      for (int row=0; row < nRow; row++)
	for (int col=0; col < nCol/2; col++) {
	  int p0 = row * nCol + col, p1 = (row+1) * nCol - col - 1;
	  T  tmp = data[p0];  data[p0] = data[p1];  data[p1] = tmp;
	}
    } else {
      T* tmp = (T*)malloc( nCol * sizeof(T) );
      for (int row=0; tmp && row < nRow/2; row++) {
	int p0 = row * nCol, p1 = (nRow-row-1) * nCol;
	memcpy( tmp,     data+p0, nCol*sizeof(T) );
	memcpy( data+p0, data+p1, nCol*sizeof(T) );
	memcpy( data+p1, tmp,     nCol*sizeof(T) );
      }
      if (tmp) free(tmp);
    }
  }
  void rotateElements(void) {
    int i, total=nRow*nCol, half=nRow*nCol/2;
    for (i = 0; i < half; i++) { 
      T tmp = data[i]; data[i] = data[total-i-1]; data[total-i-1] = tmp; 
    }
  }

  void clearNearZeroValues(T epsilon = SMALL_VALUE) {
    for (int i = 0; i < nRow * nCol; i++)
      if (ABSV(data[i]) < epsilon) data[i] = 0;
  }
  
  T Norm(char order = 'f') {
    T norm = 0, temp;
    int i, j;
    switch (order) {
    case 'f':	// Frobenius norm - square root of the sum of squares of elements
      for (i = 0; i < nRow; i++)
	for (j = 0; j < nCol; j++) norm += data[ i + nCol + j ] * data[ i + nCol + j ];
      norm = (T)sqrt((double)norm);
      break;
    case '1': 	// maximum absolute column sum norm
      for (j = 0; j < nCol; j++) {
	temp = 0;
	for (i = 0; i < nRow; i++) temp += fabs(data[ i + nCol + j ]);
	if (temp > norm) norm = temp;
      }
      break;
    case '2':	// spectral norm - square root of the max eigenvalue of A* A
      norm = 0.0;
      break;
    case 'i':	// maximum absolute column sum norm
      for (i = 0; i < nRow; i++) {
	temp = 0;
	for (j = 0; j < nCol; j++) temp += fabs(data[ i + nCol + j ]);
	if (temp > norm) norm = temp;
      }
      break;
    default:  norm = 0.0;  break;
    }
    return norm;
  }

  void normalize(void) {
    int i, size = nRow * nCol;
    T   sum = 0, len;
    for (i = 0; i < size; i++)  sum += data[i] * data[i];
    len = (T) sqrt(sum);
    if (len == 0) return;  
    for (i = 0; i < size; i++)  data[i] /= len;
  }
  
  // copy -----------------------------------------------------------
public:

  inline void copyValues(T *values) {
    memcpy ( data, values, sizeof(T) * nRow * nCol );
  }
  inline void copyValues(Matrix<T> &m) {
    memcpy ( data, m.data, sizeof(T) * nRow * nCol );
  }
  
  inline Matrix<T> &operator= (const Matrix<T> &m) {
    set(m.nRow, m.nCol);
    memcpy ( data, m.data, sizeof(T) * nRow * nCol );
    return *this;
  }
  inline void operator= (double value) {
    int i, total = nRow * nCol;
    for (i = 0; i < total; i++) data[i] = (T)value;
  }
  
  // boolean test ---------------------------------------------------
public:
  
  bool isZero(bool exact = true) {
    if (exact) {
      for (int i = 0; i < nRow; i++)
	for (int j = 0; j < nCol; j++)
	  if ( data [ i * nCol + j ] != 0 ) return false;
    } else {
      T v;
      for (int i = 0; i < nRow; i++)
	for (int j = 0; j < nCol; j++) {
	  v = data[ i * nCol + j ];
	  if ( v > 0 && v > +SMALL_VALUE ||
	       v < 0 && v < -SMALL_VALUE ) return false;
	}
    }
    return true;
  }
  bool isIdentity(bool exact = true) {
    if (nRow != nCol) return false;
    if (exact) {
      for (int i = 0; i < nRow; i++)
	for (int j = 0; j < nCol; j++)
	  if ( data [ i * nCol + j ] != ((i == j) ? 1 : 0) )
	    return false;
    } else {
      T v;
      for (int i = 0; i < nRow; i++)
	for (int j = 0; j < nCol; j++) {
	  v = data[ i * nCol + j ] - ((i == j) ? 1 : 0);
	  if ( v > 0 && v > +SMALL_VALUE ||
	       v < 0 && v < -SMALL_VALUE ) return false;
	}
    }
    return true;
  }
  bool operator== (Matrix<T> &m) {
    if (nRow != m.nRow || nCol != m.nCol) return false;
    for (int i = 0; i < nRow * nCol; i++) 
      if ( data[ i ] != m.data[ i ] ) return false;
    return true;
  }
  bool operator== (T value) {
    for (int i = 0; i < nRow * nCol; i++)
      if ( data[ i ] != value ) return false;
    return true;
  }
  inline bool operator!= (Matrix<T> &m) { return (!(*this == m)); }
  inline bool operator!= (T value)   { return (!(*this == value)); }
  void enforceSymmetry(void) {
    if (nRow != nCol) return;
    int i, j, pos0, pos1;
    for (i = 0; i < nRow; i++)
      for (j = i+1; j < nCol; j++) {
	pos0 = i * nCol + j;	pos1 = j * nCol + i;
	T value = (data[pos0] + data[pos1]) * 0.5;
	data[pos0] = data[pos1] = value;
      }
  }
  bool isSymmetric(void) {
    if (nRow != nCol || nRow == 0) return false;
    for (int i = 0; i < nRow; i++)
      for (int j = i; j < nCol; j++)
	if (data[i * nCol + j] != data[j * nCol + i]) return false;
    return true;
  }
  bool isPositiveDefinite(void) {
    // A matrix is positive definite when
    //   xt A x > 0  for all nonzero vectors x.
    //   All the eigenvalues of A are bigger than 0.
    //   All the submatrices of A have positive determinants.
    //   All the pivots (without row exchanges) are bigger than 0.
    //   There exists a nonsingular matrix W, such that A = Wt W.
    //   THIS FUNCTION IS INEFFICIENT! USE CHOLESKY DECOMPOSITION INSTEAD.
    int        i, j, k;
    if (nRow != nCol || nRow == 0) return false;
    Matrix<T> SM;
    if (nRow != nCol) return false;
    for (k = 0; k < nRow; k++) {
      SM.set( k + 1, k + 1 );
      for (i = 0; i < k; i++) for (j = 0; j < k; j++) SM(i,j) = (*this)(i,j);
      float det = SM.getDeterminant();
      if (det <= 0) return false;
    }
    return true;
  }
  bool isNegativeDefinite(void) { 
    multValue(-1);  
    bool result = isPositiveDefinite();
    multValue(-1);  
    return result;
  }
  
  // determinant -----------------------------------------------------
public:  
  inline T getDeterminant(void) { return getDeterminant(nRow, data); }
private:
  T getDeterminant(int N, T A[]) {
    // calculate the determinant of NxN matrix 'A'.
    if      (N == 1) return A[0];
    else if (N == 2)     // ad - bc
      return ( A[0] * A[3] - A[1] * A[2] );
    else if (N == 3)
      // +a11a22a33+a12a23a31+a13a21a32-a11a23a32-a12a21a33-a13a22a31
      return ( + A[0] * A[4] * A[8]
	       + A[1] * A[5] * A[6]
	       + A[2] * A[3] * A[7]
	       - A[0] * A[5] * A[7]
	       - A[1] * A[3] * A[8]
	       - A[2] * A[4] * A[6] );
    else return getDeterminantByRecursion(N, A);
  }
  T getDeterminantByRecursion(int N, T A[]) {
    // calculate the determinant of NxN matrix 'A'.
    if (N < 4) return getDeterminant(N, A);
    else {
      // Expansions in cofactors:
      // det(A) = ai1Ai1 + ai2Ai2 + ... + ainAin
      // Aij = ((-1)^(i+j)) det(Mij)
      // The minor Mij is formed by deleting row i and column j of A.
      int k, i, j, sign;
      T   sum = 0;
      Matrix<T> m(N-1, N-1);
      for (k = 0; k < N; k++) {
	if (A[k] == 0) continue;
	sign = ((k % 2) == 0 ? +1 : -1);
	for (i = 0; i < m.nRow; i++)
	  for (j = 0; j < m.nCol; j++)
	    m(i,j) = A[ (i+1) * N + j + ((j>=k) ? 1 : 0) ];
	sum += A[k] * sign * getDeterminantByRecursion(N-1, m.data);
      }
      return sum;
    }
  }

  // basic operations -----------------------------------------------
private:  
  void moveDataFrom(Matrix<T> &m) {
    nRow = m.nRow;
    nCol = m.nCol;
    if (bMemFixed || m.bMemFixed) {	// copy the data from 'm'
      memcpy(data, m.data, nRow * nCol * sizeof(T));
    } else {				// take the data from 'm'
      if (data) free(data);
      data = m.data;	
      m.data = NULL;  
      m.nRow = m.nCol = 0;
    }
  }
  
  // ================================================================
  // Algebraic Calculation
  // ================================================================
public:

  Matrix<T> &add(const Matrix<T> &m) {
    if (nRow != m.nRow || nCol != m.nCol) return *this;
    int i, total = nRow * nCol;
    for (i = 0; i < total; i++) data[i] += m.data[i];
    return *this;
  }
  Matrix<T> &add(const Matrix<T> &m1, const Matrix<T> &m2) {
    // check size of matrices;
    if (m1.nRow != m2.nRow || m1.nCol != m2.nCol) {
      fprintf(stderr, "Error : Matrix add() : (%dx%d) + (%dx%d)\n", m1.nRow, m1.nCol, m2.nRow, m2.nCol);
      return *this;
    }
    if (nRow != m1.nRow || nCol != m1.nCol) set(m1.nRow, m1.nCol);
    // add all the values;
    int i, total = nRow * nCol;
    for (i = 0; i < total; i++)  data[i] = m1.data[i] + m2.data[i];
    return *this;
  }
  Matrix<T> &addWithScale(const Matrix<T> &m1, T s, const Matrix<T> &m2) {
    // check size of matrices;
    if (m1.nRow != m2.nRow || m1.nCol != m2.nCol) {
      fprintf(stderr, "Error : Matrix addWithScale() : (%dx%d) + s*(%dx%d)\n", m1.nRow, m1.nCol, m2.nRow, m2.nCol);
      return *this;
    }
    if (nRow != m1.nRow || nCol != m1.nCol) set(m1.nRow, m1.nCol);
    // add all the values;
    int i, total = nRow * nCol;
    for (i = 0; i < total; i++)  data[i] = m1.data[i] + s * m2.data[i];
    return *this;
  }
  Matrix<T> &add(const Matrix<T> &m1, const Matrix<T> &m2, const Matrix<T> &m3) {
    // check size of matrices;
    if (m1.nRow != m2.nRow || m1.nCol != m2.nCol ||
        m1.nRow != m3.nRow || m1.nCol != m3.nCol ) {
      fprintf(stderr, "Error : Matrix add() : (%dx%d) + (%dx%d) + (%dx%d)\n",
	      m1.nRow, m1.nCol, m2.nRow, m2.nCol, m3.nRow, m3.nCol);
      return *this;
    }
    if (nRow != m1.nRow || nCol != m1.nCol) set(m1.nRow, m1.nCol);
    // add all the values;
    int i, total = nRow * nCol;
    for (i = 0; i < total; i++)  
      data[i] = m1.data[i] + m2.data[i] + m3.data[i];
    return *this;
  }
  Matrix<T> &add(const Matrix<T> &m1, const Matrix<T> &m2,
		 const Matrix<T> &m3, const Matrix<T> &m4) {
    // check size of matrices;
    if (m1.nRow != m2.nRow || m1.nCol != m2.nCol ||
	m1.nRow != m3.nRow || m1.nCol != m3.nCol ||
        m1.nRow != m4.nRow || m1.nCol != m4.nCol ) {
      fprintf(stderr, "Error : Matrix add() : (%dx%d) + (%dx%d) + (%dx%d) + (%dx%d)\n", 
	      m1.nRow, m1.nCol, m2.nRow, m2.nCol, m3.nRow, m3.nCol, m4.nRow, m4.nCol);
      return *this;
    }
    if (nRow != m1.nRow || nCol != m1.nCol) set(m1.nRow, m1.nCol);
    // add all the values;
    int i, total = nRow * nCol;
    for (i = 0; i < total; i++)  
      data[i] = m1.data[i] + m2.data[i] + m3.data[i] + m4.data[i];
    return *this;
  }
  Matrix<T> &add(const Matrix<T> &m1, const Matrix<T> &m2,
		 const Matrix<T> &m3, const Matrix<T> &m4,
		 const Matrix<T> &m5) {
    // check size of matrices;
    if (m1.nRow != m2.nRow || m1.nCol != m2.nCol ||
	m1.nRow != m3.nRow || m1.nCol != m3.nCol ||
        m1.nRow != m4.nRow || m1.nCol != m4.nCol ||
        m1.nRow != m5.nRow || m1.nCol != m5.nCol ) {
      fprintf(stderr, "Error : Matrix add() : (%dx%d) + (%dx%d) + (%dx%d) + (%dx%d) + (%dx%d)\n", 
	      m1.nRow, m1.nCol, m2.nRow, m2.nCol, m3.nRow, m3.nCol, m4.nRow, m4.nCol, m5.nRow, m5.nCol);
      return *this;
    }
    if (nRow != m1.nRow || nCol != m1.nCol) set(m1.nRow, m1.nCol);
    // add all the values;
    int i, total = nRow * nCol;
    for (i = 0; i < total; i++)
      data[i] = ( m1.data[i] + m2.data[i] + m3.data[i] + m4.data[i] + m5.data[i] );
    return *this;
  }
  Matrix<T> &sub(const Matrix<T> &m) {
    if (nRow != m.nRow || nCol != m.nCol) return *this;
    int i, total = nRow * nCol;
    for (i = 0; i < total; i++)  data[i] -= m.data[i];
    return *this;
  }
  Matrix<T> &sub(const Matrix<T> &m1, const Matrix<T> &m2) {
    // check size of matrices;
    if (m1.nRow != m2.nRow || m1.nCol != m2.nCol) {
      fprintf(stderr, "Error : Matrix sub() : (%dx%d) - (%dx%d)\n", m1.nRow, m1.nCol, m2.nRow, m2.nCol);
      return *this;
    }
    if (nRow != m1.nRow || nCol != m1.nCol) set(m1.nRow, m1.nCol);
    int i, total = nRow * nCol;
    for (i = 0; i < total; i++)  data[i] = m1.data[i] - m2.data[i];
    return *this;
  }
  Matrix<T> &multValue(T v) {
    int i, total = nRow * nCol;
    for (i = 0; i < total; i++)  data[i] *= v;
    return *this;
  }
  Matrix<T> &divValue(T v) {
    int i, total = nRow * nCol;
    for (i = 0; i < total; i++)  data[i] /= v;
    return *this;
  }
  Matrix<T> &addWithWeights(T v1, Matrix<T> &m1, T v2, Matrix<T> &m2) {
    // check size of matrices;
    if (m1.nRow != m2.nRow || m1.nCol != m2.nCol) {
      fprintf(stderr, "Error : Matrix addWithWeights() : a*(%dx%d) + b*(%dx%d)\n", m1.nRow, m1.nCol, m2.nRow, m2.nCol);
      return *this;
    }
    if (nRow != m1.nRow || nCol != m1.nCol) set(m1.nRow, m1.nCol);
    int i, total = nRow * nCol;
    for (i = 0; i < total; i++)  data[i] = v1 * m1.data[i] + v2 * m2.data[i];
    return *this;
  }
  
  // matrix multiplication ------------------------------------------
  
  Matrix<T> &mult(const Matrix<T> &m1, const Matrix<T> &m2 ) {
    // Note that " m1.mult(m1, m2); " is possible
    Matrix<T> temp;
    if (m1.nCol == m2.nRow) {
      temp.MultPrivate(m1, m2);
      this->moveDataFrom(temp);
    } else {
      fprintf(stderr, "Error : Matrix mult() : (%dx%d) * (%dx%d)\n", m1.nRow, m1.nCol, m2.nRow, m2.nCol);
    }
    return *this;
  }
  Matrix<T> &mult(const Matrix<T> &m1,
		  const Matrix<T> &m2,
		  const Matrix<T> &m3 ) {
    // Note that " m.mult(m1, m, m2); " is possible
    Matrix<T> temp;
    if (m1.nCol == m2.nRow && m2.nCol == m3.nRow) {
      temp.MultPrivate(m1, m2);
      temp.MultPrivate(m3);
      this->moveDataFrom(temp);
    } else {
      fprintf(stderr, "Error : Matrix mult() : (%dx%d) * (%dx%d) * (%dx%d)\n", 
	      m1.nRow, m1.nCol, m2.nRow, m2.nCol, m3.nRow, m3.nCol);
    }
    return *this;
  }
  Matrix<T> &mult(const Matrix<T> &m1, const Matrix<T> &m2,
		  const Matrix<T> &m3, const Matrix<T> &m4 ) {
    // Note that " m.mult(m1, m2, m, m3); " is possible
    Matrix<T> temp;
    if (m1.nCol == m2.nRow &&
	m2.nCol == m3.nRow && m3.nCol == m4.nRow) {
      temp.MultPrivate(m1, m2);
      temp.MultPrivate(m3);
      temp.MultPrivate(m4);
      this->moveDataFrom(temp);
    } else {
      fprintf(stderr, "Error : Matrix mult() : (%dx%d) * (%dx%d) * (%dx%d) * (%dx%d)\n", 
	      m1.nRow, m1.nCol, m2.nRow, m2.nCol, m3.nRow, m3.nCol, m4.nRow, m4.nCol);
    }
    return *this;
  }
  Matrix<T> &mult(const Matrix<T> &m1, const Matrix<T> &m2,
		  const Matrix<T> &m3, const Matrix<T> &m4,
		  const Matrix<T> &m5 ) {
    // Note that " m.mult(m1, m2, m, m3, m4); " is possible
    Matrix<T> temp;
    if (m1.nCol == m2.nRow && m2.nCol == m3.nRow &&
	m3.nCol == m4.nRow && m4.nCol == m5.nRow) {
      temp.MultPrivate(m1, m2);
      temp.MultPrivate(m3);
      temp.MultPrivate(m4);
      temp.MultPrivate(m5);
      this->moveDataFrom(temp);
    } else {
      fprintf(stderr, "Error : Matrix mult() : (%dx%d) * (%dx%d) * (%dx%d) * (%dx%d) * (%dx%d)\n", 
	      m1.nRow, m1.nCol, m2.nRow, m2.nCol, m3.nRow, m3.nCol, m4.nRow, m4.nCol, m5.nRow, m5.nCol);
    }
    return *this;
  }
  Matrix<T> &multABt(const Matrix<T> &A, const Matrix<T> &B) {
    // A * B' = [ Mik ] = [ Aij B'jk ] = [ Aij Bkj ] (with Einstein summation)
    set( A.nRow, B.nRow );
    int i, j, k;
    T   *r, *a, *b, sum;
    for (i = 0; i < A.nRow; i++) {
      r = this->data + this->nCol * i;		// i-th row of result
      a = A.data + A.nCol * i;			// i-th row of A
      for (k = 0; k < B.nRow; k++) {
	b = B.data + B.nCol * k;		// k-th row of B
	for (sum = j = 0; j < A.nCol; j++) sum += a[j] * b[j];
	r[k] = sum;				// (*this)(i,k) = sum
      }
    }
    return *this;
  }
  Matrix<T> &multAtB(const Matrix<T> &A, const Matrix<T> &B) {
    // A' * B = [ Mik ] = [ A'ij Bjk ] = [ Aji Bjk ] (with Einstein summation)
    set( A.nCol, B.nCol );
    int i, j, k;
    T   *r, *a, *b, sum;
    for (i = 0; i < A.nCol; i++) {
      r = this->data + this->nCol * i;		// i-th row of result
      for (k = 0; k < B.nCol; k++) {
	a = A.data + i;				// i-th col of A
	b = B.data + k;				// k-th col of B
	for (sum=j=0; j < A.nRow; j++, a+=A.nCol, b+=B.nCol) sum += (*a) * (*b);
	r[k] = sum;				// (*this)(i,k) = sum
      }
    }
    return *this;
  }
  Matrix<T> &multABAt(const Matrix<T> &A, const Matrix<T> &B) {
    // A * B * A' = [ Mil] = [ Aij Bjk A'kl ] = [ Aij Bjk Alk ] (with Einstein summation)
    set( A.nRow, A.nRow );
    int i, j, k, l;
    T   *r, *a1, *a2, *b, sum;
    for (i = 0; i < A.nRow; i++) {
      r  = this->data + this->nCol * i;		// i-th row of result
      a1 = A.data + A.nCol * i;			// i-th row of A
      for (l = 0; l < A.nRow; l++) {
	a2 = A.data + A.nCol * l;		// l-th row of A
	for (sum = j = 0; j < A.nCol; j++) {
	  b = B.data + B.nCol * j;		// j-th row of B
	  for (k = 0; k < A.nCol; k++) sum += a1[j] * b[k] * a2[k];
	}
	r[l] = sum;				// (*this)(i,l) = sum
      }
    }
    return *this;
  }
  T multVMV(const Matrix<T> &v) {
    // calculates and return  vt * M * v
    if (v.nCol != 1) { cerr << "Invalid vector for multVMV" << endl; return 0; }
    int i, j;   T result, temp;
    for (result = j = 0; j < nCol; j++) {
      for (temp = i = 0; i < nRow; i++)	temp += v.data[i] * data[i*nCol+j];
      result += temp * v.data[j];
    }
    return result;
  }

private:
  void MultPrivate(const Matrix<T> &m1, const Matrix<T> &m2) {
    // this private function assumes colume/row agreement
    set(m1.nRow, m2.nCol);  // set new matrix
    // calculate new values
    T   sum;
    for (int i = 0; i < nRow; i++) {
      for (int j = 0; j < nCol; j++) {
	sum = 0;
	for (int k = 0; k < m1.nCol; k++)
	  sum += m1.data[i * m1.nCol + k] * m2.data[k * m2.nCol + j];
	data[i * nCol + j] = sum;
      }
    }
  }
  
  void MultPrivate(const Matrix<T> &m) {
    // this private function assumes colume/row agreement
    int row, col, cnt, i, j, k;
    T   *temp, sum;
    // set new matrix
    row = this->nRow;
    col = m.nCol;
    cnt = this->nCol;
    temp = (T *) calloc (row * col, sizeof(T));
    if (temp == NULL) return;
    // calculate new values
    for (i = 0; i < row; i++) {
      for (j = 0; j < col; j++) {
	sum = 0;
	for (k = 0; k < cnt; k++)
	  sum += this->data[i * cnt + k] * m.data[k * col + j];
	temp[i * col + j] = sum;
      }
    }
    // replace with new row,col,values.
    if (!bMemFixed && data) free(data);
    this->nRow = row;
    this->nCol = col;
    this->data = temp;
  }
  
  // ================================================================
  // submatrix
  // ================================================================
  
public:
  // Be careful to use appropriate size of 'm', 'row' or 'col'.
  //   The caller matrix should be bigger than (or equal to) the callee.
  // Otherwise, you will overflow the matrix and may experience glibc errors.
  void addSubMatrixFrom(int row, int col, Matrix<T> &m) { addSubMatrixFrom(row, col, m.nRow, m.nCol, m.data); }
  void addSubMatrixFrom(int row, int col, int rs, int cs, T v[]) {
    // add the values of the array 'v' at specified location
    int i, j, k=0, pos=0;
    for (i = 0; i < rs; i++) {
      pos = (row + i) * this->nCol + col;
      for (j = 0; j < cs; j++)  data[pos++] += v[k++];
    }
  }
  void copySubMatrixTo(int row, int col, Matrix<T> &m) { copySubMatrixTo(row, col, m.nRow, m.nCol, m.data); }
  void copySubMatrixTo(int row, int col, int rs, int cs, T v[]) {
    // copy the submatrix at (row,col) with size (rs, cs) into an array
    int i, j, k=0, rend, cend;
    rend = (row+rs >= nRow ? nRow : row+rs);
    cend = (col+cs >= nCol ? nCol : col+cs);
    for (i = row; i < rend; i++)
      for (j = col; j < cend; j++)  v[k++] = (*this)( i, j );
  }
  void copySubMatrixFrom(int row, int col, Matrix<T> &m) { copySubMatrixFrom(row, col, m.nRow, m.nCol, m.data); }
  void copySubMatrixFrom(int row, int col, int rs, int cs, T v[]) {
    // copy the values of the array 'v(nRow x nCol)' at this->(row,col)
    T   *dst = this->data + row * this->nCol + col;
    T   *src = v;
    for (int i = 0; i < rs; i++) {
      memcpy( dst, src, cs * sizeof(T) );
      dst += this->nCol;
      src += cs;
    }
  }
  void copySubMatrixFromTr(int row, int col, int rs, int cs, T v[]) {
    // copy the values of the transpose of array 'v(nRow x nCol)' at (row,col)
    int i, j, k=0, pos=0;
    for (j = 0; j < cs; j++) {
      pos = row * this->nCol + col + j;	// (col+j)-th column
      for (i = 0; i < rs; i++, pos+=this->nCol)  data[pos] = v[k++];
    }
  }
  void getSubMatrixFrom(Matrix<T> &m, int row, int col, int rs, int cs) {
    // copy the submatrix at (row,col) with size (nRow, nCol) of the big matrix 'm'
    set(rs, cs);
    T *src = m.data + row * m.nCol + col;
    T *dst = this->data;
    for (int i = 0; i < nRow; i++) {
      memcpy( dst, src, nCol*sizeof(T) );
      src += m.nCol;
      dst += nCol;
    }
  }
  
  // ================================================================
  // miscellaneous
  // ================================================================
  
  void appendMatrix(Matrix<T> &m) {
    if (nCol != m.nCol || bMemFixed) return;
    data = (T*)realloc(data, (nRow+m.nRow)*nCol * sizeof(T));
    memcpy (data + nRow*nCol, m.data, m.nRow*m.nCol * sizeof(T));
    nRow = nRow + m.nRow;
  }
  
  void changeColumnSize(int ncols) {
    int total = nRow * nCol;
    if (total <= 0) return;
    if (ncols < 1) ncols = 1;
    if (ncols > total) ncols = total;
    nRow = total / ncols;
    nCol = ncols;  // Note it may lose some entries at the last row.
  }
  
  int  findMaxColInRow(int row, T *value, int except=-1, int except2=-1) {
    // Find and return the column with maximum value in the 'row'.
    if (row < 0 || row >= nRow) return -1;
    int j, max_j=-1;
    T   val, max_val=data[row*nCol]-1;
    for (j = 0; j < nCol; j++) {
      if (j == except || j == except2) continue;
      val = data[row * nCol + j];
      if (val > max_val) { max_val = val; max_j = j; }
    }
    if (max_j >= 0 && value) *value = max_val;
    return max_j;
  }
  int  findMaxRowInCol(int col, T *value, int except=-1, int except2=-1) {
    // Find and return the row with maximum value in the column 'col'.
    if (col < 0 || col >= nCol) return -1;
    int i, max_i=-1;
    T   val, max_val=data[col]-1;
    for (i = 0; i < nRow; i++) {
      if (i == except || i == except2) continue;
      val = data[i * nCol + col];
      if (val > max_val) { max_val = val; max_i = i; }
    }
    if (max_i >= 0 && value) *value = max_val;
    return max_i;
  }
  T  getSubMatrixSum(int row, int col, int rsize, int csize) {
    T sum = 0;
    int rend = row + rsize;  if (rend >= nRow) rend = nRow;
    int cend = col + csize;  if (cend >= nCol) cend = nCol;
    for (int r = row; r < rend; r++)
      for (int c = col; c < cend; c++)
	sum += data[r * nCol + c];
    return sum;
  }
  T  getRowSum(int row) {
    T sum = 0, *pp=data+row*nCol;
    for (int i=0; i<nCol; i++) sum += pp[i];
    return sum;
  }
  void sumOverCol(Matrix<T> &m) {
    // Create a (Mx1) matrix from (MxN), by summing up all the values in the same row.
    set( m.nRow, 1, true);
    for (int i=0,k=0; i < m.nRow; i++)
      for (int j=0; j < m.nCol; j++, k++) data[i] += m.data[k];
  }
  void sumOverRow(Matrix<T> &m) {
    // Create a (1xN) matrix from (MxN), by summing up all the values in the same column.
    set( 1, m.nCol, true );
    for (int i=0,k=0; i < m.nRow; i++)
      for (int j=0; j < m.nCol; j++, k++) data[j] += m.data[k];
  }
  
  // ================================================================
  // vector operation
  // ================================================================
  
 public:
  T dot(Matrix<T> &m) {
    // assuming both matrix is actually nRow x 1 vectors
    assert( nCol == 1 && m.nCol == 1 && nRow == m.nRow );
    T temp = 0;
    for (int i = 0; i < nRow; i++)  temp += data[i] * m.data[i];
    return temp;
  }
  
  T QuadraticMult(Matrix<T> &m) {
    assert ( nCol == 1 && m.nCol == nRow );
    // calculate [v]T [M] [v], assumming 'this' is m.nCol x 1 vector
    Matrix<T> temp; 
    temp.mult(m, *this);
    return this->dot(temp);
  }
  
  // ================================================================
  // file I/O
  // ================================================================
public:
  void writeFile(char *filename) {
    ofstream ofs(filename);
    ofs << "Matrix  " << nRow << " rows  " << nCol << " cols" << endl;
    int i, j, pos;
    for (i = pos = 0; i < nRow; i++) {
      for (j = 0; j < nCol; j++, pos++) {
	if (data[pos] == 0.0) ofs << "\t0\t";
	else ofs << setw(12) << data[ pos ] << "\t";
      }
      ofs << "\n";
    }
  }
  bool checkFile(char *filename) {
    FILE *fp = fopen(filename, "r");
    if (fp == NULL) return false;
    char str[256];  int row, col;
    int  n = fscanf(fp, "%s %d rows %d cols\n", str, &row, &col);
    fclose(fp);
    return (n == 3 && strcmp(str, "Matrix") == 0);
  }
  bool readFile(char *filename) {
    char *ext = filename+strlen(filename)-3;
    if (strcmp(ext, "txt")==0 || strcmp(ext, "dat")==0)
      return readTextFile(filename);
    clear();
    FILE *fp = fopen(filename, "r");
    if (fp == NULL) return false;
    int    row, col, n, i, j, pos;   
    char   str[256];    
    double value;
    n = fscanf(fp, "%s  %d rows  %d cols\n", str, &row, &col);
    if (n != 3 || strcmp(str, "Matrix") != 0) return false;
    set( row, col );
    for (i = pos = 0; i < nRow; i++)
      for (j = 0; j < nCol; j++) {
	fscanf(fp, "%s", str);
	value = ( (str[0] == '-' && str[1] == '\0' || str[0] == 'N' ) ? 0.0 : atof(str) );
	data[pos++] = (T) value;
      }
    fclose(fp);
    return true;
  }
  bool readTextFile(char *filename, int columns=1) {
    clear();
    FILE *fp = fopen(filename, "r");
    if (fp == NULL) return false;
    int    count, i, j, pos;
    double value;
    // count total number of values
    for (count = 0; fscanf(fp, "%lf", &value) > 0; count++);
    if (count <= 0) return false;
    fseek(fp, 0L, SEEK_SET);
    // read the values into the matrix
    set( count/columns, columns );
    for (i = pos = 0; i < nRow; i++)
      for (j = 0; j < nCol; j++) {
	fscanf(fp, "%lf", &value);
	data[pos++] = (T) value;
      }
    fclose(fp);
    return true;
  }

  // ================================================================
  // printing the content
  // ================================================================
  
public:
  void printInfo(char *cmmt=NULL) {
    printf("%s : ( %d x %d )\n", (cmmt ? cmmt : "Matrix<>"), nRow, nCol);
  }
  void printMatrix(char *format="%12g ", char *cmmt=NULL) {
    int i, j;  char indent[31]; 
    if (!cmmt) indent[0]='\0';
    else { for (i=0; i<30 && cmmt[i]==' '; i++) indent[i]=' '; indent[i]='\0'; }
    printf("%s : ( %d x %d )\n", (cmmt ? cmmt : "Matrix<>"), nRow, nCol);
    for (i = 0; i < nRow; i++) {
      printf("%s  [ ", indent);
      for (j = 0; j < nCol; j++) printf(format, data[i*nCol+j]);
      printf("]\n");
    }
  }
  void printSubMatrix(int rs, int cs, int rows, int cols, char *format="%12g ", char *cmmt=NULL) {
    int i, j;  char indent[30]; 
    if (!cmmt) indent[0]='\0';
    else { for (i=0; i<20 && cmmt[i]==' '; i++) indent[i]=' '; indent[i]='\0'; }
    int re = (rows<=0 ? nRow : rs+rows);  if (re>nRow) re = nRow;
    int ce = (cols<=0 ? nCol : cs+cols);  if (ce>nCol) ce = nCol;
    printf("%s : [ %d:%d, %d:%d ] from (%d x %d)\n", 
	   (cmmt ? cmmt : "Matrix<>"), rs, re-1, cs, ce-1, nRow, nCol);
    for (i = rs; i < re; i++) {
      printf("%s  [ ", indent);
      for (j = cs; j < ce; j++) printf(format, data[i*nCol+j]);
      printf("]\n");
    }
  }
  void printRow(int row, char *format="%12g ", char *cmmt=NULL) {
    printf("%s (%2d, 0:%2d): ", (cmmt ? cmmt : "MtxRow"), row, nCol-1);
    for (int col=0; col < nCol; col++)
      printf(format, data[row*nCol+col]);
    printf("\n");
  }
  void printCol(int col, char *format="%12g ", char *cmmt=NULL) {
    printf("%s (0:%2d, %2d): ", (cmmt ? cmmt : "MtxCol"), nRow-1, col);
    for (int row=0; row < nRow; row++)
      printf(format, data[row*nCol+col]);
    printf("\n");
  }
  void printAlongWith(Matrix<T> &A, char *format="%12g ", char *cmmt=NULL) {
    int  i;  char indent[31], eformat[10]; 
    if (!cmmt) indent[0]='\0';
    else { for (i=0; i<30 && cmmt[i]==' '; i++) indent[i]=' '; indent[i]='\0'; }
    int  mRow = (nRow > A.nRow ? nRow : A.nRow);
    if (! getEmptyFormatForPAW(format,eformat)) return;
    printf("%s : ( %d x %d )  ( %d x %d )\n", (cmmt ? cmmt : "Matrix<>"), nRow, nCol, A.nRow, A.nCol);
    for (i = 0; i < mRow; i++) {
      printf("%s  ", indent);
      printRowForPAW( *this, i, format, eformat );
      printRowForPAW( A,     i, format, eformat );
      printf("\n");
    }
  }
  void printAlongWith(Matrix<T> &A, Matrix<T> &B, char *format="%12g ", char *cmmt=NULL) {
    int  i;  char indent[31], eformat[10]; 
    if (!cmmt) indent[0]='\0';
    else { for (i=0; i<30 && cmmt[i]==' '; i++) indent[i]=' '; indent[i]='\0'; }
    int  mRow = (nRow > A.nRow ? (nRow > B.nRow ? nRow : B.nRow) : (A.nRow > B.nRow ? A.nRow : B.nRow));
    if (! getEmptyFormatForPAW(format,eformat)) return;
    printf("%s : ( %d x %d )  ( %d x %d )  ( %d x %d )\n", (cmmt ? cmmt : "Matrix<>"), nRow, nCol, A.nRow, A.nCol, B.nRow, B.nCol);
    for (i = 0; i < mRow; i++) {
      printf("%s  ", indent);
      printRowForPAW( *this, i, format, eformat );
      printRowForPAW( A,     i, format, eformat );
      printRowForPAW( B,     i, format, eformat );
      printf("\n");
    }
  }
  void printAlongWith(Matrix<T> &A, Matrix<T> &B, Matrix<T> &C, char *format="%12g ", char *cmmt=NULL) {
    int  i;  char indent[31], eformat[10]; 
    if (!cmmt) indent[0]='\0';
    else { for (i=0; i<30 && cmmt[i]==' '; i++) indent[i]=' '; indent[i]='\0'; }
    int  mRow = (nRow > A.nRow ? (nRow > B.nRow ? nRow : B.nRow) : (A.nRow > B.nRow ? A.nRow : B.nRow));
    mRow = (mRow > C.nRow ? mRow : C.nRow);
    if (! getEmptyFormatForPAW(format,eformat)) return;
    printf("%s : ( %d x %d )  ( %d x %d )  ( %d x %d )  ( %d x %d )\n", (cmmt ? cmmt : "Matrix<>"), nRow, nCol, A.nRow, A.nCol, B.nRow, B.nCol, C.nRow, C.nCol);
    for (i = 0; i < mRow; i++) {
      printf("%s  ", indent);
      printRowForPAW( *this, i, format, eformat );
      printRowForPAW( A,     i, format, eformat );
      printRowForPAW( B,     i, format, eformat );
      printRowForPAW( C,     i, format, eformat );
      printf("\n");
    }
  }
private:
  bool getEmptyFormatForPAW(char *format, char *eformat) {
    // This function is meant to be called by 'printAlongWith()'.
    //   return the format string for blank space
    if (!format) { printf("Error (printAlongWith): NULL format string\n"); return false; }
    int i, space=0, len=strlen(format);
    for (i=0; i<len; i++) if (format[i]==' ') space++;
    for (i=0; i<len; i++) if (format[i]=='%') break;  
    if (i==len) { printf("Error (printAlongWith): invalid format '%s'\n", format); return false; }
    sprintf(eformat, "%%%ds", (space + atoi(format+i+1)));
    return true;
  }
  void printRowForPAW(Matrix<T> &M, int i, char *format, char *eformat) {
    // This function is meant to be called by 'printAlongWith()'.
    if (i < M.nRow) {
      printf("[ ");  for (int j=0; j<M.nCol; j++) printf(format, M.data[i*M.nCol+j]);  printf("] ");
    } else {
      printf("  ");  for (int j=0; j<M.nCol; j++) printf(eformat, "");  printf("  ");
    }
  }
  
public:
  void printStatInfo(char *cmmt=NULL) {
    if (nRow==0 || nCol==0) return;
    printf("%s (%d x %d)\n", (cmmt ? cmmt : "MTX::Matrix statistical information"), nRow, nCol );
    T *minv = (T*)calloc( nRow, sizeof(T) );
    T *maxv = (T*)calloc( nRow, sizeof(T) );
    T *sumv = (T*)calloc( nRow, sizeof(T) );
    int  i, j, k;
    for (i=k=0; i< nRow; i++) {
      for (j=0; j< nCol; j++, k++) {
	sumv[j] += data[k];
	if (i==0 || data[k] < minv[j]) minv[j] = data[k];
	if (i==0 || data[k] > maxv[j]) maxv[j] = data[k];
      }
    }
    printf("  Avg: "); for (j=0; j<nCol; j++) printf("%12g ", sumv[j]/nRow); printf("\n");
    printf("  Min: "); for (j=0; j<nCol; j++) printf("%12g ", minv[j]); printf("\n");
    printf("  Max: "); for (j=0; j<nCol; j++) printf("%12g ", maxv[j]); printf("\n");
    free( minv );  free( maxv );  free( sumv );
  }
  
};


}	// end of namespace MTX


#endif  // MTX_MATRIX_HPP
