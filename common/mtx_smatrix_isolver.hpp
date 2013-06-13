
//
// MTX::SMatrixISolver<> class template
//
// Jaeil Choi
// last modified in Nov, 2006
//
// This code has implementations of following iterative solvers:
//  - Steepest Descent		SteepestDescent(A, x, b)
//  - Conjugate Gradients	CG(A, x, b)
//
// The algorithms of Steepest Descent and CG are originally from
// "An Introduction to the Conjugate Gradient Method Without Agonizing Pain"
//         by Jonathan Richard Shewchuk, 1994, CMU.
//


#ifndef MTX_SMATRIX_ITERATIVE_SOLVER_HPP
#define MTX_SMATRIX_ITERATIVE_SOLVER_HPP


#include <iostream>
#include <iomanip>
#include <cfloat>
#include "mtx_matrix.hpp"
#include "mtx_smatrix.hpp"

namespace MTX {
  
using namespace std;


// ===================================================================
// SMatrixISolver for  Ax = b
// ===================================================================

template <class T>
class SMatrixISolver {
  
public:
  SMatrixISolver()  { }
  ~SMatrixISolver() { }

  // ----------------------------------------------------------------
  // Steepest Descent
  // ----------------------------------------------------------------

public:
  int SteepestDescent(SMatrix<T> &A, T x[], T b[], 
		      int max_iter = 0, T epsilon = 0.0001) {
    if (A.nRow != A.nCol) return -1;
    int i = 0, size = A.nRow;
    T alpha, epsilon2, delta, delta_0;
    epsilon2 = epsilon * epsilon;
    max_iter = ((max_iter < 1) ? (10 * size) : max_iter);
    T *r  = (T*)malloc(size * sizeof(T));	// residual vector
    T *Ax = (T*)malloc(size * sizeof(T));	// temporary vector
    T *q  = (T*)malloc(size * sizeof(T));	// temporary vector
    A.MultVector(x, Ax);
    sub(r, size, b, Ax);			// r(0) = b - A x(0)
    delta_0 = delta = dot(size, r, r);
    while (i < max_iter && delta > epsilon2 * delta_0) {
      A.MultVector(r, q);
      alpha = delta / dot(r, q);		// alpha  = (r(i).r(i)) / (r(i) A r(i))
      if (!(alpha >= -FLT_MAX && alpha <= FLT_MAX)) return i;
      addWithScale(x, size, x, alpha, r);		// x(i+1) = x + alpha * r(i)
      if (i % 50 == 49) {
	A.MultVector(x, Ax);
	sub(r, size, b, Ax);			// restart
      } else {
	addWithScale(r, size, r, -alpha, q);	// r(i+1) = r - alpha * A * r
      }
      delta = dot(size, r, r);
      //cout << x << " x(" << i+1 << ")" << endl;
      //cout << r << " r(" << i+1 << ")" << endl;
      i++;
    }
    free(r);  free(Ax);  free(q);
    return i;
  }
  
  
  int SteepestDescent(SMatrix<T> &A, SMatrix<T> &x, SMatrix<T> &b,
		      int max_iter = 0, T epsilon = 0.0001) {
    if (!(A.nRow == A.nCol && A.nCol == x.nRow && x.nCol == 1 && x.nRow == b.nRow && b.nCol == 1)) return -1;
    T alpha, epsilon2, delta, delta_0;
    epsilon2 = epsilon * epsilon;
    max_iter = ((max_iter < 1) ? (10 * x.nRow) : max_iter);
    int i = 0;		// iteration
    SMatrix<T> r;			// residual vector ( size: x.nRow x 1 )
    SMatrix<T> Ax, q;			// temporary vectors
    Ax.mult( A, x );
    r.sub( b, Ax );			// r(0) = b - A x(0)
    delta_0 = delta = r.dot(r);
    //cout << x << " x(0)" << endl;
    //cout << r << " r(0)" << endl;
    while (i < max_iter && delta > epsilon2 * delta_0) {
      q.mult( A, r );
      alpha = delta / r.dot(q);		// alpha  = (r(i) . r(i)) / (r(i) A r(i))
      if (!(alpha >= -FLT_MAX && alpha <= FLT_MAX)) return i;
      x.addWithScale( x, alpha, r );	// x(i+1) = x + alpha * r(i)
      if (i % 50 == 49) {
	Ax.mult( A, x );
	r.sub( b, Ax );			// restart
      } else {
	r.addWithScale( r, -alpha, q );	// r(i+1) = r - alpha * A * r
      }
      delta = r.dot(r);
      //cout << x << " x(" << i+1 << ")" << endl;
      //cout << r << " r(" << i+1 << ")" << endl;
      i++;
    }
    return i;
  }
  
  
  // ----------------------------------------------------------------
  // Conjugate Gradients
  // ----------------------------------------------------------------
public:
  
  int CG(SMatrix<T> &A, T x[], T b[]) {
    if (A.nRow != A.nCol) return -1;
    int i, size = A.nRow;
    T   rho0, rho1, rho2=0, alpha, epsilon;
    T *r  = (T*)malloc(size * sizeof(T));	// residual vector
    T *Ax = (T*)malloc(size * sizeof(T));	// temporary vector
    T *p  = (T*)malloc(size * sizeof(T));	// temporary vector
    T *q  = (T*)malloc(size * sizeof(T));	// temporary vector
    
    A.MultVector(x, Ax);
    sub(r, size, b, Ax);			// residual r(0) = b - A x(0)
    epsilon = dot(size, r, r) * 0.000001;
    for (i = 0; i <= size; i++) {
      rho1 = dot(size, r, r);			// size of residual
      if (i==0) { rho0 = rho1; memcpy(p, r, size*sizeof(T)); }	// search direction vector
      else      addWithScale(p, size, r, rho1/rho2, p);
      A.MultVector(p, q);			// 
      alpha = rho1 / dot(size, p, q);		// optimal update
      addWithScale(x, size, x, +alpha, p );	// update the solution
      addWithScale(r, size, r, -alpha, q );	// update the residual
      //if (i % 50 == 49) { A.MultVector(x, Ax); sub(r, size, b, Ax); } // restart residual
      if (dot(size, r, r) < epsilon) break;	// check if it converged
      rho2 = rho1;
    }
    return i;
  }
  
  int CG(SMatrix<T> &A, SMatrix<T> &x, SMatrix<T> &b, 
	 void (*preconditioner)(SMatrix<T> &z, SMatrix<T> &r) = NULL) {
    // Conjugate Gradient Method - Use array version above instead, if possible.
    if (!(A.nRow == A.nCol && A.nCol == x.nRow && x.nCol == 1 && x.nRow == b.nRow && b.nCol == 1)) return -1;
    int        i;
    T          rho0=0, rho1=0, rho2=0, alpha, epsilon;
    SMatrix<T> r, Ax, p, q, z;
    
    r.sub( b, Ax.mult(A, x) );		// residual r(0) = b - A x(0)
    epsilon = r.dot(r) * 0.000001;
    for (i = 0; i <= x.nRow; i++) {
      if (preconditioner) {		// preconditioned CG
	preconditioner(z, r);			// solve preconditioner M * z = r
	rho1 = r.dot(z);			// size of residual
	if (i==0) { rho0 = rho1; p = z; }	// search direction vector
	else      p.addWithScale(z, rho1/rho2, p);
      } else {				// unpreconditioned CG
	rho1 = r.dot(r);			// size of residual
	if (i==0) { rho0 = rho1; p = r; }	// search direction vector
	else      p.addWithScale(r, rho1/rho2, p);
      }
      q.mult( A, p );			// 
      alpha = rho1 / p.dot(q);		// optimal update
      x.addWithScale( x, +alpha, p );	// update the solution
      r.addWithScale( r, -alpha, q );	// update the residual
      //if (i%50 == 49) r.sub( b, Ax.mult(A,x) );	// restart residual
      if (r.dot(r) < epsilon) break;	// check if it converged
      rho2 = rho1;
    }
    return i;
  }
  
  // -----------------------------------------------------------------
  // private subroutines
  // -----------------------------------------------------------------
private:
  T* add(T r[], int n, T a[], T b[]) { for (int i = 0; i < n; i++) r[i] = a[i] + b[i]; return r; }
  T* sub(T r[], int n, T a[], T b[]) { for (int i = 0; i < n; i++) r[i] = a[i] - b[i]; return r; }
  T* addWithScale(T r[], int n, T a[], T s, T b[]) { for (int i = 0; i < n; i++) r[i] = a[i] + s * b[i]; return r; }
  T  dot(int n, T a[], T b[]) { T sum = 0;  for (int i = 0; i < n; i++) sum += a[i] * b[i]; return sum; }
  
};

}	// end of namespace MTX

#endif  // MTX_SMATRIX_ITERATIVE_SOLVER_HPP

