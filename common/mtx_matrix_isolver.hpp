
//
// MatrixISolver<> class template
//
// Jaeil Choi
// last modified in July, 2003
//
// This code has implementations of various methods for solving 'Ax = b'
//  - Jacobi's method			Jacobi(A, x, b)
//  - Gauss-Seidel method		GaussSeidel(A, x, b)
//  - Successive Overrelaxation 	SOR(A, x, b)
//  - Steepest Descent			SteepestDescent(A, x, b)
//  - Conjugate Gradients		ConjugateGradients(A, x, b)
//
// The algorithms of Steepest Descent and CG are originally from
// "An Introduction to the Conjugate Gradient Method Without Agonizing Pain"
//         by Jonathan Richard Shewchuk, 1994, CMU.
//


#ifndef MTX_MATRIX_ITERATIVE_SOLVER_HPP
#define MTX_MATRIX_ITERATIVE_SOLVER_HPP


#include <iostream>
#include <iomanip>
#include <cfloat>
#include "mtx_matrix.hpp"

namespace MTX {
  
using namespace std;


// ===================================================================
// MatrixSolver for  Ax = b
// ===================================================================

#define VALID_EQUATION3(A, x, b) \
  (A.nRow == A.nCol && A.nCol == x.nRow && x.nCol == 1 && x.nRow == b.nRow && b.nCol == 1)
#define VALID_EQUATION2(A, b) \
  (A.nRow == A.nCol && A.nRow == b.nRow && b.nCol == 1)

template <class T>
class MatrixISolver {
 public:
  MatrixISolver()  { }
  ~MatrixISolver() { }

  // ----------------------------------------------------------------
  // Jacobi's method
  // ----------------------------------------------------------------
  
  int Jacobi(Matrix<T> &A, Matrix<T> &x, Matrix<T> &b, int max_iter = 20) {
    // Iterative relaxation method of dividing A into D and [L+U]
    // D x(k+1) = - [L + U] x(k) + b	

    if (!VALID_EQUATION3(A, x, b)) return -1;
    int i, j, k;
    T   *tmp, sum;
    T   *xnew = (T*) malloc( x.nRow * sizeof(T) );
    T   *xold = (T*) malloc( x.nRow * sizeof(T) );
    
    //cout << "Jacobi's relaxation method for 'A x = b' " << endl;
    //cout << x << " x(" << 0 << ")" << endl;
    memcpy( xnew, x.data, x.nRow * sizeof(T) );
    for (int k = 0; k < max_iter; k++) {
      tmp = xnew;  xnew = xold;  xold = tmp;	// swap
      for (i = 0; i < A.nRow; i++) {
	// aii * xi(k+1) = - SUM(j != i) [ aij * xj(k) ] + bi
	for (sum = j = 0; j < A.nCol; j++)  if (j != i) sum += A(i,j) * x.data[j];
	xnew[i] = (- sum + b.data[i]) / A(i,i);
      }
      //Matrix<T> xm;  xm.set(x.nRow,1, xnew);
      //cout << xm << " x(" << k+1 << ")" << endl;
    }
    memcpy( x.data, xnew, x.nRow * sizeof(T) );
    return 0;
  }
  
  // ----------------------------------------------------------------
  // Gauss-Seidel method
  // ----------------------------------------------------------------
  
  int GaussSeidel(Matrix<T> &A, Matrix<T> &x, Matrix<T> &b, int max_iter = 20) {
    // Iterative relaxation method of dividing A into [L+D] and U
    // [D + L] x(k+1) = - U x(k) + b
    // D x(k+1) = - (L x(k+1) + U x(k)) + b   (x is updated in-place)
    // if all aii > 0, Gauss-Seidel converges iff A is positive definite.
    //
    // For sparse matrices, it is much more efficient to apply this method
    // without constructing the matrix 'A' explicitly.
    
    if (!VALID_EQUATION3(A, x, b)) return -1;
    int i, j, k;
    T   sum;
    
    //cout << "GaussSeidel relaxation method for 'A x = b' " << endl;
    //cout << x << " x(" << 0 << ")" << endl;
    for (int k = 0; k < max_iter; k++) {
      for (i = 0; i < A.nRow; i++) {
	// aii * xi(k+1) = - SUM(j != i) [ aij * xj(k) ] + bi
	for (sum = j = 0; j < A.nCol; j++)  if (j != i) sum += A(i,j) * x.data[j];
	x.data[i] = (- sum + b.data[i]) / A(i,i);
      }
      //cout << x << " x(" << k+1 << ")" << endl;
    }
    return 0;
  }
  
  // ----------------------------------------------------------------
  // SOR : Successive Overrelaxation
  // ----------------------------------------------------------------
  
  int SOR(Matrix<T> &A, Matrix<T> &x, Matrix<T> &b) {
    if (!VALID_EQUATION3(A, x, b)) return -1;
    //// Under Construction
    return 0;
  }
  
  // ----------------------------------------------------------------
  // Steepest Descent
  // ----------------------------------------------------------------
  
  int SteepestDescent(Matrix<T> &A, Matrix<T> &x, Matrix<T> &b,
		      int max_iter = 0, T epsilon = 0.0001) {
    if (!VALID_EQUATION3(A, x, b)) return -1;
    T alpha, epsilon2, delta, delta_0;
    epsilon2 = epsilon * epsilon;
    max_iter = ((max_iter < 1) ? (10 * x.nRow) : max_iter);
    int i = 0;		// iteration
    Matrix<T> r;			// residual vector ( size: x.nRow x 1 )
    Matrix<T> Ax, q;			// temporary vectors
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
  
  int ConjugateGradients(Matrix<T> &A, Matrix<T> &x, Matrix<T> &b,
			 int max_iter = 0, T epsilon = 0.0001) {
    // Iterative method for solving  A x = b .
    // only for symmetric positive definite matrix A.
    // Use x as starting point, and update x with final result.
    
    if (!VALID_EQUATION3(A, x, b)) return -1;
    T alpha, epsilon2, delta_old, delta_new, delta_0, beta;
    epsilon2 = epsilon * epsilon;
    max_iter = ((max_iter < 1) ? x.nRow : max_iter);
    int i = 0;		// iteration
    Matrix<T> r;			// residual vector         ( size: x.nRow x 1 )
    Matrix<T> d;			// search direction vector ( size: x.nRow x 1 )
    Matrix<T> Ax, q;			// temporary vectors
    Ax.mult( A, x );
    r.sub( b, Ax );			// r(0) = b - A x(0)
    d = r;				// d(0) = r(0)
    delta_0 = delta_new = r.dot(r);
    //cout << x << " x(0)" << endl;
    //cout << r << " r(0)" << endl;
    //cout << d << " d(0)" << endl;
    while (i < max_iter && delta_new > epsilon2 * delta_0) {
      q.mult( A, d );
      alpha = delta_new / d.dot(q);   // alpha  = (r(i) . r(i)) / (d(i) A d(i))
      if (!(alpha >= -FLT_MAX && alpha <= FLT_MAX)) return i;
      x.addWithScale( x, alpha, d );	// x(i+1) = x + alpha * d(i)
      if (i % 50 == 49) {
	Ax.mult( A, x );
	r.sub( b, Ax );		// r(i+1) = b - A x(i+1)   (restart)
      } else {
	r.addWithScale( r, -alpha, q );	// r(i+1) = r(i) - alpha A d(i)
      }
      delta_old = delta_new;
      delta_new = r.dot(r);		// delta_new = r(i+1) . r(i+1)
      beta = delta_new / delta_old;	// beta = (r(i+1) . r(i+1)) / (r(i) . r(i))
      d.addWithScale( r, beta, d );	// d(i+1) = r(i+1) + beta d(i)
      //cout << x << " x(" << i+1 << ")" << endl;
      //cout << r << " r(" << i+1 << ")" << endl;
      //cout << d << " d(" << i+1 << ")" << endl;
      i++;
    }
    return i;
  }  // end of ConjugateGradient()  Matrix<> version
  
  int CG(Matrix<T> &A, Matrix<T> &x, Matrix<T> &b, 
	 void (*preconditioner)(Matrix<T> &z, Matrix<T> &r) = NULL) {
    // Conjugate Gradient Method
    if (!VALID_EQUATION3(A, x, b)) return -1;
    int        i;
    T          rho0, rho1, rho2, alpha, epsilon;
    Matrix<T> r, Ax, p, q, z;
    
    r.sub( b, Ax.mult(A, x) );		// residual r(0) = b - A x(0)
    epsilon = r.dot(r) * 0.00000001;
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
  
  int BiCG(Matrix<T> &A, Matrix<T> &x, Matrix<T> &b,
	   void (*preconditioner)(Matrix<T> &z, Matrix<T> &r, bool bT) = NULL) {
    // Bi-Conjugate Gradient Method for non-symmetric matrix
    if (!VALID_EQUATION3(A, x, b)) return -1;
    int        i;
    T          rho0, rho1, rho2, alpha, epsilon;
    Matrix<T> r1, r2, Ax, At, p1, p2, q1, q2, z1, z2;
    
    At.transpose(A);
    r1.sub( b, Ax.mult(A, x) );		// residual r(0) = b - A x(0)
    r2 = r1;
    epsilon = r1.dot(r1) * 0.000001;
    for (i = 0; i <= 2*x.nRow; i++) {
      if (preconditioner) {		// preconditioned BiCG
	preconditioner(z1, r1, false);		// solve preconditioner M  z1 = r1
	preconditioner(z2, r2, true);		// solve preconditioner Mt z2 = r2
	rho1 = z1.dot(r2);			// size of residual
	if (rho1 == 0) return -1;
	if (i==0) { 
	  rho0 = rho1; p1 = z1; p2 = z2;	// search direction vector
	} else {
	  p1.addWithScale(z1, rho1/rho2, p1);
	  p2.addWithScale(z2, rho1/rho2, p2);
	}
      } else {				// unpreconditioned BiCG
	rho1 = r1.dot(r2);			// size of residual
	if (rho1 == 0) return -1;
	if (i==0) { 
	  rho0 = rho1; p1 = r1; p2 = r2;	// search direction vector
	} else {
	  p1.addWithScale(r1, rho1/rho2, p1);
	  p2.addWithScale(r2, rho1/rho2, p2);
	}
      }
      q1.mult( A,  p1 );
      q2.mult( At, p2 );
      alpha = rho1 / p2.dot(q1);	// optimal update
      x.addWithScale( x, +alpha, p1 );	// update the solution
      r1.addWithScale( r1, -alpha, q2 );	// update the residual
      r2.addWithScale( r2, -alpha, q2 );	// update the residual
      if (r1.dot(r1) < epsilon) break;	// check if it converged
      rho2 = rho1;
    }
    return i;
  }
  
};
  
  
}	// end of namespace MTX


#endif  // MTX_MATRIX_ITERATIVE_SOLVER_HPP



/* ================================================================ */
#if 0	// coding example

#include "Matrix.hpp"
#include "MatrixISolver.hpp"

int main(void)
{
  Matrix<float> A, x, b;
  MatrixISolver<float> solver;
  int iter;
  
  A.assign( 2, 2,  3.0, 2.0,  2.0, 6.0 );
  b.assign( 2, 1,  2.0, -8.0 );
  x.assign( 2, 1,  -2.0, -2.0 );
  PrintMatrices( A, x, b );  cout << endl << endl;
  
  x.assign( 2, 1,  -2.0, -2.0 );
  solver.Jacobi( A, x, b );
  cout << x << "Result of Jacobi's method" << endl << endl;
  
  x.assign( 2, 1,  -2.0, -2.0 );
  solver.GaussSeidel( A, x, b );
  cout << x << "Result of Gauss-Seidel method" << endl << endl;
  
  x.assign( 2, 1,  -2.0, -2.0 );
  iter = solver.SteepestDescent( A, x, b );
  cout << x << "Result of Steepest Descent (iter=" << iter << ")" << endl << endl;
  
  x.assign( 2, 1,  -2.0, -2.0 );
  iter = solver.ConjugateGradients( A, x, b );
  cout << x << "Result of Conjugate Gradients (iter=" << iter << ")" << endl << endl;
  
  x.assign( 2, 1,  -2.0, -2.0 );
  iter = solver.CG( A, x, b );
  cout << x << "Result of new CG (iter=" << iter << ")" << endl << endl;
  
  x.assign( 2, 1,  -2.0, -2.0 );
  iter = solver.BiCG( A, x, b );
  cout << x << "Result of BiCG (iter=" << iter << ")" << endl << endl;
  
  return EXIT_SUCCESS;
}
  
#endif
/* ================================================================ */
