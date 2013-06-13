
//
// MTX::MatrixOp<> class template
// for Exponential, Root, Logarithm, and Linear Combination of matrix.
//
// Jaeil Choi
// last modified in Aug, 2005
//

#ifndef MTX_MATRIX_ADV_HPP
#define MTX_MATRIX_ADV_HPP

#include <iostream>
#include <cmath>
#include "mtx_matrix.hpp"
#include "mtx_matrix_solver.hpp"

namespace MTX {
  
// ===================================================================
// MTX::MatrixOp
// ===================================================================

template <class T>
class MatrixOp {
private:
  MTX::MatrixSolver<T> solver;
public:
  MatrixOp()  {}
  ~MatrixOp() {}

  // ================================================================
  // Exponential, Root, Logarithm, and Linear Combination
  // ================================================================
  
  void Exp( MTX::Matrix<T> &A, int q = 6) { Exp(A, A, q); }
  void Exp( MTX::Matrix<T> &A, MTX::Matrix<T> &X, int q = 6) {
    // 'A' is input matrix, and 'X' is the result.
    MTX::Matrix<T> A2, tmp, I, D, N;
    double j, c, k;
    // j = max ( 0, 1 + ceil(log2(||A||)) )
    j = 1 + floor( log(A.Norm()) / log(2.0) );
    j = ( j > 0 ? j : 0 );
    A2 = A;  A2.multValue( pow(2, -j) );
    I.setIdentity(A2.nRow);
    D = I;  N = I;  X = I;  c = 1;
    for (k = 1; k <= q; k++) {
      c = c * (q - k + 1) / (k * (2*q - k + 1));
      X.mult(A2, X);			// X = A2 * X
      tmp = X;  tmp.multValue(c);
      N.add(N, tmp);			// N = N + c * X
      tmp = X;  tmp.multValue(pow(-1,k) * c);
      D.add(D, tmp);			// D = D + (-1)^k * c * X
    }
    solver.inverseByGaussJordan(D, tmp);  
    tmp.mult(tmp, N);			// X = Dinv * N
    X = I;  c = pow(2, j);
    for (k = 0; k < c; k++) X.mult(X, tmp);	// X = X^{2j}
  }
  
  bool SquareRoot( MTX::Matrix<T> &A, double eps = 1.0e-11, int max_iteration = 20) {
    return SquareRoot( A, A, eps, max_iteration );
  }
  bool SquareRoot( MTX::Matrix<T> &A, MTX::Matrix<T> &X, double eps = 1.0e-11, int max_iteration = 20) {
    // 'A' is input matrix, and 'X' is the result.
    MTX::Matrix<T> Y, XX, iX, iY, A2;
    X = A2 = A;
    Y.setIdentity(X.nRow);
    //if (debug_matrix) cout << A << " calculating SquareRoot   iter : " ;
    for (int iter=0; iter <= max_iteration; iter++) {
      XX.mult(X, X);  XX.sub(A2);
      if (XX.Norm() <= eps) break;	// break if  || XX - A || <= eps
      if ( solver.inverseByGaussJordan(X) == false ||
	   solver.inverseByGaussJordan(Y) == false ) return false;
      //if (debug_matrix) cout << "X+iY" << iter << " ";
      X.add(X, iY);  X.multValue(0.5);	// X = (X + iY) / 2.0
      //if (debug_matrix) cout << "Y+iX" << iter << " ";
      Y.add(Y, iX);  Y.multValue(0.5);	// Y = (Y + iX) / 2.0
    }
    //if (debug_matrix) cout << endl;
    return true;
  }
  
  bool Log( MTX::Matrix<T> &M, double eps = 1.0e-11, int max_iteration = 20 ) {
    return Log( M, M, eps, max_iteration );
  }
  bool Log( MTX::Matrix<T> &M, MTX::Matrix<T> &X, double eps = 1.0e-11, int max_iteration = 20 ) {
    // C.Kenney, A.J.Laub, Condition estimates for matrix functions
    // SIAM Journal on Matrix Analysis and Applications
    // Basically, it's approximation by Taylor's expansion.
    // This function results in : theta * [  0  -nz  ny ]
    //                                    [  nz  0  -nx ]
    //                                    [ -ny  nx  0  ], where 
    // theta and n = [ nx, ny, nz ] is the axis-angle form of the rotation.
    // 'A' is input matrix, and 'X' is the result. Both can be the same matrix.
    MTX::Matrix<T> tmp, A, A2, I, Z;
    I.setIdentity(M.nRow);
    int k = 0;
    A2 = M;
    //cout << "Log start" << endl;
    for (k = 0; k <= max_iteration; k++) {
      tmp = A2;  tmp.sub(I);
      if (tmp.Norm() <= 0.5) break;	// break if  || A - I || <= 0.5
      if (SquareRoot(A2) == false) return false;
    }
    //cout <<  M << " (log) Input matrix" << endl;
    //cout << A2 << " (log) A^{1/2}^k  with k=" << k << endl;
    A = I;  A.sub(A2);			// A = I - A2
    Z = A;  X = A;  X.multValue(-1.0);
    int i = 1;
    while (i <= max_iteration) {
      //cout << Z << " (log) Z   at iteration=" << i << endl;
      //cout << X << " (log) X   at iteration=" << i << endl;
      if (Z.Norm() <= eps) break;	// break if || Z || <= eps
      Z.mult(Z, A);  i++;
      tmp = Z;  tmp.multValue(-1.0/i);
      X.add(X, tmp);			// X = X + Z/i
    }
    X.multValue( pow(2.0, k) );		// X = 2^k X
    return true;
  }
  
  void Pow( MTX::Matrix<T> &A, double r, MTX::Matrix<T> &X ) {
    // calcluate A^r
    // X = A^r = Exp(Log(A^r)) = Exp(r * Log(A))
    MTX::Matrix<T> M;
    Log(A, M);
    M.multValue(r);
    Exp(M, X);
  }
  
  void interpolate ( MTX::Matrix<T> &A, double r, MTX::Matrix<T> &X ) { Pow(A,r,X); }
    
  void LinearCombination ( MTX::Matrix<T> &A, double ra, MTX::Matrix<T> &B, double rb, MTX::Matrix<T> &X ) {
    // Marc Alexa, 'Linear Combination of Transformations', SIGGRAPH 2002
    // 'A' and 'B' is input matrices, and 'X' is the result.
    MTX::Matrix<T> AX, BX, M;
    Log(A, AX);    AX.multValue(ra);
    Log(B, BX);    BX.multValue(rb);
    M.add(AX, BX);
    Exp(M, X);
  }
  
};
  
  
}	// end of namespace MTX


#endif  // MTX_MATRIX_ADV_HPP
