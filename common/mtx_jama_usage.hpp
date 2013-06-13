
#if 0

//
// THE CONTENTS OF THIS FILE IS NOT LEGITIMATE C++ CODE
// IT'S JUST A COLLECTION OF USAGE EXAMPLES.
// ( based on JAMA 1.2.4,  http://math.nist.gov/tnt/ )
//


// -------------------------------------------------------------------

class JAMA::LU <>
{
  //   For an mxn matrix A with m >= n, 
  // the LU decomposition is an mxn unit lower triangular matrix L, 
  // an nxn upper triangular matrix U, and a permutation vector piv of length m 
  // so that A(piv,:) = L*U. If m < n, then L is mxm and U is mxn.
  //   The LU decompostion with pivoting ALWAYS EXISTS, 
  // even if the matrix is singular, so the constructor will never fail. 
  // The primary use of the LU decomposition is in the solution of 
  // square systems of simultaneous linear equations. 
  // This will fail if isNonsingular() returns false.   
  
  TNT::Array2D<double> A(M,N);  // set, and assign values of A
  
  JAMA::LU<double> jLU( A );	// LU decomposition
  if ( !jLU.isNonsingular() ) return FAILURE;
  
  TNT::Array2D<double> L, U;	// A(piv,:) = L*U
  TNT::Array1D<double> piv;
  L = jLU.getL();		// MxN unit lower triangular
  U = jLU.getU();		// NxN upper triangular
  piv = jLU.getPivot();		// M permutation vector
  
  double det = jLU.det();	// determinant
  
  TNT::Array1D<double> x(N), b(N);
  // assign values of b vector (It won't be touched)
  x = jLU.solve( b );		// solve a linear system
  
  TNT::Array1D<double> X(N,K), b(N,K);
  // assign values of B matrix (It won't be touched)
  X = jLU.solve( B );		// solve multiple (K) linear systems
  
};  
  

// -------------------------------------------------------------------

class JAMA::QR <>
{
  //   Calculate A = Q*R for an mxn matrix A, where
  // Q is mxn orthogonal, and R is nxn upper triangular matrix.
  //   The QR decompostion ALWAYS EXISTS, even if the matrix does not 
  // have full rank, so the constructor will never fail. 
  // The primary use of the QR decomposition is in the least squares 
  // solution of nonsquare systems of simultaneous linear equations. 
  // This will fail if isFullRank() returns 0 (false).
  //   The Q and R factors can be retrived via the getQ() and getR() 
  // methods. Furthermore, a solve() method is provided to find 
  // the least squares solution of Ax=b using the QR factors.

  TNT::Array2D<double> A(M,N);  // set, and assign values of A
  
  JAMA::QR<double> jQR( A );	// QR decomposition
  (jQR.isFullRank() ? : );	// TRUE/FALSE
  
  TNT::Array2D<double> Q, R, h;	// A = Q*R
  Q = jQR.getQ();		// MxN orthogonal
  R = jQR.getR();		// NxN upper triangular
  h = jQR.getHouseholder();	// 
  
  TNT::Array1D<double> x(N), b(N);
  // assign values of b vector
  x = jQR.solve( b );		// solve a linear system
  
  TNT::Array2D<double> X(N,K), b(N,K);
  // assign values of B matrix
  X = jQR.solve( B );		// solve multiple linear systems
	  
};


// -------------------------------------------------------------------

class JAMA::SVD <>
{
  //   Calculate A = U*S*V' for an mxn matrix A with m >= n, where
  // U is mxn orthogonal, S is nxn diagonal, and V is nxn orthogonal matrix.
  //   The singular values, sigma[k] = S[k][k], are ordered 
  // so that sigma[0] >= sigma[1] >= ... >= sigma[n-1].
  //   The singular value decompostion ALWAYS EXISTS, so the constructor 
  // will never fail. The matrix condition number and the effective 
  // numerical rank can be computed from this decomposition.

  TNT::Array2D<double> A(M,N);  // set, and assign values of A
  
  JAMA::SVD<double> jSVD( A );	// SVD decomposition
  
  TNT::Array2D<double> U, V, S;	// A = U*S*V'
  TNT::Array1D<double> s;
  jSVD.getU( U );		// MxN orthogonal, copied from jSVD.U
  jSVD.getV( V );		// NxN orthogonal, referred to jSVD.V
  jSVD.getS( S );		// NxN diagonal,   copied from dig( jSVD.s )
  jSVD.getSingularValues( s );	// N	referred to jSVD.s
  
  double n2 = jSVD.norm2();	// max(S) = s[0]
  double cn = jSVD.cond();	// max(S) / min(S)
  int    rk = jSVD.rank();	// number of nonnegligible singular values
  
};


// -------------------------------------------------------------------

class JAMA::Cholesky <>
{
  //   For a symmetric, positive definite matrix A, 
  // this function computes the Cholesky factorization, i.e. 
  // it computes a lower triangular matrix L such that A = L*L'. 
  // If the matrix is not symmetric or positive definite, 
  // the function computes only a partial decomposition. 
  // This can be tested with the is_spd() flag.
  
  TNT::Array2D<double> A(N,N);  // set, and assign values of A
  
  JAMA::Cholesky<double> jCholesky( A );
  
  TNT::Array2D<double> L;	// A = Q*R
  L = jCholesky.getL();		// lower triangular
  int v = jCholesky.is_spd();	// 1 if A is symmetric positive-definite
  
  if (v == 1) {
    TNT::Array1D<double> x(N), b(N);
    // assign values of b vector
    x = jCholesky.solve( b );	// solve a linear system
  
    TNT::Array2D<double> X(N,K), b(N,K);
    // assign values of B matrix
    X = jCholesky.solve( B );	// solve multiple linear systems
  }
};


// -------------------------------------------------------------------

class JAMA::Eigenvalue <>
{
  // Computes eigenvalues and eigenvectors of a real (non-complex) matrix.
  //   If A is symmetric, then A = V*D*V' where the eigenvalue matrix D 
  // is diagonal and the eigenvector matrix V is orthogonal. That is, 
  // the diagonal values of D are the eigenvalues, and V*V' = I, where 
  // I is the identity matrix. The columns of V represent the eigenvectors 
  // in the sense that A*V = V*D.
  //   If A is not symmetric, then the eigenvalue matrix D is block diagonal
  // with the real eigenvalues in 1x1 blocks and any complex eigenvalues,
  // a + i*b, in 2x2 blocks, [a, b; -b, a]. This keeps V a real matrix 
  // in both symmetric and non-symmetric cases, and A*V = V*D.

  TNT::Array2D<double> A(N,N);  // set, and assign values of A
  
  JAMA::Eigenvalue<double> jEig( A );	// SVD decomposition
  
  TNT::Array2D<double> V, D;	// A = V*D*V'
  jEig.getV( V );		// NxN orthogonal, referred to jEig.V
  jEig.getD( D );		// NxN diagonal,   copied from jEig.U
  
  TNT::Array1D<double> er, ei;
  jEig.getRealEigenvalues( er );
  jEig.getImagEigenvalues( ei );
  
};

#endif
