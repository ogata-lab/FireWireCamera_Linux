
#if 0

//
// THE CONTENTS OF THIS FILE IS NOT LEGITIMATE C++ CODE
// IT'S JUST A COLLECTION OF USAGE EXAMPLES.
// ( based on TNT 1.2.6,  http://math.nist.gov/tnt/ )
//


// -------------------------------------------------------------------

class TNT::Array1D <>
{
  double values[] = {1.0, 2.0, 3.0};
  TNT::Array1D<double> va, vc(3), vc(3, 1.5);
  TNT::Array1D<double> vb(va), vd(3,values);  // same array, no copying
  
  //// operator T**();
  vd = 3.0;			// [3 3 3]
  vd[1] = 4.0;			// [3 4 3]
  vd.copy();			// create a new instance (by copying)
  vc.inject(vd);		// copy from 'vd', if the size agrees.
  va = vd;   va.ref(vd);	// set 'va' to share the entries of 'vd'
  ( va.dim1() == va.dim() );	// true
  
  vd.subarray(1,2)		// [4 3]  (same array, no copying)
  int ref_cnt = vd.ref_count();
  
};


// -------------------------------------------------------------------

class TNT::Array2D <>
{
  double values[] = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0};
  TNT::Array2D<double> ma, mc(2,3), mc(2,3, 1.5);
  TNT::Array2D<double> mb(ma), md(3,3,values);  // same array, no copying
  
  //// operator T**();
  mc = 3.0;			// [3 3 3; 3 3 3]
  md[2][0] = 0.0;		// [1 2 3; 4 5 6; 0 8 9];
  md.copy();			// create a new instance (by copying)
  mc.inject(md);		// copy from 'md', if the size agrees.
  ma = md;   ma.ref(md);	// set 'ma' to share the entries of 'md'
  ( ma.dim1() == 2 && ma.dim2() == 3 );	// true
  
  md.subarray(1,2,0,1)		// [4 5; 0 8]  (same array, no copying)
  int ref_cnt = md.ref_count();
  int ref_cnt_dim1 = md.ref_count_dim1();
  
};


// -------------------------------------------------------------------
// [Note]: The Vector class has been superceded by the Array1D class
//         but it's retained for backward compatibility.

class TNT::Vector <>
{
  TNT::Vector<double> v, v1, v2, v3;
  double values[] = {1.0, 2.0, 3.0};
  
protected:
  va.initialize(3);		// vector of size 3 (by 'new' - value unset)
  va.copy( values );		// copy values
  va.set( 2.0 );		// set all entries to 2.0
  va.destroy();
  
public:
  TNT::Vector<double> va, vb(3), vc(v1), ve(4, "1, 2, 3, 4");
  TNT::Vector<double> vd(3, values);		// values are copied
  double *p0 = v1.begin();
  double *p9 = v1.end();
  vb.newsize(3);
  ve = vd;			// [1 2 3]   (copied, not sharing)
  ve = 3.0;			// [3 3 3]
  ( ve.dim() == ve.size() );	// true
  ve(1) = 4.0;			// [4 3 3]   (indexing: 1,2,...,N)
  ve[0] = 5.0;			// [5 3 3]   (indexing: 0,1,...,N-1)
  
  cout << vc;
  cin  >> vc;
  
  double d = dot_prod(vd, ve);
  vc = (vd - ve) + (vd * ve);	// ([1 2 3] - [5 3 3]) + ([1 2 3] .* [5 3 3])
  
};


// -------------------------------------------------------------------
// [Note]: The Matrix class has been superceded by the Array2D class,
//         but it's retained for backward compatibility.

class TNT::Matrix <>
{
  TNT::Matrix<double> m, m1, m2, m3;
  double values[] = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0};
  
protected:
  ma.initialize(2,3);		// 2x3 matrix (by 'new' - value unset)
  ma.copy( values );		// copy values
  ma.set( 2.0 );		// set all entries to 2.0
  ma.destroy();
  
public:
  TNT::Matrix<double> ma, mb(ma), mc(2,3), me(2,3, 3.0), mf(4, "1 2 3 4");
  TNT::Matrix<double> md(2,3, values);
  double *row = *ma;
  int    size  = ma.size();
  ma.newsize(2,3);
  me = md;			// [1 2 3; 4 5 6]  (copied, not sharing)
  me = 3.0;			// [3 3 3; 3 3 3]
  ( me.dim(1) == me.num_rows() );	// true
  ( me.dim(2) == me.num_cols() );	// true
  me(4) = 4.0;			// [3 3 3; 4 3 3]   (indexing: 1,2,...,N)
  me(2,1) = 5.0;		// [3 3 3; 5 3 3]   (indexing: 1,2,...,N)
  me[1][0] = 6.0;		// [3 3 3; 6 3 3]   (indexing: 0,1,...,N-1)
  
  cout << mc;
  cin  >> mc;
  
public: // basic matrix-matrix operators
  double d = dot_prod(md, me);
  mf = (md - me) + mult_element(md, me);  // ((2x3)-(2x3)) + ((2x3).*(2x3))
  mb = transpose(md);
  mc = (ma * mb);		// (2x3) * (3x2)
  mc = matmult(ma, mb);		// (2x3) * (3x2)
  int s = matmult(mc, ma, mb);	// (2x3) * (3x2)
  
public: // basic matrix-vector operators
  TNT:Vector<double> va(3);
  vb = matmult( ma, va );	// (2x3) * (3x1)
  vb = ma * va;			// (2x3) * (3x1)

};


#endif
