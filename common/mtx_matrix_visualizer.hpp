
//
// MTX::MatrixVisualizer<> class template 
// for the visualization of matrix using OpenGL
//
// Jaeil Choi
// last modified in Sep, 2005
//

#ifndef MTX_MATRIX_VISUALIZER_HPP
#define MTX_MATRIX_VISUALIZER_HPP

#include <iostream>
#include <cfloat>
#include "mtx_matrix.hpp"

namespace MTX {
using namespace std;
  
  
template <class T>
class MatrixVisualizer {

public:
  MatrixVisualizer() {}
  ~MatrixVisualizer() {}

  T getMinMax(Matrix<T> &m, T *minp=NULL, T *maxp=NULL) {
    int i, total = m.nRow * m.nCol;
    if (total <= 0) { if (minp) *minp = 0;  if (maxp) *maxp = 0; return 0; }
    T   minv, maxv;
    if (minp == NULL) minp = &minv;
    if (maxp == NULL) maxp = &maxv;
    *minp = +FLT_MAX, *maxp = -FLT_MAX;
    for (i = 0; i < total; i++) {
      if (m.data[i] > *maxp) *maxp = m.data[i];
      if (m.data[i] < *minp) *minp = m.data[i];
    }
    return (fabs(*minp) > fabs(*maxp) ? fabs(*minp) : fabs(*maxp));
  }
  
  void render(Matrix<T> &m, T max_range, bool bBorder=false, float pt_size=1.0) {
    int   i, j, size;
    float rgb[3];
    if (pt_size < 1) pt_size = 1.0;
    glPointSize(pt_size);
    glBegin(GL_POINTS);
    if (m.nRow == 1 || m.nCol == 1) {	// vector
      size = (m.nRow > 1 ? m.nRow : m.nCol);
      for (i = 0; i < size; i++) {
	DecideColor( m.data[i], rgb, max_range );
	glColor3fv( rgb );
	if (m.nRow == 1) glVertex2d( +i, -100 );
	else             glVertex2d( +100, -i );
      }
    } else {				// matrix
      for (i = 0; i < m.nRow; i++)
	for (j = 0; j < m.nCol; j++) {
	  DecideColor( m.data[i*m.nCol+j], rgb, max_range );
	  glColor3fv( rgb );
	  glVertex2d( j, -i );
	}
    }
    glEnd();
    glPointSize(1.0);
    if (bBorder) {
      glColor3ub( 100, 100, 100 );
      glBegin(GL_LINE_STRIP);
      glVertex2d(-1, +1);
      glVertex2d(-1, -m.nRow);
      glVertex2d(m.nCol, -m.nRow);
      glVertex2d(m.nCol, +1);
      glVertex2d(-1, +1);
      glEnd();
    }
    glFlush();
  }
  
  void DecideColor(T value, float rgb[3], float max_range) {
    float f = fabs(value / max_range);   if (f > 1.0) f = 1.0;
    if (value == 0.0) {
      rgb[0] = rgb[1] = rgb[2] = 0.2;
    } else if (value > 0) {
      rgb[0] = 1.0f;
      rgb[1] = 1.0f - f;
      rgb[2] = 1.0f - f;
    } else if (value < 0) {
      rgb[0] = 1.0f - f;
      rgb[1] = 1.0f - f;
      rgb[2] = 1.0f;
    }
  }
  
  void checkSymmetricity(Matrix<T> &m)
  {
    cout << endl;
    int i, j;
    // column header
    cout << "Matrix  :  nRow = " << m.nRow << "  nCol = " << m.nCol << endl;
    cout << "  ";
    for (j = 0; j < m.nCol; j++) 
      if (j%10 == 0) cout << (j/10)%10 << " ";
      else cout << "  ";
    cout << endl;
    for (i = 0; i < m.nRow; i++) {
      // row header
      if (i%10 == 0) cout << (i/10)%10 << " ";
      else cout << "  ";
      // check symmetry for all the values
      for (j = 0; j < m.nCol; j++) {
	T value = (T)fabs(m(i,j)-m(j,i));
	if (value == 0)            cout << ". ";
	else if (value < 0.000001) cout << "* ";
	else if (value < 0.001)    cout << "x ";
	else                       cout << "X ";
      }
      cout << endl;
    }
    cout << endl;
  }
  
  void Compare(Matrix<T> &m, Matrix<T> &m2)
  {
    cout << endl;
    int i, j;
    // column header
    cout << "Matrix Comparison :  nRow = " << m.nRow << "  nCol = " << m.nCol << endl;
    if (m.nRow != m2.nRow || m.nCol != m2.nCol) {
      cout << "  different size :  nRow = " << m2.nRow << "  nCol = " << m2.nCol << endl;
      return;
    }
    cout << "  ";
    for (j = 0; j < m.nCol; j++) 
      if (j%10 == 0) cout << (j/10)%10 << " ";
      else cout << "  ";
    cout << endl;
    for (i = 0; i < m.nRow; i++) {
      // row header
      if (i%10 == 0) cout << (i/10)%10 << " ";
      else cout << "  ";
      // check symmetry for all the values
      for (j = 0; j < m.nCol; j++) {
	T value = (T)fabs(m(i,j)-m2(i,j));
	if (value == 0)            cout << ". ";
	else if (value < 0.000001) cout << "* ";
	else if (value < 0.001)    cout << "x ";
	else                       cout << "X ";
      }
      cout << endl;
    }
    cout << endl;
  }
  
};

}	// end of namespace MTX

#endif  // MTX_MATRIX_VISUALIZER_HPP
