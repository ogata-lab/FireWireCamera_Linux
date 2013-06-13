
//
// MTX::SMatrixVisualizer<> class template 
// for the visualization of matrix using OpenGL
//
// Jaeil Choi
// last modified in Nov, 2005
//

#ifndef MTX_SMATRIX_VISUALIZER_HPP
#define MTX_SMATRIX_VISUALIZER_HPP

#include <iostream>
#include <cfloat>
#include "mtx_smatrix.hpp"

namespace MTX {
  
using namespace std;

template <class T>
class SMatrixVisualizer {

public:
  SMatrixVisualizer() {}
  ~SMatrixVisualizer() {}

  T getMinMax(SMatrix<T> &m, T *minp=NULL, T *maxp=NULL) {
    T   minv, maxv;
    if (minp == NULL) minp = &minv;
    if (maxp == NULL) maxp = &maxv;
    *minp = +FLT_MAX, *maxp = -FLT_MAX;
    SMatrixEntry<T> *ep;
    for (int row = 0; row < m.nRow; row++)
      for (ep = m.Row[row]; ep; ep = ep->next) {
	if (ep->value > *maxp) *maxp = ep->value;
	if (ep->value < *minp) *minp = ep->value;
      }
    return (fabs(*minp) > fabs(*maxp) ? fabs(*minp) : fabs(*maxp));
  }
  
  void render(SMatrix<T> &m, T max_range, bool bBorder=false, float pt_size=1.0) {
    int   row;
    float rgb[3];
    SMatrixEntry<T> *ep;
    if (pt_size < 1) pt_size = 1.0;
    glPointSize(pt_size);
    glBegin(GL_POINTS);
    if (m.nCol == 1) {		// column vector
      for (row = 0; row < m.nRow; row++) {
	if ((ep = m.Row[row]) == NULL) continue;
	DecideColor( ep->value, rgb, max_range );
	glColor3fv( rgb );
	glVertex2d( +100, -row );
      }
    } else if (m.nRow == 1) {	// row vector
      for (ep = m.Row[0]; ep; ep = ep->next) {
	DecideColor( ep->value, rgb, max_range );
	glColor3fv( rgb );
	glVertex2d( ep->col, -100 );
      }
    } else {			// matrix
      for (row = 0; row < m.nRow; row++)
	for (ep = m.Row[row]; ep; ep = ep->next) {
	  DecideColor( ep->value, rgb, max_range );
	  glColor3fv( rgb );
	  glVertex2d( ep->col, -ep->row );
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
  
  void checkSymmetricity(SMatrix<T> &m)
  {
    cout << "Not implemented yet" << endl;
  }
  
  void Compare(SMatrix<T> &m, SMatrix<T> &m2)
  {
    cout << "Not implemented yet" << endl;
  }
  
};

}	// end of namespace MTX

#endif  // MTX_SMATRIX_VISUALIZER_HPP
