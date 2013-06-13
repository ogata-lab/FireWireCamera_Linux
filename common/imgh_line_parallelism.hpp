
//
// IMGH::LineDetector 
//
// Jaeil Choi
// last modified in Dec, 2006
//
// -------------------------------------------------------------------
// 
// This code is a part of IMGH library (IMaGe library in Header files).
// IMGH is a very easy-to-use, light-weight and cross-platform image library,
//   which is defined entirely in header files. The purpose of this library
//   is to eliminate the installation and other setup, while providing
//   most functions of any image processing libraries.
// IMGH library has several advantages:
//   - No installation required. Just copy, and '#include' them.
//   - Operating System independent (Windows/Linux/Mac)
//   - Conversion from/to different image format can be easily done without copying.
//   - Provides a tool to draw lines, circles and polygons on the image.
//   - Supports pixel formats of 'integer', 'float', or even 'void* pointer',
//     which are useful for scalar field, or a map of pointers to complex data structures.
// IMGH library consists of:
//   IMGH::Image		in imgh_common.hpp	for creation & basic handling
//   IMGH::ImageFileIO		in imgh_fileio.hpp	for file I/O (JPG,PNG,PGM, and more)
//   IMGH::ImageConverter	in imgh_converter.hpp	for pixel format conversion
//   IMGH::ImageEditor		in imgh_editor.hpp	for image edition (cut,resize,merge,etc.)
//   IMGH::ImageDrawer		in imgh_drawer.hpp	for drawing on the image
//   IMGH::ImageFilter		in imgh_filter.hpp	for Gaussian filtering
//   IMGH::ImageGradient	in imgh_gradient.hpp	for gradient images
//   IMGH::CornerDetector	in imgh_corner_detector.hpp
//   IMGH::EdgeDetector		in imgh_edge_detector.hpp
//   IMGH::SuperPixels		in imgh_super_pixels.hpp
//   IMGH::IntegralImage	in imgh_integral.hpp
//   and some more. (The rest are under development, and have not been tested enough)
// 
// --------------------------------------------------------------------------
//
//   This program is free software; you can redistribute it and/or modify it
//  under the terms of the GNU General Public License as published by the
//  Free Software Foundation; either version 2 of the License or later.
//  This program is distributed in the hope that it will be useful, but 
//  WITHOUT ANY WARRANTY. See the GNU General Public License for more details.
//
// -------------------------------------------------------------------
// 
// This class is designed for:
//   - extracting parallel models from a set of lines
// This class requires:
//   - IMGH::Image		in 'imgh_common.hpp'
//   - IMGH::EdgeDetector	in 'imgh_edge_detector.hpp'
//   - IMGH::LineDetector	in 'imgh_line_detector.hpp'
//   - IMGH::LineMapper		in 'imgh_line_mapper.hpp'
//

#ifndef IMGH_LINE_PARALLELISM_HPP
#define IMGH_LINE_PARALLELISM_HPP

#include <iostream>
#include <cmath>
#include "imgh_common.hpp"
#include "imgh_drawer.hpp"
#include "imgh_line_detector.hpp"
#include "algh_linked_list.hpp"
#include "algh_sorter.hpp"
#include "mtx_matrix_solver.hpp"

#ifndef M_PI
#define M_PI (3.14159265358979323846)
#endif

namespace IMGH {
  
class PLGroup {	// group of parallel lines
public:
  int		w, h;		// width and height of the image
  bool		vertical;	// whether or not the lines are vertical
  bool		projective;	// whether or not the lines were undergone projection (with vanishing point)
  double	vp[2];		// (projective ? VanishingPoint[2] : AvgNormalVector[2])
  double	error[2];	// error per unit length  ([0]:average and [1]:max)
  ALGH::LinkedList<LDLine*>	lines;	// list of lines, sorted by their locations
public:
  PLGroup() {}
  ~PLGroup() {}
public:
  void   getLineEquation(double x, double y, double eq[3]) {
    if (projective) { eq[0] = -(vp[1]-y); eq[1] = +(vp[0]-x); eq[2] = -(eq[0]*x+eq[1]*y); }
    else            { eq[0] =   vp[0];    eq[1] =   vp[1];    eq[2] = -(eq[0]*x+eq[1]*y); }
    double len = sqrt(eq[0] * eq[0] + eq[1] * eq[1]);
    if (len>0) { eq[0] /= len;  eq[1] /= len;  eq[2] /= len; }
  }
  double getLocation(double px, double py) {
    // The return value of this function is used to determine the order 
    //   of points, sorted by the given parallel model.
    if (projective) {	// eq of line through (px,py): (py-vp[1])*x + (vp[0]-px)*y + c = 0
      if (vertical) 	// projected point at y=0 : x = - c / (py-vp[1]);
	return ( px + (vp[0] - px) * py / (py - vp[1]) );
      else if (vp[0]<0)	// projected point at x=0 : y = - c / (vp[0]-px);
	return ( py + (py - vp[1]) * px / (vp[0] - px) );
      else if (vp[0]>w)	// projected point at x=w : y = - c / (vp[0]-px) - (py-vp[1])*w/(vp[0]-px);
	return ( py + (py - vp[1]) * (px - w) / (vp[0] - px) );
      else {		// use angle instead (degree from downward(0,1), [-180,+180])
	double dir[2], len;
	IMGH2V_SET( dir, (px - vp[0]), (py - vp[1]) );
	len = IMGH2V_LEN( dir );  if (len > 0) IMGH2V_DIV( dir, len ); 
	return ( acos( dir[1] ) * 180 / M_PI );  // (right +, left -)
      }
    } else {		// eq of line through (px,py): vp[0]*x + vp[1]*y + c = 0
      if (vertical) 	// projected point at y=0 : x = - c / vp[0]
	return (px + vp[1] * py / vp[0]);
      else		// projected point at x=0 : y = - c / vp[1]
	return ( vp[0] * px / vp[1] + py );
    }
  }
  void  getLineLocation(LDLine *lp, double loc[2]) {
    loc[0] = getLocation( lp->p[0][0], lp->p[0][1] );
    loc[1] = getLocation( lp->p[1][0], lp->p[1][1] );
  }
  void  updateLineLocation(LDLine *lp) {
    float avg[2];
    IMGH2V_AVG( avg, lp->p[0], lp->p[1] );
    lp->loc = this->getLocation( avg[0], avg[1] );
  }
  void  updateLineLocations(void) {
    // The value 'lp->loc' is used to sort lines in the same parallel group.
    LDLine *lp;
    for (lines.goFirst(); (lp = lines.getCurr()); lines.goNext()) 
      updateLineLocation( lp );
  }
};

#define MAX_PLG	6
  

class LineParallelism
{
public:
  int			w, h, total;	// frequently used values
  int			nplg;		// number of groups of parallel lines
  IMGH::PLGroup		plg[MAX_PLG];	// groups
  IMGH::LineDetector	*ldt;
private:
  IMGH::LineDetector	ldt_dummy;
  bool			debug;
  int	count_m_old, count_m_new;
    
public:
  LineParallelism() : nplg(0), ldt(&ldt_dummy), debug(false) {}
  ~LineParallelism() {}
  
  void clear(void) { for (int i=0; i<nplg; i++) plg[i].lines.clear();  nplg = 0; }
  int  countTotalLines(void) {
    int i, count;
    for (i = count = 0; i < nplg; i++) count += plg[i].lines.count;
    return count;
  }
  
  // -----------------------------------------------------------------
  // 
  // -----------------------------------------------------------------
public:  
  int create(IMGH::LineDetector *ldt, bool improve=true) {
    if (ldt) this->ldt = ldt;
    else     this->ldt = ldt = &ldt_dummy;
    this->w = ldt->w;  this->h = ldt->h;  this->total = w * h;
    int     it, i, j, acount, bcount, count;
    double  aerror[2], berror[2], avp[2], bvp[2], *vp, *error, lerror;
    bool    projective;
    LDLine *lp, **alist, **blist, **llist;
    alist = (LDLine**) malloc( ldt->nlines * sizeof(LDLine*) );
    blist = (LDLine**) malloc( ldt->nlines * sizeof(LDLine*) );
    debug = true;
    clear();
    // clear the parallel line group membership
    for (i=0; i < ldt->nlines; i++) ldt->lines[i]->plg = -1;
    
    for (it = 0; true; it++) {
      switch (it) {
      case 0:	// find vertical lines -------------------------------
	for (i=0, acount=bcount=0; i < ldt->nlines; i++) {
	  lp = ldt->lines[i];
	  if (lp->plg != -1) continue;
	  if (fabs(lp->eq[0]) < 0.98 || lp->tlen < ldt->h/10) continue;
	  if (lp->nplines>0 && lp->tlen < lp->plines[0]->tlen) continue;
	  alist[acount++] = lp;   blist[bcount++] = lp;
	}
	fitParallelModel( alist, &acount, true,  0.04, avp, aerror );
	fitParallelModel( blist, &bcount, false, 0.04, bvp, berror );
	projective = (acount > bcount);
	if (projective) { llist = alist; count = acount; vp = avp; error = aerror; }
	else            { llist = blist; count = bcount; vp = bvp; error = berror; }
	plg[nplg].vertical   = true;
	break;
      case 1:	// find parallel horizontal lines --------------------
	for (i=0, acount=bcount=0; i < ldt->nlines; i++) {
	  lp = ldt->lines[i];
	  if (lp->plg != -1) continue;
	  if (fabs(lp->eq[1]) < 0.99 || lp->tlen < ldt->w/10) continue;
	  if (lp->nplines>0 && lp->tlen < lp->plines[0]->tlen) continue;
	  alist[acount++] = lp;   blist[bcount++] = lp;
	}
	fitParallelModel( alist, &acount, true,  0.03, avp, aerror );
	fitParallelModel( blist, &bcount, false, 0.03, bvp, berror );
	projective = (acount >= bcount);
	if (projective) { llist = alist; count = acount; vp = avp; error = aerror; }
	else            { llist = blist; count = bcount; vp = bvp; error = berror; }
	plg[nplg].vertical   = false;
	break;
      default:	// find projected horizontal lines -------------------
	for (i=0, acount=bcount=0; i < ldt->nlines; i++) {
	  lp = ldt->lines[i];
	  if (lp->plg != -1) continue;
	  if (fabs(lp->eq[1]) < 0.5 || lp->tlen < ldt->w/10) continue;
	  if (lp->nplines>0 && lp->tlen < lp->plines[0]->tlen) continue;
	  alist[acount++] = lp;
	}
	fitParallelModel( alist, &acount, true, 0.025, avp, aerror );
	{ llist = alist; count = acount; vp = avp; error = aerror; }
	plg[nplg].vertical   = false;
	projective = true;
	break;
      }
      if (it != 0 && count < 3) { if (it > 1) break; else continue; }
      // save the parallelism model
      plg[nplg].projective = projective;
      plg[nplg].error[0] = error[0];  plg[nplg].error[1] = error[1];
      plg[nplg].vp[0]    = vp[0];     plg[nplg].vp[1]    = vp[1];  
      // add all the similar lines to the list
      for (i=0, count=0; i < ldt->nlines; i++) {
	lp = ldt->lines[i];
	lerror = getError( lp, plg[nplg].vp, plg[nplg].projective );
	if (lerror < plg[nplg].error[1]) { llist[count++] = lp; lp->plg = nplg; }
      }
      // make sure all the parallel lines are in a same group
      for (i=0, acount=count; i < acount; i++) {
	lp = llist[i];
	for (j=0; j < lp->nplines; j++)
	  if (lp->plines[j]->plg != nplg) { 
	    lp->plines[j]->plg = nplg;  llist[count++] = lp->plines[j];
	  }
      }
      // sort the list in the ascending order of projected locations
      for (i=0; i<count; i++) plg[nplg].updateLineLocation( llist[i] );
      ALGH::Sorter<LDLine*> sorter;
      sorter.QuickSort( llist, count, compare_lines_by_loc );
      // save the lines in the new list
      for (i=0; i<count; i++) plg[nplg].lines.append( llist[i] );
      if (++nplg >= MAX_PLG) break;
    }
    free(alist);  free(blist);
//     if (improve) {
//       mergeConnectedLines();
//       removeInvalidatedLines();
//     }
    return nplg;
  }
  
private:
  void fitParallelModel(LDLine *llist[], int *nlines, bool projective, 
			double threshold, double vp[2], double error[2]) {
    // Try to find a optimal vanishing point that agrees with as many lines as possible
    int    i, count=*nlines, worst_idx=-1;
    double *A = (double*) malloc( 2 * count * sizeof(double) );
    double *B = (double*) malloc(     count * sizeof(double) );
    if (A == NULL) { std::cerr << "Error (IMGH::LineParallelism): malloc failed" << std::endl; return; }
    double len, lerror, worst_error=0, esum=0, wsum=0;
    MTX::MatrixSolver<double> solver;
    LDLine *lp;
    while (count >= 3) {	
      if (projective) {	// try to find a vanishing point
	for (i = 0; i < count; i++) {	// setup the equation A x = b
	  lp = llist[i];  A[ 2*i+0 ] = lp->eq[0];  A[ 2*i+1 ] = lp->eq[1];  B[i] = -lp->eq[2];
	}
	if (!solver.solveByQR( count, 2, A, vp, B )) { *nlines = 0;  error[0] = error[1] = -1;  break; }
      } else {		// calculate the weighted average of line normals
	for (vp[0]=vp[1]=i=0; i < count; i++) {	
	  lp = llist[i];
	  vp[0] += lp->tlen * fabs(lp->eq[0]);  
	  vp[1] += lp->tlen * fabs(lp->eq[1]);
	}
	len = sqrt(vp[0] * vp[0] + vp[1] * vp[1]);
	vp[0] /= len;  vp[1] /= len;
      }
      for (worst_error=esum=wsum=i=0; i < count; i++) {	// find the worst
	lerror = getError( llist[i], vp, projective );
	esum += llist[i]->tlen * lerror;
	wsum += llist[i]->tlen;
	if (lerror > worst_error) { worst_error = lerror; worst_idx = i; }
      }
      if (worst_error >= threshold) {	// remove the worst
	esum -= llist[worst_idx]->tlen * worst_error;
	wsum -= llist[worst_idx]->tlen;
	for (i = worst_idx; i < count-1; i++) llist[i] = llist[i+1];
	count--;
      } else break;
    }
    free(A);  free(B);
    if (count > 2) {
      error[0] = esum / wsum;
      error[1] = worst_error;
    } else {
      error[0] = error[1] = -1;
    }
//     if (projective) printf("  Vanishing Point (%7.1f %7.1f)  (%2d -> %2d)  error[]=(%.4f %.4f) \n", vp[0], vp[1], *nlines, count, error[0], error[1]);
//     else            printf("  Parallel lines  (%7.4f %7.4f)  (%2d -> %2d)  error[]=(%.4f %.4f) \n", vp[0], vp[1], *nlines, count, error[0], error[1]);
    *nlines = count;
  }
  
  double getError(LDLine *lp, double vp[2], bool projective) {
    double *normal, dummy[2], len, cosv;
    if (!projective) normal = vp;
    else {
      normal = dummy;
      normal[0] = -(vp[1] - (lp->p[0][1] + lp->p[1][1])/2);
      normal[1] = +(vp[0] - (lp->p[0][0] + lp->p[1][0])/2);
      len = IMGH2V_LEN( normal );  IMGH2V_DIV( normal, len );
    }
    cosv   = fabs( IMGH2V_DOT( normal, lp->eq ) );
    return (cosv >= 1.0 ? 0 : sqrt( 1 - cosv * cosv ) / cosv);	// error (per unit length) for the line
  }
  
  void mergeConnectedLines(void) {
    int    i, j, pidx, merge_count;
    LDLine *la, *lb;
    count_m_old = count_m_new = 0;
    for (i = 0; i < nplg; i++) {
      for (plg[i].lines.goFirst(); (la = plg[i].lines.getCurr()); plg[i].lines.goNext()) {
	if (la->isInvalid()) continue;  // skip invalidated lines
	for (j = merge_count = 0; j < la->nlinks; j++) {
	  pidx = la->links[j].pidx;
	  lb   = la->links[j].lp;
	  if (la->plg != lb->plg || IMGH2V_DOT(la->eq, lb->eq) < 0.9) continue;
	  double d[5];
	  la->getRelativePosition( lb, d, true );
	  if (d[2] < 1 && d[3] < 1 || d[2] > d[4]-1 && d[3] > d[4]-1) {
	    // merge 'la' and 'lb' into 'la', and invalidate 'lb'
	    ldt->mergeTwoLines( la, lb, false );  merge_count++;
	  }
	}
	if (merge_count) { count_m_old += merge_count;  count_m_new++; }
      }
    }
  }

  void removeInvalidatedLines(void) {
    LDLine *la;
    for (int i = 0; i < nplg; i++) {
      for (plg[i].lines.goFirst(); (la = plg[i].lines.getCurr()); ) {
	if (la->isInvalid()) plg[i].lines.remCurr();
	else plg[i].lines.goNext();
      }
    }
    ldt->removeInvalidatedLines();
  }

public:
  void printInfo(void) {
    int  i, count, pcount=0;
    if (ldt == NULL) return;
    if (ldt == &ldt_dummy) ldt->printInfo();
    printf("IMGH::LineParallelism \n");
    for (i=0; i < MAX_PLG && i < nplg; i++) {
      if (plg[i].projective) {
	printf("  parallel line group %d : %3d  %s  %s (%7.1f %7.1f)  error=(%.4f %.4f) \n", 
	       i, plg[i].lines.count, (plg[i].vertical ? "Vert" : "Hori"), 
	       (plg[i].projective ? "VPoint" : "PLines"), 
	       plg[i].vp[0], plg[i].vp[1], plg[i].error[0], plg[i].error[1]);
      } else {
	printf("  parallel line group %d : %3d  %s  %s (%7.4f %7.4f)  error=(%.4f %.4f) \n", 
	       i, plg[i].lines.count, (plg[i].vertical ? "Vert" : "Hori"),
	       (plg[i].projective ? "VPoint" : "PLines"), 
	       plg[i].vp[0], plg[i].vp[1], plg[i].error[0], plg[i].error[1]);
      }
      pcount += plg[i].lines.count;
    }
    for (i=count=0; i < ldt->nlines; i++) if (ldt->lines[i]->plg == -1) count++;
    printf("  unmarked lines (%3d)  : %3d  \n", ldt->nlines - pcount, count);
    printf("  p merge: %3d lines replaced with %3d lines  (%4d)\n", count_m_old, count_m_new, count_m_new - count_m_old);
  }
  
  void showLines(void *image_buffer, pixel_t type, int count=0) {
    if (ldt == NULL) return;
    Drawer<double> ip( w, h, type, image_buffer );
    LDLine *lp;
    double *p0, *p1;
    int    i;
    if (count <= 0) {
      for (i=0; i < ldt->nlines; i++) {
	lp = ldt->lines[i];
	p0 = lp->p[0];  p1 = lp->p[1];  
	if (p0[1] < 0 || p1[1] < 0) { lp->printInfo(); continue; }
	ip.drawLine( p0, p1, 80, 80, 80 );
      }
    }
    int  r[] = { 255, 255,   0, 100, 200, 200,   0 };
    int  g[] = { 255,   0, 255, 100, 200,   0,   0 };
    int  b[] = { 255,   0,   0, 200,   0, 200, 255 };
    for (i = 0; i < nplg && (count <= 0 || i < count); i++)
      for (plg[i].lines.goFirst(); (lp = plg[i].lines.getCurr()); plg[i].lines.goNext()) {
	p0 = lp->p[0];  p1 = lp->p[1];  
	ip.drawLine( p0, p1, r[i], g[i], b[i] );
      }
  }
  
  LDLine* findNearestLine( double x, double y, bool line_segment=true ) {
    LDLine *lp, *lp_best = NULL;
    double   dist, min_dist = 9999;
    int      i;
    for (i=0; i < ldt->nlines; i++) {
      lp = ldt->lines[i];
      dist = lp->getDistanceToPoint( x, y, line_segment );
      if (dist < min_dist) { min_dist = dist; lp_best = lp; }
    }
    for (i = 0; i < nplg; i++) {
      for (plg[i].lines.goFirst(); (lp = plg[i].lines.getCurr()); plg[i].lines.goNext()) {
	dist = lp->getDistanceToPoint( x, y, line_segment );
	if (dist < min_dist) { min_dist = dist; lp_best = lp; }
      }
    }
    return lp_best;
  }
  
};

  
}	// namespace IMGH


#endif	// IMGH_LINE_PARALLELISM_HPP


