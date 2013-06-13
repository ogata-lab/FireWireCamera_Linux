
//
// IMGH::RectDetector 
//
// Jaeil Choi
// last modified in Dec, 2006
// GVU / RIM (Robotics & Intelligent Machines)
// College of Computing, Georgia Tech
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
//   - extracting rectangles from the extracted lines of a image
// This class requires:
//   - IMGH::Image		in 'imgh_common.hpp'
//   - IMGH::EdgeDetector	in 'imgh_edge_detector.hpp'
//   - IMGH::LineDetector	in 'imgh_line_detector.hpp'
//   - IMGH::LineParallelism	in 'imgh_line_parallelism.hpp'
//

#ifndef IMGH_RECT_DETECTOR_HPP
#define IMGH_RECT_DETECTOR_HPP

#include <iostream>
#include <cmath>
#include "imgh_common.hpp"
#include "imgh_line_detector.hpp"
#include "imgh_line_parallelism.hpp"
#include "algh_linked_list.hpp"
#include "algh_sorter.hpp"

namespace IMGH {
  
class RDRect {
public:
  double	ip[4][2];	// four corners 
  LDLine*	llp[4];		// left side, bottom, right side, top
  int		pllidx;		// index of the parallel model that it belongs to
  char		type;
public:
  RDRect() { 
    memset( ip,  0, 8*sizeof(double) );
    memset( llp, 0, 4*sizeof(LDLine*) );
    pllidx = -1;  type = ' ';
  }
  ~RDRect() {}
public:
  void setup(LDLine *ll, LDLine *lb, LDLine *lr, LDLine *lt, int pllidx) {
    IMGH4V_SET( llp,  ll, lb, lr, lt );
    this->pllidx = pllidx;
    getCorners( ip );
  }
  void getCorners(double ip[4][2], double center[2]=NULL) {
    // top left corner
    if (llp[0] && llp[3]) llp[0]->intersect( llp[3], ip[0] );
    else if (llp[0]) memcpy( ip[0], llp[0]->p[0], 2*sizeof(double) );
    else if (llp[3]) memcpy( ip[0], llp[3]->p[0], 2*sizeof(double) );
    else memset( ip[0], 0, 2*sizeof(double) );
    // bottom left corner
    if (llp[0] && llp[1]) llp[0]->intersect( llp[1], ip[1] );
    else if (llp[0]) memcpy( ip[1], llp[0]->p[1], 2*sizeof(double) );
    else if (llp[1]) memcpy( ip[1], llp[1]->p[0], 2*sizeof(double) );
    else memset( ip[1], 0, 2*sizeof(double) );
    // bottom right corner
    if (llp[1] && llp[2]) llp[1]->intersect( llp[2], ip[2] );
    else if (llp[1]) memcpy( ip[2], llp[1]->p[1], 2*sizeof(double) );
    else if (llp[2]) memcpy( ip[2], llp[2]->p[1], 2*sizeof(double) );
    else memset( ip[2], 0, 2*sizeof(double) );
    // bottom left corner
    if (llp[2] && llp[3]) llp[2]->intersect( llp[3], ip[3] );
    else if (llp[2]) memcpy( ip[3], llp[2]->p[0], 2*sizeof(double) );
    else if (llp[3]) memcpy( ip[3], llp[3]->p[1], 2*sizeof(double) );
    else memset( ip[3], 0, 2*sizeof(double) );
    if (center) {
      int  i, count;
      for (i=count=0, center[0]=center[1]=0; i < 4; i++)
	if (ip[i][0]!=0 || ip[i][1]!=0) { 
	  IMGH2V_ADD( center, center, ip[i] );
	  count++;
	}
      IMGH2V_DIV( center, count );
    }
  }
  void getCenter(double center[2]) {
    int  i, count;
    for (i=count=0, center[0]=center[1]=0; i < 4; i++)
      if (ip[i][0]!=0 || ip[i][1]!=0) { 
	IMGH2V_ADD( center, center, ip[i] );
	count++;
      }
    IMGH2V_DIV( center, count );
  }
  void printInfo(char *comment=NULL, int w=320) {
    printf("%s %x: \n", (comment ? comment : "Rect"), (unsigned int)this);
    printf("  (%3.0f,%3.0f)  %5d  (%3.0f,%3.0f)\n", 
	   ip[0][0], ip[0][1], (llp[3] ? llp[3]->head : 0), ip[3][0], ip[3][1]);
    printf("      %5d     %5s     %5d  \n", 
	   (llp[0] ? llp[0]->head : 0), "", (llp[2] ? llp[2]->head : 0));
    printf("  (%3.0f,%3.0f)  %5d  (%3.0f,%3.0f)\n",
	   ip[1][0], ip[1][1], (llp[1] ? llp[1]->head : 0), ip[2][0], ip[2][1]);
  }
};
  
class RectDetector
{
public:
  int			w, h, total;	// frequently used values
  ALGH::LinkedList<RDRect*>	rects;
  IMGH::LineDetector	*ldt;		// frequently used pointer
  IMGH::LineParallelism	*lpl;		// frequently used pointer
  bool		debug;
private:
  IMGH::LineDetector	ldt_dummy;
  IMGH::LineParallelism	lpl_dummy;
  
public:
  RectDetector() : w(0), h(0), total(0), 
		   ldt(&ldt_dummy), lpl(&lpl_dummy), debug(false) {}
  ~RectDetector() {}
  
  // -----------------------------------------------------------------
  // -----------------------------------------------------------------
  
public:  
  int findRects(int w, int h, pixel_t type, uchar_t *udata ) {
    ldt = &ldt_dummy;  ldt->findLines( w, h, type, udata );
    return findRects( ldt );
  }
  int findRects(LineDetector *ldt, LineParallelism *lpl=NULL) {
    // set up preliminary line information
    if (ldt == NULL) return 0;
    this->ldt = ldt;
    if (lpl) this->lpl = lpl;
    else { this->lpl = &lpl_dummy;  this->lpl->create(ldt); }
    this->w = ldt->w;  this->h = ldt->h;  this->total = w * h;
    rects.clearAndFree( );
    // get the list of vertical lines
    int      nvlines = lpl->plg[0].lines.count;
    LDLine  **vlines = (LDLine**)malloc( nvlines * sizeof(LDLine*) );
    lpl->plg[0].lines.getListInArray( vlines );
    // find rectangles
    int     i, j, k, pllidx, n, n2;
    LDLine  *vla, *vlb, *vlp=NULL, *hlines[40], *hlines2[40], *rlines[4];
    double  aloc[2], bloc[2], alen, blen, eclipse[2];
    ALGH::Sorter<LDLine*> sorter;
    RDRect  *rp=NULL;
    for (pllidx=1; pllidx < this->lpl->nplg; pllidx++) {
      for (i = 0; i < nvlines-1; i++) {
	vla = vlines[i];
	lpl->plg[pllidx].getLineLocation( vla, aloc );
	IMGH2V_SET( eclipse, +1e9, -1e9 );
	for (k = 0, vlp = NULL; k < vla->nplines; k++) 
	  if (vla->plines[k]->loc > vla->loc && 
	      (vlp == NULL || vla->plines[k]->tlen > vlp->tlen)) vlp = vla->plines[k];
	if (vlp && vlp->tlen < vla->tlen/2) vlp = NULL;
	for (j = i; j < nvlines; j++) {
	  vlb = vlines[j];
	  lpl->plg[pllidx].getLineLocation( vlb, bloc );
	  if (vlb->loc - vla->loc < 5) continue;  // too small
	  if (IMGH2V_DOT( vla->eq, vlb->eq ) > 0) continue;  // same direction
	  if (bloc[0] < eclipse[1]-1 && bloc[1] > eclipse[0] + 1 || 
	      bloc[1] > eclipse[0]+1 && bloc[0] < eclipse[1] - 1 ) continue; // eclipsed
	  // find a horizontal line that is shared by both 'vla' and 'vlb'
	  n = ldt->findLinesConnectingTwo( vla, vlb, hlines, 40 );
	  // get the list of lines that belongs to the same parallel model
	  for (k = n2 = 0; k < n; k++) 
	    if (hlines[k]->plg == pllidx) hlines2[n2++] = hlines[k];
	  // sort the horizontal lines
	  sorter.QuickSort( hlines2, n2, compare_lines_by_loc );
	  for (k = 0; k < n2; k++) {
	    if (vlp && vlp->findLinkTo( hlines2[k] )) continue;
	    if (k == 0) {	// look for a rectangle like ( |_| )
	      if ((alen = hlines2[0]->loc - aloc[0]) > h/10 &&
		  (blen = hlines2[0]->loc - bloc[0]) > h/10 &&
		  (vla->eq[0] * hlines2[0]->eq[1] < 0)) {
		IMGH4V_SET( rlines, vla, hlines2[0], vlb, NULL );
		tryToCompleteRectangle( rlines, aloc, bloc, pllidx );
		rp = new RDRect;
		rp->setup( rlines[0], rlines[1], rlines[2], rlines[3], pllidx );
		rects.append( rp );
		if (bloc[0] < eclipse[0]) eclipse[0] = bloc[0];
		if (hlines2[0]->loc > eclipse[1]) eclipse[1] = hlines2[0]->loc;
	      }
	    }
	    if (k == n2-1) {	// look for a rectangle like ( |^| )
	      if ((alen = aloc[1] - hlines2[k]->loc) > h/10 &&
		  (blen = bloc[1] - hlines2[k]->loc) > h/10 &&
		  (vla->eq[0] * hlines2[k]->eq[1] > 0) ) {
		IMGH4V_SET( rlines, vla, NULL, vlb, hlines2[k] );
		tryToCompleteRectangle( rlines, aloc, bloc, pllidx );
		rp = new RDRect;
		rp->setup( rlines[0], rlines[1], rlines[2], rlines[3], pllidx );
		rects.append( rp );
		if (bloc[1] > eclipse[1]) eclipse[1] = bloc[1];
		if (hlines2[k]->loc < eclipse[0]) eclipse[0] = hlines2[k]->loc;
	      }
	    } else {		// there's a rectangle with four sides
	      if (hlines2[k+1]->loc - hlines2[k]->loc > 5 &&
		  (vla->eq[0] * hlines2[k]->eq[1] > 0) && 
		  (vla->eq[0] * hlines2[k+1]->eq[1] < 0) ) {
		IMGH4V_SET( rlines, vla, hlines2[k+1], vlb, hlines2[k] );
		tryToCompleteRectangle( rlines, aloc, bloc, pllidx );
		rp = new RDRect;
		rp->setup( rlines[0], rlines[1], rlines[2], rlines[3], pllidx );
		rects.append( rp );
		if (hlines2[k+0]->loc < eclipse[0]) eclipse[0] = hlines2[k+0]->loc;
		if (hlines2[k+1]->loc > eclipse[1]) eclipse[1] = hlines2[k+1]->loc;
	      }
	    }
	  }
	}	// end of the loop for 'vlb'
      }		// end of the loop for 'vla'
    }		// end of the loop for parallel models
    free( vlines );
    return rects.count;
  }
  
  void getRelativeLinePosition(LDLine *vlp, LDLine *hlp, double vll[2], double hll[2]) {
    double *hp0, *hp1;
    if (hlp->p[0][0] < hlp->p[1][0]) {  hp0 = hlp->p[0];  hp1 = hlp->p[1]; }
    else {  hp0 = hlp->p[1];  hp1 = hlp->p[0]; }
    vll[0] = fabs( hlp->getDistanceToPoint( vlp->p[0][0], vlp->p[0][1] ));
    vll[1] = fabs( hlp->getDistanceToPoint( vlp->p[1][0], vlp->p[1][1] ));
    hll[0] = fabs( vlp->getDistanceToPoint( hp0[0], hp0[1] ));
    hll[1] = fabs( vlp->getDistanceToPoint( hp1[0], hp1[1] ));
  }
  
  bool tryToCompleteRectangle(LDLine *lines[4], double aloc[2], double bloc[2], int pllidx) {
    // Complete a rectangle with three lines (two sides and either top or bottom line)
    LDLine *hl, *hl2 = NULL;
    double hloc[2];
    int    i;
    if        (lines[1] == NULL) {	// rectangle like ( |^| )
      if (aloc[1] > bloc[1]) {		// find bottom line from left line
	for (i=0; i < lines[0]->nlinks; i++) {
	  hl = lines[0]->links[i].lp;
	  if (hl->plg != pllidx || hl->loc < bloc[1]) continue;
	  lpl->plg[0].getLineLocation( hl, hloc );
	  if (hloc[1] >= lines[2]->loc-2 &&
	      (hl2 == NULL || hl->loc < hl2->loc) ) hl2 = hl;
	}
      } else {                		// find bottom line from right line
	for (i=0; i < lines[2]->nlinks; i++) {
	  hl = lines[2]->links[i].lp;
	  if (hl->plg != pllidx || hl->loc < aloc[1]) continue;
	  lpl->plg[0].getLineLocation( hl, hloc );
	  if (hloc[0] <= lines[0]->loc+2 &&
	      (hl2 == NULL || hl->loc < hl2->loc) ) hl2 = hl;
	}
      }
      lines[1] = hl2;
    } else if (lines[3] == NULL) {	// rectangle like ( |_| )
      if (aloc[0] < bloc[0]) {		// find top line from left line
	for (i=0; i < lines[0]->nlinks; i++) {
	  hl = lines[0]->links[i].lp;
	  if (hl->plg != pllidx || hl->loc > bloc[0]) continue;
	  lpl->plg[0].getLineLocation( hl, hloc );
	  if (hloc[1] >= lines[2]->loc-2 &&
	      (hl2 == NULL || hl->loc > hl2->loc) ) hl2 = hl;
	}
      } else {                		// find top line from right line
	for (i=0; i < lines[2]->nlinks; i++) {
	  hl = lines[2]->links[i].lp;
	  if (hl->plg != pllidx || hl->loc > aloc[0]) continue;
	  lpl->plg[0].getLineLocation( hl, hloc );
	  if (hloc[0] <= lines[0]->loc+2 &&
	      (hl2 == NULL || hl->loc > hl2->loc) ) hl2 = hl;
	}
      }
      lines[3] = hl2;
    }
    return (hl2 != NULL);
  }
  
public:
  void printInfo(void) {
    if (ldt == &ldt_dummy) ldt->printInfo();
    if (lpl == &lpl_dummy) lpl->printInfo();
    printf("IMGH::RectDetector \n");
    printf("  rects  (linked list)  : %3d\n", rects.count);
  }
  
  void showRects(void *image_buffer, pixel_t type) {
    // Visualize the rectangles on given image,
    //   assumming the image buffer has the same size, only with its own pixel type.
    Drawer<double> imgd( w, h, type, image_buffer );
    double  ip[4][2];
    int     r, g, b;
    RDRect *rp;
    for (rects.goFirst(); (rp = rects.getCurr()); rects.goNext()) {
      //rp->getCorners( w, ip[0], ip[1], ip[2], ip[3] );
      rp->getCorners( ip );
      srand( (unsigned int)rp );
      r = 50+rand()%206;  g = 50+rand()%206;  b = 50+rand()%206;
      for (int i = 0; i < 4; i++) {
	if (rp->llp[i]) imgd.drawLine( ip[i], ip[(i+1)%4], r, g, b );
	else            imgd.drawLine( ip[i], ip[(i+1)%4], r/3, g/3, b/3 ); 
      }
    }
  }
  
  void showRect(void *image_buffer, pixel_t type, RDRect *rp) {
    // Visualize the rectangles on given image,
    //   assumming the image buffer has the same size, only with its own pixel type.
    Drawer<double> imgd( w, h, type, image_buffer );
    double  ip[4][2];
    int     r, g, b;
    //rp->getCorners( w, ip[0], ip[1], ip[2], ip[3] );
    rp->getCorners( ip );
    srand( (unsigned int)rp );
    r = 50+rand()%206;  g = 50+rand()%206;  b = 50+rand()%206;
    for (int i = 0; i < 4; i++) {
      if (rp->llp[i]) imgd.drawLine( ip[i], ip[(i+1)%4], r, g, b );
      else            imgd.drawLine( ip[i], ip[(i+1)%4], r/3, g/3, b/3 ); 
    }
  }
  
  RDRect* findNearestRect(int x, int y, double *cdist) {
    RDRect *rp, *min_rp=NULL;
    double dist, min_dist=1e9, xy[2]={x,y}, ip[4][2], center[2];
    for (rects.goFirst(); (rp = rects.getCurr()); rects.goNext()) {
      rp->getCorners( ip, center );
      dist = IMGH2V_DIST( xy, center );
      if (dist < min_dist) { min_dist = dist; min_rp = rp; }
    }
    if (cdist) *cdist = min_dist;
    return min_rp;
  }
  
};

  
}	// namespace IMGH


#endif	// IMGH_RECT_DETECTOR_HPP


