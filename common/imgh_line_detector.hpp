
//
// IMGH::LineDetector
//   Line Detector with pixel map and link information
//
// Jaeil Choi
// last modified in Sep, 2009
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
//   - extracting straight lines from edges
//   - capturing intersections and connections between lines
//   - maintaining a pixel map for lines
// This class requires:
//   - IMGH::Image		in 'imgh_common.hpp'
//   - IMGH::EdgeDetector	in 'imgh_edge_detector.hpp'
//

#ifndef IMGH_LINE_DETECTOR_HPP
#define IMGH_LINE_DETECTOR_HPP

#include <iostream>
#include <cmath>
#include "imgh_common.hpp"
#include "imgh_drawer.hpp"
#include "imgh_edge_detector.hpp"
#include "algh_sorter.hpp"
#include "mtx_matrix_solver.hpp"
#ifndef ROUND
#define ROUND(a) ((int)((a)+0.5f))
#endif

#ifndef M_PI
#define M_PI (3.14159265358979323846)
#endif

namespace IMGH {
  
// ===================================================================
// Links between lines; a line has a list of links
// ===================================================================
class LDLine;
class LDLink {
public:
  int		pidx;	// the pixel of the intersection
  LDLine	*lp;	// 
public:
  LDLink() : pidx(0), lp(NULL) {}
  LDLink(int pidx, LDLine *lp) { this->pidx = pidx; this->lp = lp; }
  ~LDLink() {}
public:
  void setup(int pidx, LDLine *lp) { this->pidx = pidx; this->lp = lp; }
  void operator=(LDLink &lk) { this->pidx = lk.pidx; this->lp = lk.lp; }
};
  
// ===================================================================
// Line; 
// ===================================================================
class LDLine {
public:
  int    head, tail;	// the first and the last pixel of the line (edge)
  int    tlen;		// pixel count of the thread
  float  pp[2][2];	// starting/ending positions (may not agree with 'head'/'tail')
  float  eq[3];		// line equation  ax+by+c=0  (in image coordinates, (0,0) at TL)
  float  extra[2];	// extra information for the use of application
  
  int		nplines;	// the number of parallel lines
  LDLine	**plines;	// parallel lines
  int		nlinks;		// the number of links
  LDLink	*links;		// links on the line
  float		loc;	// projected location calculated by PLGroup::getLocation()
  int		tag;	// temporary information (set and use it only in a single function)
  char		plg;	// PLGroup index for IMGH::LineParallelism class
  bool		merged;

public:
  LDLine() 
    : head(0), tail(0), tlen(0), 
      nplines(0), plines(NULL), nlinks(0), links(NULL), loc(0), plg(-1), merged(false) {}
  LDLine(double x0, double y0, double x1, double y1) 
    : head(0), tail(0), tlen(0), 
      nplines(0), plines(NULL), nlinks(0), links(NULL), loc(0), plg(-1), merged(false) { 
    setEndPoints(x0, y0, x1, y1); 
  }
  LDLine(LDLine *lp) 
    : head(0), tail(0), tlen(0), 
      nplines(0), plines(NULL), nlinks(0), links(NULL), loc(0), plg(-1), merged(false) { 
    head = lp->head;  tail = lp->tail;  tlen = lp->tlen;
    memcpy( pp, lp->pp, 4*sizeof(float) );
    memcpy( eq, lp->eq, 3*sizeof(float) );
  }
  ~LDLine() { if (plines) free(plines);  if (links) free(links); }
  
  bool isInvalid(void) { return (tlen <= 1); }
  void invalidate(void) { setEndPoints(0,0,0,0, false); }
  
  void setEndPoints(double x0, double y0, double x1, double y1, bool update_eq=true) { 
    pp[0][0] = (float)x0; pp[0][1] = (float)y0;
    pp[1][0] = (float)x1; pp[1][1] = (float)y1;
    double len, dir[2] = {(x1-x0), (y1-y0)};
    if (dir[0] > dir[1]) { tlen = (int)dir[0]+1; }
    else                 { tlen = (int)dir[1]+1; }
    if (tlen == 0) IMGH3V_SET( eq, 0, 0, 0 );
    else if (update_eq) {
      eq[0] = (float)-dir[1];  eq[1] = (float)dir[0];  len = sqrt(eq[0]*eq[0] + eq[1]*eq[1]);
      if (len > 0) { eq[0] /= (float)len;  eq[1] /= (float)len; }
      eq[2] = (float) -( eq[0] * x0 + eq[1] * y0 );
    }
  }
  
  // -----------------------------------------------------------------

  double getDistanceToPoint(double x, double y, bool line_segment=false) {
    if (!line_segment) {	// infinite line
      // positive if the point is in the left of the line direction (positive rotation)
      return ( eq[0] * x + eq[1] * y + eq[2] );
    } else {			// line segment (distance is non-negative)
      double v01[2], v10[2], v0p[2], v1p[2];
      IMGH2V_SUB( v01, pp[1], pp[0] );  IMGH2V_SUB( v10, pp[0], pp[1] );
      v0p[0] = x-pp[0][0];  v0p[1] = y-pp[0][1];
      v1p[0] = x-pp[1][0];  v1p[1] = y-pp[1][1];
      if        ( IMGH2V_DOT( v0p, v01 ) < 0 ) {
	return sqrt((x-pp[0][0])*(x-pp[0][0]) + (y-pp[0][1])*(y-pp[0][1]));
      } else if ( IMGH2V_DOT( v1p, v10 ) < 0 ) {
	return sqrt((x-pp[1][0])*(x-pp[1][0]) + (y-pp[1][1])*(y-pp[1][1]));
      } else {
	return fabs( eq[0] * x + eq[1] * y + eq[2] );
      }
    }
  }
  void getPointOnLine(double x, double y, double nxy[2], int wh[2]=NULL) {
    // find the closest point (double) on the line, given the line equation
    double  v[2], vo[2], len;
    double  vl[2] = { +eq[1], -eq[0] };	// vl[2] : unit vector of the line direction
    // get a point 'v[2]' on the line
    if (fabs(eq[0]) > fabs(eq[1])) { v[1] = y;  v[0] = -(eq[1]*v[1]+eq[2])/eq[0]; }  // ax = - by - c
    else                           { v[0] = x;  v[1] = -(eq[0]*v[0]+eq[2])/eq[1]; }  // by = - ax - c
    // project the input point to line using 'v[2]'
    vo[0] = x - v[0];  vo[1] = y - v[1];		// vo[] = xy[] - v[]
    len = IMGH2V_DOT( vo, vl );				// len  = (vo . vl)
    IMGH2V_SADD( nxy, v, len, vl );			// nxy[] = v[] + len * vl[]
    if (wh) {
      if      (nxy[0]<    0 ) { nxy[0] = 0;        nxy[1] = -(eq[0]*nxy[0]+eq[2])/eq[1]; }
      else if (nxy[0]>=wh[0]) { nxy[0] = wh[0]-1;  nxy[1] = -(eq[0]*nxy[0]+eq[2])/eq[1]; }
      else if (nxy[1]<    0 ) { nxy[1] = 0;        nxy[0] = -(eq[1]*nxy[1]+eq[2])/eq[0]; }
      else if (nxy[1]>=wh[1]) { nxy[1] = wh[1]-1;  nxy[0] = -(eq[1]*nxy[1]+eq[2])/eq[0]; }
    }
  }
  int samplePointsAlongLine(double points[]) {
    // Sample points along the line in 'points[]' array,
    //   and return the number of samples
    // Note that 'points[]' will have the size of 2*(ReturnValue).
    double xy[2], dx, dy;
    if (fabs(eq[0]) > fabs(eq[1])) { dy = 1.0;  dx = - (eq[1] / eq[0]); }
    else                           { dx = 1.0;  dy = - (eq[0] / eq[1]); }
    getPointOnLine( pp[0][0], pp[0][1], xy );
    for (int i = 0; i < tlen && i < 1024; i++, points+=2, xy[0]+=dx, xy[1]+=dy) {
      points[0] = xy[0];  points[1] = xy[1];
    }
    return tlen;
  }
  bool intersect(LDLine *lp, double xp[2], int tolerance=0, bool debug=false) {
    xp[0] = xp[1] = 0;
    double det = ( eq[0] * lp->eq[1] - eq[1] * lp->eq[0] );
    if (fabs(det) < 0.0001) return false;
    xp[0] = - (+ lp->eq[1] * eq[2] - eq[1] * lp->eq[2]) / det;
    xp[1] = - (- lp->eq[0] * eq[2] + eq[0] * lp->eq[2]) / det;
    if (tolerance <= 0) return true;	// intersection of infinite lines
    double dist, d0, d1;
    dist = fabs(pp[0][0] - pp[1][0]) + fabs(pp[0][1] - pp[1][1]) + tolerance;
    d0   = fabs(pp[0][0] - xp[0]) + fabs(pp[0][1] - xp[1]);
    d1   = fabs(pp[1][0] - xp[0]) + fabs(pp[1][1] - xp[1]);
    if (d0 > dist || d1 > dist) return false;
    dist = fabs(lp->pp[0][0] - lp->pp[1][0]) + fabs(lp->pp[0][1] - lp->pp[1][1]) + tolerance;
    d0   = fabs(lp->pp[0][0] - xp[0]) + fabs(lp->pp[0][1] - xp[1]);
    d1   = fabs(lp->pp[1][0] - xp[0]) + fabs(lp->pp[1][1] - xp[1]);
    if (d0 > dist || d1 > dist) return false;
    return true;
  }
  void getRelativePosition(LDLine *lp, double d[5], bool more) {
    // calculate relative position of itself against 'lp'
    d[0] = lp->getDistanceToPoint(pp[0][0], pp[0][1]);  // dist from pp[0] to infinite lp
    d[1] = lp->getDistanceToPoint(pp[1][0], pp[1][1]);  // dist from pp[1] to infinite lp
    if (!more) return;
    double dir[2], a0[2], a1[2];
    IMGH2V_SUB( dir, lp->pp[1], lp->pp[0] );  
    d[4] = IMGH2V_LEN( dir );  if (d[4]>0) IMGH2V_DIV( dir, d[4] );  // length of lp
    IMGH2V_SUB( a0,  pp[0], lp->pp[0] );
    IMGH2V_SUB( a1,  pp[1], lp->pp[0] );
    d[2] = IMGH2V_DOT( dir, a0 );  // dist of pp[0] from lp->pp[0] in the direction of lp
    d[3] = IMGH2V_DOT( dir, a1 );  // dist of pp[1] from lp->pp[0] in the direction of lp
  }
  
  // -----------------------------------------------------------------

  LDLink* findLinkTo(LDLine *lp) {
    for (int i=0; i<nlinks; i++) if (links[i].lp == lp) return links + i;
    return NULL;
  }
  bool removeLinkTo(LDLine *lp) {
    int  i, j = -1;
    for (i = 0; i < nlinks; )
      if (links[i].lp == lp) {
	for (j = i; j < nlinks-1; j++) links[j] = links[j+1];
	nlinks--;
      } else i++;
    return (j >= 0 ? true : false);
  }
  bool replaceLinkTo(LDLine *lp_old, LDLine *lp_new) {
    bool ret = false;
    for (int i = 0; i < nlinks; i++)
      if (links[i].lp == lp_old) { ret=true; links[i].lp = lp_new; }
    return ret;
  }
  void removeLinks(bool cautious) {
    if (cautious) for (int i=0; i<nlinks; i++) links[i].lp->removeLinkTo(this);
    if (nlinks && links) { free(links); links = NULL; }  nlinks = 0; 
  }
  bool isParallelWith(LDLine *lp) { 
    for (int i=0; i<nplines; i++) if (plines[i]==lp) return true;
    return false;
  }
  void removePLines(void) { if (plines) free(plines);  plines = NULL;  nplines = 0; }
  void removePLineTo(LDLine *lp) {
    for (int i = 0; i < nplines; i++)
      if (plines[i] == lp) { plines[i] = plines[nplines-1];  nplines--;  break; }
  }
  void replacePLineTo(LDLine *lp_old, LDLine *lp_new) {
    for (int i = 0; i < nplines; i++)
      if (plines[i] == lp_old) plines[i] = lp_new;
  }
  void takeLinksFrom(LDLine *lp) {
    // Take all the links of 'lp', leaving it with 0 links. (for merging lines)
    int i, total = this->nlinks + lp->nlinks;
    if (total < 1) return;
    LDLink *llinks = (LDLink*) malloc( total * sizeof(LDLink) );
    memcpy( llinks,  this->links, this->nlinks * sizeof(LDLink) );
    memcpy( llinks + this->nlinks, lp->links, lp->nlinks * sizeof(LDLink) );
    for (i=nlinks; i<total; i++) llinks[i].lp->replaceLinkTo(lp, this);
    this->removeLinks(false);  lp->removeLinks(false);
    this->links = llinks;   this->nlinks = total;
  }
  
  // -----------------------------------------------------------------

  void printInfo(char* comment=NULL, char debug_type=' ') { 
    int  i, j, k;
    switch (debug_type) {
    case 'k':	// show link information
      printf("%s (%2d%2d) plg(%2d %6.1f) ", (comment ? comment : "LDLine:"), 
	     nplines, nlinks, plg, loc);  fflush(stdout);
      printf("head=%5d  t=%d  links(%x)[%d]: ", head, tag, (unsigned int)links, nlinks);  fflush(stdout);
      for (j=0; j < nlinks; j++) { printf("%d ", links[j].lp->head); fflush(stdout); }
      printf("\n");
      break;
    case 'p':	// show parallelism information
      printf("%s (%2d%2d) plg(%2d %6.1f) ", (comment ? comment : "LDLine:"), 
	     nplines, nlinks, plg, loc);  fflush(stdout);
      printf("head=%5d  t=%d  plines(%x)[%d]: ", head, tag, (unsigned int)plines, nplines);  fflush(stdout);
      for (j=0; j < nplines && plines; j++) {
	if (!plines[j]) { printf("NULL ");  fflush(stdout); }
	else {
	  for (k=0; k < plines[j]->nplines; k++) if (plines[j]->plines[k] == this) break;
	  printf("%d%s ", plines[j]->head, (k == plines[j]->nplines ? "X" : ""));
	}
      }
      printf("\n");
      break;
    default:	// 
//       printf("%s head=%5d (%5.1f,%5.1f)-(%5.1f,%5.1f) len=%3d eq(%7.4f %7.4f %9.4f) (pl=%d lk=%d)\n", 
// 	     (comment ? comment : "Line:"), head, pp[0][0], pp[0][1], pp[1][0], pp[1][1], tlen, eq[0], eq[1], eq[2], nplines, nlinks ); 
      printf("%s: len=%3d  (%3.0f,%3.0f)-(%3.0f,%3.0f) eq(%7.4f %7.4f %9.4f) (pl=%d lk=%d) (h=%d,t=%d)\n", 
	     (comment ? comment : "Line"), tlen, pp[0][0], pp[0][1], pp[1][0], pp[1][1], eq[0], eq[1], eq[2], nplines, nlinks, head, tail ); 
      for (i = 0; i < nlinks; i++) 
	if (!links[i].lp->findLinkTo(this)) { printf(" (LinkERR %d,%x)", i, (unsigned int)links[i].lp); break; }
      break;
    }
  }
  
};
  
// static int imgline_compare_by_len(LDLine *a, LDLine *b) { return (int) -(a->tlen - b->tlen); }  // descending
inline int compare_lines_by_len(LDLine *a, LDLine *b) { return (a->tlen - b->tlen); }  // ascending
inline int compare_lines_by_loc(LDLine *a, LDLine *b) { return (int) +(a->loc - b->loc); }    // ascending
inline int compare_lines_by_tag(LDLine *a, LDLine *b) { return (int) +(a->tag - b->tag); }    // ascending
  
}	// namespace IMGH


// ===================================================================
// ===================================================================
  
namespace IMGH {
  
class LineDetector
{
public:
  int		w, h, total;		// frequently used values
  IMGH::LDLine		**lines;	// the list   of LDLine instances
  int			nlines;		// the number of LDLine instances
  IMGH::Image		lmap;		// map of the pixels all the lines (PIXEL_VOIDP)
  IMGH::EdgeDetector	*edt;
private:
  IMGH::EdgeDetector	edt_dummy;
  IMGH::Image	*gx, *gy, *gz;		// frequently used pointers (PIXEL_FLOAT)
  int	count_init, count_removed;
  int	count_im_old, count_im_new;
  int	count_pm_old, count_pm_new;
public:
  int   min_length;
  bool	debug;
  bool	verbose;
    
public:
  LineDetector() 
    : lines(NULL), nlines(0), edt(NULL), gx(NULL), gy(NULL), gz(NULL), debug(false), verbose(false) {}
  LineDetector(Image *img) 
    : lines(NULL), nlines(0), edt(NULL), gx(NULL), gy(NULL), gz(NULL), debug(false), verbose(false) { 
    findLines(img->w, img->h, img->type, img->data); 
  }
  ~LineDetector() {}
  void  clear(void) { 
    for (int i=0; i < nlines; i++) delete( lines[i] );
    free( lines );  lines = NULL;
    nlines = 0;
  }
  
  // -----------------------------------------------------------------
  // 
  // -----------------------------------------------------------------
public:  
  void findLines(int w, int h, pixel_t type, uchar_t *udata, int min_length=10,
		 bool get_parallelism=false, bool get_links=false) {
    edt_dummy.findEdges( w, h, type, udata );
    findLines( &edt_dummy, min_length, get_parallelism, get_links );
  }
  void findLines(IMGH::EdgeDetector *edt, int min_length=10, 
		 bool get_parallelism=false, bool get_links=false) {
    clear();
    // check required data (gx, gy and gz)
    if (edt == NULL || !edt->emap.sameSize( &edt->grd->gz )) return;
    this->edt = edt;
    this->gx = &(edt->grd->gx);  this->gy = &(edt->grd->gy);  this->gz = &(edt->grd->gz);
    this->w  = gx->w;  this->h = gx->h;  this->total = w * h;
    this->min_length = min_length;
    
    // create lines using the result of edge detector
    if (verbose) printf("[LineD] creating lines from edges\n");
    createLines( min_length );
    
    // merge lines that are connected sequentially in a big straight line
    if (verbose) printf("[LineD] merging connected straight lines\n");
    mergeLinesConnectedStraight();
    
    // sort the list of LDLine*
    if (verbose) printf("[LineD] sorting lines in the order of descending length\n");
    ALGH::Sorter<LDLine*> sorter;  // sort them in the ascending order of length
    sorter.QuickSort( lines, nlines, compare_lines_by_len );
    
    // find parallel neighbors for each line
    if (get_parallelism) {
      if (verbose) printf("[LineD] finding parallel lines\n");
      findParallelLines();
      //mergeUsingParallelism();
    }
    
    // remove invalidated lines;
    if (verbose) printf("[LineD] removing invalidated lines\n");
    removeInvalidatedLines();
    
    // find links to other lines (only after 'removeInvalidatedLines()' called)
    if (get_links) {
      if (verbose) printf("[LineD] finding links between lines\n");
      findLineLinks();
    }
  }
  
  // -----------------------------------------------------------------
  // converting edge points to lines
  // -----------------------------------------------------------------

private:
  int  createLines(int min_length=10) {
    lines = (LDLine**)malloc( edt->edge_count * sizeof(LDLine*) );
    int    i, idx, count=0, np;
    int    *plist = (int*)malloc( (w+h) * sizeof(int) );
    double *xy, *xylist = (double*) malloc( (w+h) * 2 * sizeof(double) );
    lmap.setImage( w, h, PIXEL_VOIDP, true );
    int    *epp = (int*)(edt->emap.data);
    for (idx = 0; idx < total; idx++) {
      //debug = (idx == 915);
      if (idx == 0 || idx != epp[idx] || lmap.getPointer(idx)) continue;
      np = edt->getEdgePoints(idx, plist, (w+h));
      //if (debug) printf("  creating a line for head=%d : len=%d\n", idx, np);
      if (np < min_length) continue;
      for (i=0, xy=xylist; i<np; i++, xy+=2) { xy[0] = plist[i]%w; xy[1] = plist[i]/w; }
      count = createLineRecursively( plist, xylist, np );
    }
    free( plist );  free( xylist );
    for (i=0; i<total; i++)		// complete line pointer map 
      if (epp[i] && epp[i] != i) lmap.setPointer( i, lmap.getPointer(epp[i]) );
    count_removed = 0;
    count_init = nlines;
    return count;
  }
  
  int createLineRecursively(int *plist, double *xylist, int np) {
    LDLine line, *lp;
    //if (debug) printf("  creating a line for head=%d with %d points\n", plist[0], np);
    if (fitLine( &line, plist, xylist, np )) {
      //if (debug) { printf("    1: "); line.printInfo(); }
      lp = new LDLine( &line );		// create the line
      lines[nlines++] = lp;  lmap.setPointer(lp->head, lp);
      return 1;
    } else if ( breakIntoTwoLines( &line, plist, xylist, np, false ) ||
		breakIntoTwoLines( &line, plist, xylist, np, true  ) ) {
      //if (debug) { printf("    2: "); line.printInfo(); }
      lp = new LDLine( &line );		// create the line
      lines[nlines++] = lp;  lmap.setPointer(lp->head, lp);
      int  i, j;			// update edge map
      for (i=0; i<np; i++) if (plist[i]==line.head) break;  
      for (j=0; j<line.tlen; j++) edt->emap.setInt( plist[i+j], plist[i] );
      return 1 + createLineRecursively( plist, xylist, np - line.tlen );
    } else {
      //if (debug) { printf("    3: fitting failed - a curve?\n"); }
      return 0;
    }
  }
  
  bool breakIntoTwoLines(LDLine *lp, int *plist, double *xylist, int np, bool brutal) {
    // find a line for the later part of the list, save the information in 'lp',
    if (brutal == false) {	// fast heuristic --------------------
      int    i, trial, idx, straight_count=0;
      float  ogv[2], cgv[2], straight_gv[2]={0,0};
      bool   straight = false;
      edt->grd->getGradientVector( np-1, ogv, true );
      float  c_same[2] = { 0.985f, 0.939f };   // (less than) 10 or 20 degree
      float  c_diff[2] = { 0.906f, 0.707f };   // (more than) 25 or 45 degree
      for (trial = 0; trial < 2; trial++) {
	for (i = np-2; i >= 0; i--) {	// scan backward
	  idx = plist[i];
	  edt->grd->getGradientVector( idx, cgv, true );
	  if (IMGH2V_DOT(ogv, cgv) > c_same[trial]) {
	    if (++straight_count > 5) { 
	      straight = true;
	      IMGH2V_COPY( straight_gv, ogv );
	      straight_count = 0;
	    }
	  } else {
	    straight_count = 0;
	    IMGH2V_COPY( ogv, cgv );
	  }
	  if (straight && IMGH2V_DOT(straight_gv, cgv) < c_diff[trial]) break;
	}
	//if (debug) printf("    break (quick ) %3d points starting at %3d => %s  (break at %d (%d,%d) with %d points)\n", np, plist[0], (i>0 ? "Yes" : "No"), plist[i+1], plist[i+1]%w, plist[i+1]/w, np-i-1); 
	if      (i > 0) break;	// found the breaking point
	else if (straight && trial==0 && np>100) { i = np/2; break; }
      }
      if (i <= 0) return false;
      bool ret = fitLine( lp, plist+(i+1), xylist+2*(i+1), np-(i+1) );
      //if (debug) printf("    break (quick ) fitLine => %s  np=%d  i=%d\n", (ret ? "Y":"n"), np, i);
      if        (!ret && np > 120 && i < np/2) {
	i = np*2/3; ret = fitLine( lp, plist+(i+1), xylist+2*(i+1), np-(i+1) );
      } else if (!ret && np >  80 && i < np/2) {
	i = np/2;   ret = fitLine( lp, plist+(i+1), xylist+2*(i+1), np-(i+1) );
      }
      return ret;
    } else {			// brutal search ---------------------
      for (int len = 4; len < np-8; len++) {
	if (!fitLine( lp, plist+(np-len), xylist+2*(np-len), len )) return false;
	double *p, d[4];
	p = xylist+2*(np-len-2);  d[0] = lp->getDistanceToPoint( p[0], p[1] );
	p = xylist+2*(np-len-3);  d[1] = lp->getDistanceToPoint( p[0], p[1] );
	p = xylist+2*(np-len-4);  d[2] = lp->getDistanceToPoint( p[0], p[1] );
// 	p = xylist+2*(np-len-4);  d[4] = lp->getDistanceToPoint( p[0], p[1] );
	bool ret = ( d[0] > +1 && d[1] > d[0] && d[2] > d[1] ||
		     d[0] < -1 && d[1] < d[0] && d[2] < d[1] );
	//if (debug) printf("    break (brutal) %3d points starting at %3d => %c  (break at %d with %d points d[]=%.1f %.1f %.1f)\n", np, plist[0], (ret ? 'Y':'n'), plist[np-len], len, d[0], d[1], d[2]); 
	if (ret) return true;
      }
      return false;
    }
  }
  
public:
  bool fitLine(LDLine *lp, int *plist, double *xylist, int np) {
    //// calculate error and eq[]
    // Optimize a line equation for sequential list of points,
    //   (XY coordinates), and return the error measure of the optimization.
    if (!xylist || np == 0) return false;
    double mx=0, my=0, *xy, *a, *A = (double*)malloc(2*np*sizeof(double)), ratio;
    int    i;
    for (i=0, xy=xylist; i < np; i++, xy+=2) { mx += xy[0]; my += xy[1]; }
    mx /= np;   my /= np;
    for (i=0, xy=xylist, a=A; i < np; i++, xy+=2, a+=2) {
      a[0] = xy[0] - mx;  a[1] = xy[1] - my;
    }
    // solve the coefficients of the line equation using least square
    MTX::MatrixSolver<double> solver;  // min(A*eq[2]) (least square)
    double eq[3];
    solver.solveNullBySVD( np, 2, A, eq, NULL, &ratio );
    eq[2] = - ( eq[0] * mx + eq[1] * my );
    IMGH3V_SET( lp->eq, (float)eq[0], (float)eq[1], (float)eq[2] );
    free(A);
    //if (debug) printf("a ratio=%.4f  3.5/np=%.4f\n", ratio, 3.5/np);
    if (ratio > 3.0/np) return false;	// not a straight line
    // normalize the line equation
    float len = (float)IMGH2V_LEN( lp->eq );
    if (len > 0) IMGH3V_DIV( lp->eq, len );
    // set up the line
    int    stt_pos=0, end_pos=2*(np-1), wh[2]={w,h};
    double stt[2], end[2];
    lp->getPointOnLine( xylist[stt_pos+0], xylist[stt_pos+1], stt, wh );
    lp->getPointOnLine( xylist[end_pos+0], xylist[end_pos+1], end, wh );
    lp->setEndPoints( stt[0], stt[1], end[0], end[1], false );
    lp->head = plist[0];  lp->tail = plist[np-1];  lp->tlen = np;
    // check end points
    double d0 = fabs(IMGH2V_LPDIST( lp->eq, xylist+stt_pos ));
    double d1 = fabs(IMGH2V_LPDIST( lp->eq, xylist+end_pos ));
    double dmax = 2.5 * np;  //(np>100 ? np/100.0 : 1.0);
    //if (debug) printf("a  d0=%.2f d1=%.2f  dmax=%.2f\n", d0, d1, dmax);
    if (d0 > dmax || d1 > dmax) return false;  // end points are too far away
    // decide the direction
    double gv[2] = {0,0};
    for (i = 0; i < np; i++) {
      gv[0] += gx->getFloat( plist[i] );
      gv[1] += gy->getFloat( plist[i] );
    }
    if (IMGH2V_DOT(lp->eq, gv) < 0) IMGH3V_MUL( lp->eq, -1 );
    return true;
  }
  
  
  // -----------------------------------------------------------------
  // -----------------------------------------------------------------

  void mergeLinesConnectedStraight(void) {
    // Merge lines that are connected sequentially in a big straight line
    // Note that this function should be called before finding
    //   parallel neighbor lines or links to other lines.
    int    checked_already=746291660;  // magic number to mark checked lines
    int    i, j, k, nl, np, cnt, total, pidx;
    int    *tlist = (int*)   malloc(     (w+h) * sizeof(int)   );
    int    *plist = (int*)   malloc(     (w+h) * sizeof(int)   );
    double *clist = (double*)malloc( 2 * (w+h) * sizeof(double));
    LDLine *lp, *cllist[40];
    LDLine tmp;
    count_im_old = count_im_new = 0;
    for (i = 0; i < nlines; i++) lines[i]->tag = 0;
    for (i = 0; i < nlines; i++) {
      lp = lines[i];
      if (lp->isInvalid() || lp->tag == checked_already) continue;
      //if (lp->tail == 111715) printf("A line(%d-%d) \n", lp->head, lp->tail);
      // find a list of lines that are connected sequentially
      nl = findLinesConnectedHeadToTail( lp, cllist, 40 );
      if (nl < 2) continue;
      int stt=0, end=(nl-1);
      bool good=false;
      while (!good && stt<end) {
	// get a big new line that is optimized for all the lines
	for (j = stt, cnt = 0; j <= end; j++) {
	  //if (cllist[j]->tail == 111715) { printf("B connected with "); lp->printInfo(); }
	  np = getLinePoints( cllist[j]->head, tlist, (w+h) );
	  for (k = 0; k < np; k++, cnt+=2) {
	    plist[cnt/2] = tlist[k];		// pixel index
	    clist[cnt+0] = tlist[k] % w;	// pixel X position
	    clist[cnt+1] = tlist[k] / w;	// pixel Y position
	  }
	}
	total = cnt / 2;
	fitLine( &tmp, plist, clist, total );
	// check whether or not the new line is good
	for (j = stt; j <= end; j++)
	  if ( fabs(IMGH2V_LPDIST( tmp.eq, cllist[j]->pp[0] )) > 1 ||
	       fabs(IMGH2V_LPDIST( tmp.eq, cllist[j]->pp[1] )) > 1 ) break;
	if      (j >end) { good = true; break; }
	else if (j==stt && (end-stt)>=2) { stt++; continue; }
	else if (j==end && (end-stt)>=2) { end--; continue; }
	else { good = false; break; }
      }
      if (!good) continue;
      for (j = stt; j <= end; j++) cllist[j]->tag = checked_already;
      // replace the first line with the new line, and set the rest invalid
      *(cllist[stt]) = tmp;  cllist[stt]->merged = true;
      for (j = stt+1; j <= end; j++) cllist[j]->invalidate();
      for (k = cnt = 0; k < total; k++, cnt+=2) {
	pidx = (int)(clist[cnt+1] * w + clist[cnt+0]);
	lmap.setPointer( pidx, cllist[stt] );
      }
      count_im_old += nl;
      count_im_new ++;
    }
    free(tlist);  free(plist);  free(clist);
  }
  
  int  findLinesConnectedHeadToTail(LDLine *lp, LDLine *cllist[], int max_count) {
    LDLine *blist[40], *flist[40];
    int i, nb, nf, count=0;
    nb = traceLinesHeadToTail( lp, blist, false, 40 );
    nf = traceLinesHeadToTail( lp, flist, true , 40 );
    for (i=0; i<nb && count<max_count; i++) cllist[count++] = blist[nb-1-i];
    if (count < max_count)                  cllist[count++] = lp;
    for (i=0; i<nf && count<max_count; i++) cllist[count++] = flist[i];
    return count;
  }
  int  traceLinesHeadToTail(LDLine *lp, LDLine *llist[], bool forward, int max_count) {
    LDLine  *curr = lp, *next, *max_lp=NULL;
    int     i, cidx, pu, nb[20], count = 0, w2 = 2*w, tw2 = total-2*w;
    double  dot, max_dot;
    //if (curr->tail == 111715 && forward) printf("line(%d-%d) \n", curr->head, curr->tail);
    while (curr && count < max_count) {
      cidx = (forward ? curr->tail : curr->head);
      if (cidx < w2 || cidx >= tw2) break;
      pu = cidx % w;  if (pu < 2 || pu >= w-2) break;
      IMGH4V_SET( nb+0,  cidx-w2-1, cidx-w2, cidx-w2+1, cidx-w-2 );
      IMGH4V_SET( nb+4,  cidx-w-1, cidx-w, cidx-w+1, cidx-w+2 );
      IMGH4V_SET( nb+8,  cidx-2, cidx-1, cidx+1, cidx+2 );
      IMGH4V_SET( nb+12, cidx+w-2, cidx+w-1, cidx+w, cidx+w+1 );
      IMGH4V_SET( nb+16, cidx+w+2, cidx+w2-1, cidx+w2, cidx+w2+1 );
      for (i = 0, max_dot = 0; i < 20; i++) {
	next = ((LDLine**)(lmap.data))[ nb[i] ];
	if (!next || next == curr || 
	    nb[i] != (forward ? next->head : next->tail)) continue;
	dot = IMGH2V_DOT(curr->eq, next->eq);
	if (dot > max_dot) { max_dot = dot;  max_lp = next; }
      }
      //if (curr->tail == 111715 && forward) printf("line(%d-%d) and line(%d-%d) dot=%.2f \n", curr->head, curr->tail, max_lp->head, max_lp->tail, max_dot);
      if (max_dot < 0.9) break;
      llist[ count++ ] = curr = max_lp;
    }
    return count;
  }
  
  // -----------------------------------------------------------------
  // Parallel lines (lines that are very close with opposite direction)
  // -----------------------------------------------------------------
private:
  void findParallelLines(void) {
    LDLine *la, *lb;
    int     i, j;
    // count the number of parallel lines for each line
    for (i = 0; i < nlines; i++) {
      la = lines[i];  
      if (la->isInvalid()) continue;
      for (j = i+1; j < nlines; j++) {
	lb = lines[j];
	if (lb->isInvalid() || lb->tlen < 10) continue;
	if (!checkLineParallel( la, lb )) continue;
	la->nplines++;  lb->nplines++;
      }
    }
    // allocate memory
    for (i = 0; i < nlines; i++) {
      lines[i]->plines = (LDLine**)malloc(lines[i]->nplines * sizeof(LDLine*));
      lines[i]->nplines = 0;
    }
    // find parallel lines for each line
    for (i = 0; i < nlines; i++) {
      la = lines[i];
      if (la->isInvalid()) continue;
//       debug = (la->head == 68086 || la->head == 66805 || la->head == 66670);
//       if (debug) la->printInfo(" pline 0: ");
      for (j = i+1; j < nlines; j++) {
	lb = lines[j];
	if (lb->isInvalid() || lb->tlen < 10) continue;
	if (!checkLineParallel( la, lb )) continue;
// 	if (debug) lb->printInfo("       1: ");
	if (addLineAsParallel(la, lb, true) && addLineAsParallel(lb, la, true)) {
	  addLineAsParallel(la, lb, false);
	  addLineAsParallel(lb, la, false);
	}
// 	if (debug) lb->printInfo("       2: ", 'p');
      }
//       if (debug) la->printInfo("       9: ", 'p');
    }
  }
  
  bool checkLineParallel(LDLine *la, LDLine *lb) {
    double d[5];
    // check if two lines are parallel
    if (IMGH2V_DOT(la->eq, lb->eq) > -0.9 || fabs(la->eq[2]+lb->eq[2]) > 10) return false;
    la->getRelativePosition( lb, d, true );
    // ignore lines more than 5 pixels far away
    if (fabs(d[0]) > 5 || fabs(d[1]) > 5) return false;
    // ignore lines that intersect
    if (d[0] * d[1] < 0) return false;
    // ignore lines that do not overlap in the direction of the line
    if ((d[2] < 3 && d[3] < 3) || (d[2] > d[4]-3 && d[3] > d[4]-3)) return false;
    return true;
  }
  
  bool addLineAsParallel(LDLine *la, LDLine *lb, bool test_only) {
    // Try to add 'la' to the list of parallel lines of 'lb'
    // If 'test_only == true', it only tests without making any changes.
    LDLine *lc;  
    double  da[5], dc[5], ad[2], adist, cd[2], cdist;
    la->getRelativePosition( lb, da, true );
    ad[0] = (da[2] >= 0 ?     da[0] : (IMGH2V_DOT(la->eq, lb->eq)>0 ? -1 : +1) * IMGH2V_LPDIST(la->eq, lb->pp[0]));
    ad[1] = (da[3] <= da[4] ? da[1] : (IMGH2V_DOT(la->eq, lb->eq)>0 ? -1 : +1) * IMGH2V_LPDIST(la->eq, lb->pp[1]));
    adist = (fabs(ad[0]) > fabs(ad[1]) ? ad[0] : ad[1]);  // distance from 'lb' to 'la'
    // check 'la' against existing parallel lines of 'lb'
    for (int i = 0; i < lb->nplines; i++) {
      lc = lb->plines[i];
      lc->getRelativePosition( lb, dc, true );
      cd[0] = (dc[2] >= 0 ?     dc[0] : (IMGH2V_DOT(lc->eq, lb->eq)>0 ? -1 : +1) * IMGH2V_LPDIST(lc->eq, lb->pp[0]));
      cd[1] = (dc[3] <= dc[4] ? dc[1] : (IMGH2V_DOT(lc->eq, lb->eq)>0 ? -1 : +1) * IMGH2V_LPDIST(lc->eq, lb->pp[1]));
      cdist = (fabs(cd[0]) > fabs(cd[1]) ? cd[0] : cd[1]);  // distance from 'lb' to 'lc'
      if (adist * cdist < 0 ||  // lines on opposite side of 'lb'
	  da[3] < dc[2]+1 || da[2] > dc[3]-1) continue;  // does not overlap
      if (fabs(adist) > fabs(cdist)) return false;  // 'la' is worse than overlapped 'lc'
      if (!test_only) {
	// remove 'lc' from the list of parallel lines of 'lb'
	lc->removePLineTo( lb );
	lb->plines[i] = lb->plines[ lb->nplines - 1 ];  
	lb->nplines--;  i--;
      }
    }
    if (!test_only) lb->plines[ lb->nplines++ ] = la;
    return true;
  }
  
  void mergeUsingParallelism(void) {
    ALGH::Sorter<LDLine*> sorter;
    LDLine *la, *lb, *lc, *lp;
    double d[5], dir[2], tmp[2], lb_end, lc_end;
    int    i, j, merge_count;
    count_pm_old = count_pm_new = 0;
    for (i = nlines-1; i >= 0; i--) {
      la = lines[i];
//       debug = true; // (la->head == 4485);
      if (la->isInvalid() || la->nplines < 2) continue;
      // sort the list of parallel neighbors in the direction of 'la'
      IMGH2V_SUB( dir, la->pp[1], la->pp[0] );
      for (j=0; j < la->nplines; j++) {
	IMGH2V_SUB( tmp, la->plines[j]->pp[0], la->pp[0] );
	la->plines[j]->tag = (int)IMGH2V_DOT( dir, tmp );
      }
      sorter.QuickSort( la->plines, la->nplines, compare_lines_by_tag );
//       if (debug) la->printInfo("Merge la: ", 'p');
      // merge parallel neighbor lines
      for (j=merge_count=0, lb = lc = NULL; j < la->nplines; j++) {
	lp = la->plines[j];
	if (lp->isInvalid()) continue;
	lp->getRelativePosition( la, d, true );
	if (d[3] < 5 || d[2] > d[4]-5) continue;
	if (d[0] > 0) {  // merge all the parallel lines above 'la'
	  if (lb == NULL) { lb = la->plines[j];  lb_end = d[3]; }
	  else if (d[2] < lb_end -  1) ;
	  else if (d[2] - lb_end > h/10) lb = NULL;
	  else {
// 	    if (debug) { lb->printInfo("  prev:",'p'); lp->printInfo("  next:",'p'); }
	    mergeTwoLines( lb, la->plines[j], true );  merge_count++;
	    lb_end = d[3];
	  }
	} else {          // merge all the parallel lines below 'la'
	  if (lc == NULL) { lc = la->plines[j];  lc_end = d[3]; }
	  else if (d[2] < lc_end -  1) ;
	  else if (d[2] - lc_end > h/10) lc = NULL;
	  else {
// 	    if (debug) { lc->printInfo("  prev:",'p'); lp->printInfo("  next:",'p'); }
	    mergeTwoLines( lc, la->plines[j], true );  merge_count++;
	    lc_end = d[3];
	  }
	}
      }
      if (merge_count > 0) { count_pm_old += merge_count+1; count_pm_new++; }
//       if (debug) la->printInfo("      la: ", 'p');
      // remove the pointers to invalidated parallel lines from the list
      for (j=0; j < la->nplines; j++) {
	if (la->plines[j]->isInvalid()) { 
	  la->plines[j] = la->plines[ la->nplines-1 ];
	  la->nplines--;  j--;
	}
      }
//       if (debug) la->printInfo("      la: ", 'p');
    }
  }
  
  // -----------------------------------------------------------------
  // Line links (connections between two lines)
  // -----------------------------------------------------------------
private:
  void findLineLinks(void) {
    // Look for connections between lines
    // Note that this function should be called only after 'removeInvalidatedLines()'.
    int    i;
    LDLine **lmapp = (LDLine**)(lmap.data);
    // count the number of links for each line
    for (i = 0; i < nlines; i++) setupLineLinks( lines[i], lmapp, true );
    // allocate memory for link information
    for (i = 0; i < nlines; i++) {
      lines[i]->links = (LDLink*)malloc( lines[i]->nlinks * sizeof(LDLink) );
      lines[i]->nlinks = 0;
    }
    // save all the links for each line
    for (i = 0; i < nlines; i++) setupLineLinks( lines[i], lmapp, false );
    
    // search for cross of two lines using parallel lines
    int    j, k;
    LDLine *la, *lb, *lc;
    double d[2], x, y;
    for (i = 0; i < nlines; i++) {
      la = lines[i];
      for (j = 0; j < la->nplines; j++) {
	lb = la->plines[j];
	if (lb->tlen < la->tlen/2) continue;
	for (k = 0; k < lb->nlinks; k++) {
	  lc = lb->links[k].lp;
	  if (la->findLinkTo( lc )) continue;
	  lc->getRelativePosition( la, d, false );
	  x = lb->links[k].pidx % w;
	  y = lb->links[k].pidx / w;
	  if (d[0] * d[1] > 0 || la->getDistanceToPoint(x,y,true) >= 5) continue;
	  // add a new link to 'lc' into the link list of 'la'
	  la->links = (LDLink*)realloc( la->links, (la->nlinks+1)*sizeof(LDLink) );
	  la->getPointOnLine( x, y, d );
	  la->links[la->nlinks].pidx = (int)(ROUND(d[1])) * w + (int)ROUND(d[0]);
	  la->links[la->nlinks].lp   = lc;
	  la->nlinks++;
	  // add a new link to 'la' into the link list of 'lc'
	  if (lc->findLinkTo( la )) continue;
	  lc->links = (LDLink*)realloc( lc->links, (lc->nlinks+1)*sizeof(LDLink) );
	  lc->getPointOnLine( x, y, d );
	  lc->links[lc->nlinks].pidx = (int)(ROUND(d[1])) * w + (int)ROUND(d[0]);
	  lc->links[lc->nlinks].lp   = la;
	  lc->nlinks++;
	}
      }
    }
  }
  
  void setupLineLinks(LDLine *la, LDLine **lmapp, bool count_only) {
    LDLine  *lb;  LDLink lks[40];
    int     i, j, k, ht[2], nb[20], nlinks=0, w2 = 2*w, cpidx, lks_pidx[40];
    if (la->isInvalid()) return;
    ht[0] = la->head;  ht[1] = la->tail;
    for (k = 0; k < 2; k++) {
      // search the 1st ring of the head/tail of the line
      IMGH4V_SET( nb+0,  ht[k]-w  , ht[k]+w  , ht[k]  -1, ht[k]  +1 );
      IMGH4V_SET( nb+4,  ht[k]-w-1, ht[k]-w+1, ht[k]+w-1, ht[k]+w+1 );
      for (i = 0; i < 8; i++) {
	lb = lmapp[ nb[i] ];
	if (!lb || lb == la || lb->isInvalid() || 
	    la->isParallelWith(lb)) continue;
	for (j = 0; j < nlinks; j++) if (lks[j].lp == lb) break;
	if (j < nlinks) continue;
	if ((nb[i] == lb->head || nb[i] == lb->tail) && la < lb) continue;
	lks_pidx[nlinks] = ht[k];
	lks[nlinks++].setup( nb[i], lb );	// save the link
      }
      // search the 2nd ring of the head/tail of the line
      IMGH4V_SET( nb+0,  ht[k]-w2, ht[k]+w2, ht[k]-2, ht[k]+2 );
      IMGH4V_SET( nb+4,  ht[k]-w2-1, ht[k]-w2+1, ht[k]+w2-1, ht[k]+w2+1 );
      IMGH4V_SET( nb+8,  ht[k]-w-2,  ht[k]+w-2,  ht[k]-w+2,  ht[k]+w+2 );
      IMGH4V_SET( nb+12, ht[k]-w2-2, ht[k]-w2+2, ht[k]+w2-2, ht[k]+w2+2 );
      for (i = 0; i < 16; i++) {
	lb = lmapp[ nb[i] ];
	if (!lb || lb == la || lb->isInvalid() || 
	    fabs(IMGH2V_DOT(la->eq, lb->eq)) > 0.86) continue;
	for (j = 0; j < nlinks; j++) if (lks[j].lp == lb) break;
	if (j < nlinks) continue;
	if ((nb[i] == lb->head || nb[i] == lb->tail) && la < lb) continue;
	// check the connectivity through the pixels between them
	cpidx = ht[k] + (nb[i]-ht[k])/2;
	if (lmapp[ cpidx ]) continue;
// 	gz[0] = ldt->edt->gz->getFloat( ht[k] );
// 	gz[1] = ldt->edt->gz->getFloat( nb[i] );
// 	gz[2] = ldt->edt->gz->getFloat( cpidx );
// 	if (gz[2] < (gz[0]+gz[1])/4) continue;
	lks_pidx[nlinks] = ht[k];
        lks[nlinks++].setup( nb[i], lb );	// save the link
      }
    }
    if (count_only) {
      for (i=0; i<nlinks; i++) { la->nlinks++; lks[i].lp->nlinks++; }
    } else {
      for (i=0; i<nlinks; i++) { 
	lb = lks[i].lp;
	LDLink *lk = la->findLinkTo(lb);
	if (lk == NULL) {
	  la->links[ la->nlinks++ ].setup( lks_pidx[i], lks[i].lp );
	  lb->links[ lb->nlinks++ ].setup( lks[i].pidx, la );
	} else continue; // connected between endpoints (visited twice)
      }
    }
  }
  
  // -----------------------------------------------------------------
  // functions to be used later on
  // -----------------------------------------------------------------
public:
  int findLinesConnectingTwo(LDLine *la, LDLine *lb, LDLine **clines, int max_count=20) {
    // Find lines that connect 'la' and 'lb'
    int k, count=0;
    for (k = 0; k < la->nlinks && count < max_count; k++) {
      LDLine *cl = la->links[k].lp;
      if (lb->findLinkTo( cl )) clines[count++] = cl;
    }
    return count;
  }

  bool mergeTwoLines(LDLine *la, LDLine *lb, bool connect) {
    // Merge two lines 'la' and 'lb' into 'la'
    //   'lb' will lose all the links, and its 'lp' will be set invalid.
    if (!la || !lb || la == lb) return false;
    int i, k, n, *plist;
    int na, *alist = (int*) malloc( (w+h) * sizeof(int) );
    int nb, *blist = (int*) malloc( (w+h) * sizeof(int) );
    int nc, *clist = (int*) malloc( (w+h) * sizeof(int) );
    // find all the pixels that belongs to two lines
    na = getLinePoints( la->head, alist, (w+h) );
    nb = getLinePoints( lb->head, blist, (w+h) );
    if (na < 2 || nb < 2) return false;
    int     a0 = alist[0], a1 = alist[na-1], b0 = blist[0], b1 = blist[nb-1];
    double  dx, dy, dist0, dist1;
    dx = (a1 % w - b0 % w);  dy = (a1 / w - b0 / w);  dist0 = sqrt(dx*dx + dy*dy);
    dx = (a0 % w - b1 % w);  dy = (a0 / w - b1 / w);  dist1 = sqrt(dx*dx + dy*dy);
    if (dist0 <= dist1) {	// 'la' followed by 'lb'
      if (dist0 > h/5) {
	fprintf(stderr, "Warning (IMGH::LineMapper::mergeTwoLines) %.1f pixels away (%d and %d pixels) \n", dist0, na, nb);
	la->printInfo("  merge a:");  lb->printInfo("  merge b:");
	printf("  a: %d %d %d ... %d %d %d\n", alist[0], alist[1], alist[2], alist[na-3], alist[na-2], alist[na-1]);
	printf("  b: %d %d %d ... %d %d %d\n", blist[0], blist[1], blist[2], blist[nb-3], blist[nb-2], blist[nb-1]);
      }
      nc = (connect ? findPixelsBetween( a1, b0, clist ) : 0);
      memcpy( alist + na, clist, nc * sizeof(int) );
      memcpy( alist + na + nc, blist, nb * sizeof(int) );
      plist = alist;  free(blist);  free(clist);
      n = na + nc + nb;
    } else {			// 'lb' followed by 'la'
      if (dist1 > h/5) {
	fprintf(stderr, "Warning (IMGH::LineMapper::mergeTwoLines) %.1f pixels away (%d and %d pixels) \n", dist1, na, nb);
	la->printInfo("  merge a:");  lb->printInfo("  merge b:");
      }
      nc = (connect ? findPixelsBetween( b1, a0, clist ) : 0);
      memcpy( blist + nb, clist, nc * sizeof(int) );
      memcpy( blist + nb + nc, alist, na * sizeof(int) );
      plist = blist;  free(alist);  free(clist);
      n = nb + nc + na;
    }
    // find the optimum line equation of 'la' for the list of pixels
    double *dlist = (double*)malloc( n * 2 * sizeof(double) );
    for (i = 0; i < n; i++) { dlist[2*i+0] = plist[i]%w; dlist[2*i+1] = plist[i]/w; }
    fitLine( la, plist, dlist, n );
    free(dlist);
    // update the line map 'lmap'
    for (i = 0; i < n; i++) lmap.setPointer( plist[i], la );
    // create new list of the links from 'la' and 'lb'
    LDLink *linklist = (LDLink*) malloc( (la->nlinks + lb->nlinks) * sizeof(LDLink) );
    int linkcount = 0;
    for (i = 0; i < la->nlinks; i++)	// copy links of 'la'
      if (la->links[i].lp != lb) linklist[linkcount++] = la->links[i];
    for (i = 0; i < lb->nlinks; i++)	// copy links of 'lb'
      if (lb->links[i].lp != la) linklist[linkcount++] = lb->links[i];
    for (i = 0; i < linkcount; i++) linklist[i].lp->replaceLinkTo( lb, la );
    // create new list of parallel lines from 'la' and 'lb'
    LDLine **pllist = (LDLine**) malloc((la->nplines+lb->nplines) * sizeof(LDLine*));
    int nplines = 0;
    for (i = 0; i < la->nplines; i++) { // copy plines of 'la'
      for (k=0; k<nplines; k++) if (la->plines[i] == pllist[k]) break;
      if (k == nplines) pllist[nplines++] = la->plines[i];
    }
    for (i = 0; i < lb->nplines; i++) { // copy plines of 'lb'
      for (k=0; k<nplines; k++) if (lb->plines[i] == pllist[k]) break;
      if (k == nplines) {
	lb->plines[i]->replacePLineTo( lb, la );
	pllist[nplines++] = lb->plines[i];
      }
    }
    // update 'la' with the new links, plines, head and tail
    la->removeLinks(false);  la->removePLines();
    la->links = linklist;    la->nlinks = linkcount;  
    la->plines = pllist;     la->nplines = nplines;
    la->head = plist[0]; la->tail = plist[n-1];
    la->merged = true;
    // clear 'lb' 
    lb->removeLinks(false);  lb->removePLines();  lb->invalidate();
    free(plist);
    return true;
  }
  
  int  getLinePoints(int idx, int array[], int maxlen=0) {
    // Assuming all the pixels are mapped to lines in 'lmap' and
    //   'idx' is the index of the first pixel of the line,
    //   return the list of pixels for the line.
    if (maxlen == 0) maxlen = w+h;
    int   i, cidx=0, pu, pv, nidx=0, nu, nv, count=0;
    IMGH::LDLine  **lmapp = (LDLine**)(lmap.data);
    IMGH::LDLine   *lp = lmapp[idx];
    if (lp == NULL) return 0;
    array[ count++ ] = cidx = idx;
    lmapp[ idx ] = NULL;        // change the pixel value
    int du[24] = {  0,  0, -1, +1,  -1, +1, -1, +1,   0,  0, -2, +2,  -1, +1, -1, +1,  -2, -2, +2, +2,  -2, +2, -2, +2 };
    int dv[24] = { -1, +1,  0,  0,  -1, -1, +1, +1,  -2, +2,  0,  0,  -2, -2, +2, +2,  -1, +1, -1, +1,  -2, -2, +2, +2 };
    while (count < maxlen) {	// search 1st ring of the current pixel
      pu = cidx % w;  pv = cidx / w;
      for (i = 0; i < 24; i++) {  // search neighborhood (5x5-1) for connected edge pixel
	nu = pu + du[i];  nv = pv + dv[i];
	if (nu < 0 || nu >= w || nv < 0 || nv >= h) continue;
	nidx = nv * w + nu;
	if (lmapp[nidx] == lp) break;
      }
      if (i==24) break;
      array[ count++ ] = cidx = nidx;              // save the next pixel
      lmapp[ nidx ] = NULL;                         // change the pixel value
    }
    for (i=0; i<count; i++) lmapp[array[i]] = lp;  // restore the pixel values
    return count;
  }
  
  int findPixelsBetween(int pa, int pb, int plist[], bool overwrite=true) {
    int   a[2] = { pa%w, pa/w }, b[2] = { pb%w, pb/w };
    int   dx = b[0] - a[0],  dy = b[1] - a[1];
    int   i, count, pidx, n = (abs(dx) > abs(dy) ? abs(dx) : abs(dy));
    float x, y, sx = dx/(float)n,  sy = dy/(float)n;
    for (i = 1, count=0, x = a[0]+sx, y = a[1]+sy; i < n; i++, x+=sx, y+=sy) {
      pidx = ( int(ROUND(y)) * w + int(ROUND(x)) );
      LDLine *lp = (LDLine*)lmap.getPointer(pidx);
      if (overwrite) {
	if (!lp || pidx != lp->head) plist[count++] = pidx;
      } else {
	if (!lp) plist[count++] = pidx;
      }
    }
    return count;
  }
  
  // -----------------------------------------------------------------
  // -----------------------------------------------------------------
public:
  void removeInvalidatedLines(void) {
    IMGH::Image tmp( w, h, PIXEL_GRAY, true );
    LDLine *lp;
    // remove invalid lines from the 'lines[]'
    int    i, j, n, *plist=(int*)malloc((w+h)*sizeof(int));
    for (i = 0; i < nlines; i++) {
      lp = lines[i];
      if (lp->isInvalid()) {	// remove invalid lines
	lines[i] = lines[nlines-1];  nlines--;  i--;
	count_removed++;
      } else {			// mark all the pixels on the line
	n = getLinePoints( lp->head, plist, (w+h) );
	for (j = 0; j < n; j++) tmp.setChar( plist[j], 's' );
      }
    }
    // remove invalid pixels from 'lmap'
    char    *src = (char*)(tmp.data);
    LDLine **dst = (LDLine**)(lmap.data);
    for (i = 0; i < total; i++) 
      if (dst[i] != NULL && src[i] == 0) dst[i] = NULL;
    free( plist );
  }
  
  // -----------------------------------------------------------------
  // 
  // -----------------------------------------------------------------
public:
  void printInfo(char *cmmt=NULL) {
    if (edt == NULL) return;
    if (edt == &edt_dummy) edt->printInfo();
    printf("%s\n", (cmmt ? cmmt : "IMGH::LineDetector"));
    printf("  lines  : %3d => %3d  \n", count_init, nlines);
    lmap.printInfo("  lmap");
    printf("  i merge: %3d lines replaced with %3d lines  (%4d)\n", count_im_old, count_im_new, count_im_new - count_im_old );
    printf("  p merge: %3d lines replaced with %3d lines  (%4d)\n", count_pm_old, count_pm_new, count_pm_new - count_pm_old);
    LDLine *la, *lb;
    int    lidx, i, j, pair_count=0, link_count=0;
    for (lidx = 0; lidx < nlines; lidx++) {
      la = lines[lidx];
      if (la->head <= 0 || la->head > total) { la->printInfo("  Error (IMGH::LineMapper::printInfo;HEAD)", true); continue; }
      if (la->nplines > 0) {			// check parallel lines
	for (i=0; i < la->nplines; i++) {
	  lb = la->plines[i];
	  for (j=0; j < lb->nplines; j++) if (lb->plines[j] == la) break;
	  if (j == lb->nplines) {
	    printf("  Error (IMGH::LineMapper::printInfo): invalid parallel line\n");
	    la->printInfo("    LDLine", 'p');
	  }
	}
	pair_count += la->nplines;
      }
      if (la->nlinks < 0 || la->nlinks > 100) { la->printInfo("  Error (IMGH::LineMapper::printInfo;NLINKS)", true); continue; }
      for (i = 0; i < la->nlinks; i++) {	// check link information
	if (la->links[i].pidx <= 0 || la->links[i].pidx >= total) {
	  printf("  Error (IMGH::LineMapper::printInfo): invalid link at pidx=%d\n", la->links[i].pidx);
	  la->printInfo("  LDLine", 'k');
	}
	lb = la->links[i].lp;
	if (lb == NULL || lb->nlinks < 0 || lb->nlinks > 100) { la->printInfo("  Error (IMGH::LineMapper::printInfo;LINK)", true); continue; }
	if (lb->findLinkTo( la ) == NULL)
	  printf("  Error (IMGH::LineMapper::printInfo): invalid line link (%d -> %d, not %d <- %d)\n", la->head, lb->head, la->head, lb->head);
      }
      if (la->nlinks > 0) link_count += la->nlinks;
    }
    printf("  final  : %d lines (%d pairs, %d links)  (%d removed from %d)\n",
	   nlines, pair_count/2, link_count/2, count_removed, count_init);
  }
  
  void printLine(int pidx) {
    for (int i=0; i<nlines; i++)
      if (lines[i]->head == pidx) { lines[i]->printInfo(); break; }
  }
  
  void showLines(Image *vimg, bool clear=true, bool random_color=true) {
    if (!vimg || nlines<=0) return;
    vimg->setImage( w, h, PIXEL_RGB, clear);
    showLines( vimg->data, PIXEL_RGB, clear, random_color );
  }
  void showLines(void *image_buffer, pixel_t type, bool clear=true, bool random_color=true) {
    if (!edt || !image_buffer || nlines<=0) return;
    ImageDrawer<float> ip( w, h, type, image_buffer );
    if (clear) ip.clearImage();
    for (int i=0; i<nlines; i++) {
      LDLine *lp = lines[i];
      if (lp->pp[0][1] < 0 || lp->pp[1][1] < 0) { lp->printInfo("Invalid Line"); continue; }
      if (random_color) {
	srand( lp->head );
	ip.drawLine( lp->pp[0], lp->pp[1], 50+rand()%206, 50+rand()%206, 50+rand()%206 );
      } else {
	float prob = lp->extra[1] - 0.5f;  // [ -0.5 ~ +0.5 ]
	float  wh[3]={200,200,200}, rd[3]={255,0,0}, bl[3]={50,50,255}, color[3];
	if (prob >= 0) IMGH3V_WADD( color, +prob*2, rd, (1-prob*2), wh );
	else           IMGH3V_WADD( color, -prob*2, bl, (1+prob*2), wh );
	ip.drawLine( lp->pp[0], lp->pp[1], (int)color[0], (int)color[1], (int)color[2] );
      }
    }
  }
  
  void showInfo(void *image_buffer, pixel_t type, char data_type='l') {
    int  lidx, i, pa[2], pb[2];
    ImageDrawer<float> ip( w, h, type, image_buffer );
    LDLine *lp;  LDLine **src = (LDLine**)(lmap.data);
    switch (data_type) {
    case 'l':	// visualize lines with random colors
      for (lidx = 0; lidx < nlines; lidx++) {
	lp = lines[lidx];
	srand( lp->head );
	ip.drawLine( lp->pp[0], lp->pp[1], 50+rand()%206, 50+rand()%206, 50+rand()%206 );
      }
      break;
    case 'm': case 'M':	// visualize line map with random colors
      for (i = 0; i < total; i++) {
	if (src[i]) srand( src[i]->head ); else continue;
	if (data_type == 'M' && !src[i]->merged) continue;
	ip.setPixelWithRandomColor( i );
      }
      break;
    case 'k':	// visualize line links
      for (i = 0; i < total; i++)
	if (src[i]) ip.setPixel( i, 30, 30, 30 );
      for (lidx = 0; lidx < nlines; lidx++) {
	lp = lines[lidx];
	for (i = 0; i < lp->nlinks; i++) {
	  LDLine *lb = lp->links[i].lp;
	  LDLink *lka = lp->links + i;
	  LDLink *lkb = lb->findLinkTo( lp );  if (!lkb) continue;
	  if (lka->pidx == 712 || lkb->pidx == 712) printLinkInfo(lka);
	  srand( (lka->pidx < lkb->pidx ? lka->pidx : lkb->pidx) );
	  IMGH2V_SET( pa, lka->pidx%w, lka->pidx/w );
	  IMGH2V_SET( pb, lkb->pidx%w, lkb->pidx/w );
	  //if (IMGH2V_DIST(pa, pb) > 5) printf("LINK between %d and %d : %d - %d : dist=%g\n",
	  //				     lp->head, lb->head, lka->pidx, lkb->pidx, IMGH2V_DIST(pa, pb));
	  ip.drawLine( (float)pa[0], (float)pa[1], (float)pb[0], (float)pb[1], 
		       50+rand()%206, 50+rand()%206, 50+rand()%206 );
	}
      }
      break;
    case 'p':	// visualize parallel lines with random colors
      for (lidx = 0; lidx < nlines; lidx++) {
	lp = lines[lidx];
	if (lp->nplines < 1) continue;
	srand( lp->head );
	ip.drawLine( lp->pp[0][0], lp->pp[0][1], 
		     lp->pp[1][0], lp->pp[1][1], 
		     50+rand()%206, 50+rand()%206, 50+rand()%206 );
      }
      break;
    }
  }
  
  LDLine* findNearestLine(double x, double y) {
    LDLine *lp, *lp_best = NULL;
    double   dist, min_dist = 9999;
    for (int lidx = 0; lidx < nlines; lidx++) {
      lp = lines[lidx];
      if (lp->isInvalid()) continue;
      dist = fabs(lp->getDistanceToPoint( x, y, true ));
      if (dist < min_dist) { min_dist = dist; lp_best = lp; }
    }
    return lp_best;
  }
  void showLineInfo(void *image_buffer, pixel_t type, LDLine *lp) {
    ImageDrawer<double> ip( w, h, type, image_buffer );
    LDLine *llp2;
    int    i, k, np, plist[800], r, g, b;
    // draw the line
    np = getLinePoints( lp->head, plist, 800 );
    for (k = 0; k < np; k++) ip.setPixel( plist[k], 255, 255, 255 );
//     printf("Points of %d: ", lp->head);
//     for (k = 0; k < np && k < 10; k++) printf("%d ", plist[k]);
//     printf("\n");
    // draw its links
    for (i = 0; i < lp->nlinks; i++) {
      llp2 = lp->links[i].lp;
      np = getLinePoints( llp2->head, plist, 800 );
      srand( llp2->head );  r = 50+rand()%206; g = 50+rand()%206; b = 50+rand()%206;
      for (k = 0; k < np; k++) ip.setPixel( plist[k], 255, 0, 0 );
    }
    // draw its parallel lines
    for (i = 0; i < lp->nplines; i++) {
      llp2 = lp->plines[i];
      np = getLinePoints( llp2->head, plist, 800 );
      srand( llp2->head );  r = 50+rand()%206; g = 50+rand()%206; b = 50+rand()%206;
      for (k = 0; k < np; k++) ip.setPixel( plist[k], 0, 255, 0 );
    }
  }
  void printLinkInfo(LDLink *lk, char *cmmt=NULL) {
    int i, pa=lk->pidx, pb;
    LDLine *la, *lb = lk->lp;
    if (lb == NULL) return;
    for (i=0; i < lb->nlinks; i++)
      if (lb->links[i].lp && lb->links[i].lp->findLinkTo( lb )) break;
    if (i < lb->nlinks) {  la = lb->links[i].lp;  pb = lb->links[i].pidx; }
    else { la = NULL;  pb = 0; }
    printf("%s : from pixel %d on line %d  to  pixel %d on line %d\n", 
	   (cmmt ? cmmt : "Link"), 
	   pa, (la ? la->head : 0), pb, (lb ? lb->head : 0));
  }
};

  
}	// namespace IMGH


#endif	// IMGH_LINE_DETECTOR_HPP


