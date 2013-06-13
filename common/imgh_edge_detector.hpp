
//
// IMGH::EdgeDetector 
//
// Jaeil Choi
// last modified in Nov, 2006
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
//   - extracting edges from a given image
// This class requires:
//   - IMGH::Image		in 'imgh_common.hpp'
//   - IMGH::ImageGradient	in 'imgh_gradient.hpp'
//

#ifndef IMGH_EDGE_DETECTOR_HPP
#define IMGH_EDGE_DETECTOR_HPP

#include <iostream>
#include <cmath>
#include "imgh_common.hpp"
#include "imgh_gradient.hpp"

namespace IMGH {
  
class EdgeDetector
{
  // Every member variable is ReadOnly.
public:
  int        w, h, total;	// frequently used values
  IMGH::Image emap;		// edge map (image pixel type IMGH::PIXEL_INT)
  float      low_val, high_val;	// low/high gradient strength
  int        edge_count;
  int        max_edge_length;
  IMGH::ImageGradient	*grd;
  IMGH::Image *gx, *gy, *gz;	// frequently used pointers
  float     *vedge_count;
  
private:
  IMGH::ImageGradient	grd_dummy;
  
public:
  bool       edge_following;
  
public:
  EdgeDetector() : grd(NULL), gx(NULL), gy(NULL), gz(NULL), vedge_count(NULL), edge_following(true) {}
  EdgeDetector(IMGH::Image *img, float low_thr=0.05, float high_thr=0.1, int filter_size=0) 
    : grd(NULL), gx(NULL), gy(NULL), gz(NULL), vedge_count(NULL), edge_following(true) {
    findEdges( img, low_thr, high_thr, filter_size );
  }
  ~EdgeDetector() { clear(); }
  void clear(void) { 
    emap.clear(); grd_dummy.clear(); 
    if (vedge_count) { free(vedge_count); vedge_count = NULL; }
  }
  
  // -----------------------------------------------------------------
  // 
  // -----------------------------------------------------------------
  
public:
  void findEdges(IMGH::Image *img, 
		 float low_thr=0.05, float high_thr=0.1, int filter_size=0 ) {
    // Find edges with low/high thresholds w.r.t. the maximum size of gradient.
    if (filter_size<0 || filter_size>500) {  printf("Error: filter size (%d) is too big\n", filter_size); return; }
    grd_dummy.getGradientImage( img, filter_size, true );
    findEdges( &grd_dummy, low_thr, high_thr );
  }
  void findEdges(int w, int h, pixel_t type, void *udata, 
		 float low_thr=0.05, float high_thr=0.1, int filter_size=0 ) {
    // Find edges with low/high thresholds w.r.t. the maximum size of gradient.
    if (filter_size<0 || filter_size>500) {  printf("Error: filter size (%d) is too big\n", filter_size); return; }
    grd_dummy.getGradientImage( w, h, type, udata, filter_size, true );
    findEdges( &grd_dummy, low_thr, high_thr );
  }
  void findEdges(IMGH::ImageGradient *grd, float low_thr=0.05, float high_thr=0.1) {
    // check required data (gx, gy and gz)
    if (!grd || grd->gx.w == 0 || grd->gx.h == 0 || !grd->gx.sameSize( &grd->gy )) return;
    if (!grd->gz.sameSize( &grd->gx )) grd->getGradientStrength();
    // save frequently used values
    this->grd = grd;
    this->gx = &(grd->gx);  this->gy = &(grd->gy);  this->gz = &(grd->gz);
    this->w  = gx->w;  this->h = gx->h;  this->total = w * h;
    // set up the edge map image
    emap.setImage( w, h, PIXEL_INT, true );
    edge_count = max_edge_length = 0;
    // decide threshold values
    if (low_thr <= 0.01) low_thr = 0.01f;   if (low_thr >= 0.99) low_thr = 0.99f;
    if (high_thr<= 0.01) high_thr = 0.01f;  if (high_thr>= 0.99) high_thr = 0.99f;
    low_val  = grd->max_strength * low_thr;
    high_val = grd->max_strength * high_thr;
    //printf("E  low=%.4f high=%.4f (%.4f %.4f max=%.4f) edge_following=%d\n", low_val, high_val, low_thr, high_thr, grd->max_strength, edge_following);
    float strength, gv[2];
    int   x, y, i, idx, na, nb, len, *epp=(int*)(emap.data), *hist=NULL, w1=w-1, h1=h-1;
    if (edge_following) hist = (int*) malloc( 3 * w * sizeof(int) );
    for (y = idx = 0; y < h; y++) {
      for (x = 0; x < w; x++, idx++) {
	//bool debug = (idx==498 || idx==498+640-1);
	strength = gz->getFloat(idx);
	//if (debug) printf("E %d  str=%.4f low=%.4f high=%.4f  epp[]=%d\n", idx, strength, low_val, high_val, epp[idx]);
	if (strength < high_val) continue;	// skip pixels with weak strength
	if (epp[idx] != 0) continue;		// skip pixels that were visited
	grd->getGradientVector( idx, gv, true );
	// non-maximum suppression (search for peaks in edge strength)
	if        (fabs(gv[0]) > fabs(gv[1]) * 1.414213) {	// vertical edge   (|)
	  if (x>0  && strength < gz->getFloat( idx-1 )) continue;
	  if (x<w1 && strength < gz->getFloat( idx+1 )) continue;
	} else if (fabs(gv[1]) > fabs(gv[0]) * 1.414213) {	// horizontal edge (-)
	  if (y>0  && strength < gz->getFloat( idx-w )) continue;
	  if (y<h1 && strength < gz->getFloat( idx+w )) continue;
	} else if (gv[0] * gv[1] < 0) {				// diagonal edge   (/)
	  if (x>0  && y>0  && strength < gz->getFloat( idx-w-1 )) continue;
	  if (x<w1 && y<h1 && strength < gz->getFloat( idx+w+1 )) continue;
	} else {						// diagonal edge   (\)
	  if (x<w1 && y>0  && strength < gz->getFloat( idx-w+1 )) continue;
	  if (x>0  && y<h1 && strength < gz->getFloat( idx+w-1 )) continue;
	}
	if (edge_following) {
	  // start searching for neighbor edge points from this pixel (backward and forward),
	  //   and set all the points on the edge to the index of the first pixel.
	  na = findEdgePointsByLoop( idx, idx, gv, false, hist );  // backward
	  // if (debug) { printf("Edge %d : na=%d (", idx, na); for (i=0;i<na;i++) printf("%d ", hist[i]); printf(")\n"); }
	  if (na>1)  { hist[0] = hist[na-1]; for (i=1; i<na; i++) epp[ hist[i] ] = hist[0]; }
	  nb = findEdgePointsByLoop( (na>1 ? hist[0]:idx), idx, gv, true, hist+na-1 ); // forward
	  // if (debug) { printf("Edge %d : nb=%d (", idx, nb); for (i=0;i<nb;i++) printf("%d ", hist[na-1+i]); printf(")\n"); }
	  len = 1 + (na-1) + (nb-1);
	  if (len > max_edge_length) max_edge_length = len;
	  // if (debug) { printf("Edge %d : len=%d (", idx, len); for (i=0;i<len;i++) printf("%d ", hist[i]); printf(")\n"); }
	  if (len < 2) for (i=0; i<len; i++) epp[ hist[i] ] = 0;
	} else {
	  epp[ idx ] = 1;
	}
	edge_count++;
      }
    }
    if (hist) free(hist);
    int plist[1024];
    getEdgePoints( 943, plist, 1024 );
  }
  
  float countVerticalEdgePixels(bool as_ratio=true) {
    // Count the number of vertical edge pixels for each column in 'vedge_count[]'.
    if (vedge_count) { free(vedge_count); vedge_count = NULL; }
    if (emap.w <= 0 || emap.h <= 0) return 0;
    vedge_count = (float*)calloc( emap.w, sizeof(float) );
    int  i, total=emap.w*emap.h, w=emap.w;  float cmax=0;
    for (i = w; i < total-w; i++)
      if (emap.getInt(i)>0 && emap.getInt(i-w)>0 && emap.getInt(i+w)>0) vedge_count[i%w]+=1;
    for (i = 0; i < w; i++) {
      if (as_ratio) vedge_count[i] /= emap.h;
      if (vedge_count[i]>cmax) cmax = vedge_count[i];
    }
    return cmax;
  }
  
private:
  int findEdgePointsByLoop(int oidx, int cidx, float gv[2], bool forward, int hist[]) {
    float strength, stmax, cgv[2], ngv[3][2], ogv[2], once_straight_gv[2]={0,0};
    int   i, idx, prev=-1, nb[3], midx=0, len=0;
    int   straight_count=0, x, y, *epp=(int*)(emap.data);
    bool  once_straight = false;
    IMGH2V_COPY( cgv, gv );  IMGH2V_COPY( ogv, gv );
    epp[ cidx ] = oidx;				// set this pixel to 'oidx'
    if (hist) hist[len++] = cidx;
    // bool debug = (cidx == 91646 || cidx == 99316);
    while (true) {
      if        (fabs(cgv[0]) > fabs(cgv[1]) * 1.414213) {	// vertical edge   (|)
	// if (debug) printf("%d| ", cidx);
	if (forward) IMGH3V_SET( nb, cidx+w-1, cidx+w, cidx+w+1 );
	else         IMGH3V_SET( nb, cidx-w-1, cidx-w, cidx-w+1 );
	if (prev==nb[0]||prev==nb[1]||prev==nb[2]) {  
	  // resolve ambiguity because forward/backward is decided by X axis
	  if (forward) IMGH3V_SET( nb, cidx-w-1, cidx-w, cidx-w+1 );
	  else         IMGH3V_SET( nb, cidx+w-1, cidx+w, cidx+w+1 );
	}
      } else if (fabs(cgv[1]) > fabs(cgv[0]) * 1.414213) {	// horizontal edge (-)
	// if (debug) printf("%d- ", cidx);
	if (forward) IMGH3V_SET( nb, cidx-w+1, cidx+1, cidx+w+1 );
	else         IMGH3V_SET( nb, cidx-w-1, cidx-1, cidx+w-1 );
      } else if (cgv[0] * cgv[1] > 0) {				// diagonal edge   (/)
	// if (debug) printf("%d/ ", cidx);
	if (forward) IMGH3V_SET( nb, cidx-w, cidx-w+1, cidx+1 );
	else         IMGH3V_SET( nb, cidx-1, cidx+w-1, cidx+w );
      } else {							// diagonal edge   (\)
	// if (debug) printf("%d\\ ", cidx);
	if (forward) IMGH3V_SET( nb, cidx+1, cidx+w, cidx+w+1 );
	else         IMGH3V_SET( nb, cidx-1, cidx-w-1, cidx-w );
      }
      for (stmax = 0, i = 0; i < 3; i++) {
	idx = nb[i];  x = idx % w;  y = idx / w;
	if (x<1 || x>=w-1 || y<1 || y>=h-1) continue;	// terminate at boundary
	strength = gz->getFloat( idx );
	if (strength < low_val) continue;		// skip pixels with weak strength
	grd->getGradientVector( idx, ngv[i], true );
	if (IMGH2V_DOT( cgv, ngv[i] ) < 0.7) continue;	// skip pixels with different gradient vector
	if (strength > stmax) { stmax = strength; midx = i; }
      }
      if (stmax <= 0 || epp[ nb[midx] ] != 0) break;
      if (IMGH2V_DOT(ogv, ngv[midx]) > 0.80) {	// a heuristic to check the straightness
	if (++straight_count > 5) { 		//   of the edge to break at a corner
	  once_straight = true;
	  IMGH2V_COPY( once_straight_gv, ogv );
	  straight_count = 0;
	}
      } else {
	straight_count = 0;
	IMGH2V_COPY( ogv, ngv[midx] );
      }
      if (once_straight &&  // break at a corner turning by more than 60 degrees
	  IMGH2V_DOT(once_straight_gv, ngv[midx]) < 0.5) break;
      // check jittering
      if (len>=3 && abs(nb[midx]-hist[len-2])==2 && abs(cidx-hist[len-3])==w) {
	int pidx = hist[len-2]; 
	epp[pidx] = 0; 	hist[len-2] = cidx;  len--;
      }
      // set values to proceed
      prev = cidx;
      cidx = nb[midx];  IMGH2V_COPY( cgv, ngv[midx] );
      epp[ cidx ] = oidx;				// set this pixel to 'oidx'
      if (hist) hist[len] = cidx;
      len++;
    }
    // if (debug) printf("\n");
    return len;
  }
  
  // -----------------------------------------------------------------
  // 
  // -----------------------------------------------------------------
  
public:
  int  getEdgePoints(int idx, int array[], int maxlen=0) {
    // Assuming all the edges were found by 'findEdges()' and 
    //   'idx' is the index of the first pixel of the edge,
    //   return the list of pixels for the edge. 
    if (maxlen == 0) maxlen = max_edge_length;
    int   i, cidx, nb[16], count=0, w2 = 2*w, x=0, y=0, nbc=0;
    int   *emapp=(int*)(emap.data);
    idx = emapp[idx];	// start from the head of the edge
    if (idx == 0) return 0;
    array[ count++ ] = cidx = idx;
    emapp[ idx ] = 0;                             // change the pixel value
    //bool debug = (idx == 943);
    //if (debug) printf("getEdgePoints from %d : ", idx);
    while (count < maxlen) {	// search 1st ring of the current pixel
      y = cidx / w;  x = cidx % w;
      if        (y==0) {
	if      (x==0)   { nbc=3; IMGH3V_SET( nb, cidx+w, cidx+1, cidx+w+1 ); }
	else if (x==w-1) { nbc=3; IMGH3V_SET( nb, cidx+w, cidx-1, cidx+w-1 ); }
	else             { nbc=5; IMGH5V_SET( nb, cidx+w, cidx-1, cidx+1, cidx+w-1, cidx+w+1 ); }
      } else if (y==h-1) {
	if      (x==0)   { nbc=3; IMGH3V_SET( nb, cidx-w, cidx+1, cidx-w+1 ); }
	else if (x==w-1) { nbc=3; IMGH3V_SET( nb, cidx-w, cidx-1, cidx-w-1 ); }
	else             { nbc=5; IMGH5V_SET( nb, cidx-w, cidx-1, cidx+1, cidx-w-1, cidx-w+1 ); }
      } else if (x==0) {
	nbc=5; IMGH5V_SET( nb, cidx+1, cidx-w, cidx+w, cidx+1-w, cidx+1+w );
      } else if (x==w-1) {
	nbc=5; IMGH5V_SET( nb, cidx-1, cidx-w, cidx+w, cidx-1-w, cidx-1+w );
      } else {				// search 1st ring of the current pixel
	nbc=8;
	IMGH4V_SET( nb+0,  cidx-w  , cidx+w  , cidx  -1, cidx  +1 );
	IMGH4V_SET( nb+4,  cidx-w-1, cidx-w+1, cidx+w-1, cidx+w+1 );
      }
      //if (debug) printf("(%d,%d)%c ", x,y, (emapp[nb[3]]==idx ? 'Y':'N'));
      for (i = 0; i < nbc; i++) if (emapp[nb[i]] == idx) break;
      if (i < nbc) {			// found
	array[ count++ ] = cidx = nb[i];              // save the next pixel
	emapp[nb[i]] = 0;                             // change the pixel value
      } else if (y > 1 && y < h-2) {	// search 2nd ring of the current pixel
	IMGH4V_SET( nb+0,  cidx-w2, cidx+w2, cidx-2, cidx+2 );
	IMGH4V_SET( nb+4,  cidx-w2-1, cidx-w2+1, cidx+w2-1, cidx+w2+1 );
	IMGH4V_SET( nb+8,  cidx-w-2,  cidx+w-2,  cidx-w+2,  cidx+w+2 );
	IMGH4V_SET( nb+12, cidx-w2-2, cidx-w2+2, cidx+w2-2, cidx+w2+2 );
//       IMGH4V_SET( nb+0, cidx-w3, cidx+w3, cidx-3, cidx+3 );
//       IMGH4V_SET( nb+4, cidx-w3-1, cidx-w3+1, cidx+w3-1, cidx+w3+1 );
//       IMGH4V_SET( nb+8, cidx-w-3, cidx+w-3, cidx-w+3, cidx+w+3 );
	for (i = 0; i < 16; i++) if (emapp[nb[i]] == idx) break;
	if (i == 16) break;  // not found
	array[ count++ ] = cidx = nb[i];              // save the next pixel
	emapp[nb[i]] = 0;                             // change the pixel value
      } else break;
    }
    for (i=0; i<count; i++) emapp[array[i]] = idx;  // restore the pixel values
    //if (debug) printf("\n");
    return count;
  }
  
public:
  int breakIntoTwoLines(int plist[], int n) {
    // Break a line into two, and return the size of the first.
    float  ogv[2], cgv[2], once_straight_gv[2];
    int    i, j, trial, idx, straight_count=0;
    bool   once_straight = false;
    grd->getGradientVector( n-1, ogv, true );
    float  c_same[2] = { 0.985f, 0.939f };   // (less than) 10 or 20 degree
    float  c_diff[2] = { 0.906f, 0.707f };   // (more than) 25 or 45 degree
    for (trial = 0; trial < 2; trial++) {
      for (i = n-2; i >= 0; i--) {	// scan backward
	idx = plist[i];
	grd->getGradientVector( idx, cgv, true );
	if (IMGH2V_DOT(ogv, cgv) > c_same[trial]) {
	  if (++straight_count > 5) { 
	    once_straight = true;
	    IMGH2V_COPY( once_straight_gv, ogv );
	    straight_count = 0;
	  }
	} else {
	  straight_count = 0;
	  IMGH2V_COPY( ogv, cgv );
	}
	if (once_straight && IMGH2V_DOT(once_straight_gv, cgv) < c_diff[trial]) break;
      }
      if (i > 0) break;
    }
    if (i > 0) {
      // update the edge map image 'emap'
      for (j = i+1; j < n; j++)  emap.setInt( plist[j], plist[i+1] );
      return (i+1);
    } else return 0;
  }
  
  // -----------------------------------------------------------------
  // Test 
  // -----------------------------------------------------------------
public:
  int findEdgePathUsingStraightLine(int p0, int p1, int array[], double tol=1.5) {
    // Find a straight edge that connects pixel 'p0' and 'p1',
    //   with a tolerance 'tol' (distance from the straight line).
    int    x0=p0%w, y0=p0/w, xy1[2]={p1%w, p1/w}, w2 = 2*w;
    int    i, min_i, pos=0, cidx, nb[9], *emapp=(int*)(emap.data);
    double eq[3], len, cxy[2], ldist, min_ldist, ddist, last_ddist=9999;
    IMGH2V_SET( eq, -(xy1[1]-y0), (xy1[0]-x0) );
    len = IMGH2V_LEN( eq );  IMGH2V_DIV( eq, len );
    eq[2] = - (eq[0] * x0 + eq[1] * y0);
    // find the starting edge point from (cx,cy)
    cidx = y0 * w + x0;
    if (!emapp[cidx]) {
      IMGH4V_SET( nb+0,  cidx-w, cidx+w, cidx-1, cidx+1 );
      float gv[2], cosv, max_cosv=-1;
      for (i=0; i<4; i++) {
	grd->getGradientVector( nb[i], gv, true );
	cosv = (float)IMGH2V_DOT( eq, gv );
	if (cosv > max_cosv) { cidx = nb[i]; max_cosv = cosv; }
      }
      if (max_cosv < 0.5) return -1;
    }
    array[ pos++ ] = cidx;
    // follow the edges while maintaining distance from the straight line
    while (true) {
      IMGH4V_SET( nb+0,  cidx-w, cidx+w, cidx-1, cidx+1 );
      IMGH4V_SET( nb+4,  cidx-w-1, cidx-w+1, cidx+w-1, cidx+w+1 );
      for (i = 0, min_ldist = 2*tol; i < 8; i++) {
	if (!emapp[nb[i]]) continue;
	IMGH2V_SET( cxy, nb[i]%w, nb[i]/w );
	ddist = IMGH2V_DIST( xy1, cxy );  if (ddist<last_ddist) continue;
	ldist = IMGH2V_LPDIST( eq, cxy );  if (ldist>tol) continue;
	if (ldist < min_ldist) { min_ldist = ldist; min_i = i; }
      }
      if (min_ldist > tol) {	// search 2nd ring of the current pixel
	IMGH4V_SET( nb+0,  cidx-w2, cidx+w2, cidx-2, cidx+2 );
	IMGH4V_SET( nb+4,  cidx-w2-1, cidx-w2+1, cidx+w2-1, cidx+w2+1 );
	IMGH4V_SET( nb+8,  cidx-w-2,  cidx+w-2,  cidx-w+2,  cidx+w+2 );
	IMGH4V_SET( nb+12, cidx-w2-2, cidx-w2+2, cidx+w2-2, cidx+w2+2 );
	for (i = 0; i < 16; i++) {
	  if (!emapp[nb[i]]) continue;
	  IMGH2V_SET( cxy, nb[i]%w, nb[i]/w );
	  ddist = IMGH2V_DIST( xy1, cxy );  if (ddist<last_ddist) continue;
	  ldist = IMGH2V_LPDIST( eq, cxy );  if (ldist>tol) continue;
	  if (ldist < min_ldist) { min_ldist = ldist; min_i = i; }
	}
      }
      if (min_ldist > tol) break;  // not found
      array[ pos++ ] = cidx = nb[min_i];
      last_ddist = ddist;
    }
    return (last_ddist < tol ? pos : -1);
  }
  void setEdgePath(int array[], int n) {
    // Update the edge map 'emap' using a new path in 'array[]'.
    int  i, j, k, cidx, pidx, ne, eidx;
    int  *emapp=(int*)(emap.data);
    int  *oldedge = (int*)malloc((w+h)*sizeof(int));
    if (!oldedge) return;
    for (i = 0; i < n; i++) {
      cidx = array[i];   pidx = emapp[cidx];
      if (pidx <= 0) continue;
      ne = getEdgePoints( pidx, oldedge, w+h );	// get old edge path
      for (j = 0; j < ne; j++) {  // break old edge path into pieces
	eidx = oldedge[j];
	for (k = 0; k < n; k++) if (array[k]==eidx) break;
	if (k<n) { emapp[eidx] = -1; pidx = 0; }
	else if (pidx == 0) emapp[eidx] = pidx = eidx;
	else emapp[eidx] = pidx;
      }
    }
    for (i = 0; i < n; i++) emapp[ array[i] ] = array[0];
    free(oldedge);
  }
  
  
  // -----------------------------------------------------------------
  // 
  // -----------------------------------------------------------------
  
public:
  void showEdges(Image *vimg, int eidx=0) {
    if (w <= 0 || h <= 0 || !emap.data || !vimg) return;
    vimg->setImage( w, h, PIXEL_RGB, true);
    showEdges( vimg->data, PIXEL_RGB, eidx );
  }
  void showEdges(void* image_buffer, pixel_t type, int eidx=0) {
    // Visualize edges on 'image_buffer'.
    // If 'eidx' <= 0, visualize all the edges.
    if (w <= 0 || h <= 0 || !emap.data || !image_buffer) return;
    int  i, *src = (int*)emap.data;
    Image  dst( w, h, type, image_buffer );
    switch (type) {
    case PIXEL_GRAY:  case PIXEL_GRAYA:
      for (i = 0; i < total; i++) {
	if (src[i]) srand( src[i] ); else continue;
	if (eidx > 0 && src[i] != eidx) continue;
	dst.setGray( i, 100 + rand()%156 );
      }
      break;
    case PIXEL_RGB:  case PIXEL_RGBA:  case PIXEL_BGR: 
      for (i = 0; i < total; i++) {
	if (src[i]) srand( src[i] ); else continue;
	if (eidx > 0 && src[i] != eidx) continue;
	dst.setRGB( i, 50+rand()%206, 50+rand()%206, 50+rand()%206 );
      }
      break;
    default: break;
    }
  }
  
  void showEdgeStrength(void* image_buffer, pixel_t type) {
    // Save the normalized edge strength in 'image_buffer'.
    // Note that image 'gz' should not be used directly, because
    //   it has edge strength values before non-maximum suppression.
    if (w <= 0 || h <= 0 || emap.data == NULL) return;
    int  i, *mpp = (int*)emap.data;
    float  *spp = (float*)(gz->data), minv, maxv;
    Image  dst( w, h, type, image_buffer );
    gz->getPixelStat( &minv, &maxv );
    if (type == PIXEL_FLOAT) {
      for (i = 0; i < total; i++) {
	if (mpp[i] == 0) continue;
	float *pp = (float*)dst.getPixel( i );
	pp[0] = spp[i] / maxv;
      }
    } else {
      for (i = 0; i < total; i++) {
	if (mpp[i] == 0) continue;
	unsigned char *pp = (unsigned char*)dst.getPixel( i );
	pp[0] = (unsigned char)(255 * spp[i]/maxv);
      }
    }
  }
  
  void getEdgeInfo(char buf[]) {
    // Get one-line information of the edges
    sprintf(buf, "emap=(%d x %d)  threshold=(%.2f %.2f)  edge_count=%d  max_edge_length=%d", 
	    emap.w, emap.h, low_val, high_val, edge_count, max_edge_length );
  }
  
  void printInfo(void) {
    if (grd == &grd_dummy) grd->printInfo();
    printf("IMGH::EdgeDetector \n");
    emap.printInfo("  emap");
    printf("  threshold=(%.2f %.2f)  edge_count=%d  max_edge_length=%d \n", 
	   low_val, high_val, edge_count, max_edge_length);
    int *epp = (int*)(emap.data);
    for (int i = 0; i < total; i++) {
      if (epp[i] == 0 || epp[i] == i) continue;
      if (epp[ epp[i] ] != epp[i]) {
	std::cerr << "Warning (IMGH::EdgeDetector): edge map 'emap' is invalid at " << epp[i] << std::endl;
      }
    }
  }
  
  int  findNearestEdgePoint(int x, int y) {
    int du[25] = { 0,  0,  0, -1, +1,  -1, +1, -1, +1,   0,  0, -2, +2,  -1, +1, -1, +1,  -2, -2, +2, +2,  -2, +2, -2, +2 };
    int dv[25] = { 0, -1, +1,  0,  0,  -1, -1, +1, +1,  -2, +2,  0,  0,  -2, -2, +2, +2,  -1, +1, -1, +1,  -2, -2, +2, +2 };
    int  i, xx, yy, nidx=0;
    for (i = 0; i < 25; i++) {
      xx = x + du[i];  yy = y + dv[i];
      if (xx < 0 || xx >= w || yy < 0 || yy >= h) continue;
      nidx = yy * w + xx;
      if (emap.getInt(nidx) > 0) break;
    }
    return (i==25 ? -1 : nidx);
  }
  void printNearestEdgeInfo(int x, int y) {
    int nidx, n, plist[4096];
    if ((nidx = findNearestEdgePoint( x, y )) < 0) return;
    n = getEdgePoints( emap.getInt(nidx), plist, 4096 );
    printf("NearestEdge %6d : len=%3d  p[%d]:%d(%3d,%3d) -> p[%d]:%d(%3d,%3d)\n", nidx, n, 0, plist[0], plist[0]%w, plist[0]/w, n-1, plist[n-1], plist[n-1]%w, plist[n-1]/w);
  }
  
  void printEdge(int idx, char *cmmt=NULL) {
    int i, n, plist[4096];
    if ((n = getEdgePoints( idx, plist, 4096 )) <= 0) return;
    printf("%s at %d (%d) : ", (cmmt ? cmmt:"  edge"), idx, n);
    for (i=0; i<n; i++) printf("%d ", plist[i]);
    printf("\n");
  }
};

  
}	// namespace IMGH


#endif	// IMGH_EDGE_DETECTOR_HPP


// ===================================================================
#if 0	// start of the example code
// ===================================================================
#include <iostream>
#include "imgh_fileio.hpp"
#include "imgh_edge_detector.hpp"
using namespace std;
int main(int argc, char **argv)
{
  // read the image from a file
  IMGH::Image  img;
  IMGH::ImageFileIO ifile;
  if (! ifile.readFile( argv[1], &img ) ) return EXIT_FAILURE;
  // find edges in the image
  IMGH::EdgeDetector edt;
  if      (1) edt.findEdges( &img );
  else if (0) edt.findEdges( img.w, img.h, IMGH::PIXEL_RGB, img.data );
  else        edt.findEdges( &img, 0.1, 0.5, 3 );  // with gaussian smoothing ('img' changed)
  // The result is in the edge map image 'emap', which has the same size
  // as the input image, with integer pixel type 'IMGH::PIXEL_INT'.
  edt.emap.printInfo();
  cout << "Number of edges : " << edt.edge_count << endl;
  cout << "Max edge length : " << edt.max_edge_length << endl;
  // access edges (series of edge points) one by one
  int i, npixels, total = img.w * img.h;
  int plist[2048], *emap = (int*)edt.emap.data;
  for (i=0; i < total; i++) {
    if (emap[i] == 0) continue;	// the pixel is not on edge
    if (emap[i] != i) continue;	// it's on edge, but not the head of the edge
    npixels = edt.getEdgePoints( i, plist );
    // Do whatever you want ... like fitting a line to the edge.
  }
  return EXIT_SUCCESS;
}
// ===================================================================
#endif	// end of the example code
// ===================================================================
