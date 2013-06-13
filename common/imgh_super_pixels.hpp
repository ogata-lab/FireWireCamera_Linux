
//
// IMGH::SuperPixels - super pixel segmentation
//   based on the original code from "Efficient Graph-Based Image Segmentation"
//   by Pedro F. Felzenszwalb and Daniel P. Huttenlocher. IJCV, V59, N2, 2004.
//
// Jaeil Choi
// last modified in Apr, 2007
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

#ifndef IMGH_SUPER_PIXELS_HPP
#define IMGH_SUPER_PIXELS_HPP

#include <iostream>
#include <cmath>
#include "vm_macros.h"
#include "imgh_common.hpp"
#include "imgh_fileio.hpp"
#include "imgh_edge_detector.hpp"
#include "algh_linked_list.hpp"
#include "algh_sorter.hpp"

namespace IMGH {
  
#define PxL(p,img)	(p - img->pixel_size)
#define PxR(p,img)	(p + img->pixel_size)
#define PxU(p,img)	(p - img->row_size)
#define PxD(p,img)	(p + img->row_size)
#define PxUL(p,img)	(p - img->row_size - img->pixel_size)
#define PxUR(p,img)	(p - img->row_size + img->pixel_size)
#define PxDL(p,img)	(p + img->row_size - img->pixel_size)
#define PxDR(p,img)	(p + img->row_size + img->pixel_size)

#define PiL(idx,img)	(pidx - 1)
#define PiR(idx,img)	(pidx + 1)
#define PiU(idx,img)	(pidx - img->w)
#define PiD(idx,img)	(pidx + img->w)
#define PiUL(idx,img)	(pidx - img->w - 1)
#define PiUR(idx,img)	(pidx - img->w + 1)
#define PiDL(idx,img)	(pidx + img->w - 1)
#define PiDR(idx,img)	(pidx + img->w + 1)

typedef struct { int rank, p, size; unsigned char minv[3], maxv[3]; } spinfo;
typedef struct { float w; int a, b; } spedge;
//static bool operator<(const spedge &a, const spedge &b) { return a.w < b.w; }
static int compare_edges(spedge a, spedge b) { return (a.w == b.w ? 0 : (a.w < b.w ? -1 : +1)); }
  
// Superpixel information
class Spp {
public:
  // variables available after calling 'segmentIntoSuperpixels()'
  int	p;		// spp ID (index of the seed pixel)
  int	size;		// spp size (number of all the pixels in the superpixel)
  // variables available after calling 'createSegmentList()'
  int	idx;		// index in the sequential list 'spplist' or array 
  int	avg[3];		// average color
  int	bbox[4];	// bounding box of the spp
public:
  Spp() : p(-1), size(0) {}
  ~Spp() {}
};
  
class SuperPixels {
private:
  // Accessible after calling 'segmentIntoSuperpixels()'
  spinfo	*elts;		// trace information for each pixel
  int		num;		// number of superpixels
  spedge	*spedges;	// pixel connections 
  int		num_spedges;	// number of pixel connections (w*h*4)
  int		w, h;		// input image size
  
public:
  // Accessible only after 'createSegmentList()'
  IMGH::Image	simg;		// (PIXEL_VOID*) map of the pointers to segments
  ALGH::LinkedList<Spp*> spplist;
  
public:
  bool		verbose;
  
  // If EdgeDetector is given, the extent of superpixels is also restricted by edges
  IMGH::EdgeDetector *edt;
  
public:
  SuperPixels() : elts(NULL), num(0), spedges(NULL), num_spedges(0), 
		  verbose(false), edt(NULL) { }
  ~SuperPixels() { clear();  }
  void clear(void) { 
    if (elts) delete [] elts;  elts = NULL;  num = 0; 
    if (spedges) delete [] spedges;  spedges = NULL;  num_spedges = 0;
    spplist.clearAndFree();
    simg.clear();
  }
  inline bool isReady(void)    { return (elts!=NULL); }
  // available after 'sementImage()'
  inline int getSegCount(void) { return (spplist.count>0 ? spplist.count : num);  }
  inline int getSegID(int pidx)      { return find(pidx);  }
  inline int getSegID(int x, int y)  { return find(y*w+x); }
  // available after 'createSegmentList()'
  Spp* findSegInfo(int pidx)     { return (Spp*)simg.getPointer(pidx);  }
  Spp* findSegInfo(int x, int y) { return (Spp*)simg.getPointer(y*w+x); }
  
  // -------------------------------------------------------------------
  // Oversegmentation (superpixels)
  // -------------------------------------------------------------------

public:
  int  segmentIntoSuperpixels(Image *img, int gfsize=0, float thr=200, int min_size=50,
			      bool use_diagonal=false, bool merge_same_color=false) {
    // Segment a graph, and 
    //   save the disjoint-set forest representing the segmentation.
    //   img     : image to segment
    //   gfsize  : size of Gaussian filter
    //   thr     : threshold constant (The bigger, the more likely to merge)
    //   min_size: minimum component size in pixels 
    clear();
    if (!img) return 0;
    if (verbose) printf("SuperPixels segmentation \n");
  
    //   int x, y, pidx, width = img.w, height = img.h;
    if (img->type != IMGH::PIXEL_RGB) {
      if (verbose) printf("  Error (SuperPixels::segment): not uchar RGB\n");
      return 0;
    }
    
    // smoothe the input image
    IMGH::Image *iimg=NULL, smoothed_img;
    if (gfsize < 3) iimg = img;
    else {
      IMGH::ImageFilter gfilter;
      smoothed_img.setImage( img->w, img->h, PIXEL_RGB );
      gfilter.convoluteWithGaussian( img, &smoothed_img, gfsize );
      iimg = &smoothed_img;
    }
    // build the graph
    spedges = new spedge[iimg->w * iimg->h * (use_diagonal ? 4:2)];
    num_spedges = buildGraph( iimg, spedges, use_diagonal );
  
    if (verbose) printf("  greedy merging with threshold_const = %g\n", thr);
    //std::sort(spedges, spedges + num_spedges);	// sort spedges by weight
    ALGH::Sorter<spedge> sorter;
    sorter.QuickSort( spedges, num_spedges, compare_edges );

    // make a disjoint-set forest
    this->w = iimg->w;
    this->h = iimg->h;
    int i, npixels = iimg->w * iimg->h;
    init(npixels);

    // init thresholds
    float *threshold = new float[npixels];
    for (i = 0; i < npixels; i++)  threshold[i] = (thr / 1);  // THRESHOLD(1,thr)
    // for each edge, in non-decreasing weight order...
    for (i = 0; i < num_spedges; i++) {
      spedge *pedge = &spedges[i];
      // components conected by this spedge
      int a = find(pedge->a);
      int b = find(pedge->b);
      if (a == b || a < 0 || b < 0) continue;
//       if ((a==19250 || b==19250) && pedge->w > 14) continue;
      if ((pedge->w <= threshold[a]) && (pedge->w <= threshold[b])) {
// 	if (a==19250 || b==19250) printf("edge %d join %d and %d at (%d and %d): rank(%d,%d) size(%d,%d) edgew=%.2f  thr(%.2f,%.2f)\n", i, a, b, pedge->a, pedge->b, elts[a].rank, elts[b].rank, elts[a].size, elts[b].size, pedge->w, threshold[a], threshold[b]);
	join(a, b);
	a = find(a);
	threshold[a] = pedge->w + (thr / size(a));  // THRESHOLD(size(a),thr)
      }
    }
  
    // post process small components  //// waste in 'for' loop ?
    if (verbose) printf("  removing small components with minsize = %d  from %d segmentations \n", min_size, num);
    for (i = 0; i < num_spedges; i++) {
      int a = find(spedges[i].a);
      int b = find(spedges[i].b);
      if ( a != b && a >= 0 && b >= 0 &&
	  (size(a) < min_size || size(b) < min_size)) join(a, b);
    }
    
    // merge segments with similar color distribution [default: false]
    if (merge_same_color && iimg->type==PIXEL_RGB) {
      for (i = 0; i < npixels; i++) {	// initialize color value
	uchar_t *p = (uchar_t*)iimg->getPixel(i);
	G3V_COPY( elts[i].minv, p );  G3V_COPY( elts[i].maxv, p );
      }
      for (i = 0; i < npixels; i++) {
	uchar_t *ip = (uchar_t*)iimg->getPixel(i);
	int s = find(i);
	G3V_MAKE_BELOW(elts[s].minv, ip);  G3V_MAKE_ABOVE(elts[s].maxv, ip);  
      }
      for (i = 0; i < num_spedges; i++) {
	spedge *pedge = &spedges[i];
	int a = find(pedge->a);
	int b = find(pedge->b);
	if (a == b || a < 0 || b < 0) continue;
	uchar_t *aminv = elts[a].minv, *amaxv = elts[a].maxv;
	uchar_t *bminv = elts[b].minv, *bmaxv = elts[b].maxv;
	if ( (pedge->w <= threshold[a]*2) && (pedge->w <= threshold[b]*2) &&
	     (G3V_TEST_ABOVE( bminv, aminv ) && G3V_TEST_BELOW( bmaxv, amaxv ) ||
	      G3V_TEST_ABOVE( aminv, bminv ) && G3V_TEST_BELOW( amaxv, bmaxv )) ) join(a, b);
      }
    }
    
    // merge small disconnected components between two edges
    if (edt && !use_diagonal) mergePixelsBetweenTwoEdges( img, threshold );
   
    delete threshold;
    if (verbose)  printf("  result : %d segmentations\n", num);
    simg.setImage( w, h, PIXEL_VOIDP, true );
    return num;
  }
  
  void showOversegmentation(Image *vimg) {
    if (!vimg) return;
    vimg->setImage( w, h, PIXEL_RGB, true);
    showOversegmentation( vimg->data, PIXEL_RGB );
  }
  void showOversegmentation(void *buffer, pixel_t type) {
    if (!elts || !buffer) return;
    // pick random colors for each component
    IMGH::Image dest( w, h, type, buffer );
    dest.clearImage();
    int i, total = w * h * dest.pixel_size;
    unsigned char *colors = new unsigned char[total];
    srand(0);
    for (i = 0; i < total; i++)  colors[i] = (uchar_t)(rand()%256);
  
    int x, y, cidx, pidx = 0, ppidx;
    switch (type) {
    case PIXEL_INT:
      for (y = pidx = 0; y < dest.h; y++) 
	for (x = 0; x < dest.w; x++, pidx++) {
	  cidx = find(pidx);
	  dest.setInt( pidx, (cidx<0 ? -1 : cidx) );
	}
      break;
    case PIXEL_GRAY:  case PIXEL_GRAYA:
      for (y = pidx = 0; y < dest.h; y++) 
	for (x = 0; x < dest.w; x++, pidx++) {
	  if ((cidx = find(pidx)) < 0) continue;
	  dest.setGray( pidx, colors[cidx] );
	}
      break;
    case PIXEL_RGB:  case PIXEL_RGBA:  case PIXEL_BGR: 
      for (y = pidx = 0; y < dest.h; y++) 
	for (x = 0; x < dest.w; x++, pidx++) {
	  if ((ppidx = find(pidx)) < 0) continue;
	  cidx = ppidx * 3;
	  dest.setRGB( pidx, colors[cidx+0], colors[cidx+1], colors[cidx+2] );
	}
      break;
    default: break;
    }
    delete(colors);
  }
  
  // -----------------------------------------------------------------
  // private functions
  // -----------------------------------------------------------------
  
 private:
  void init(int npixels) {
    if (elts) { delete [] elts; }
    elts = new spinfo[npixels];
    num = npixels;
    for (int i = 0; i < npixels; i++) {
      elts[i].rank = 0;  elts[i].size = 1;  elts[i].p = i;
      if (edt && edt->emap.getInt(i) != 0) { num--; elts[i].p = -1; }
    }
  }
  int find(int x) {
    int y = x;
    while (y >= 0 && y != elts[y].p) y = elts[y].p;
    elts[x].p = y;
    return y;
  }
  int join(int x, int y) {
    if (x < 0 || y < 0) printf("strange join %d %d\n", x, y);
    if (elts[x].rank > elts[y].rank) {
      elts[y].p = x;   elts[x].size += elts[y].size;
      num--;  return x;
    } else {
      elts[x].p = y;   elts[y].size += elts[x].size;
      if (elts[x].rank == elts[y].rank)	elts[y].rank++;
      num--;  return y;
    }
  }
  int size(int x) const { return elts[x].size; }
  int num_sets()  const { return num; }
  
  float diff(uchar_t *p1, uchar_t *p2) {
#if 1
    // Euclidean distance in its 3-channel color space
    return sqrt( (p1[0] - (float)p2[0]) * (p1[0] - (float)p2[0]) +
		 (p1[1] - (float)p2[1]) * (p1[1] - (float)p2[1]) +
		 (p1[2] - (float)p2[2]) * (p1[2] - (float)p2[2]) );
#else
    float d0 = fabs(p1[0]-(float)p2[0]);
    float d1 = fabs(p1[1]-(float)p2[1]);
    float d2 = fabs(p1[2]-(float)p2[2]);
    return (d0 + d1 + d2);
//     return (d0>d1 ? (d0>d2 ? d0 : d2) : (d1>d2 ? d1 : d2));
#endif
  }
  
#define ADD_SPEDGE(aa, bb, ww) \
  do { if (!emp || emp[bb]<=0) { spedges[num_spedges].a = aa; spedges[num_spedges].b = bb; spedges[num_spedges].w = ww; num_spedges++; } } while(0)
    
  int  buildGraph(Image *img, spedge *spedges, bool use_diagonal) {
    int x, y, pidx, w = img->w, h = img->h;
    int num_spedges = 0, *emp=(edt ? (int*)(edt->emap.data):NULL);
    if (img->pixel_size == 3) {
      for (y = pidx = 0; y < h; y++) {
	uchar_t *p = (uchar_t*)img->data + y * img->row_size;
	for (x = 0; x < w; x++, pidx++, p+=3) {
	  if (emp && emp[pidx]>0) continue;
	  if (x < w-1)			// spedge to R
	    ADD_SPEDGE( pidx, PiR(pidx,img), diff(p,PxR(p,img)) );
	  if (y < h-1)			// spedge to D
	    ADD_SPEDGE( pidx, PiD(pidx,img), diff(p,PxD(p,img)) );
	  if (use_diagonal) {
	    if ((x < w-1) && (y < h-1))	// spedge to DR
	      ADD_SPEDGE( pidx, PiDR(pidx,img), diff(p,PxDR(p,img)) );
	    if ((x < w-1) && (y > 0))	// spedge to UR
	      ADD_SPEDGE( pidx, PiUR(pidx,img), diff(p,PxUR(p,img)) );
	  }
	}
      }
    } else if (img->pixel_size == 1) {
      for (y = pidx = 0; y < h; y++) {
	uchar_t *p = (uchar_t*)img->data + y * img->row_size;
	for (x = 0; x < w; x++, pidx++, p++) {
	  if (emp && emp[pidx]>0) continue;
	  if (x < w-1)			// spedge to R
	    ADD_SPEDGE( pidx, PiR(pidx,img), fabs((float)*p - *PxR(p, img)) );
	  if (y < h-1)			// spedge to D
	    ADD_SPEDGE( pidx, PiD(pidx,img), fabs((float)*p - *PxD(p, img)) );
	  if (use_diagonal) {
	    if ((x < w-1) && (y < h-1))	// spedge to DR
	      ADD_SPEDGE( pidx, PiDR(pidx,img), fabs((float)*p - *PxDR(p, img)) );
	    if ((x < w-1) && (y > 0))	// spedge to UR
	      ADD_SPEDGE( pidx, PiUR(pidx,img), fabs((float)*p - *PxUR(p, img)) );
	  }
	}
      }
    }
    return num_spedges;
  }
  
  void mergePixelsBetweenTwoEdges(Image *img, float threshold[]) {
    int  i, j, w2=w*2, pmin=w*3+3, pmax=w*h-w*3-3, a, am, ap, b, bm, bp;
    for (i = pmin; i < pmax; i++) {
      a = find(i);  if (a < 0) continue;
      am = edt->emap.getInt(i-1);
      ap = edt->emap.getInt(i+1);
      j = i;
      if (am>0 && ap>0 && elts[a].size < (w>>1)) {  // pixels between two vertical edges
	if      ((bm = edt->emap.getInt(i+w+0)) == am &&
		 (bp = edt->emap.getInt(i+w+2)) == ap)   j = (i+w+1);
	else if ((b = find(i+w-1)) >= 0 && b != a &&		// DL
		 (bm = edt->emap.getInt(i+w-2)) == am &&
		 (bp = edt->emap.getInt(i+w-0)) == ap)   j = i+w-1;
	else if ((bm = edt->emap.getInt(i+w+0)) == am &&	// D.R
		 (bp = edt->emap.getInt(i+w+1)) == ap &&
		 (bm = edt->emap.getInt(i+w2+0)) == am &&
		 (bp = edt->emap.getInt(i+w2+2)) == ap)  j = i+w2+1;
	else if ((bm = edt->emap.getInt(i+w-1)) == am &&	// D.L
		 (bp = edt->emap.getInt(i+w-0)) == ap &&
		 (bm = edt->emap.getInt(i+w2-2)) == am &&
		 (bp = edt->emap.getInt(i+w2-0)) == ap)  j = i+w2-1;
	if (i != j && (b=find(j)) >= 0 && a != b && elts[b].size<(w>>1)) {
	  uchar_t *ip = (uchar_t*)img->getPixel(i);
	  uchar_t *jp = (uchar_t*)img->getPixel(j);
	  float ww = (img->type == PIXEL_GRAY ? abs(ip[0]-jp[0]) : diff( ip, jp ));
	  if (ww <= threshold[a]*2 && ww <= threshold[b]*2) join(a,b);
	}
      }
      am = edt->emap.getInt(i-w);
      ap = edt->emap.getInt(i+w);
      j = i;
      if (am>0 && ap>0 && elts[a].size<(w>>1)) {  // pixels between two horizontal edges
	if      ((bm = edt->emap.getInt(i-w2+1)) == am &&
		 (bp = edt->emap.getInt(i   +1)) == ap)  j = i-w+1;
	else if ((bm = edt->emap.getInt(i   +1)) == am &&
		 (bp = edt->emap.getInt(i+w2+1)) == ap)  j = i+w+1;
	else if ((bm = edt->emap.getInt(i-w+1)) == am &&	// R.U
		 (bp = edt->emap.getInt(i  +1)) == ap &&
		 (bm = edt->emap.getInt(i-w2+2)) == am &&
		 (bp = edt->emap.getInt(i   +2)) == ap)  j = i-w+2;
	else if ((bm = edt->emap.getInt(i  +1)) == am &&	// R.D
		 (bp = edt->emap.getInt(i+w+1)) == ap &&
		 (bm = edt->emap.getInt(i   +2)) == am &&
		 (bp = edt->emap.getInt(i+w2+2)) == ap)  j = i+w+2;
	if (i != j && (b=find(j)) >= 0 && a != b && elts[b].size<(w>>1)) {
	  uchar_t *ip = (uchar_t*)img->getPixel(i);
	  uchar_t *jp = (uchar_t*)img->getPixel(j);
	  float ww = (img->type == PIXEL_GRAY ? abs(ip[0]-jp[0]) : diff( ip, jp ));
	  if (ww <= threshold[a] && ww <= threshold[b]) join(a,b);
	}
      }
    }
  }
    
  
  // -----------------------------------------------------------------
  // Segmentation List
  // -----------------------------------------------------------------
public:  
  void createSegmentList(Image *img, int min_size=0) {
    if (!img) return;
    int i, ppidx, x, y, total=img->w * img->h;
    unsigned char *pixel;
    Spp *spp;
    // clear old segment information
    spplist.clearAndFree();
    // create new segment information
    simg.setImage( img->w, img->h, PIXEL_VOIDP, true );
    for (i=0; i<total; i++) {
      x = i % img->w;   y = i / img->w;
      if ((ppidx = find(i)) < 0) continue;
      pixel = (unsigned char*)img->getPixel(i);
      if (simg.getPointer(ppidx) == NULL) {
	spp = new Spp;
	spp->p = ppidx;
	spp->size = 1;
	if (img->pixel_size == 3) {
	  spp->avg[0] = (int)pixel[0];
	  spp->avg[1] = (int)pixel[1];
	  spp->avg[2] = (int)pixel[2];
	} else if (img->pixel_size == 1) {
	  spp->avg[0] = (int)pixel[0];
	}
	spp->bbox[0] = spp->bbox[2] = x;
	spp->bbox[1] = spp->bbox[3] = y;
	simg.setPointer( ppidx, spp );
	spplist.append( spp );
      } else {
	spp = (Spp*)simg.getPointer( ppidx );
	spp->size++;
	if (img->pixel_size == 3) {
	  spp->avg[0] += (int)pixel[0];
	  spp->avg[1] += (int)pixel[1];
	  spp->avg[2] += (int)pixel[2];
	} else if (img->pixel_size == 1) {
	  spp->avg[0] += (int)pixel[0];
	}
	if (x < spp->bbox[0]) spp->bbox[0] = x;
	if (y < spp->bbox[1]) spp->bbox[1] = y;
	if (x > spp->bbox[2]) spp->bbox[2] = x;
	if (y > spp->bbox[3]) spp->bbox[3] = y;
      }
      simg.setPointer( i, spp );
    }
    for (spplist.goFirst(); (spp = spplist.getCurr()); spplist.goNext()) {
      spp->avg[0] /= spp->size;  spp->avg[1] /= spp->size;  spp->avg[2] /= spp->size;
    }
    // remove too small (less than 'min_size') superpixels between edges
    if (edt && min_size > 1) {
      for (i = 0; i < total; i++) {
	spp = (Spp*)simg.getPointer( i );
	if (spp && spp->size < min_size) {
	  simg.setPointer( i, NULL );  elts[i].p = -1;
	}
      }
      for (spplist.goFirst(); (spp = spplist.getCurr()); )
	if (spp->size < min_size) spplist.remCurr(); else spplist.goNext();
    }
  }
  
  void updateSegmentAvgColor(Image *img) {
  }
  
  void mergeSegmentByAvgColor(float min_diff) {
    if (min_diff <= 0) return;
    int oldnum = num, joined=0, joinedseg=0;
    for (int i = 0; i < num_spedges; i++) {
      if (spedges[i].w >= min_diff) continue;
      int a = find(spedges[i].a);
      int b = find(spedges[i].b);
      if (a == b || a < 0 || b < 0) continue;
      Spp *sa = (Spp*)simg.getPointer(a);
      Spp *sb = (Spp*)simg.getPointer(b);
      float diff = sqrt( (sa->avg[0] - (float)sb->avg[0]) * (sa->avg[0] - (float)sb->avg[0]) +
			 (sa->avg[1] - (float)sb->avg[1]) * (sa->avg[1] - (float)sb->avg[1]) +
			 (sa->avg[2] - (float)sb->avg[2]) * (sa->avg[2] - (float)sb->avg[2]) );
      if (diff >= min_diff) continue;
      int winner = join(a, b);
      joined++;
      if (sa && sb) {
	int avg[3], size, bbox[4];
	size = (sa->size + sb->size);
	avg[0] = (sa->avg[0] * sa->size + sb->avg[0] * sb->size) / size;
	avg[1] = (sa->avg[1] * sa->size + sb->avg[1] * sb->size) / size;
	avg[2] = (sa->avg[2] * sa->size + sb->avg[2] * sb->size) / size;
	bbox[0] = (sa->bbox[0] < sb->bbox[0] ? sa->bbox[0] : sb->bbox[0]); // xmin
	bbox[1] = (sa->bbox[1] < sb->bbox[1] ? sa->bbox[1] : sb->bbox[1]); // ymin
	bbox[2] = (sa->bbox[2] > sb->bbox[2] ? sa->bbox[2] : sb->bbox[2]); // xmax
	bbox[3] = (sa->bbox[3] > sb->bbox[3] ? sa->bbox[3] : sb->bbox[3]); // ymax
	if (winner == a) {
	  sa->p = a;
	  sa->size = size;
	  memcpy(sa->avg,  avg,  3*sizeof(int));
	  memcpy(sa->bbox, bbox, 4*sizeof(int));
	  simg.setPointer( b, NULL );
	  spplist.findAndRemove( sb );
	} else {
	  sb->p = b;
	  sb->size = size;
	  memcpy(sb->avg,  avg,  3*sizeof(int));
	  memcpy(sb->bbox, bbox, 4*sizeof(int));
	  simg.setPointer( a, NULL );
	  spplist.findAndRemove( sa );
	}
	joinedseg++;
      }
    }
    if (verbose) printf("  segments merged : %d => %d   (joined=%d  jseg=%d)\n", 
			oldnum, num, joined, joinedseg);
  }
  int findSegment(Spp **slist, int max_count, bool absolute_color, 
		  int r, int rt, int g, int gt, int b, int bt, 
		  int bbox[4]=NULL, int whmin[2]=NULL, int whmax[2]=NULL, 
		  float bbMinFill=0.0) {
    // Find segments with similar color within the region 'bbox[4]', 
    //   and save them in 'slist' in the order of decreasing segment size.
    int i, count = 0, bbw, bbh;
    IMGH::Spp *spp=NULL, *tmp=NULL;
    for (spplist.goFirst(); (spp = spplist.getCurr()); spplist.goNext()) {
      // check its color
      if (absolute_color) {
	if (rt!=0 && abs(spp->avg[0] - r) > rt) continue;
	if (gt!=0 && abs(spp->avg[1] - g) > gt) continue;
	if (bt!=0 && abs(spp->avg[2] - b) > bt) continue;
      } else {
	if (rt!=0 && abs(spp->avg[0] - (spp->avg[1]+r)) > rt) continue;
	if (gt!=0 && abs(spp->avg[1] - (spp->avg[2]+g)) > gt) continue;
	if (bt!=0 && abs(spp->avg[2] - (spp->avg[0]+b)) > bt) continue;
      }
      // check its bounding box (position and size)
      bbw = (spp->bbox[2] - spp->bbox[0] + 1);
      bbh = (spp->bbox[3] - spp->bbox[1] + 1);
      if ( bbox &&
	   (spp->bbox[0] < bbox[0] || spp->bbox[1] < bbox[1] ||  // out of bounding box
	    spp->bbox[2] > bbox[2] || spp->bbox[3] > bbox[3]) ) continue;
      if ( whmin && (bbw < whmin[0] || bbh < whmin[1]) ) continue;
      if ( whmax && (bbw > whmax[0] || bbh > whmax[1]) ) continue;
      if ( bbMinFill>0 && spp->size < (bbw * bbh * bbMinFill) ) continue;
      // put it in 'slist[]'
      for (i = 0; i < count; i++) if (slist[i]->size < spp->size) break;
      if (i >= max_count) {	// 'spp' should be the last, but 'max_count' reached.
	continue;
      } else if (i >= count) {	// 'spp' will be the last
	slist[i] = spp;  count++;
      } else {			// 'spp' will be inserted in the middle
	for (; i <= count && i < max_count; i++) {
	  tmp = slist[ i ];
	  slist[i] = spp;
	  spp = tmp;
	}
	if (count < max_count) count++;
      }
    }
    return count;
  }
  
  int findSegmentExtremePixel(Spp *spp, float dx, float dy, bool inside) {
    // Find a pixel at maximum distance in the direction of (dx, dy) in the segment.
    if (spp == NULL || elts == NULL) return -1;
    int  i, checked, ppidx, pbest = -1, total=w*h;
    int  x, y, x2, y2, nb;
    for (i=checked=0; i<total; i++) {
      ppidx = find(i);
      if (ppidx != spp->p && (inside || ppidx >= 0)) continue;
      x  = i % w;   y = i / w;
      x2 = (int)(x - dx + 0.5);
      y2 = (int)(y - dy + 0.5);
      nb = y2 * w + x2;
      if (nb < 0 || nb >= total || find(nb) != spp->p) continue;
      if (pbest < 0 || x*dx + y*dy > (pbest%w)*dx + (pbest/w)*dy) pbest = i;
      checked++;
    }
    //printf("pbest=%d  %d pixesl checked\n", pbest, checked);
    return pbest;
  }
  
  void printSegmentArray(char* cmmt, Spp **slist, int count) {
    char buf[256];
    printf("%s (%d)\n", (cmmt ? cmmt : "SegmentArray"), count);
    for (int i = 0; i < count; i++) {
      if (getSegmentInfo(slist[i], buf) == NULL) continue;
      printf("  %s\n", buf);
    }
  }
  
  // -----------------------------------------------------------------
  // 
  // -----------------------------------------------------------------
public:
  char* getSegmentInfo(Spp *s, char *buf) {
    if (!s || !buf) return NULL;
    sprintf(buf, "Seg: size=%5d avg(%3d %3d %3d) bbox(%3d %3d)-(%3d %3d):(%d %d)",
	    s->size, s->avg[0], s->avg[1], s->avg[2], 
	    s->bbox[0], s->bbox[1], s->bbox[2], s->bbox[3], 
	    (s->bbox[2] - s->bbox[0]+1), (s->bbox[3] - s->bbox[1]+1) );
    return buf;
  }
  void printInfo(void) {
    printf("Oversegmentation using SuperPixels\n");
    printf("  input : (%d x %d) image  \n", w, h);
    printf("  output: num_elts = %d  num_edges = %d \n", num, num_spedges);
    printf("  segments: total=%d \n", spplist.count);
    Spp *spp=NULL;
    char     buf[256];
    for (spplist.goFirst(); (spp = spplist.getCurr()); spplist.goNext()) {
      if (spp->size != elts[spp->p].size) {
	printf("    Error : %s\n", getSegmentInfo(spp, buf));
	printf("            elts[%d]: size=%d\n", spp->p, elts[spp->p].size);
      }
      if (spp->size < 20) continue;
      getSegmentInfo( spp, buf );
      printf("    %s\n", buf);
    }
  }

};
  
}	// namespace IMGH


#endif	// IMGH_SUPER_PIXELS_HPP


