
//
// IMGH::CornerDetector 
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

#ifndef IMGH_CORNER_DETECTOR_HPP
#define IMGH_CORNER_DETECTOR_HPP

#include <iostream>
#include "imgh_common.hpp"
#include "imgh_drawer.hpp"
#include "imgh_converter.hpp"

namespace IMGH {
  
class CornerPoint {
 public:
  float  xy[2];		// location (potentially subpixel)
  float  strength;	// 
 public:
  CornerPoint() {}
  CornerPoint(float x, float y, float strength) { set(x, y, strength); }
  ~CornerPoint() {}
  void set(float x, float y, float strength) { 
    this->xy[0] = x;   this->xy[1] = y;
    this->strength = strength;
  }
};

  
class CornerDetector
{
public:
  int        w, h, total;	// frequently used values
  IMGH::Image   smap;		// signal strength image (PIXEL_FLOAT)
  IMGH::Image   cmap;		// corner point image    (PIXEL_VOIDP)
  float        max_strength;	// maximum of the corner strength
  float        low_strength;	// low-limit for corner strength
  int          corner_count;
  IMGH::ImageGradient *grd;
private:
  IMGH::ImageGradient grd_dummy;
  
public:
  CornerDetector() : corner_count(0) { }
  CornerDetector(IMGH::Image *img, float thr=0.1) : corner_count(0) { findCorners(img, thr); }
  ~CornerDetector() { }
  
  void clearMap(void) {
    if (cmap.data) {
      int i, total = cmap.w * cmap.h;
      CornerPoint **mp = (CornerPoint**)(cmap.data);
      for (i=0; i<total; i++) if (mp[i]) { free(mp[i]); mp[i] = NULL; }
    }
  }
  
  // -----------------------------------------------------------------
  // -----------------------------------------------------------------
public:
  void findCorners(IMGH::Image *img, float thr=0.1) {
    grd_dummy.getGradientImage( img );
    findCorners( &grd_dummy, thr );
  }
  void findCorners(int w, int h, pixel_t type, void *udata, float thr=0.1) {
    grd_dummy.getGradientImage( w, h, type, udata );
    findCorners( &grd_dummy, thr );
  }
  int  findCorners(ImageGradient *grd, float thr=0.1) {
    // Find corners using Harris corner detector, assumming 
    //   gradient images 'grd->gx' and 'grd->gy' are ready (by getGradient()).
    if (!grd || grd->gx.w == 0 || grd->gx.h == 0 || !grd->gx.sameSize( &grd->gy )) return 0;
    this->grd = grd;
    this->w  = grd->gx.w;  this->h = grd->gx.h;  this->total = w * h;
    int   x, y, idx, margin=4;
    clearMap();
    cmap.setImage( w, h, PIXEL_VOIDP );
    smap.setImage( w, h, PIXEL_FLOAT );
    // Average intensity change in direction (u,v) can be expressed as a linear form:
    //   E(u,v) = [u v] M [u],  where M =  SUM  W(x,y) [ Ix^2 IxIy ]
    //                    [v]             (x,y)        [ IxIy Iy^2 ]
    // Describe a point in terms of eigenvalues of M (measurement of corner response):
    //   R = Lambda1 * Lambda2 - kappa * (Lambda1 + Lambda2),  where kappa=0.04~0.06
    //     = det(M) - kappa * tr(M)^2               (OpenCV : icvCalcHarris() in icvCornerEigenValsVecs())
    //     = (Ix2 Iy2 - IxIy IxIy) - k (Ix2+Iy2)^2  (Gandalf: calc_strength() in gan_harris_corner_q())
    //   R = Ix2 Iyy + Iy2 Ixx - 2 Ix Iy IxIy       (OpenCV : PreCornerDetect())
    // A good corner point should have a large intensity change in all directions,
    //   i.e. R should be large positive. (Lambda_min should be large.)
    float *gxp = (float*)(grd->gx.data), *gyp = (float*)(grd->gy.data);
    float gx, gy, gx2, gy2, gxy, value, *curr, *prev, *next;
    // calculate corner strength for each pixel
    for (max_strength = idx = 0; idx < total; idx++) {
      gx  = gxp[idx];  gy  = gyp[idx];  
      gx2 = gx * gx;   gy2 = gy * gy;   gxy = gx * gy;
      value = (float)(fabs((gx2 * gy2 - gxy * gxy) - 0.04 * (gx2 + gy2) * (gx2 + gy2)));  //////
      if (value < 0.0) value = 0.0;
      smap.setFloat( idx, value );
      if (value > max_strength) max_strength = value;
    }
    // decide threshold values (default : 10% of maximum strength)
    if (thr <= 0.01) thr = 0.01;
    if (thr >= 0.99) thr = 0.99;
    // low_strength = (float)pow(pow((double)max_strength,0.25) * thr, 2.0);
    low_strength = max_strength * thr;
    // search for peaks in corner strength
    for (y = margin, corner_count=0; y < h-margin; y++) {
      curr = (float*)(smap.data) + y * w + margin;
      prev = curr - w;
      next = curr + w;
      for (x = margin; x < w-margin; x++, curr++,  prev++,  next++) {
	if (curr[0] < low_strength) continue;
	// non-maximum suppression
	if ( curr[0] < prev[-1] || curr[0] < prev[0] || curr[0] < prev[+1] ||
	     curr[0] < curr[-1] ||                      curr[0] < curr[+1] ||
	     curr[0] < next[-1] || curr[0] < next[0] || curr[0] < next[+1] ) continue;
	//// quadratic fitting for subpixel accuracy
	cmap.setPointer( x, y, new CornerPoint(x, y, curr[0]) );
	corner_count++;
      }
    }
    return corner_count;
  }
  
  // -----------------------------------------------------------------
  // -----------------------------------------------------------------
public:  
  void printInfo(void) {
    printf("IMGH::CornerDetector \n");
    smap.printInfo("  smap ");
    cmap.printInfo("  cmap ");
    std::cout << "  strength : low_limit=" << low_strength << "  max=" << max_strength << std::endl;
  }
  
  void showSignalStrength(void *image_buffer, pixel_t type) {
    IMGH::ImageConverter conv;
    conv.convertImage( smap.data, image_buffer, smap.w, smap.h, PIXEL_FLOAT, type);
  }
  void showCorners(Image *vimg) {
    vimg->setImage( cmap.w, cmap.h, PIXEL_RGB );
    showCorners( vimg->data, vimg->type );
  }
  void showCorners(void* image_buffer, pixel_t type) {
    // render a small red cross for each corner on the given image.
    int  idx, x, y, w = cmap.w, h = cmap.h, total = cmap.w * cmap.h;
    if (w <= 0 || h <= 0 || cmap.data == NULL) return;
    CornerPoint **mp = (CornerPoint**)(cmap.data);
    ImageDrawer<float> drawing( w, h, type, image_buffer );
    for (idx = 0; idx < total; idx++) {
      if (mp[idx] == NULL) continue;
      x = idx % w;  y = idx / w;
      if (y < 3 || y >= h-3) continue;
      drawing.drawCross( x, y, 2, '+', 255, 0, 0 );
    }
  }

};

  
}	// namespace IMGH


#endif	// IMGH_CORNER_DETECTOR_HPP

