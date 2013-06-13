
//
// IMGH::IntegralImage 
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

#ifndef IMGH_INTEGRAL_IMAGE_HPP
#define IMGH_INTEGRAL_IMAGE_HPP

#include <iostream>
#include "imgh_common.hpp"

namespace IMGH {
  
  
class IntegralImage
{
public:
  int	w, h;	// image width and height
  int	*data[3];
public:
  IntegralImage() { data[0]=data[1]=data[2]=NULL; }
  IntegralImage(Image *img) { data[0]=data[1]=data[2]=NULL; createIntegralImage(img); }
  IntegralImage(int w, int h, pixel_t type, uchar_t *src) { data[0]=data[1]=data[2]=NULL; createIntegralImage(w,h,type,src); }
  ~IntegralImage() {}
public:
  inline void clear(void) { 
    w = h = 0;
    if (data[0]) free(data[0]);  data[0] = NULL;
    if (data[1]) free(data[1]);  data[1] = NULL;
    if (data[2]) free(data[2]);  data[2] = NULL;
  }
  inline bool isReady(void) { return (data[0]!=NULL); }
  
  // -----------------------------------------------------------------
  // Integral Image
  // -----------------------------------------------------------------

  bool createIntegralImage(Image *img, bool need_to_clear=true) { 
    return create(img->w, img->h, img->type, img->data, need_to_clear); 
  }
  bool createIntegralImage(int w, int h, pixel_t type, uchar_t *src, bool need_to_clear=true) {
    // create an integral image of the input image 'src'
    if (w <= 0 || h <= 0 || src == NULL) return false;
    if (need_to_clear) clear();
    int i, x, y, total = w * h, *da=NULL, *db=NULL, *dc=NULL;
    switch(type) {
    case PIXEL_GRAY:  case PIXEL_GRAYA: 
      this->w = w;  this->h = h;
      if (data[0]) da = data[0];			// allocate buffer
      else         da = data[0] = (int*)malloc( total * sizeof(int) );
      if (type == PIXEL_GRAY) {
	for (i=0; i<total; i++, da++,db++,dc++) {	// vertical sum
	  *da = (i<w ? *src : *src + da[-w]);  src++;
	}
      } else {
	for (i=0; i<total; i++, da++,db++,dc++) {	// vertical sum
	  *da = (i<w ? *src : *src + da[-w]);  src += 2;
	}
      }
      for (y = 0; y < h; y++) {				// horizontal sum
	for (x=1, da=this->data[0]+y*w+1; x<w; x++, da++)  *da += da[-1];
      }
      break;
    case PIXEL_RGB:  case PIXEL_RGBA:  case PIXEL_BGR:  case PIXEL_YUV:  
      this->w = w;  this->h = h;
      if (data[0]) da = data[0];			// allocate buffer
      else         da = data[0] = (int*)malloc( total * sizeof(int) );
      if (data[1]) db = data[1];
      else         db = data[1] = (int*)malloc( total * sizeof(int) );
      if (data[2]) dc = data[2];
      else         dc = data[2] = (int*)malloc( total * sizeof(int) );
      if (type != PIXEL_RGBA) {
	for (i=0; i<total; i++, da++,db++,dc++) {	// vertical sum
	  *da = (i<w ? *src : *src + da[-w]);  src++;
	  *db = (i<w ? *src : *src + da[-w]);  src++;
	  *dc = (i<w ? *src : *src + da[-w]);  src++;
	}
      } else {
	for (i=0; i<total; i++, da++,db++,dc++) {	// vertical sum
	  *da = (i<w ? *src : *src + da[-w]);  src++;
	  *db = (i<w ? *src : *src + da[-w]);  src++;
	  *dc = (i<w ? *src : *src + da[-w]);  src += 2;
	}
      }
      for (y = 0; y < h; y++) {				// horizontal sum
	for (x=1, da=this->data[0]+y*w+1; x<w; x++, da++)  *da += da[-1];
	for (x=1, db=this->data[1]+y*w+1; x<w; x++, db++)  *db += db[-1];
	for (x=1, dc=this->data[2]+y*w+1; x<w; x++, dc++)  *dc += dc[-1];
      }
      break;
    default:
      std::cerr << "Error (IMGH::IntegralImage::create): Unsupported pixel type" << std::endl;
      break;
    }
    return (this->data != NULL);
  }
  
  // rectangular sum from (x0,y0) to (x1, y1), inclusive, assumming   x1 > x0  &&  y1 > y0
  inline int getSum(int x,  int y) { 
    // Get the integral value of the given pixel (ONLY FOR SINGLE-CHANNEL IMAGE!)
    // Note this function do NOT check, for speed, whether or not the data[] is ready.
    //   Use appropiately, otherwise it will crash.
    //// if (x<0 || y<0) return 0; //// commented out for speed
    return data[0][ (w * y + x) ];
  }
  inline int getSum(int x, int y, int w, int h) { 
    // Get the sum of pixel values in the window (ONLY FOR SINGLE-CHANNEL IMAGE!)
    // Note this function do NOT check, for speed, whether or not the data[] is ready.
    //   Use appropiately, otherwise it will crash.
    int x0 = x-1, y0 = y-1, x1 = x+w-1, y1 = y+h-1;
    return ( getSum(x1,y1) - getSum(x1,y0) - getSum(x0,y1) + getSum(x0,y0) ); 
  }
  void getSum(int x, int y, int val[3]) {
    // Get the integral value of the given pixel (ONLY FOR THREE-CHANNEL IMAGE!)
    // Note this function do NOT check, for speed, whether or not the data[] is ready.
    //   Use appropiately with valid window range (x,y>0), otherwise it will crash.
    //// if (x<0 || y<0) { val[0] = val[1] = val[2] = 0; return; } //// commented out for speed
    int idx = y * w + x;
    val[0] = data[0][idx];  val[1] = data[1][idx];  val[2] = data[2][idx];
  }
  void getSum(int x, int y, int w, int h, int val[3]) { 
    // Get the sum of pixel values in the window (ONLY FOR THREE-CHANNEL IMAGE!)
    // Note this function do NOT check, for speed, whether or not the data[] is ready.
    //   Use appropiately with valid window range (x,y>0), otherwise it will crash.
    int x0 = x-1, y0 = y-1, x1 = x+w-1, y1 = y+h-1;
    x0 -= 1;  y0 -= 1;
    getSum( x0, y0, v00 );   getSum( x1, y0, v10 );
    getSum( x0, y1, v01 );   getSum( x1, y1, v11 );
    val[0] = v11[0] - v10[0] - v01[0] + v00[0];
    val[1] = v11[1] - v10[1] - v01[1] + v00[1];
    val[2] = v11[2] - v10[2] - v01[2] + v00[2];
  }
  
  void printInfo(const char *cmmt=NULL) {
    int nch = 0;
    for (i = 0; i < 3; i++) if (data[i]!=NULL) nch++;
    printf("%s (%3dx%3d) with %d channels of integrated pixel values\n", (cmmt ? cmmt:"IntegralImage"), nch);
  }
  
};

  
}	// namespace IMGH


#endif	// IMGH_INTEGRAL_IMAGE_HPP
