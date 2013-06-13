
//
// IMGH::ImageGradient 
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
//   IMGH::ImageFileIO 		in imgh_fileio.hpp	for file I/O (JPG,PNG,PGM, and more)
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
//   - creating gradient images (gx, gy, gz) from a given image
// This class requires:
//   - IMGH::Image		in 'imgh_common.hpp'
//   - IMGH::ImageFilter	in 'imgh_filter.hpp'
//

#ifndef IMGH_IMAGE_GRADIENT_HPP
#define IMGH_IMAGE_GRADIENT_HPP
#ifndef M_PI
#define M_PI		3.14159265358979323846
#endif

#include <iostream>
#include <cmath>
#include "imgh_common.hpp"
#include "imgh_converter.hpp"
#include "imgh_filter.hpp"

namespace IMGH {
  
  
class ImageGradient
{
public:
  IMGH::Image gx;	// gradient image for ( 1,0) (Ix)  type:PIXEL_FLOAT
  IMGH::Image gy;	// gradient image for ( 0,1) (Iy)  type:PIXEL_FLOAT
  IMGH::Image gz;	// gradient strength image         type:PIXEL_FLOAT
  float max_strength;
  
public:
  ImageGradient() {}
  ~ImageGradient() { clear(); }
  void clear(void) { gx.clear(); gy.clear(); gz.clear(); }
  
public:
  void getGradientImage(int w, int h, pixel_t type, void *udata, 
			int filter_size=0, bool strength_image=true) {
    Image img(w, h, type, udata);  getGradientImage( &img, filter_size );
  }
  void getGradientImage(Image *img, int filter_size=0, bool strength_image=true) {
    int		i, j, w = img->w, h = img->h;
    uchar_t	*pp, *cp, *np;
    gx.setImage( w, h, PIXEL_FLOAT, true );
    gy.setImage( w, h, PIXEL_FLOAT, true );
    float *xp = ((float*)gx.data)+w, *yp = ((float*)gy.data)+w, gxy[2], xpv[3], ypv[3];
    if (filter_size > 0) { 
      ImageFilter filter;
      filter.convoluteWithGaussian(img, NULL, filter_size);
    }
    switch (img->type) {
    case PIXEL_GRAY:
      for (i = 1; i < h-1; i++) {				// all the valid rows
	xp++; yp++;
	cp = img->data + i*w + 1;  pp = cp - w;  np = cp + w;	// all the valid columns
	if (filter_size > 0) {	// central difference
	  for (j = 1; j < w-1; j++, pp++, cp++, np++) {
	    *xp++ = (float)(cp[+1] - cp[-1]);  *yp++ = (float)(np[+0] - pp[+0]);
	  }
	} else {		// Sobel operator
	  for (j = 1; j < w-1; j++, pp++, cp++, np++) {
	    getSobelGray( pp, cp, np, gxy );
	    *xp++ = gxy[0]/8;  *yp++ = gxy[1]/8;
	  }
	}
	xp++; yp++;
      }
      break;
    case PIXEL_RGB: case PIXEL_BGR:
      for (i = 1; i < h-1; i++) {				// all the valid rows
	xp++; yp++;
	cp = img->data + 3*(i*w+1);  pp = cp-3*w;  np = cp+3*w;	// all the valid columns
	if (filter_size > 0) {	// central difference
	  for (j = 1; j < w-1; j++, pp+=3, cp+=3, np+=3) {
	    xpv[0] = (float)(cp[+3]-cp[-3]);  xpv[1] = (float)(cp[+4]-cp[-2]);  xpv[2] = (float)(cp[+5]-cp[-1]);
	    ypv[0] = (float)(np[+0]-pp[+0]);  ypv[1] = (float)(np[+1]-pp[+1]);  ypv[2] = (float)(np[+2]-pp[+2]);
	    *xp++ = IMGH3V_MAX( xpv[0], xpv[1], xpv[2] );
	    *yp++ = IMGH3V_MAX( ypv[0], ypv[1], ypv[2] );
	  }
	} else {		// Sobel operator
	  for (j = 1; j < w-1; j++, pp+=3, cp+=3, np+=3) {
	    getSobelRGB( pp, cp, np, gxy );
	    *xp++ = gxy[0]/8;  *yp++ = gxy[1]/8;
	  }
	}
	xp++; yp++;
      }
      break;
    case PIXEL_YUV:
      for (i = 1; i < h-1; i++) {				// all the valid rows
	xp++; yp++;
	cp = img->data + 3*(i*w+1);  pp = cp-3*w;  np = cp+3*w;	// all the valid columns
	if (filter_size > 0) {	// central difference
	  for (j = 1; j < w-1; j++, pp+=3, cp+=3, np+=3) {
	    xpv[1] = (float)(cp[+4]-cp[-2]);  xpv[2] = (float)(cp[+5]-cp[-1]);  //xpv[0] = cp[+3]-cp[-3];
	    ypv[1] = (float)(np[+1]-pp[+1]);  ypv[2] = (float)(np[+2]-pp[+2]);  //ypv[0] = np[+0]-pp[+0];
	    *xp++ = IMGH2V_MAX( xpv[1], xpv[2] );  //IMGH3V_MAX(xpv[0], xpv[1], xpv[2]);
	    *yp++ = IMGH2V_MAX( ypv[1], ypv[2] );  //IMGH3V_MAX(ypv[1], ypv[1], ypv[2]);
	  }
	} else {		// Sobel operator
	  for (j = 1; j < w-1; j++, pp+=3, cp+=3, np+=3) {
	    getSobelRGB( pp, cp, np, gxy );
	    *xp++ = gxy[0]/8;  *yp++ = gxy[1]/8;
	  }
	}
	xp++; yp++;
      }
      break;
    default:
      std::cerr << "Error (IMGH::ImageGradient::getGradientImage) invalid pixel format" << std::endl;
    } 
    getGradientImageBoundary(img);
    if (strength_image) getGradientStrength();
  }
  
  inline float getGrad3(uchar_t *cp,int a, int b) {
    uchar_t *ca = (cp + a*3), *cb = (cp + b*3);  float t[3];
    IMGH3V_SET( t, ca[0] - (float)cb[0], ca[1] - (float)cb[1], ca[2] - (float)cb[2] );
    return IMGH3V_MAX( t[0], t[1], t[2] );
  }
  void getGradientImageBoundary(Image *img) {
    if (!img || !img->sameSize(&gx) || !img->sameSize(&gy)) return;
    int		i, j, w = img->w, h = img->h;
    uchar_t	*cp = NULL;
    float *xp = NULL, *yp = NULL;
    switch (img->type) {
    case PIXEL_GRAY:
      // Y
      yp = (float*)(gy.data);  cp = (uchar_t*)img->data;
      for (i=0; i<w  ; i++) { *yp = ((float)*(cp+w) - (float)*(cp)) / 2;  yp++; cp++; }
      for (j=1; j<h-1; j++) { 
	*yp = ((float)*(cp+w) - (float)*(cp-w));  yp += w-1;  cp += w-1;
	*yp = ((float)*(cp+w) - (float)*(cp-w));  yp += 1;    cp += 1;
      }
      for (i=0; i<w  ; i++) { *yp = ((float)*(cp) - (float)*(cp-w)) / 2;  yp++; cp++; }
      // X
      xp = (float*)(gx.data);  cp = (uchar_t*)img->data;
      for (j=0; j<1; j++) {
	*xp = ((float)*(cp+1) - (float)*(cp)) / 2;  xp += 1;  cp += 1;
	for (i=1;i<w-1;i++) { *xp = ((float)*(cp+1) - (float)*(cp-1));  xp += 1;  cp += 1; }
	*xp = ((float)*(cp) - (float)*(cp-1)) / 2;  xp += 1;  cp += 1;
      }
      for (  ; j<h-1; j++) {
	*xp = ((float)*(cp+1) - (float)*(cp)) / 2;  xp += w-1;  cp += w-1;
	*xp = ((float)*(cp) - (float)*(cp-1)) / 2;  xp += 1;  cp += 1;
      }
      for (  ; j<h; j++) {
	*xp = ((float)*(cp+1) - (float)*(cp)) / 2;  xp += 1;  cp += 1;
	for (i=1;i<w-1;i++) { *xp = ((float)*(cp+1) - (float)*(cp-1));  xp += 1;  cp += 1; }
	*xp = ((float)*(cp) - (float)*(cp-1)) / 2;  xp += 1;  cp += 1;
      }
      break;
    case PIXEL_RGB: case PIXEL_BGR:
      // Y
      yp = (float*)(gy.data);  cp = (uchar_t*)img->data;
      for (i=0; i<w  ; i++) { *yp = getGrad3(cp, w, 0)/2;  yp++; cp+=3; }
      for (j=1; j<h-1; j++) { 
	*yp = getGrad3(cp, +w, -w);  yp+=(w-1);  cp += 3*(w-1);
	*yp = getGrad3(cp, +w, -w);  yp++;    cp += 3;
      }
      for (i=0; i<w  ; i++) { *yp = getGrad3(cp, 0, -w)/2;  yp++; cp+=3; }
      // X
      xp = (float*)(gx.data);  cp = (uchar_t*)img->data;
      for (j=0; j<1; j++) {
	*xp = getGrad3(cp, +1, 0)/2;  xp++;  cp += 3;
	for (i=1;i<w-1;i++) { *xp = getGrad3(cp, +1, -1);  xp++;  cp += 3; }
	*xp = getGrad3(cp, 0, -1)/2;  xp++;  cp += 3;
      }
      for (  ; j<h-1; j++) {
	*xp = getGrad3(cp, +1, 0)/2;  xp+=(w-1);  cp += 3*(w-1);
	*xp = getGrad3(cp, 0, -1)/2;  xp++;  cp += 3;
      }
      for (  ; j<h; j++) {
	*xp = getGrad3(cp, +1, 0)/2;  xp++;  cp += 3;
	for (i=1;i<w-1;i++) { *xp = getGrad3(cp, +1, -1);  xp++;  cp += 3; }
	*xp = getGrad3(cp, 0, -1)/2;  xp++;  cp += 3;
      }
      break;
    default:  break;
    } 
  }
  
  void getGradientStrength(int type=0) {
    // Calculate gradient strength for each pixel in 'gz',
    //   assumming 'gx' and 'gy' is ready (by getGradient()).
    if (gx.w == 0 || gx.h == 0 || !gx.sameSize(&gy)) return;
    int  i, w = gx.w, h = gx.h, total = gx.w * gx.h;
    int  stt = w + 1, end = total - w - 1;
    if (!gz.sameSize(w,h)) gz.setImage( w, h, PIXEL_FLOAT );
    else                   gz.clearImage();
    float *xp, *yp, *zp;
    float xv, yv, v1, v2;
    max_strength = 0;
    switch (type) {
    case 0:	// L2 norm of (gx, gy) 
      xp = ((float*)gx.data);  yp = ((float*)gy.data);  zp = ((float*)gz.data);
      for (i=0; i<total; i++, xp++, yp++, zp++) {
	*zp = (float)sqrt((*xp)*(*xp)+(*yp)*(*yp));
	if (*zp > max_strength) max_strength = *zp;
      }
      break;
    case 1:	// L1 norm of (gx, gy)   (the gain is not strong enough)
      xp = ((float*)gx.data);  yp = ((float*)gy.data);  zp = ((float*)gz.data);
      for (i=0; i<total; i++, xp++, yp++, zp++) {
	*zp = (float)(fabs(*xp)+fabs(*yp));  
	if (*zp > max_strength) max_strength = *zp;
      }
      break;
    case 2:	// ratio of orthogonal directions (gx vs. gy)
      xp = ((float*)gx.data)+stt;  yp = ((float*)gy.data)+stt;  zp = ((float*)gz.data)+stt;
      for (i = stt; i < end; i++, xp++, yp++, zp++) {
	xv = fabs(*xp);  yv = fabs(*yp);
	if (xv > yv) {
	  v1 = fabs(*(yp-1));  v2 = fabs(*(yp+1));  yv = IMGH3V_MIN( yv, v1, v2 );
	  *zp = ((yv > 2) ? (xv/yv) : xv/2);
	} else {
	  v1 = fabs(*(xp-w));  v2 = fabs(*(xp+w));  xv = IMGH3V_MIN( xv, v1, v2 );
	  *zp = ((xv > 2) ? (yv/xv) : yv/2);
	}
	if (*zp > max_strength) max_strength = *zp;
      }
      break;
    default: std::cerr << "Error (IMGH::ImageGradient::getGradientStrengthImage) invalid type" << std::endl;
    }
  }
  
  void getGradientCovariance(void) {
    // Gradient covariance for corner detection.
    // Calculate gradient covariance for each pixel in 'gx' 'gy' and 'gz',
    //   assumming 'gx' and 'gy' is ready (by getGradient()).
    // Note that 'gx', 'gy' and 'gz' is updated with         [ Ix^2 IxIy ] 
    //   gradient covariance values (Ix*Ix, Iy*Iy, Ix*Iy)    [ IxIy Iy^2 ]
    int   i, w = gx.w, h = gx.h, total = gx.w * gy.h;
    gz.setImage( w, h, PIXEL_FLOAT );
    float *xp=(float *)(gx.data), *yp=(float *)(gy.data), *zp=(float *)(gz.data);
    for (i = 0; i < total; i++, xp++, yp++, zp++) {
      *zp = (*xp) * (*yp);	// Ix*Iy
      *xp = (*xp) * (*xp);	// Ix^2
      *yp = (*yp) * (*yp);	// Iy^2
    }
  }
  
  void getGradientVector(int idx, float gv[2], bool normalize=false) {
    IMGH2V_SET( gv, gx.getFloat( idx ), gy.getFloat( idx ) );
    if (normalize) { 
      float len = (gz.w == gx.w ? gz.getFloat(idx) : IMGH2V_LEN( gv ));
      if (len>0) IMGH2V_DIV( gv, len ); 
    }
  }
  
private:
  char getGradientAngle( float gx, float gy ) {
    // get the discretized gradient angle [0 - 35]  with 0 for (1,0) in CCW.
    int angle;
    if (gx == 0) angle = (gx >= 0 ? 0 : 180);
    else         angle = (int)( atanf( gy / gx ) * 180 / M_PI );
    return ( angle >= 0 ? angle/10 : 35 + angle/10 );
  }
  void getSobelGray(uchar_t *pp, uchar_t *cp, uchar_t *np, float gxy[2]) {
    gxy[0] = 2.0f * ((int)cp[+1] - (int)cp[-1]) + (int)pp[+1] - (int)pp[-1] + (int)np[+1] - (int)np[-1];
    gxy[1] = 2.0f * ((int)np[+0] - (int)pp[+0]) + (int)np[-1] - (int)pp[-1] + (int)np[+1] - (int)pp[+1];
  }
  void getSobelRGB(uchar_t *pp, uchar_t *cp, uchar_t *np, float gxy[2]) {
    float tv[3][2];
    // Red channel
    tv[0][0] = 2.0f * ((int)cp[+3] - (int)cp[-3]) + (int)pp[+3] - (int)pp[-3] + (int)np[+3] - (int)np[-3];
    tv[0][1] = 2.0f * ((int)np[+0] - (int)pp[+0]) + (int)np[-3] - (int)pp[-3] + (int)np[+3] - (int)pp[+3];
    // Green channel
    tv[1][0] = 2.0f * ((int)cp[+4] - (int)cp[-2]) + (int)pp[+4] - (int)pp[-2] + (int)np[+4] - (int)np[-2];
    tv[1][1] = 2.0f * ((int)np[+1] - (int)pp[+1]) + (int)np[-2] - (int)pp[-2] + (int)np[+4] - (int)pp[+4];
    // Blue channel
    tv[2][0] = 2.0f * ((int)cp[+5] - (int)cp[-1]) + (int)pp[+5] - (int)pp[-1] + (int)np[+5] - (int)np[-1];
    tv[2][1] = 2.0f * ((int)np[+2] - (int)pp[+2]) + (int)np[-1] - (int)pp[-1] + (int)np[+5] - (int)pp[+5];
    // get maximum
    gxy[0] = tv[0][0];  gxy[1] = tv[0][1];
    if (fabs(tv[1][0]) > fabs(gxy[0])) gxy[0] = tv[1][0];  if (fabs(tv[1][1]) > fabs(gxy[1])) gxy[1] = tv[1][1];
    if (fabs(tv[2][0]) > fabs(gxy[0])) gxy[0] = tv[2][0];  if (fabs(tv[2][1]) > fabs(gxy[1])) gxy[1] = tv[2][1];
  }
  
public:
  void printInfo(void) {
    printf("IMGH::ImageGradient \n");
    gx.printInfo("  gx  ");
    gy.printInfo("  gy  ");
    gz.printInfo("  gz  ");
  }
  
  void showGradient(void* image_buffer, pixel_t type, char name='x') {
    IMGH::ImageConverter conv;
    switch (name) {
    case 'x': conv.convertImage(gx.data, image_buffer, gx.w, gx.h, PIXEL_FLOAT, type);  break;
    case 'y': conv.convertImage(gy.data, image_buffer, gy.w, gy.h, PIXEL_FLOAT, type);  break;
    case 'z': conv.convertImage(gz.data, image_buffer, gz.w, gz.h, PIXEL_FLOAT, type);  break;
    }
  }

};

  
}	// namespace IMGH


#endif	// IMGH_IMAGE_GRADIENT_HPP


