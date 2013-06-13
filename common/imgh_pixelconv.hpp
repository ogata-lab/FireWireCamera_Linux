
//
// IMGH::PixelConverter 
//
// Jaeil Choi
// last modified in Dec, 2009
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
//   IMGH::PixelConverter	in imgh_pixelconv.hpp	for pixel color space conversion
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

#ifndef IMGH_PIXEL_CONVERTER_HPP
#define IMGH_PIXEL_CONVERTER_HPP

#include <iostream>

namespace IMGH {
  
#define PCONV_MINV(a) ( (a)[0]<(a)[1] ? ((a)[0]<(a)[2] ? (a)[0] : (a)[2]) : ((a)[1]<(a)[2] ? (a)[1] : (a)[2]) )
#define PCONV_MAXV(a) ( (a)[0]>(a)[1] ? ((a)[0]>(a)[2] ? (a)[0] : (a)[2]) : ((a)[1]>(a)[2] ? (a)[1] : (a)[2]) )
  
template <class T>
class PixelConverter
{
public:
  PixelConverter() {}
  ~PixelConverter() {}
  
  // -----------------------------------------------------------------
  // RGB, BGR, and GRAY
  //   RGB : [0,255],[0,255],[0,255]
  //   Gray: [0,255]
  // -----------------------------------------------------------------
public:
  inline void convRGB2Gray(T rgb[3], T gray[1]) {	// Craig's formula
    gray[0] = ((30*rgb[0] + 59*rgb[1] + 11*rgb[2])/100);
  }
  inline void convBGR2Gray(T bgr[3], T gray[1]) {	// Craig's formula
    gray[0] = ((30*bgr[2] + 59*bgr[1] + 11*bgr[0])/100);
  }
  inline void convRGB2Gray2(T rgb[3], T gray[1]) {
    gray[0] = ((6969*(int)rgb[0] + 23434*(int)rgb[1] + 2365*(int)rgb[2])/32768);
  }
  inline void convBGR2Gray2(T bgr[3], T gray[1]) {
    gray[0] = ((6969*(int)bgr[2] + 23434*(int)bgr[1] + 2365*(int)bgr[0])/32768);
  }
  
  // -----------------------------------------------------------------
  // YUV
  //   RGB : [0,255],[0,255],[0,255]
  //   YUV : [0,255],[0,255],[0,255]
  //   RGB to YUV Conversion:
  //     Y  =      (0.257 * R) + (0.504 * G) + (0.098 * B) + 16
  //     Cr = V =  (0.439 * R) - (0.368 * G) - (0.071 * B) + 128
  //     Cb = U = -(0.148 * R) - (0.291 * G) + (0.439 * B) + 128
  //   YUV to RGB Conversion:
  //     R = 1.164(Y - 16) + 1.596(V - 128)
  //     G = 1.164(Y - 16) - 0.813(V - 128) - 0.391(U - 128)
  //     B = 1.164(Y - 16)                  + 2.018(U - 128)
  // -----------------------------------------------------------------
#define PCONV_RGB_YUV(r,g,b,y,u,v) do { \
    y = (int)( 0.257 * (r) + 0.504 * (g) + 0.098 * (b)) +  16; \
    v = (int)( 0.439 * (r) - 0.368 * (g) - 0.071 * (b)) + 128; \
    u = (int)(-0.148 * (r) - 0.291 * (g) + 0.439 * (b)) + 128; } while(0)
#define PCONV_YUV_RGB(y,u,v,r,g,b) do { \
    r = (int)(1.164*((y)-16) + 1.596*((v)-128)); \
    g = (int)(1.164*((y)-16) - 0.813*((v)-128) - 0.391*((u) - 128)); \
    b = (int)(1.164*((y)-16)                   + 2.018*((u) - 128)); } while(0)
#define PCONV_MIN_MAX(a,b,c) do { \
    if (a < 0) a = 0;  if (a > 255) a = 255; \
    if (b < 0) b = 0;  if (b > 255) b = 255; \
    if (c < 0) c = 0;  if (c > 255) c = 255; } while(0)
  
public:
  void convRGB2YUV(T rgb[3], T yuv[3]) {
    int y, u, v;
    PCONV_RGB_YUV( rgb[0], rgb[1], rgb[2], y, u, v );  PCONV_MIN_MAX(y,u,v);
    yuv[0] = (T)y;  yuv[1] = (T)u;  yuv[2] = (T)v;
  }
  void convBGR2YUV(T bgr[3], T yuv[3]) {
    int y, u, v;
    PCONV_RGB_YUV( bgr[2], bgr[1], bgr[0], y, u, v );  PCONV_MIN_MAX(y,u,v);
    yuv[0] = (T)y;  yuv[1] = (T)u;  yuv[2] = (T)v;
  }
  
  void convYUV2RGB(T yuv[3], T rgb[3]) {
    int r, g, b;
    PCONV_YUV_RGB( yuv[0], yuv[1], yuv[2], r, g, b );  PCONV_MIN_MAX(r,g,b);
    rgb[0] = (T)r;  rgb[1] = (T)g;  rgb[2] = (T)b;
  }
  void convYUV2BGR(T yuv[3], T bgr[3]) {
    int r, g, b;
    PCONV_YUV_RGB( yuv[0], yuv[1], yuv[2], r, g, b );  PCONV_MIN_MAX(r,g,b);
    bgr[0] = (T)b;  bgr[1] = (T)g;  bgr[2] = (T)r;
  }
  void convUYV2RGB(T uyv[3], T rgb[3]) {
    int r, g, b;
    PCONV_YUV_RGB( uyv[1], uyv[0], uyv[2], r, g, b );  PCONV_MIN_MAX(r,g,b);
    rgb[0] = (T)r;  rgb[1] = (T)g;  rgb[2] = (T)b;
  }
  
  // -----------------------------------------------------------------
  // YUYV    (only for <unsigned char>)
  //   RGB  : [0,255],[0,255],[0,255]
  //   YUYV : [0,255],[0,255],[0,255]
  // -----------------------------------------------------------------
public:
  void convYUYV2RGBRGB(T yuyv[4], T rgb[6]) {
    int r, g, b;
    PCONV_YUV_RGB( yuyv[0], yuyv[1], yuyv[3], r, g, b );  PCONV_MIN_MAX(r,g,b);
    rgb[0] = (T)r;  rgb[1] = (T)g;  rgb[2] = (T)b;
    PCONV_YUV_RGB( yuyv[2], yuyv[1], yuyv[3], r, g, b );  PCONV_MIN_MAX(r,g,b);
    rgb[3] = (T)r;  rgb[4] = (T)g;  rgb[5] = (T)b;
  }
  void convYUYV2BGRBGR(T yuyv[4], T bgr[6]) {
    convYUYV2RGBRGB( yuyv, bgr );
    T a = bgr[0];  bgr[0] = bgr[2];  bgr[2] = a;
    T b = bgr[3];  bgr[3] = bgr[5];  bgr[5] = b;
  }
  void convYUYV2YUVYUV(T yuyv[4], T yuv[6]) {
    yuv[0] = yuyv[0];  yuv[1] = yuyv[1];  yuv[2] = yuyv[2];
    yuv[3] = yuyv[2];  yuv[4] = yuyv[1];  yuv[5] = yuyv[2];
  }
  void convUYVY2RGBRGB(T uyvy[4], T rgb[6]) {
    int r, g, b;
    PCONV_YUV_RGB( uyvy[1], uyvy[0], uyvy[2], r, g, b );  PCONV_MIN_MAX(r,g,b);
    rgb[0] = (T)r;  rgb[1] = (T)g;  rgb[2] = (T)b;
    PCONV_YUV_RGB( uyvy[3], uyvy[0], uyvy[2], r, g, b );  PCONV_MIN_MAX(r,g,b);
    rgb[3] = (T)r;  rgb[4] = (T)g;  rgb[5] = (T)b;
  }
  void convUYYVYY2RGB4(T uyyvyy[6], T rgb[12]) {
    int r, g, b;
    PCONV_YUV_RGB( uyyvyy[1], uyyvyy[0], uyyvyy[3], r, g, b );  PCONV_MIN_MAX(r,g,b);
    rgb[0] = (T)r;  rgb[1] = (T)g;  rgb[2] = (T)b;
    PCONV_YUV_RGB( uyyvyy[2], uyyvyy[0], uyyvyy[3], r, g, b );  PCONV_MIN_MAX(r,g,b);
    rgb[3] = (T)r;  rgb[4] = (T)g;  rgb[5] = (T)b;
    PCONV_YUV_RGB( uyyvyy[4], uyyvyy[0], uyyvyy[3], r, g, b );  PCONV_MIN_MAX(r,g,b);
    rgb[6] = (T)r;  rgb[7] = (T)g;  rgb[8] = (T)b;
    PCONV_YUV_RGB( uyyvyy[5], uyyvyy[0], uyyvyy[3], r, g, b );  PCONV_MIN_MAX(r,g,b);
    rgb[9] = (T)r;  rgb[10] = (T)g;  rgb[11] = (T)b;
  }
  
  // -----------------------------------------------------------------
  // HSV       (only for <float> and <double>)
  //   HSV : hue, saturation, value [0,359],[0,1],[0,1]
  //   RGB : normalized RGB values  [0,1],[0,1],[0,1]
  // -----------------------------------------------------------------
public:
  void convRGB2HSV(T rgb[3], T hsv[3]) {
    // Convert normalized RGB [0~1] to HSV color space.
    T minv=PCONV_MINV(rgb), maxv=PCONV_MAXV(rgb);
    // Hue [0 ~ 360]
    if      (maxv==minv)   hsv[0] = 0;
    else if (maxv==rgb[0]) hsv[0] = ((int)(60 * (rgb[1]-rgb[2])/(maxv-minv) + 360)) % 360;
    else if (maxv==rgb[1]) hsv[0] = ((int)(60 * (rgb[2]-rgb[0])/(maxv-minv) + 120));
    else if (maxv==rgb[2]) hsv[0] = ((int)(60 * (rgb[0]-rgb[1])/(maxv-minv) + 240));
    // Saturation [0 ~ 1]
    if (maxv==minv) hsv[1] = 0;
    else            hsv[1] = 1 - minv / maxv;
    // Value [0 ~ 1]
    hsv[2] = maxv;
  }
  void convHSV2RGB(T hsv[3], T rgb[3]) {
    // Convert HSV to normalized RGB [0~1]
    double h60 = (hsv[0]/60);
    double f = h60 - (int)h60;
    double p = hsv[2] * ( 1 - hsv[1] );
    double q = hsv[2] * ( 1 - f * hsv[1] );
    double t = hsv[2] * ( 1 - (1-f) * hsv[1] );
    switch ((int)(h60) % 6) {
    case 0:  rgb[0] = hsv[2];  rgb[1] = t;       rgb[2] = p;       break;
    case 1:  rgb[0] = q;       rgb[1] = hsv[2];  rgb[2] = p;       break;
    case 2:  rgb[0] = p;       rgb[1] = hsv[2];  rgb[2] = t;       break;
    case 3:  rgb[0] = p;       rgb[1] = q;       rgb[2] = hsv[2];  break;
    case 4:  rgb[0] = t;       rgb[1] = p;       rgb[2] = hsv[2];  break;
    case 5:  rgb[0] = hsv[2];  rgb[1] = p;       rgb[2] = q;       break;
    }
  }
  
  // -----------------------------------------------------------------
  // HSL       (only for <float> and <double>)
  //   HSL : hue, saturation, lightness [0,359],[0,1],[0,1]
  //   RGB : normalized RGB values      [0,1],[0,1],[0,1]
  // -----------------------------------------------------------------
public:
  void convRGB2HSL(T rgb[3], T hsl[3]) {
    // Convert normalized RGB [0~1] to HSL color space.
    T minv=PCONV_MINV(rgb), maxv=PCONV_MAXV(rgb);
    // Hue [0 ~ 360]
    if      (maxv==minv)   hsl[0] = (T)0;
    else if (maxv==rgb[0]) hsl[0] = (T)((int)(60 * (rgb[1]-rgb[2])/(maxv-minv) + 360)) % 360;
    else if (maxv==rgb[1]) hsl[0] = (T)((int)(60 * (rgb[2]-rgb[0])/(maxv-minv) + 120));
    else if (maxv==rgb[2]) hsl[0] = (T)((int)(60 * (rgb[0]-rgb[1])/(maxv-minv) + 240));
    // Saturation [0 ~ 1]
    double l = (maxv + minv) * 0.5;
    if  (maxv==minv) hsl[1] = (T)0;
    else if (l<=0.5) hsl[1] = (T)((maxv - minv) / (2*l));
    else if (l> 0.5) hsl[1] = (T)((maxv - minv) / (2 - 2*l));
    // lightness [0 ~ 1]
    hsl[2] = (T) l;
  }
  
  
};

  
}	// namespace IMGH


#endif	// IMGH_PIXEL_CONVERTER_HPP


