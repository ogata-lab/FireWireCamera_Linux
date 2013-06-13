
//
// ImageIconEditor 
//
// Jaeil Choi
// last modified in July, 2008
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
//

#ifndef IMGH_IMAGE_ICON_EDITOR_HPP
#define IMGH_IMAGE_ICON_EDITOR_HPP

#include <iostream>
#include "imgh_common.hpp"
#include "imgh_editor.hpp"
#include "imgh_drawer.hpp"
#include "imgh_pixelconv.hpp"

namespace IMGH {
  
#define IMGH3V_INTERPOLATE(sp, ep, r, dp)	\
  do { dp[0] = (uchar_t)(sp[0] * (1-r) + ep[0] * r); \
       dp[1] = (uchar_t)(sp[1] * (1-r) + ep[1] * r); \
       dp[2] = (uchar_t)(sp[2] * (1-r) + ep[2] * r); } while(0)
#define IMGH4V_INTERPOLATE(sp, ep, r, dp)	\
  do { dp[0] = (uchar_t)(sp[0] * (1-r) + ep[0] * r); \
       dp[1] = (uchar_t)(sp[1] * (1-r) + ep[1] * r); \
       dp[2] = (uchar_t)(sp[2] * (1-r) + ep[2] * r); \
       dp[3] = (uchar_t)(sp[3] * (1-r) + ep[3] * r); } while(0)
  
class ImageIconEditor
{
private:
  ImageEditor	edt;
  Image		dst;
  PixelConverter<unsigned char> pc;
  
public:
  ImageIconEditor() {}
  ~ImageIconEditor() {}
  
public:  
  bool processImage(IMGH::Image *img, char *cmd, char *msg=NULL) {
    // Process the editing command.
    if (!img || !cmd) return false;
    int    n, idx, x, y, w, h, c, r, g, b;
    unsigned char rgb[3];
    float darker;  char shades[40];  shades[0]='\0';
    if          (strncmp (cmd, "copyToAlpha:", 12)==0) {
      n = sscanf(cmd, "copyToAlpha:%d", &idx);
      if (n != 1) { if (msg) sprintf(msg, "invalid args for copyToAlpha:idx"); return false; }
      setAlphaWithRGBChannel(img, idx);
    } else if   (strncmp (cmd, "fixRect", 7)==0) {
      n = sscanf(cmd, "fixRect:%d:%d:%d:%d", &x, &y, &w, &h);
      if (n != 4) { if (msg) sprintf(msg, "invalid args for fixRect:x:y:w:h"); return false; }
      return fixRect( img, x, y, w, h );
    } else if (strncmp (cmd, "reduceButton", 12)==0) {
      n = sscanf(cmd, "reduceButton:%d:%d:%d", &w, &h, &c);
      if (n != 3) { if (msg) sprintf(msg, "invalid args for reduceButton:w:h:c"); return false; }
      return reduceRoundedButton( img, w, h, c );
    } else if (strncmp (cmd, "changeColor", 11)==0) {
      n = sscanf(cmd, "changeColor:%d:%d:%d:%f", &r, &g, &b, &darker);
      if (n != 4) { if (msg) sprintf(msg, "invalid args for changeColor:r:g:b:s"); return false; }
      IMGH3V_SET( rgb, r, g, b );
      return changeColor( img, rgb, darker );
    } else if (strcmp(cmd, "drawPrevButton")==0) {
      return drawButton( img, "Prev", 255, 255, 255 );
    } else if (strcmp(cmd, "drawNextButton")==0) {
      return drawButton( img, "Next", 255, 255, 255 );
    } else if (strcmp(cmd, "drawGridButton")==0) {
      return drawButton( img, "Grid", 255, 255, 255 );
    } else if (strcmp(cmd, "drawGridButton")==0) {
      return drawButton( img, "Grid", 255, 255, 255 );
    } else if (strcmp(cmd, "drawBooksButton")==0) {
      return drawButton( img, "Books", 255, 255, 255 );
    } else if (strcmp(cmd, "drawChapButton")==0) {
      return drawButton( img, "Chap", 255, 255, 255 );
    } else if (strncmp (cmd, "drawRect", 7)==0) {
      n = sscanf(cmd, "imgRect(%d:%d)(%d:%d:%d)%s", &w, &h, &r, &g, &b, shades);
      if (n < 5) { if (msg) sprintf(msg, "invalid args for imgRect(w:h)(r:g:b)"); return false; }
      return drawControlBoundary( img, w, h, r, g, b, (n>5 ? shades:(char*)"DD") );  // "UDD"
    } else return edt.processImage( img, cmd, msg );
    return true;
  }
  
  // -----------------------------------------------------------------
  // 
  // -----------------------------------------------------------------
  
public:
  bool setAlphaWithRGBChannel(IMGH::Image *img, int channel) {
    if (!img || img->type != PIXEL_RGBA) return false;
    int      i, total=img->w*img->h;
    uchar_t *p = img->data;
    for (i=0; i<total; i++, p+=4) p[3] = p[channel];
    return true;
  }
  
//   bool copyChannelFromImage(Image *img, Image *aimg, char c='A') {
//     if (!img || !aimg || img->type != PIXEL_RGBA) return false;
//     aimg->setImage( img->w, img->h, PIXEL_RGB );
//     uchar_t *src = img->data, *dst = aimg->data;
//     int   i, total = img->w * img->h;
//     if      (c=='A') for (i=0; i<total; i++, src+=4, dst+=3) dst[0] = dst[1] = dst[2] = src[3];
//     else if (c=='R') for (i=0; i<total; i++, src+=4, dst+=3) dst[0] = dst[1] = dst[2] = src[0];
//     else if (c=='G') for (i=0; i<total; i++, src+=4, dst+=3) dst[0] = dst[1] = dst[2] = src[1];
//     else if (c=='B') for (i=0; i<total; i++, src+=4, dst+=3) dst[0] = dst[1] = dst[2] = src[2];
//     else             for (i=0; i<total; i++, src+=4, dst+=3) IMGH3V_COPY( dst, src );
//     return true;
//   }
  
  bool fixRect(IMGH::Image *img, int x, int y, int w, int h) {
    if (img->type != PIXEL_RGBA) return false;
    if (!adjustRect( img, x, y, w, h, 1 )) return false;
    for (int j = y; j < y+h; j++) {
      uchar_t *pl = img->pixel( (x-1), j );
      uchar_t *pr = img->pixel( (x+w), j );
      for (int i = 0; i < w; i++) {
	uchar_t *pc = img->pixel( (x+i), j );
	float r = (i+0.5)/w;
	IMGH3V_INTERPOLATE( pl, pr, r, pc );
      }
    }
    return true;
  }
  
  bool reduceRoundedButton(IMGH::Image *img, int w, int h, int corner) {
    if (!img || img->type != PIXEL_RGBA) return false;
    if (w > img->w || h > img->h || w < 2*corner || h < 2*corner) return false;
    dst.setImage( w, h, PIXEL_RGBA );
    dst.clearImage( 0);
    uchar_t *sr, *dr, *r0=NULL, *r1=NULL, *sp, *ep, *dp;
    int     i, j, k, wdiff=img->w-w, wfill=w-2*corner, hfill=h-2*corner;
    for (j=0; j<corner; j++) {
      for (k=0; k<2; k++) {
	if (k<1) { sr = img->data + 4*img->w*j; dr = dst.data + 4*dst.w*j; }
	else { sr = img->data + 4*img->w*(img->h-1-j); dr = dst.data + 4*dst.w*(dst.h-1-j); }
	for (i=0; i<w; i++) {
	  if        (i<corner)  { sp = sr + 4*i; dp = dr + 4*i; IMGH4V_COPY( dp, sp ); }
	  else if (i>=w-corner) { sp = sr + 4*(i+wdiff); dp = dr + 4*i; IMGH4V_COPY( dp, sp ); }
	  else {
	    sp = sr + 4*(corner-1);  ep = sr + 4*(corner+wfill+wdiff);  dp = dr + 4*i;
	    float r = (i - corner + 0.5f) / wfill;
	    IMGH4V_INTERPOLATE( sp, ep, r, dp );
	  }
	}
	if (k<1) r0 = dr;
	else     r1 = dr;
      }
    }
    for (j=0; j < hfill; j++) {
      float r = (j + 0.5f) / hfill;
      dr = dst.data + 4 * dst.w * (corner + j);
      for (i=0; i < w; i++) {
	sp = r0 + 4 * i;  ep = r1 + 4 * i;  dp = dr + 4 * i;
	IMGH4V_INTERPOLATE( sp, ep, r, dp );
      }
    }
    img->swapImage( &dst );
    return true;
  }
  
  bool changeColor(IMGH::Image *img, unsigned char rgb[3], float darker=1.0) {
    if (!img || img->type != PIXEL_RGBA) return false;
    int      i, total=img->w*img->h;
    uchar_t *p = img->data, yuv[3], yuv_set[3];
    pc.convRGB2YUV( rgb, yuv_set );
    for (i=0; i<total; i++, p+=4) {
      if (p[3] < 0.95) continue;
      pc.convRGB2YUV( p, yuv );
      // change luminance
      if (darker!=1.0) yuv[0] = (uchar_t)( pow(yuv[0]/255.0f, darker) * 255.0f );
      // change hue
      yuv[1] = yuv_set[1];  yuv[2] = yuv_set[2];
      pc.convYUV2RGB( yuv, p );
    }
    return true;
  }
  
private:
  bool adjustRect(IMGH::Image *img, int &x, int &y, int &w, int &h, int m=0) {
    if (x < m) x = m;  if (y < m) y = m;
    if (x+w+m >= img->w) w = img->w - 1 - x - m;  if (w < 1) return false;
    if (y+h+m >= img->h) h = img->h - 1 - y - m;  if (h < 1) return false;
    //printf(" x=%d y=%d w=%d h=%d\n", x, y, w, h);
    return true;
  }
  
  // -----------------------------------------------------------------
  // 
  // -----------------------------------------------------------------
  
public:
  bool drawButton(IMGH::Image *img, char *type, int r, int g, int b, int a=-2) {
    int ww=img->w, hh=img->h;
    float pa[2], pb[2], pc[2];
    IMGH::ImageDrawer<float> drw(img);
    if        (strcmp(type, "Prev")==0) {
      IMGH2V_SET( pa, ww*0.25, hh*0.50 );
      IMGH2V_SET( pb, ww*0.75, hh*0.75 );
      IMGH2V_SET( pc, ww*0.75, hh*0.25 );
      drw.drawFilledTriangle( pa, pb, pc, r, g, b, a );
    } else if (strcmp(type, "Next")==0) {
      IMGH2V_SET( pa, ww*0.75, hh*0.50 );
      IMGH2V_SET( pb, ww*0.25, hh*0.25 );
      IMGH2V_SET( pc, ww*0.25, hh*0.75 );
      drw.drawFilledTriangle( pa, pb, pc, r, g, b, a );
    } else if (strcmp(type, "Grid")==0) {
      int i, j, mw = ww/4, mh = hh/4, sw = ww/15, sh = hh/15;
      for (i=0; i<3; i++)
	for (j=0; j<3; j++)
	  drw.drawFilledRectangle( mw+i*sw*3, mh+j*sw*3, sw*2, sh*2, r, g, b, a );
    } else if (strcmp(type, "Books")==0) {
      int mw = ww/4-1, mh = hh/4, bw = ww/2/4, bh = hh/2+1;
      drw.drawFilledRectangle( mw,   mh+0, bw,   bh-0, r, g, b, a );  
      drw.drawFilledRectangle( mw+1, mh+2, 1,  8, r/2, g/2, b/2, a/2 );  mw += bw+1;
      drw.drawFilledRectangle( mw,   mh+1, bw+1, bh-1, r, g, b, a ); 
      drw.drawFilledRectangle( mw+1, mh+4, 2,  4, r/2, g/2, b/2, a/2 );  mw += bw+2;
      drw.drawFilledRectangle( mw, mh+0, bw, bh-0, r, g, b, a );
      drw.drawFilledRectangle( mw+1, mh+1, 1, 11, r/2, g/2, b/2, a/2 );   mw += bw+1;
      drw.drawFilledRectangle( mw, mh+2, bw+1, bh-2, r, g, b, a );   
      drw.drawFilledRectangle( mw+1, mh+4, 2, 6, r/2, g/2, b/2, a/2 );   
    } else if (strcmp(type, "Chap")==0) {
      int mw=ww/4, mh = hh/4, bw = ww/8, bh = hh/2+1, cw = ww/2;
      drw.drawFilledRectangle( cw - mw +  0,   mh, bw, bh, r, g, b, a );
      drw.drawFilledRectangle( cw + mw - bw+1, mh, bw, bh, r, g, b, a );
      drw.drawFilledRectangle( cw - bw/2, cw-2*bw, bw, bw, r, g, b, a );
      drw.drawFilledRectangle( cw - bw/2, cw+1*bw, bw, bw, r, g, b, a );
    } else return false;
    return true;
  }
  
  bool drawControlBoundary(IMGH::Image *img, int w, int h, int r, int g, int b, char *shades) {
    // create an image of rectangular control with shaded boundary
    img->setImage( w, h, PIXEL_RGBA );
    img->clearImage( r, g, b, 255 );
    int i, n=strlen(shades);
    uchar_t color[3], dark[3], bright[3], white[3]={255,255,255};
    IMGH3V_SET( color, (uchar_t)r, (uchar_t)g, (uchar_t)b );	// input color
    IMGH3V_COPY( dark, color );  IMGH3V_DIV( dark, 2 );		// dark   boundary
    IMGH3V_INTERPOLATE( color, white, 0.67, bright );		// bright boundary
    uchar_t *ct=color, *cl=color, *cr=color, *cb=color;
    IMGH::ImageDrawer<float> drw(img);
    for (i = 0; i < n; i++) {
      if      (shades[i]=='U') { ct = cl = bright; cr = cb = dark;   }
      else if (shades[i]=='D') { ct = cl = dark;   cr = cb = bright; }
      drw.drawLine( i, i, w-1-i, i, ct[0], ct[1], ct[2] );		// top
      drw.drawLine( i, i, i, h-1-i, cl[0], cl[1], cl[2] );		// left
      drw.drawLine( w-1-i, i, w-1-i, h-1-i, cr[0], cr[1], cr[2] );	// right
      drw.drawLine( i, h-1-i, w-1-i, h-1-i, cr[0], cr[1], cr[2] );	// bottom
    }
    return true;
  }
};

  
}	// namespace IMGH


#endif	// IMGH_IMAGE_ICON_EDITOR_HPP


