
//
// IMGH::OpenGLRenderer : rendering an IMGH::Image on OpenGL window
//
// Jaeil Choi
// last modified in Dec. 2006
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


#ifndef IMGH_OPENGL_RENDERER_HPP
#define IMGH_OPENGL_RENDERER_HPP

#include <iostream>
#if defined(__APPLE__) || defined(MACOSX)	// beginning of OS X version
#include "OpenGL/gl.h"
#include "OpenGL/glu.h"
#else
#include <GL/gl.h>
#include <GL/glu.h>
#endif
#include "imgh_common.hpp"

namespace IMGH {

#ifndef COLOR_STR2RGB
#define COLOR_STR2RGB(str, rgb) \
  do { \
    if (str == NULL || str[0] != '#') break; \
    if      (str[1] >= 'a') rgb[0]  = (str[1] - 'a' + 10) * 16; \
    else if (str[1] >= 'A') rgb[0]  = (str[1] - 'A' + 10) * 16; \
    else                    rgb[0]  = (str[1] - '0'     ) * 16; \
    if      (str[2] >= 'a') rgb[0] += (str[2] - 'a' + 10); \
    else if (str[2] >= 'A') rgb[0] += (str[2] - 'A' + 10); \
    else                    rgb[0] += (str[2] - '0'     ); \
    if      (str[3] >= 'a') rgb[1]  = (str[3] - 'a' + 10) * 16; \
    else if (str[3] >= 'A') rgb[1]  = (str[3] - 'A' + 10) * 16; \
    else                    rgb[1]  = (str[3] - '0'     ) * 16; \
    if      (str[4] >= 'a') rgb[1] += (str[4] - 'a' + 10); \
    else if (str[4] >= 'A') rgb[1] += (str[4] - 'A' + 10); \
    else                    rgb[1] += (str[4] - '0'     ); \
    if      (str[5] >= 'a') rgb[2]  = (str[5] - 'a' + 10) * 16; \
    else if (str[5] >= 'A') rgb[2]  = (str[5] - 'A' + 10) * 16; \
    else                    rgb[2]  = (str[5] - '0'     ) * 16; \
    if      (str[6] >= 'a') rgb[2] += (str[6] - 'a' + 10); \
    else if (str[6] >= 'A') rgb[2] += (str[6] - 'A' + 10); \
    else                    rgb[2] += (str[6] - '0'     ); \
  } while(0)
#endif

class OpenGLRenderer {

public:
  int alpha_mode;	// 0: normal  1: alpha_ignored  2: alpha_only

public:
  OpenGLRenderer() { alpha_mode = 0; }
  ~OpenGLRenderer() {}

  void render(IMGH::Image *img, char *bkg_color = NULL, bool draw_rect = true) {
    unsigned char rgb[4], bkg[3]={0,0,0};
    float a0, a1;
    int   x, y, y2;  
    // (0,0) is LL corner for OpenGL, and UL corner for Img
  
    if (img == NULL || img->w <= 0 || img->h <= 0) return;
    if (bkg_color) COLOR_STR2RGB( bkg_color, bkg );
    
    switch (img->type) {
    case IMGH::PIXEL_RGB:
      if (!draw_rect) {		// draw a pixel as a point
	glBegin(GL_POINTS);
	for (y = 0; y < img->h; y++) {
	  y2 = img->h - 1 - y;
	  for (x = 0; x < img->w; x++) {
	    img->getPixel(x, y, rgb);
	    glColor3ub(rgb[0], rgb[1], rgb[2]);
	    glVertex2i(x, y2);
	  }
	}
	glEnd();
      } else {			// draw a pixel as a rectangle
	glBegin(GL_QUADS);
	for (y = 0; y < img->h; y++) {
	  y2 = img->h - 1 - y;
	  for (x = 0; x < img->w; x++) {
	    img->getPixel(x, y, rgb);
	    glColor3ub(rgb[0], rgb[1], rgb[2]);
	    glVertex2i(x  , y2+1);	glVertex2i(x+1, y2+1);
	    glVertex2i(x+1, y2);	glVertex2i(x  , y2);
	  }
	}
	glEnd();
      }
      break;
    case IMGH::PIXEL_RGBA:
      if (!draw_rect) {		// a pixel as a point
	glBegin(GL_POINTS);
	for (y = 0; y < img->h; y++) {
	  y2 = img->h - 1 - y;
	  for (x = 0; x < img->w; x++) {
	    img->getPixel(x, y, rgb);
	    if (alpha_mode == 2) { // alpha_only (render the mask)
	      rgb[0] = rgb[1] = rgb[2] = rgb[3];
	    } else if (alpha_mode == 0 && rgb[3] < 255) {
	      if (bkg_color == NULL)
		bkg[0] = bkg[1] = bkg[2] = ( ((y/10 + x/10) % 2 == 0) ? 180 : 190 );
	      a0 = rgb[3]/255.0f;  a1 = (255-rgb[3])/255.0f;
	      rgb[0] = (unsigned char)( a0 * rgb[0] + a1 * bkg[0] );
	      rgb[1] = (unsigned char)( a0 * rgb[1] + a1 * bkg[1] );
	      rgb[2] = (unsigned char)( a0 * rgb[2] + a1 * bkg[2] );
	    }
	    glColor3ub(rgb[0], rgb[1], rgb[2]);
	    glVertex2i(x, y2);
	  }
	}
	glEnd();
      } else {			// draw a pixel as a rectangle
	glBegin(GL_QUADS);
	for (y = 0; y < img->h; y++) {
	  y2 = img->h - 1 - y;
	  for (x = 0; x < img->w; x++) {
	    img->getPixel(x, y, rgb);
	    if (alpha_mode == 2) { // alpha_only (render the mask)
	      rgb[0] = rgb[1] = rgb[2] = rgb[3];
	    } else if (alpha_mode == 0 && rgb[3] < 255) {
	      if (bkg_color == NULL)
		bkg[0] = bkg[1] = bkg[2] = ( ((y/10 + x/10) % 2 == 0) ? 180 : 190 );
	      a0 = rgb[3]/255.0f;  a1 = (255-rgb[3])/255.0f;
	      rgb[0] = (unsigned char)( a0 * rgb[0] + a1 * bkg[0] );
	      rgb[1] = (unsigned char)( a0 * rgb[1] + a1 * bkg[1] );
	      rgb[2] = (unsigned char)( a0 * rgb[2] + a1 * bkg[2] );
	    }
	    glColor3ub(rgb[0], rgb[1], rgb[2]);
	    glVertex2i(x  , y2+1);	glVertex2i(x+1, y2+1);
	    glVertex2i(x+1, y2);	glVertex2i(x  , y2);
	  }
	}
	glEnd();
    
      }
      break;
    default: break;
    }
    glFlush();
  }

};

}

#endif  // IMGH_OPENGL_RENDERER_HPP
