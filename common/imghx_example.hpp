
//
// ImageExample 
//
// Jaeil Choi
// last modified in Dec, 2006
//

#ifndef IMGH_IMAGE_EXAMPLE_HPP
#define IMGH_IMAGE_EXAMPLE_HPP

#include <iostream>
#include <cmath>
#include "imgh_common.hpp"
#include "util_color_mapper.hpp"

namespace IMGH {
  
  
class ImageExample
{
public:
  ImageExample() {}
  ~ImageExample() {}
  
public:
  void drawRandomPixels(IMGH::Image *img) {
    int i, total = img->w * img->h;
    switch (img->type) {
    case PIXEL_GRAY: for (i=0; i<total; i++) img->setChar(i, (char)(rand()%256));  break;
    case PIXEL_RGB:  case PIXEL_RGBA:  case PIXEL_BGR:  case PIXEL_YUV:  
      for (i=0; i<total; i++) img->setRGB(i, rand()%256, rand()%256, rand()%256);  break;
    default: break;
    }
  }
  
  void drawMandelbrot(IMGH::Image *img, double cx, double cy, double scale, int max_iteration=0) {
    unsigned char *pp = img->data;
    bool gray = (img->type==IMGH::PIXEL_GRAY || img->type==IMGH::PIXEL_GRAYA);
    UTIL::ColorMapper<unsigned char> cmap("b9r{}g[]psmkoyGt",0,1000,true);
    unsigned char rgb[3];
    if (max_iteration <= 0) max_iteration = 1000;
    for (int row=0; row < img->h; row++) {
      for (int col=0,iter=0; col < img->w; col++, pp+=img->pixel_size) {
	double x=0, x0 = cx + (col - img->w/2) / scale;
	double y=0, y0 = cy - (row - img->h/2) / scale;
	for (iter=0; x*x+y*y <= (2*2) && iter < max_iteration; iter++) {
	  double xtemp = x*x - y*y + x0;
	  y = 2*x*y + y0;
	  x = xtemp;
	}
	if (gray) *pp = (iter==max_iteration ? 0 : 255);
	else if (iter==max_iteration) IMGH3V_SET( pp, 0, 0, 0 );
	else {
	  cmap.getColor( iter, rgb );
	  memcpy( pp, rgb, 3*sizeof(unsigned char) );
	}
      }
    }
  }
  
  // -----------------------------------------------------------------
  // Image examples for camera calibration
  // -----------------------------------------------------------------
public:
  void drawCalibPattern(IMGH::Image *img, int rrows=5, int rcols=8) {
    // Draw camera calibration pattern on the given image.
    if ( !img || img->w < 1 || img->h < 1 ) return;
    int  size, rgb[3] = {0, 0, 0};
    int  x, y, mx, my, nx, ny, i, j;
    img->clearImage( 255 );
    nx = rcols;  ny = rrows;
    if (img->w/rcols > img->h/rrows) size = img->h / (rrows+2);
    else                             size = img->w / (rcols+2);
    mx   = (img->w - rcols * size) / 2;
    my   = (img->h - rrows * size) / 2;
    for (i = 0, x = mx; i < nx; i++, x+=size)
      for (j = 0, y = my; j < ny; j++, y+=size) 
	if ((i+j)%2 == 0) drawCalibRect( img, x, y, size, rgb );
  }
private:
  void drawCalibRect  (IMGH::Image* img, int x, int y, int size, int rgb[3]) {
    int     i, j, idx;
    uchar_t pv[4] = { (uchar_t)(rgb[0]), (uchar_t)(rgb[1]), (uchar_t)(rgb[2]) };
    for (j = y; j < y + size && j > 0 && j < img->h-1; j++) 
      for (i = x; i < x + size && i > 0 && i < img->w-1; i++) {
	idx = img->getIndex( i, j );
	img->setPixel( idx, pv );
      }
  }
  
public:
  void draw4Corners(IMGH::Image *img) {
    // Draw camera calibration pattern on the given image.
    if ( !img || img->w < 1 || img->h < 1 ) return;
    int  w = img->w,  h = img->h, size, rgb[3] = {0, 0, 0};
    img->clearImage( 255 );
    size = ((int)(h/100)) * 10;
    drawCalibCorner( img, 0 + 2*size, 0 + 2*size, size, rgb, false );
    drawCalibCorner( img, w - 2*size, 0 + 2*size, size, rgb, true  );
    drawCalibCorner( img, 0 + 2*size, h - 2*size, size, rgb, true  );
    drawCalibCorner( img, w - 2*size, h - 2*size, size, rgb, false );
  }
private:  
  void drawCalibCorner(IMGH::Image* img, int x, int y, int radius, int rgb[3], bool dir) {
    // draw a corner (for camera calibration) at pixel (x,y)
    int     i, j, idx;
    float   dx, dy, dist, alpha;
    uchar_t pv[4];
    for (i = -radius; i < +radius; i++) {
      dx = i + 0.5f;
      if (x+dx < 0 || x+dx > img->w) continue;
      for (j = -radius; j < +radius; j++) {
	dy = j + 0.5f;
	if (y+dy < 0 || y+dy > img->h) continue;
	if ( dir && dx*dy > 0) continue;
	if (!dir && dx*dy < 0) continue;
	if ((dist = (float)sqrt(dx*dx + dy*dy)) > radius) continue;
	alpha = (1.0 - dist / radius);
	idx = img->getIndex( x+i, y+j );
	img->getPixel( idx, pv );
	if (img->type == PIXEL_GRAY) {
	  pv[0] = (unsigned char)( rgb[0] * alpha + pv[0] * (1-alpha));
	} else {
	  pv[0] = (unsigned char)( rgb[0] * alpha + pv[0] * (1-alpha));
	  pv[1] = (unsigned char)( rgb[1] * alpha + pv[1] * (1-alpha) );
	  pv[2] = (unsigned char)( rgb[2] * alpha + pv[2] * (1-alpha) );
	}
	img->setPixel( idx, pv );
      }
    }
  }
  
public:
  void drawTestPatternForSIFT(IMGH::Image *img) {
    // Draw a black rectangle on the image, which is used for testing SIFT.
    if ( !img || img->w < 1 || img->h < 1 ) return;
    unsigned char black[3]={0,0,0}, white[3]={255,255,255};
    int xl[2] = { img->w/3, img->w/3 + img->w/3 };
    int yl[2] = { img->h*2/5, img->h*2/5 + img->h/5 };
    for (int y = 0; y < img->h; y++) {
      for (int x = 0; x < img->w; x++) {
	bool in = (x >= xl[0] && x < xl[1] && y >= yl[0] && y < yl[1]);
	if (in) img->setPixel( x, y, black );
	else    img->setPixel( x, y, white );
      }
    }
  }
};

  
}	// namespace IMGH


#endif	// IMGH_IMAGE_EXAMPLE_HPP
