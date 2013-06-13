
//
// UTIL::Timer class object
//
// Jaeil Choi
// last modified in Jan, 2004
//

#ifndef UTIL_PIXFMT_CONVERTER_HPP
#define UTIL_PIXFMT_CONVERTER_HPP

#include "util_color.h"
#include "stdlib.h"
#include "string.h" 
// LIMIT: convert a 16.16 fixed-point value to a byte, with clipping.
#define LIMIT(x) ((x)>0xffffff?0xff: ((x)<=0xffff?0:((x)>>16)))

namespace UTIL {
  
class PixFmtConverter
{
public:
  PixFmtConverter() {}
  ~PixFmtConverter() {}
  
  // =================================================================
  // =================================================================
  
public:
  void rotate180(void *buf, int w, int h, int bpp) {
    // Rotate the image 180 degrees.  (bpp : bytes per pixel)
    int  i, total = w * h;
    unsigned char *a = (unsigned char*)buf;
    unsigned char *b = (unsigned char*)buf + (total-1)*bpp;
    unsigned char *tmp = (unsigned char*)malloc(bpp);
    if (tmp) {
      for (i = 0; i < total/2; i++) {
	memcpy( tmp,   a, bpp );
	memcpy(   a,   b, bpp );
	memcpy(   b, tmp, bpp );
	a += bpp;  b -= bpp;
      }
      free(tmp);
    }
  }
  
  void flipV(void *buf, int w, int h, int bpp) {
    // Flip the image upside down.  (bpp : bytes per pixel)
    int  i, row_size=w*bpp;
    unsigned char *a = (unsigned char*)buf;
    unsigned char *b = (unsigned char*)buf + (h-1)*row_size;
    unsigned char *tmp = (unsigned char*)malloc(row_size);
    if (tmp) {
      for (i = 0; i < h/2; i++) {
	memcpy( tmp,   a, row_size );
	memcpy(   a,   b, row_size );
	memcpy(   b, tmp, row_size );
	a += row_size;  b -= row_size;
      }
      free(tmp);
    }
  }
  
  // =================================================================
  // conversion to RGB24
  // =================================================================

public:
  void convert_gray8_rgb24(void *source, void *dest, int w, int h) {
    int           i, total = w*h;
    unsigned char *src=(unsigned char*)source;
    unsigned char *dst=(unsigned char*)dest;
    for (i = 0; i < total; i++, src++, dst+=3) 
      dst[0] = dst[1] = dst[2] = src[0];
  }
  
  void convert_gray16_rgb24(void *source, void *dest, int w, int h) {
    int           i, total = w*h;
    unsigned char *src=(unsigned char*)source;
    unsigned char *dst=(unsigned char*)dest,  value;
    for (i = 0; i < total; i++, src+=2, dst+=3) {
      value = (unsigned char)( *((unsigned short int*)src) / 65535.0f * 255.0f );
      dst[0] = dst[1] = dst[2] = value;
    }
  }
  
  void convert_yuyv_rgb24(void *source, void *dest, int w, int h) {
    int           i, total=w*h;
    unsigned char *src=(unsigned char*)source;
    unsigned char *dst=(unsigned char*)dest;
    for (i=0; i<total; i+=2, src+=4, dst+=6) COLOR_YUYV_RGB( src, dst );
  }
  
  void convert_uyvy_rgb24(void *source, void *dest, int w, int h) {
    int           i, total=w*h;			// YUV422
    unsigned char *src=(unsigned char*)source;
    unsigned char *dst=(unsigned char*)dest;
    for (i=0; i<total; i+=2, src+=4, dst+=6) COLOR_UYVY_RGB( src, dst );
  }
  
  void convert_uyv_rgb24(void *source, void *dest, int w, int h) {
    int           i, total=w*h;			// YUV444
    unsigned char *src=(unsigned char*)source;
    unsigned char *dst=(unsigned char*)dest;
    for (i=0; i<total; i++, src+=3, dst+=3) COLOR_UYV_RGB( src, dst );
  }
  
  void convert_uyyvyy_rgb24(void *source, void *dest, int w, int h) {
    int           i, total=w*h;			// YUV411
    unsigned char *src=(unsigned char*)source;
    unsigned char *dst=(unsigned char*)dest;
    for (i=0; i<total; i+=4, src+=6, dst+=12) COLOR_UYYVYY_RGB( src, dst );
  }

  void convert_yuv420p_rgb24(void *pIn0, void *pOut0, int width, int height) {
    // Converts from planar YUV420P to RGB24.
    //   Consider a YUV420P image of 8x2 pixels.
    //   A plane of Y values    A B C D E F G H
    //                          I J K L M N O P
    //   A plane of U values    1   2   3   4 
    //   A plane of V values    1   2   3   4 ....
    //   The U1/V1 samples correspond to the ABIJ pixels.
    //       U2/V2 samples correspond to the CDKL pixels.
    const int numpix = width * height;
    const int bytes = 24 >> 3;
    int i, j, y00, y01, y10, y11, u, v;
    unsigned char *pY = (unsigned char*)pIn0;
    unsigned char *pU = pY + numpix;
    unsigned char *pV = pU + numpix / 4;
    unsigned char *pOut = (unsigned char*)pOut0;

    for (j = 0; j <= height - 2; j += 2) {
      for (i = 0; i <= width - 2; i += 2) {
	y00 = *pY;
	y01 = *(pY + 1);
	y10 = *(pY + width);
	y11 = *(pY + width + 1);
	u = (*pU++) - 128;
	v = (*pV++) - 128;
      
	// Note that U and V are swapped... 
	move_420_block(y00, y01, y10, y11, v, u, width, pOut);
    
	pY += 2;
	pOut += 2 * bytes;
      }
      pY += width;
      pOut += width * bytes;
    }
  }

  // Converts from interlaced YUV420 to RGB24.
  //   Consider a YUV420 image of 6x2 pixels.
  //   A B C D U1 U2
  //   I J K L V1 V2
  //   The U1/V1 samples correspond to the ABIJ pixels.
  //       U2/V2 samples correspond to the CDKL pixels.
  void convert_yuv420_rgb24(void *pIn0, void *pOut0, int width, int height) {
    const int bytes = 24 >> 3;
    int i, j, y00, y01, y10, y11, u, v;
    unsigned char *pY = (unsigned char*)pIn0;
    unsigned char *pU = pY + 4;
    unsigned char *pV = pU + width;
    unsigned char *pOut = (unsigned char*)pOut0;

    for (j = 0; j <= height - 2; j += 2) {
      for (i = 0; i <= width - 4; i += 4) {
	y00 = *pY;
	y01 = *(pY + 1);
	y10 = *(pY + width);
	y11 = *(pY + width + 1);
	u = (*pU++) - 128;
	v = (*pV++) - 128;

	move_420_block(y00, y01, y10, y11, u, v, width, pOut);
    
	pY += 2;
	pOut += 2 * bytes;

	y00 = *pY;
	y01 = *(pY + 1);
	y10 = *(pY + width);
	y11 = *(pY + width + 1);
	u = (*pU++) - 128;
	v = (*pV++) - 128;

	move_420_block(y00, y01, y10, y11, u, v, width, pOut);
    
	pY += 4; // skip UV
	pOut += 2 * bytes;
      }
      pY += width;
      pOut += width * bytes;
    }
  }
  
  // =================================================================
  // private functions
  // =================================================================
private:
  /*
   * Turn a YUV4:2:0 block into an RGB block
   *
   * Video4Linux seems to use the blue, green, red channel
   * order convention-- rgb[0] is blue, rgb[1] is green, rgb[2] is red.
   *
   * Color space conversion coefficients taken from the excellent
   * http://www.inforamp.net/~poynton/ColorFAQ.html
   * In his terminology, this is a CCIR 601.1 YCbCr -> RGB.
   * Y values are given for all 4 pixels, but the U (Pb)
   * and V (Pr) are assumed constant over the 2x2 block.
   *
   * To avoid floating point arithmetic, the color conversion
   * coefficients are scaled into 16.16 fixed-point integers.
   * They were determined as follows:
   *
   *  double brightness = 1.0;  (0->black; 1->full scale) 
   *  double saturation = 1.0;  (0->greyscale; 1->full color)
   *  double fixScale = brightness * 256 * 256;
   *  int rvScale = (int)(1.402 * saturation * fixScale);
   *  int guScale = (int)(-0.344136 * saturation * fixScale);
   *  int gvScale = (int)(-0.714136 * saturation * fixScale);
   *  int buScale = (int)(1.772 * saturation * fixScale);
   *  int yScale = (int)(fixScale);   
   */
  void move_420_block(int yTL, int yTR, int yBL, int yBR, int u, int v, 
		      int rowPixels, unsigned char * rgb) {
    const int rvScale = 91881;
    const int guScale = -22553;
    const int gvScale = -46801;
    const int buScale = 116129;
    const int yScale  = 65536;
    int r, g, b;
    g = guScale * u + gvScale * v;
    //  if (force_rgb) {
    //      r = buScale * u;
    //      b = rvScale * v;
    //  } else {
    r = rvScale * v;
    b = buScale * u;
    //  }
    yTL *= yScale; yTR *= yScale;
    yBL *= yScale; yBR *= yScale;
    /* Write out top two pixels */
    rgb[0] = LIMIT(b+yTL); rgb[1] = LIMIT(g+yTL);
    rgb[2] = LIMIT(r+yTL);
    rgb[3] = LIMIT(b+yTR); rgb[4] = LIMIT(g+yTR);
    rgb[5] = LIMIT(r+yTR);
    /* Skip down to next line to write out bottom two pixels */
    rgb += 3 * rowPixels;
    rgb[0] = LIMIT(b+yBL); rgb[1] = LIMIT(g+yBL);
    rgb[2] = LIMIT(r+yBL);
    rgb[3] = LIMIT(b+yBR); rgb[4] = LIMIT(g+yBR);
    rgb[5] = LIMIT(r+yBR);
  }
  
  // =================================================================
  // conversion to Gray8
  // =================================================================
  
public:
  void convert_rgb24_gray8(void *source, void *dest, int w, int h) {
    int      i, total = w*h;
    unsigned char *src=(unsigned char*)source;
    unsigned char *dst=(unsigned char*)dest;
    for (i=0; i<total; i++, src+=3, dst++) dst[0] = COLOR_GRAY(src);
  }
  void convert_yuyv_gray8(void *source, void *dest, int w, int h) {
    int      i, total = w*h;			// 
    unsigned char *src=(unsigned char*)source;
    unsigned char *dst=(unsigned char*)dest;
    for (i=0; i<total; i++, src+=2, dst++) dst[0] = src[0];
  }
  void convert_uyvy_gray8(void *source, void *dest, int w, int h) {
    int      i, total = w*h;			// YUV422
    unsigned char *src=(unsigned char*)source;
    unsigned char *dst=(unsigned char*)dest;
    for (i=0; i<total; i++, src+=2, dst++) dst[0] = src[1];
  }
  void convert_uyv_gray8(void *source, void *dest, int w, int h) {
    int      i, total = w*h;			// YUV444
    unsigned char *src=(unsigned char*)source;
    unsigned char *dst=(unsigned char*)dest;
    for (i=0; i<total; i++, src+=3, dst++) dst[0] = src[1];
  }
  void convert_uyyvyy_gray8(void *source, void *dest, int w, int h) {
    int      i, total = w*h;			// YUV411
    unsigned char *src=(unsigned char*)source;
    unsigned char *dst=(unsigned char*)dest;
    for (i=0; i<total; i+=4, src+=6, dst+=4) {
      dst[0] = src[1];  dst[1] = src[2];
      dst[2] = src[4];  dst[3] = src[5];
    }
  }
  void convert_gray16_gray8(void *source, void *dest, int w, int h) {
    int      i, total = w*h;
    unsigned char *src=(unsigned char*)source;
    unsigned char *dst=(unsigned char*)dest;
    for (i=0; i<total; i++, src+=2, dst++) 
      dst[0] = (unsigned char)( *((unsigned short int*)src) / 65535.0f * 255.0f );
  }

  // =================================================================
  // conversion to YUV24
  // =================================================================
  
public:
  void convert_yuyv_yuv24(void *source, void *dest, int w, int h) {
    int      row, col, dst_stride=w*3;
    unsigned char *src=(unsigned char*)source;
    unsigned char *dst=(unsigned char*)dest;
    for (row = 0; row < h; row++) {
      dst = ((unsigned char*)dest) + row * dst_stride;
      for (col = 0; col < w; col+=2) {
	dst[0] = src[0];  dst[1] = src[1];  dst[2] = src[3];
	dst[3] = src[2];  dst[4] = src[1];  dst[5] = src[3];
	src+=4;  dst+=6;
      }
    }
  }
};

} 	// end of namespace UTIL

#endif	// UTIL_PIXFMT_CONVERTER_HPP

