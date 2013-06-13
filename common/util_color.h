
//
// color macro definitions
//
// Jaeil Choi
// last modified in Dec, 2004
//
//
  

#ifndef UTIL_COLOR_HPP
#define UTIL_COLOR_HPP

#include <iostream>

typedef enum { 
  PIXEL_UNKNOWN=0, 	// 
  PIXEL_GRAY, 		// 1 x 8 bits
  PIXEL_RGB, 		// 3 x 8 bits  
  PIXEL_BGR, 		// 3 x 8 bits  
  PIXEL_YUV,		// 3 x 8 bits
  PIXEL_UYV,		// 3 x 8 bits  FourCC:IYU2  1394:YUV444 mode=0
  PIXEL_UYVY, 		// 2 x 8 bits  FourCC:UYVY  1394:YUV422 mode=1
  PIXEL_YUYV, 		// 2 x 8 bits  FourCC:YUY2
  PIXEL_UYYVYY		// 12 bits     FourCC:IYU1  1394:YUV411 mode=2
} PixelFormat_t;


// ===================================================================
// Color in String  "#ff7700"
// ===================================================================

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
    

// ===================================================================
// Gray scale
//   Charles Poynton's Color FAQ, <http://www.inforamp.net/~poynton/> 
//   Copyright (c) 1998-01-04 Charles Poynton poynton@inforamp.net
//   Y = 0.21268 * R + 0.7151 * G + 0.07217 * B
//   Y = (6969 * R + 23434 * G + 2365 * B)/32768
// ===================================================================

#define COLOR_GRAY(rgb) \
        ((6969*(int)(rgb)[0] + 23434*(int)(rgb)[1] + 2365*(int)(rgb)[2])/32768)

#define COLOR_RGB_GRAY(rgb, gray) \
        do { gray = ((6969*(int)(rgb)[0] + 23434*(int)(rgb)[1] + 2365*(int)(rgb)[2])/32768); } while(0)

#define COLOR_BGR_GRAY(bgr, gray) \
        do { gray = ((6969*(int)(bgr)[2] + 23434*(int)(bgr)[1] + 2365*(int)(bgr)[0])/32768); } while(0)

#define COLOR_GRAY_RGB(gray, rgb) \
        do { r = g = b = gray; } while(0)

#define COLOR_GRAY_BGR(gray, bgr) \
        do { r = g = b = gray; } while(0)

#define COLOR_BGR_GRAY_IMAGE(src, dst, w, h, s_stride, d_stride) \
        do { unsigned char *s, *d; \
          for (int y = 0; y < h; y++) { \
            s = (unsigned char*)src + y * s_stride; \
            d = (unsigned char*)dst + y * d_stride; \
            for (int x = 0; x < w; x++) { \
              COLOR_BGR_GRAY( s, d[0] ); \
              s += 3;  d++; \
            } \
          } \
        } while (0)
              

// ===================================================================
// YUV
//   RGB to YUV Conversion:
//     Y  =      (0.257 * R) + (0.504 * G) + (0.098 * B) + 16
//     Cr = V =  (0.439 * R) - (0.368 * G) - (0.071 * B) + 128
//     Cb = U = -(0.148 * R) - (0.291 * G) + (0.439 * B) + 128
//   YUV to RGB Conversion:
//     B = 1.164(Y - 16)                  + 2.018(U - 128)
//     G = 1.164(Y - 16) - 0.813(V - 128) - 0.391(U - 128)
//     R = 1.164(Y - 16) + 1.596(V - 128)
// ===================================================================

#define COLOR_RGB_YUV_YUV(r, g, b, y, u, v, yuv) \
   do { \
     y = (int)( 0.257 * r + 0.504 * g + 0.098 * b) +  16; \
     v = (int)( 0.439 * r - 0.368 * g - 0.071 * b) + 128; \
     u = (int)(-0.148 * r - 0.291 * g + 0.439 * b) + 128; \
     if (y < 0) y = 0;  if (y > 255) y = 255; \
     if (u < 0) u = 0;  if (u > 255) u = 255; \
     if (v < 0) v = 0;  if (v > 255) v = 255; \
     (yuv)[0] = y;  (yuv)[1] = u;  (yuv)[2] = v; \
   } while (0)

#define COLOR_YUV_RGB_RGB(y, u, v, r, g, b, rgb) \
   do { \
     b = (int)(1.164*(y - 16)                   + 2.018*(u - 128)); \
     g = (int)(1.164*(y - 16) - 0.813*(v - 128) - 0.391*(u - 128)); \
     r = (int)(1.164*(y - 16) + 1.596*(v - 128)); \
     if (r < 0) r = 0;  if (r > 255) r = 255; \
     if (g < 0) g = 0;  if (g > 255) g = 255; \
     if (b < 0) b = 0;  if (b > 255) b = 255; \
     (rgb)[0] = r;  (rgb)[1] = g;  (rgb)[2] = b; \
   } while (0)

#if 0
#define COLOR_RGB_YUV_YUV(r, g, b, y, u, v, yuv) \
   do { \
     y = (( 9798 * (int)r) + ( 19235 * (int)g) + ( 3736 * (int)b)) / 32768; \
     v = ((-4784 * (int)r) + ( -9437 * (int)g) + (14221 * (int)b)) / 32768 + 128; \
     u = ((20218 * (int)r) + (-16941 * (int)g) + (-3277 * (int)b)) / 32768 + 128; \
     if (y < 0) y = 0;  if (y > 255) y = 255; \
     if (u < 0) u = 0;  if (u > 255) u = 255; \
     if (v < 0) v = 0;  if (v > 255) v = 255; \
     (yuv)[0] = y;  (yuv)[1] = u;  (yuv)[2] = v; \
   } while (0)
#define COLOR_YUV_RGB_RGB(y, u, v, r, g, b, rgb) \
   do { \
     b = y + ((int(u) * 1814) >> 10); \
     g = y - ((int(u) * 352 + int(v) * 731) >> 10); \
     r = y + ((int(v) * 1436) >> 10); \
     if (r < 0) r = 0;  if (r > 255) r = 255; \
     if (g < 0) g = 0;  if (g > 255) g = 255; \
     if (b < 0) b = 0;  if (b > 255) b = 255; \
     (rgb)[0] = r;  (rgb)[1] = g;  (rgb)[2] = b; \
   } while (0)
#endif

#define COLOR_YUV_RGB(yuv, rgb) \
   do { int r, g, b; \
     COLOR_YUV_RGB_RGB( (yuv)[0], (yuv)[1], (yuv)[2], r, g, b, rgb ); \
   } while (0)

#define COLOR_RGB_YUV(rgb, yuv) \
   do { int y, u, v; \
     COLOR_RGB_YUV_YUV( (rgb)[0], (rgb)[1], (rgb)[2], y, u, v, yuv ); \
   } while (0)

// convert 2 pixels of YUYV (YUY2) to RGB values (YUYV -> RGBRGB)
#define COLOR_YUYV_RGB(yuyv, rgb) \
   do { int r, g, b; \
     COLOR_YUV_RGB_RGB( (yuyv)[0], (yuyv)[1], (yuyv)[3], r, g, b, rgb   );  \
     COLOR_YUV_RGB_RGB( (yuyv)[2], (yuyv)[1], (yuyv)[3], r, g, b, rgb+3 );  \
   } while (0)

// convert 2 pixels of UYVY (UYVY) to RGB values (UYVY -> RGBRGB)
#define COLOR_UYVY_RGB(uyvy, rgb) \
   do { int r, g, b; \
     COLOR_YUV_RGB_RGB( (uyvy)[1], (uyvy)[0], (uyvy)[2], r, g, b, rgb   );  \
     COLOR_YUV_RGB_RGB( (uyvy)[3], (uyvy)[0], (uyvy)[2], r, g, b, rgb+3 );  \
   } while (0)

// convert 2 pixels of RGB values to YUYV (YUY2) (RGBRGB -> YUYV)
#define COLOR_RGB_YUYV(rgb, yuyv) \
   do { int ar, ag, ab, y, u, v; \
     ar = ((int)(rgb)[0] + (int)(rgb)[3]) / 2; \
     ag = ((int)(rgb)[1] + (int)(rgb)[4]) / 2; \
     ab = ((int)(rgb)[2] + (int)(rgb)[5]) / 2; \
     v = ((-4784 * ar) + ( -9437 * ag) + (14221 * ab)) / 32768 + 128; \
     u = ((20218 * ar) + (-16941 * ag) + (-3277 * ab)) / 32768 + 128; \
     if (u < 0) u = 0;  if (u > 255) u = 255; \
     if (v < 0) v = 0;  if (v > 255) v = 255; \
     yuyv[1] = u;  yuyv[3] = v;  \
     y = (( 9798 * (int)(rgb)[0]) + ( 19235 * (int)(rgb)[1]) + ( 3736 * (int)(rgb)[2])) / 32768; \
     if (y < 0) y = 0;  if (y > 255) y = 255; \
     yuyv[0] = y; \
     y = (( 9798 * (int)(rgb)[3]) + ( 19235 * (int)(rgb)[4]) + ( 3736 * (int)(rgb)[5])) / 32768; \
     if (y < 0) y = 0;  if (y > 255) y = 255; \
     yuyv[2] = y; \
   } while (0)

// convert YUV444(IYU2) to RGB (UYV -> RGB)
#define COLOR_UYV_RGB(uyv, rgb) \
   do { int r, g, b; \
     COLOR_YUV_RGB_RGB( (uyv)[1], (uyv)[0], (uyv)[2], r, g, b, rgb ); \
   } while (0)

// convert YUV411(IYU1) to RGB (UYYVYY -> RGBRGBRGBRGB)
#define COLOR_UYYVYY_RGB(uyyvyy, rgb) \
   do { int r, g, b; \
     COLOR_YUV_RGB_RGB( (uyyvyy)[1], (uyyvyy)[0], (uyyvyy)[3], r, g, b, (rgb)   ); \
     COLOR_YUV_RGB_RGB( (uyyvyy)[2], (uyyvyy)[0], (uyyvyy)[3], r, g, b, (rgb)+3 ); \
     COLOR_YUV_RGB_RGB( (uyyvyy)[4], (uyyvyy)[0], (uyyvyy)[3], r, g, b, (rgb)+6 ); \
     COLOR_YUV_RGB_RGB( (uyyvyy)[5], (uyyvyy)[0], (uyyvyy)[3], r, g, b, (rgb)+9 ); \
   } while (0)


#endif // UTIL_COLOR_HPP

