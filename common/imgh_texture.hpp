
//
// IMGH::Texture
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

#ifndef IMGH_TEXTURE_HPP
#define IMGH_TEXTURE_HPP

#include <iostream>
#include "imgh_common.hpp"
#include "imgh_editor.hpp"

namespace IMGH {
  
class Texture
{
public:
  Image	*img;
  int    id;
  int    w, h;
  bool   texture_ready;
 private:
  char filename[256];
public:
  Texture () : w(0), h(0), img(NULL) { }
  Texture (IMGH::Image *img, int id=0) : w(0), h(0), img(NULL) { create(img, id); }
  ~Texture() {}
  void clear(void) {
    id = 0;
    texture_ready = false;
  }
  
public:
  bool create(IMGH::Image *img, int id=0) {
    if (img == NULL || img->type != PIXEL_RGB) return false;
    // check width and height
    for (int tw = img->w; tw > 1; tw = tw/2) 
      if (tw % 2 != 0) { 
	std::cerr << "Error (IMGH::Texture): invalid texture size (" << img->w << " x " << img->h << ")" << std::endl;
	return false;
      }
    for (int th = img->h; th > 1; th = th/2)
      if (th % 2 != 0) { 
	std::cerr << "Error (IMGH::Texture): invalid texture size (" << img->w << " x " << img->h << ")" << std::endl;
	return false;
      }
    IMGH::ImageEditor editor;
    editor.flip( img, true );
  
    // create the texture;
    glBindTexture(GL_TEXTURE_2D, id);
    //glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
    //glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
    //glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
    //glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR); // GL_NEAREST);
    glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_DECAL); // GL_MODULATE);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, img->w, img->h, 0, GL_RGB, GL_UNSIGNED_BYTE, img->data);
    int eno = glGetError();
    if (eno != 0) {
      cerr << "Error(Texture::loadTexture): " << gluErrorString(eno) << endl;
      return false;
    }
    this->id = id;  this->w = img->w;  this->h = img->h;
    texture_ready = true;
  
    editor.flip( img, true );
    return true;
  }
  
};

  
}	// namespace IMGH


#endif	// IMGH_TEXTURE_HPP


