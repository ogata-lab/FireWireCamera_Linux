
//
// GUIH              : GUI in Header files
// GUIH::ImageWindow : window for image visualization
//
// (c) 2006  Jaeil Choi
// last modified in Sep, 2009
//
// --------------------------------------------------------------------------
//  This program is free software; you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation; either version 2 of the License, or
//  (at your option) any later version.
//
//  This program is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
//
//  You should have received a copy of the GNU General Public License
//  along with this program; if not, write to the Free Software
//  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
// --------------------------------------------------------------------------
// 
// GUIH is a very simple, light-weight and cross-platform Graphic User
//  Interface that can visualize images and OpenGL graphics with
//  keyboard and mouse interactions on MS Windows, Apple OS X, and Linux.
// It is intended as a tool for developers and researchers, who don't want
//  to spend time in designing/creating an user interface, but who want 
//  all the features of a general user interface and the freedom to modify. 
//  If you are not satisfied, like me, with GLUT, GLUI, FLTK or HighGUI, 
//  you may find GUIH very useful and convenient.
// Key features include dialog windows, timers, default key bindings, 
//  image capture, captions, mapping of mouse click to world-space coordinates,
//  and the built-in cameras for 2D/3D OpenGL.
// GUIH doesn't support a window with controls, except simple input dialogs.
//  For complex settings/configuration of an application, I recommend using
//  a configuration file, rather than a window full of controls.
//
// GUIH package consists of eight C++ header files:
//   guih_common.hpp     : the base window   GUIH::Window
//   guih_image.hpp      : image window      GUIH::ImageWindow
//   guih_opengl2d.hpp   : 2D OpenGL Window  GUIH::OpenGL2DWindow
//   guih_opengl3d.hpp   : 3D OpenGL Window  GUIH::OpenGL3DWindow
//   guih_camera2d.hpp   : 2D camera for     GUIH::OpenGL2DWindow
//   guih_camera3d.hpp   : 3D camera for     GUIH::OpenGL3DWindow
//   guih_args.hpp       : command-line argument processing tool (GUIH::Args)
//   guih_config.hpp     : configuration file processing tool    (GUIH::Config)
// GUIH package defines four window classes:
//   Window              : the base window (Not for direct use!)
//    +- ImageWindow     : for image visualization
//    +- OpenGL2DWindow  : for 2D OpenGL rendering (with built-in 2D camera)
//    +- OpenGL3DWindow  : for 3D OpenGL rendering (with built-in 3D camera)
// GUIH package defines two extra classes:
//   Args                : command-line argument processor
//   Config              : configuration file processor
//
// To compile on Windows, 
//   - use Visual Studio 8.0 or higher (for GDI+)
//   - create an empty 'Win32 Console Application' project
//   - add 'src' directory to Project->Properties->C/C++->General->AdditionalIncludeDirectories
//   - add '_CRT_SECURE_NO_DEPRECATE' to Project->Properties->C/C++->Preprocessor->PreprocessorDefinitions
//   - add 'opengl32.lib glu32.lib gdiplus.lib' to Project->Properties->Linker->Input->AddtionalDependencies
// To compile on Apple OS X,
//   - install 'GLUT for Mac OS X' (for rendering caption on OpenGL windows), and
//   - use g++ compiler flags : '-framework Carbon -framework AGL -framework OpenGL -framework GLUT'
// To compile on Linux,
//   - install GTKGLExt, and add `pkg-config --cflag --libs gtkglext' to Makefile flags
//
// --------------------------------------------------------------------------
//
// GUIH::ImageWindow
//   This class provides a GUI for image file I/O and visualization, 
//   together with simple functions to draw lines, rectangles, etc.
//   The image can be accessed directly using the member variables
//   - imgw, imgh, imgtype, and imgdata.
//
//   Required library: nothing (Windows, OS X, and Linux)
//   Member functions:
//      void clearImage( unsigned char v=0 );
//      void removeImage( );
//      bool createImage( width, height, pixel_t type=PIXEL_RGB );
//	bool showImage( int w, int h, GUIH::pixel_t type, void *data);
//	bool showImage( char *fname );
//      bool saveImage( char *fname );
//      bool copyImage( int w, int h, unsigned char *data, GUIH::pixel_t type=PIXEL_RGB );
//
// For API example, read the example code at the bottom.
// For default key bindings, press Shift+F1 during execution.
//


#ifndef GUIH_IMAGE_HPP
#define GUIH_IMAGE_HPP

#include <iostream>
#include "guih_common.hpp"


// ===================================================================
#ifdef  WIN32   // beginning of Window version
// ===================================================================

namespace GUIH {
  
class ImageWindow : public Window
{
public:
  int		imgw, imgh;	// image width and height
  pixel_t	imgtype;	// pixel format (GRAY, RGB, BGR, YUV)
  unsigned char*imgdata;	// pixel values
  bool		img_allocated;	// whether or not the pixel buffer was allocated
  bool		ignore_alpha;	// whether or not to ignore alpha channel when draw the image
public:
  ImageWindow () 
    : imgw(0), imgh(0), imgdata(NULL), img_allocated(false), ignore_alpha(false) { }
  ImageWindow (int w, int h, char *title = NULL, 
	       void (*init) (Window *w) = NULL,
	       void (*draw) (Window *w) = NULL,
	       void (*keys) (Window *w, int k) = NULL, 
	       void (*mouse)(Window *w, int button, int state, int xy[]) = NULL,
	       void (*close)(Window *w) = NULL) 
    : imgw(0), imgh(0), imgdata(NULL), img_allocated(false), ignore_alpha(false) {
    openWindow(w, h, title, init, draw, keys, mouse, close); 
  }
  ~ImageWindow() { }
public:
  void openWindow(int w, int h, const char *title = NULL, 
		  void (*init) (Window *w) = NULL,
		  void (*draw) (Window *w) = NULL,
		  void (*keys) (Window *w, int k) = NULL, 
		  void (*mouse)(Window *w, int button, int state, int xy[]) = NULL,
		  void (*close)(Window *w) = NULL) {
    closeWindow();  removeImage();
    setCallbacks(draw, keys, mouse, close);
    createBasicObjects( w, h, title );
  }
  
  // image -----------------------------------------------------------
  int  getPixelSize(pixel_t type) {
    int size;
    switch (type) {	// Note that the unit is 'sizeof(unsigned char)'
    case PIXEL_GRAY:  size = 1;  break;
    case PIXEL_RGB:   size = 3;  break;
    case PIXEL_BGR:   size = 3;  break;
    case PIXEL_YUV:   size = 3;  break;
    case PIXEL_RGBA:  size = 4;  break;
    case PIXEL_GRAYA: size = 2;  break;
    default: 0;
    }
    return size;
  }
  void clearImage(unsigned char v=0) { if (imgdata) memset(imgdata, v, imgw * imgh * getPixelSize(imgtype)); }
  void removeImage(void) { 
    imgw = imgh = 0;  
    if (img_allocated && imgdata) { free(imgdata); img_allocated = false; }
    imgdata = NULL;
  }
  bool createImage(int w, int h, pixel_t type=PIXEL_RGB) {
    // Create a new image by allocating memory (only if it's necessary).
    if (w <= 0 || h <= 0) return false;
    if (img_allocated && w==imgw && h==imgh && type==imgtype) return true;
    if (img_allocated && imgdata) free(imgdata);
    imgw = w;  imgh = h;  imgtype = type;
    int  pixel_size = getPixelSize(type);
    imgdata = (unsigned char*)calloc( w * h * pixel_size, sizeof(unsigned char) );
    img_allocated = true;
    resize( w, h, true );
    return true;
  }
  bool showImage( int w, int h, pixel_t type, void *data) {
    if (!hwnd || w <= 0 || h <= 0 || data == NULL) { removeImage(); return false; }
    setCaption();
    if (type <= PIXEL_UNKNOWN || type > PIXEL_YUV) {
      fprintf(stderr, "Error (GUIH::ImageWindow::showImage): invalid pixel format\n"); 
      removeImage();
      return false; 
    }
    int padding = (w % 4 == 0 ? 0 : (4 - w % 4));
    if (padding != 0) {
      std::cerr << "Warning: The image has a width that is not a multiple of four." << std::endl;
      std::cerr << "         A copy of the image was created for showImage() with extra padding." << std::endl;
      createImage( w + padding, h, GUIH::PIXEL_RGB );
      int   row, srcstride = w * sizeof(unsigned char) * 3, dststride = imgw * sizeof(unsigned char) * 3;
      char *src = (char*)data, *dst = (char*)imgdata;
      for (row = 0; row < imgh; row++) {
	memcpy( dst, src, srcstride );
	src += srcstride;  dst += dststride;
      }
      img_allocated = true;
    } else {
      if (img_allocated && imgdata) free(imgdata);
      imgw = w;  imgh = h;  imgtype = type;
      imgdata = (unsigned char*)data;
      img_allocated = false;
    }
    resize( w, h, true );
    return true;
  }
  bool showImage( const char *fname ) { 
    if (!hwnd || !fname) { removeImage(); return false; }
    setCaption();
    wchar_t wfname[512];
    mbstowcs( wfname, fname, (int)strlen(fname)+1 );
    Gdiplus::Bitmap bitmap(wfname);
    if (bitmap.GetWidth() > 0) {
      if (!hwnd) openWindow( imgw, imgh, fname );
      else resize( imgw, imgh, true );
      Gdiplus::Rect rect(0, 0, bitmap.GetWidth(), bitmap.GetHeight());
      Gdiplus::BitmapData bitmapData;
      bitmap.LockBits( &rect, Gdiplus::ImageLockModeRead, PixelFormat24bppRGB, &bitmapData );
      int padding = (bitmapData.Width % 4 == 0 ? 0 : (4 - bitmapData.Width % 4));
      if (padding != 0) {
	std::cerr << "Warning: The image has a width that is not a multiple of four." << std::endl;
	std::cerr << "         A copy of the image was created for showImage() with extra padding." << std::endl;
      }
      createImage( bitmapData.Width + padding, bitmapData.Height, GUIH::PIXEL_RGB );
      int   row, rowstride = imgw * sizeof(unsigned char) * 3;
      char *src = (char*)bitmapData.Scan0, *dst;
      for (row = 0, dst = (char*)imgdata; row < imgh; row++) {
	memcpy( dst, src, rowstride );
	src += bitmapData.Stride;  dst += rowstride;
      }
      img_allocated = true;
      bitmap.UnlockBits( &bitmapData );
      setTitle(fname);
      return true;
    } else {
      char buffer[256];
      sprintf( buffer, "failed to open '%s'", fname );
      openWindow( 300, 300, buffer );
      return false;
    }
  }
  bool saveImage(char *fname) { 
    // This function saves the image buffer 'imgdata' into a file.
    // To save what you see on the window, use 'captureScreen(fname)', instead.
    if (imgw % 4 != 0) {
      std::cerr << "Error: failed to save image. Its width must be a multiple of four." << std::endl;
      std::cerr << "       This is a restriction given by Windows GDI+ library." << std::endl;
      return false;
    }
    Gdiplus::Bitmap bitmap(imgw, imgh, imgw*3*sizeof(unsigned char), PixelFormat24bppRGB, imgdata);
    wchar_t wfname[512];
    mbstowcs( wfname, fname, (int)strlen(fname)+1 );
    CLSID encoderClsid;  if (!getEncoderClsid( fname, &encoderClsid )) return false;
    Gdiplus::Status stat = bitmap.Save( wfname, &encoderClsid, NULL );
    return (stat==0 ? true : false);
  }
  bool copyImage(int w, int h, unsigned char *data, pixel_t type=PIXEL_RGB) { 
    createImage( w, h, type );
    memcpy( imgdata, data, w * h * getPixelSize(type) );
    img_allocated = true;
    return true;
  }
  
  // -----------------------------------------------------------------
  // callback functions
  // -----------------------------------------------------------------

  bool runOnInit(void) { if (cb_init) cb_init( this ); return true; }
  void runOnDestroy(void) { 
    // This function is called automatically by a event, and cannot be cancelled
    if (cb_close) cb_close(this); 
    removeImage();
  }
  bool runOnDraw(void) {
    if (!hwnd) return false;
	Gdiplus::Graphics graphics(hwnd);
    switch (imgtype) {
    case PIXEL_RGB:
	  {
	    Gdiplus::Bitmap bitmap(imgw, imgh, imgw*3*sizeof(unsigned char), PixelFormat24bppRGB, imgdata);
		graphics.DrawImage( &bitmap, 0, 0 );
	  }
      break;
    case PIXEL_RGBA:
    case PIXEL_GRAY:
    case PIXEL_GRAYA:
	default: break;
    }
    if (cb_draw) cb_draw( this );
	return true;
  }
  void runOnKeyboard(int key) {
    bool skey = (key==ESCAPE_KEY || key>200);
    // Special Keys with default behaviors
    if (!skey_by_user) {
      if (key == HOME_KEY) {
	char *ptype;   
	switch (imgtype) {
	case PIXEL_GRAY:  ptype = "GRAY ";  break;
	case PIXEL_RGB:   ptype = "RGB  ";  break;
	case PIXEL_BGR:   ptype = "BGR  ";  break;
	case PIXEL_YUV:   ptype = "YUV  ";  break;
	case PIXEL_GRAYA: ptype = "GRAYA";  break;
	case PIXEL_RGBA:  ptype = "RGBA ";  break;
	default: ptype = "---- "; break;
	}
	printf("GUIH::ImageWindow image info: (%d x %d) %s(%d)  row_size=%d \n", 
	       imgw, imgh, ptype, getPixelSize(imgtype), imgw * (imgtype == PIXEL_GRAY ? 1 : 3) );
      }
    }
    if (key == F1 && shift) {
      std::cerr << "---------------------------------------------------------------------" << std::endl;
      std::cerr << "  GUIH::ImageWindow control keys                                     " << std::endl;
      std::cerr << "    Home   (Cmd-o) : Show image information                         " << std::endl;
      std::cerr << "    Insert (Cmd-i) : capture the image of the window                " << std::endl;
      std::cerr << "    Delete         : clear caption                                  " << std::endl;
      std::cerr << "---------------------------------------------------------------------" << std::endl;
    }
    // run user key press callback function
    //printf("key=%d  skey=%d  skey_by_user=%d  ctrl=%d alt=%d\n", key, skey, skey_by_user, ctrl, alt);
    if (cb_keys && (!skey||skey_by_user)) cb_keys( this, key );
  }
  void runOnMouse(int button, int state, int xy[]) { 
    if (button == MOUSE_LEFT && state == MOUSE_CLICKED) {
      char  str[80];  // default behavior: show the clicked position
      sprintf( str, "(%d, %d)", xy[0], xy[1]);
      setCaption( str );
    }
    if (cb_mouse) cb_mouse( this, button, state, xy ); 
  }
};

}	// end of GUIH namespace

// ===================================================================
#elif defined(__APPLE__) || defined(MACOSX)	// beginning of OS X version
// ===================================================================

namespace GUIH {

class ImageWindow : public Window
{
public:
  CGImageRef	imgref;
  int		imgw, imgh;	// image width and height
  pixel_t	imgtype;	// pixel format (GRAY, RGB, BGR, YUV)
  unsigned char*imgdata;	// pixel values
  bool		img_allocated;	// whether or not the pixel buffer was allocated
  bool		ignore_alpha;	// whether or not to ignore alpha channel when draw the image
public:
  ImageWindow () 
    : imgref(NULL), imgw(0), imgh(0), imgdata(NULL), img_allocated(false), ignore_alpha(false) { }
  ImageWindow (int w, int h, char *title = NULL, 
	       void (*init) (Window *w) = NULL,
	       void (*draw) (Window *w) = NULL,
	       void (*keys) (Window *w, int k) = NULL, 
	       void (*mouse)(Window *w, int button, int state, int xy[]) = NULL,
	       void (*close)(Window *w) = NULL) 
    : imgref(NULL), imgw(0), imgh(0), imgdata(NULL), img_allocated(false), ignore_alpha(false) {
    openWindow(w, h, title, init, draw, keys, mouse, close); 
  }
  virtual ~ImageWindow() { 
    if (imgref)  CGImageRelease( imgref );
    if (img_allocated && imgdata) free(imgdata); 
  }
public:
  void openWindow(int w, int h, const char *title = NULL, 
		  void (*init) (Window *w) = NULL,
		  void (*draw) (Window *w) = NULL,
		  void (*keys) (Window *w, int k) = NULL, 
		  void (*mouse)(Window *w, int button, int state, int xy[]) = NULL,
		  void (*close)(Window *w) = NULL) {
    closeWindow();  removeImage();
    setCallbacks(draw, keys, mouse, close);
    createBasicWindow( w, h, title, false );  // create the window and other widgets
    ShowWindow( window );
    SelectWindow( window );
  }
  void redraw(void) {
    Rect rect; 
    GetWindowBounds(window, kWindowContentRgn, &rect); 
    rect.left = 0;   rect.right = rect.right - rect.left;
    rect.top = 0;    rect.bottom = rect.bottom - rect.top;
    InvalWindowRect(window, &rect); 
    if (imgw > 0 && imgh > 0 && imgdata) {
      if (imgref) CGImageRelease( imgref );   imgref = NULL;
      updateImageRef();  // set up the CGImageRef 'imgref' with 'imgdata'
    }
  }
  // image -----------------------------------------------------------
  int  getPixelSize(pixel_t type) {
    switch (type) {	// Note that the unit is 'sizeof(unsigned char)'
    case PIXEL_GRAY:  return 1;  break;
    case PIXEL_RGB:   
    case PIXEL_BGR:   
    case PIXEL_YUV:   return 3;  break;
    case PIXEL_RGBA:  return 4;  break;
    case PIXEL_GRAYA: return 2;  break;
    default: return 0;
    }
    return 0;
  }
  void clearImage(unsigned char v=0) { 
    if (imgdata) memset(imgdata, v, imgw * imgh * getPixelSize(imgtype)); 
    redraw();
  }
  void removeImage(void) { 
    if (imgref) CGImageRelease( imgref );   imgref = NULL;
    imgw = imgh = 0;  
    if (img_allocated && imgdata) free(imgdata);
    imgdata = NULL;
  }
  bool createImage(int w, int h, pixel_t type=PIXEL_RGB) {
    // Create a new image by allocating memory (only if it's necessary).
    if (w <= 0 || h <= 0) return false;
    if (type != PIXEL_GRAY && type != PIXEL_GRAYA &&
	type != PIXEL_RGB  && type != PIXEL_RGBA) {
      fprintf(stderr, "Error (GUIH::ImageWindow::createImage): invalid pixel format\n"); 
      return false; 
    }
    if (img_allocated && w==imgw && h==imgh && type==imgtype) {
    } else {
      removeImage();
      imgw = w;  imgh = h;  imgtype = type;
      int  pxsize = getPixelSize(type);
      imgdata = (unsigned char*)malloc( w * h * pxsize * sizeof(unsigned char) );
      img_allocated = true;
    }
    updateImageRef();  // set up the CGImageRef 'imgref' with 'imgdata'
    resize( w, h, true );
    return true;
  }
  bool showImage( int w, int h, pixel_t type, void *data) {
    removeImage();
    if (w <= 0 || h <= 0 || data == NULL) { return false; }
    if (type != PIXEL_GRAY && type != PIXEL_GRAYA &&
	type != PIXEL_RGB  && type != PIXEL_RGBA) {
      fprintf(stderr, "Error (GUIH::ImageWindow::showImage): invalid pixel format\n"); 
      return false; 
    }
    imgw = w;  imgh = h;  imgtype = type;
    imgdata = (unsigned char*)data;
    img_allocated = false;
    updateImageRef();  // set up the CGImageRef 'imgref' with 'imgdata'
    resize( w, h, true );
    return true;
  }
  bool showImage( const char *fname ) { 
    // Create an image source using a URL of the pathname
    CFURLRef url = CFURLCreateFromFileSystemRepresentation( NULL, (const UInt8*)fname, strlen(fname), false );
    CGImageSourceRef isrc = CGImageSourceCreateWithURL(url, NULL);
    CFRelease( url );    if (isrc == NULL) { removeImage(); return false; }
    // Create an CGImage from the first item in the image source.
    CGImageRef imgread = CGImageSourceCreateImageAtIndex(isrc, 0, NULL);
    CFRelease(isrc);  if (imgread == NULL) { removeImage(); return false; }
    removeImage();
    // allocate memory for the new iamge
    imgw = CGImageGetWidth( imgread );
    imgh = CGImageGetHeight( imgread );
    imgtype = PIXEL_RGBA;
    int pxsize = getPixelSize(imgtype);
    imgdata = (unsigned char*)malloc( imgw * imgh * pxsize );
    // copy pixel values using CGBitmapContext
    //   for possible pixel format, refer to http://developer.apple.com/qa/qa2001/qa1037.html
    CGColorSpaceRef cs = CGColorSpaceCreateWithName( kCGColorSpaceGenericRGB );
    CGContextRef bc = CGBitmapContextCreate( imgdata, imgw, imgh, 8, pxsize*imgw, cs, kCGImageAlphaNoneSkipLast );
    CGRect rect = CGRectMake( 0, 0, imgw, imgh );
    CGContextDrawImage( bc, rect, imgread );  // convert pixel format
    CGColorSpaceRelease( cs );
    CGContextRelease( bc );
    CGImageRelease( imgread );
    updateImageRef();  // set up the CGImageRef 'imgref' with 'imgdata'
    resize( imgw, imgh, true ); 
    if (!imgref) return false;
    if (!window) openWindow( imgw, imgh, (char*)fname );
    else ShowWindow( window );
    return true;
  }
  bool saveImage(char *fname) { return captureScreen( fname ); }
  bool copyImage(int w, int h, unsigned char *data, pixel_t type=PIXEL_RGB) { 
    createImage( w, h, type );
    memcpy( imgdata, data, w * h * getPixelSize(type) );
    img_allocated = true;
    return true;
  }
private:
  void updateImageRef(void) {
    int  pxsize = getPixelSize(imgtype);
    CGColorSpaceRef cs=NULL;
    if      (imgtype==PIXEL_GRAY)  cs = CGColorSpaceCreateWithName( kCGColorSpaceGenericGray );
    else if (imgtype==PIXEL_GRAYA) cs = CGColorSpaceCreateWithName( kCGColorSpaceGenericGray );
    else if (imgtype==PIXEL_RGB)   cs = CGColorSpaceCreateWithName( kCGColorSpaceGenericRGB );
    else if (imgtype==PIXEL_RGBA)  cs = CGColorSpaceCreateWithName( kCGColorSpaceGenericRGB );
    else return;
    CGDataProviderRef dp = CGDataProviderCreateWithData( NULL, imgdata, pxsize*imgw*imgh, NULL );
    imgref = CGImageCreate( imgw, imgh, 8, 8*pxsize, pxsize*imgw, cs, kCGBitmapByteOrderDefault, 
			    dp, NULL, false, kCGRenderingIntentDefault );
    CGColorSpaceRelease( cs );  CGDataProviderRelease( dp );
  }
  
  bool captureScreenInRGBBuffer(unsigned char *buffer) { 
    // Capture the screen in RGB pixel format
    int i, total, w, h;
    getSize(&w, &h, true);  total = w * h;
    unsigned char *src=imgdata, *dst=buffer;
    if (!imgref) return false;
    switch (imgtype) {
    case PIXEL_GRAY:  for (i=0; i<total; i++, src+=1, dst+=3) dst[0]=dst[1]=dst[2]=src[0];  break;
    case PIXEL_GRAYA: for (i=0; i<total; i++, src+=2, dst+=3) dst[0]=dst[1]=dst[2]=src[0];  break;
    case PIXEL_RGB:   for (i=0; i<total; i++, src+=3, dst+=3) memcpy(dst, src, 3*sizeof(unsigned char));  break;
    case PIXEL_RGBA:  for (i=0; i<total; i++, src+=4, dst+=3) memcpy(dst, src, 3*sizeof(unsigned char));  break;
    default : return false;
    }
    return true;
  }

  // -----------------------------------------------------------------
  // callback functions
  // -----------------------------------------------------------------
public:
  bool runOnInit(WindowRef w) { if (cb_init) cb_init( this ); return true; }
  void runOnDestroy(void) { 
    // This function is called automatically by a event, and cannot be cancelled
    if (cb_close) cb_close(this); 
    removeImage();
  }
  bool runOnDraw(void) {
    // Note that CGContext is given as 'windowContext'
    if (!windowContext) return false;
    if (imgref) {
      CGRect rect = CGRectMake( 0, 0, imgw, imgh );
      // CGContextSetGrayFillColor( windowContext, 0.5, 1.0 );
      CGContextSetRGBFillColor( windowContext, 0.5, 0.5, 0.0, 1.0 );
      CGContextFillRect( windowContext, rect );
      CGContextDrawImage( windowContext, rect, imgref );
    } else {
      Rect   tmp;  GetWindowBounds(window, kWindowContentRgn, &tmp); 
      CGRect rect = CGRectMake( 0, 0, (tmp.right-tmp.left), (tmp.bottom-tmp.top) );
      //CGContextClearRect( windowContext, rect );
    }
    if (cb_draw) cb_draw( this );
    return true;
  }
  void runOnKeyboard(int key) {
    bool skey = (key==ESCAPE_KEY || key>200);
    // Special Keys with default behaviors
    if (!skey_by_user) {
      if (key == HOME_KEY) {
	char *ptype;   
	switch (imgtype) {
	case PIXEL_GRAY:  ptype = "GRAY ";  break;
	case PIXEL_RGB:   ptype = "RGB  ";  break;
	case PIXEL_BGR:   ptype = "BGR  ";  break;
	case PIXEL_YUV:   ptype = "YUV  ";  break;
	case PIXEL_GRAYA: ptype = "GRAYA";  break;
	case PIXEL_RGBA:  ptype = "RGBA ";  break;
	default: ptype = "---- "; break;
	}
	printf("GUIH::ImageWindow image info: (%d x %d) %s(%d)  row_size=%d \n", 
	       imgw, imgh, ptype, getPixelSize(imgtype), imgw * (imgtype == PIXEL_GRAY ? 1 : 3) );
      }
    }
    if (key == F1 && shift) {
      std::cerr << "---------------------------------------------------------------------" << std::endl;
      std::cerr << "  GUIH::ImageWindow control keys                                     " << std::endl;
      std::cerr << "    Home   (Cmd-o) : Show image information                         " << std::endl;
      std::cerr << "    Insert (Cmd-i) : capture the image of the window                " << std::endl;
      std::cerr << "    Delete         : clear caption                                  " << std::endl;
      std::cerr << "---------------------------------------------------------------------" << std::endl;
    }
    // run user key press callback function
    //printf("key=%d  skey=%d  skey_by_user=%d  ctrl=%d alt=%d\n", key, skey, skey_by_user, ctrl, alt);
    if (cb_keys && (!skey||skey_by_user)) cb_keys( this, key );
  }
  void runOnMouse(int button, int state, int xy[]) { 
    if (button == MOUSE_LEFT && state == MOUSE_CLICKED) {
      char  str[80];  // default behavior: show the clicked position
      sprintf( str, "(%d, %d)", xy[0], xy[1]);
      setCaption( str );
    }
    if (cb_mouse) cb_mouse( this, button, state, xy ); 
  }
};

}	// end of GUIH namespace

// ===================================================================
#else		// beginning of Linux/Unix version
// ===================================================================

namespace GUIH {

class ImageWindow : public Window
{
public:
  int		imgw, imgh;	// image width and height
  pixel_t	imgtype;	// pixel format (GRAY, RGB, BGR, YUV)
  unsigned char*imgdata;	// pixel values
  bool		img_allocated;	// whether or not the pixel buffer was allocated
  bool		ignore_alpha;	// whether or not to ignore alpha channel when draw the image
public:
  ImageWindow () 
    : imgw(0), imgh(0), imgdata(NULL), img_allocated(false), ignore_alpha(false) { }
  ImageWindow (int w, int h, char *title = NULL, 
	       void (*init) (Window *w) = NULL,
	       void (*draw) (Window *w) = NULL,
	       void (*keys) (Window *w, int k) = NULL, 
	       void (*mouse)(Window *w, int button, int state, int xy[]) = NULL,
	       void (*close)(Window *w) = NULL) 
    : imgw(0), imgh(0), imgdata(NULL), img_allocated(false), ignore_alpha(false) {
    openWindow(w, h, title, init, draw, keys, mouse, close); 
  }
  virtual ~ImageWindow() { if (img_allocated && imgdata) free(imgdata); }
public:
  void openWindow(int w, int h, const char *title = NULL, 
		  void (*init) (Window *w) = NULL,
		  void (*draw) (Window *w) = NULL,
		  void (*keys) (Window *w, int k) = NULL, 
		  void (*mouse)(Window *w, int button, int state, int xy[]) = NULL,
		  void (*close)(Window *w) = NULL) {
    closeWindow();  removeImage();
    createBasicWidgets( w, h, title );	// create the window and other widgets
    RegisterBasicEvents( init );	// setup default event handlers
    setCallbacks(draw, keys, mouse, close);
    gtk_widget_show( window );		// show the window
  }
  // image -----------------------------------------------------------
  int  getPixelSize(pixel_t type) {
    switch (type) {	// Note that the unit is 'sizeof(unsigned char)'
    case PIXEL_GRAY:  return 1;  break;
    case PIXEL_RGB:   
    case PIXEL_BGR:   
    case PIXEL_YUV:   return 3;  break;
    case PIXEL_RGBA:  return 4;  break;
    case PIXEL_GRAYA: return 2;  break;
    default: return 0;
    }
    return 0;
  }
  void clearImage(unsigned char v=0) { if (imgdata) memset(imgdata, v, imgw * imgh * getPixelSize(imgtype)); }
  void removeImage(void) { 
    imgw = imgh = 0;  
    if (img_allocated && imgdata) free(imgdata);
    imgdata = NULL;
  }
  bool createImage(int w, int h, pixel_t type=PIXEL_RGB) {
    // Create a new image by allocating memory (only if it's necessary).
    if (w <= 0 || h <= 0) return false;
    if (img_allocated && w==imgw && h==imgh && type==imgtype) return true;
    if (img_allocated && imgdata) free(imgdata);
    imgw = w;  imgh = h;  imgtype = type;
    int  pixel_size = getPixelSize(type);
    imgdata = (unsigned char*)malloc( w * h * pixel_size * sizeof(unsigned char) );
    img_allocated = true;
    resize( w, h, true );
    return true;
  }
  bool showImage( int w, int h, pixel_t type, void *data) {
    if (w <= 0 || h <= 0 || data == NULL) { removeImage(); return false; }
//     if (this->imgw == w && this->imgtype == type && 
// 	this->imgh == h && this->imgdata == data) return true; 
    if (type <= PIXEL_UNKNOWN || type > PIXEL_YUV) {
      fprintf(stderr, "Error (GUIH::ImageWindow::showImage): invalid pixel format\n"); 
      removeImage();
      return false; 
    }
    if (img_allocated && imgdata) free(imgdata);
    imgw = w;  imgh = h;  imgtype = type;
    imgdata = (unsigned char*)data;
    img_allocated = false;
    resize( w, h );
    return true;
  }
  bool showImage( char *fname ) { 
    if (!fname) { removeImage(); return false; }
    GdkPixbuf *pbuf = gdk_pixbuf_new_from_file( fname, NULL );  // RGB
    if (pbuf == NULL) { removeImage(); return false; }
    if (imgdata && img_allocated) { free(imgdata); imgdata = NULL; img_allocated = false; }
    imgw = gdk_pixbuf_get_width( pbuf );
    imgh = gdk_pixbuf_get_height( pbuf );
    int   row, rowstride = gdk_pixbuf_get_rowstride( pbuf );
    char *src = (char*)gdk_pixbuf_get_pixels( pbuf ), *dst;
    imgtype = PIXEL_RGB;
    imgdata = (unsigned char*)malloc( imgw * imgh * sizeof(unsigned char) * 3 );
    for (row = 0, dst = (char*)imgdata; row < imgh; row++) {
      memcpy( dst, src, rowstride );
      src += rowstride;  dst += imgw * sizeof(unsigned char) * 3;
    }
    img_allocated = true;
    gdk_pixbuf_unref( pbuf );
    if (!window) openWindow( imgw, imgh, fname );
    else { gtk_widget_show( window );  resize( imgw, imgh, true ); }
    return true;
  }
  bool saveImage(char *fname) { return captureScreen( fname ); }
  bool copyImage(int w, int h, unsigned char *data, pixel_t type=PIXEL_RGB) { 
    createImage( w, h, type );
    memcpy( imgdata, data, w * h * getPixelSize(type) );
    img_allocated = true;
    return true;
  }
  
  // -----------------------------------------------------------------
  // callback functions
  // -----------------------------------------------------------------

  bool runOnInit(GtkWidget *widget) { if (cb_init) cb_init( this ); return true; }
  void runOnDestroy(void) { 
    // This function is called automatically by a event, and cannot be cancelled
    if (cb_close) cb_close(this); 
    if (img_allocated && imgdata) free(imgdata);  imgdata=NULL; 
  }
  bool runOnDraw(GtkWidget *widget, GdkEventExpose *event) {
    unsigned char *p = (unsigned char*)imgdata;
    if (p && imgw > 0 && imgh > 0) {
      GdkColor color; 
      switch (imgtype) {
      case PIXEL_GRAY:
	gdk_draw_gray_image( area->window, area->style->fg_gc[GTK_STATE_NORMAL],
			     0, 0, imgw, imgh, GDK_RGB_DITHER_NONE,
			     (guchar*)imgdata, imgw * sizeof(unsigned char) );
	break;
      case PIXEL_RGB:
	gdk_draw_rgb_image ( area->window, area->style->fg_gc[GTK_STATE_NORMAL],
			     0, 0, imgw, imgh, GDK_RGB_DITHER_NONE,
			     (guchar*)imgdata, imgw * 3 * sizeof(unsigned char) );
	break;
      case PIXEL_GRAYA:
	if (ignore_alpha) {
	  for (int y = 0; y < imgh; y++) 
	    for (int x = 0; x < imgw; x++, p+=2) {
	      color.red   = (guint)(65535 * p[0] / 255);
	      color.green = (guint)(65535 * p[0] / 255);
	      color.blue  = (guint)(65535 * p[0] / 255);
	      gdk_gc_set_rgb_fg_color( area->style->fg_gc[GTK_STATE_NORMAL], &color); 
	      gdk_draw_point( area->window, area->style->fg_gc[GTK_STATE_NORMAL], x, y);
	    }
	} else {
	  for (int y = 0; y < imgh; y++) 
	    for (int x = 0; x < imgw; x++, p+=2) {
	      if (p[1] == 0) 
		color.red = color.green = color.blue = (guint)(65535 / 5);
	      else
		color.red = color.green = color.blue = (guint)(65535 * p[0]/255 * p[1]/255);
	      gdk_gc_set_rgb_fg_color( area->style->fg_gc[GTK_STATE_NORMAL], &color); 
	      gdk_draw_point( area->window, area->style->fg_gc[GTK_STATE_NORMAL], x, y);
	    }
	}
	break;
      case PIXEL_RGBA:
	if (ignore_alpha) {
	  for (int y = 0; y < imgh; y++) 
	    for (int x = 0; x < imgw; x++, p+=4) {
	      color.red   = (guint)(65535 * p[0] / 255);
	      color.green = (guint)(65535 * p[1] / 255);
	      color.blue  = (guint)(65535 * p[2] / 255);
	      gdk_gc_set_rgb_fg_color( area->style->fg_gc[GTK_STATE_NORMAL], &color); 
	      gdk_draw_point( area->window, area->style->fg_gc[GTK_STATE_NORMAL], x, y);
	    }
	} else {
	  for (int y = 0; y < imgh; y++) 
	    for (int x = 0; x < imgw; x++, p+=4) {
	      if (p[3] == 0) {
		color.red = color.green = color.blue = (guint)(65535 / 5);
	      } else {
		color.red   = (guint)(65535 * p[0] / 255 * p[3] / 255);
		color.green = (guint)(65535 * p[1] / 255 * p[3] / 255);
		color.blue  = (guint)(65535 * p[2] / 255 * p[3] / 255);
	      }
	      gdk_gc_set_rgb_fg_color( area->style->fg_gc[GTK_STATE_NORMAL], &color); 
	      gdk_draw_point( area->window, area->style->fg_gc[GTK_STATE_NORMAL], x, y);
	    }
	}
	break;
      default: std::cerr << "Error (GUIH::ImageWindow::runOnDraw): invalid pixel type" << std::endl;  break;
      }
    }
    if (cb_draw) cb_draw( this );
    return true;
  }
  void runOnKeyboard(int key) {
    bool skey = (key==ESCAPE_KEY || key>200);
    // Special Keys with default behaviors
    if (!skey_by_user) {
      if (key == HOME_KEY) {
	char *ptype;   
	switch (imgtype) {
	case PIXEL_GRAY:  ptype = "GRAY ";  break;
	case PIXEL_RGB:   ptype = "RGB  ";  break;
	case PIXEL_BGR:   ptype = "BGR  ";  break;
	case PIXEL_YUV:   ptype = "YUV  ";  break;
	case PIXEL_GRAYA: ptype = "GRAYA";  break;
	case PIXEL_RGBA:  ptype = "RGBA ";  break;
	default: ptype = "---- "; break;
	}
	printf("GUIH::ImageWindow image info: (%d x %d) %s(%d)  row_size=%d \n", 
	       imgw, imgh, ptype, getPixelSize(imgtype), imgw * (imgtype == PIXEL_GRAY ? 1 : 3) );
      }
    }
    if (key == F1 && shift) {
      std::cerr << "---------------------------------------------------------------------" << std::endl;
      std::cerr << "  GUIH::ImageWindow control keys                                     " << std::endl;
      std::cerr << "    Home     : Show image information                                " << std::endl;
      std::cerr << "    Insert   : capture the image of the window                       " << std::endl;
      std::cerr << "    Delete   : clear caption                                         " << std::endl;
      std::cerr << "---------------------------------------------------------------------" << std::endl;
    }
    // run user key press callback function
    //printf("key=%d  skey=%d  skey_by_user=%d  ctrl=%d alt=%d\n", key, skey, skey_by_user, ctrl, alt);
    if (cb_keys && (!skey||skey_by_user)) cb_keys( this, key );
  }
  void runOnMouse(int button, int state, int xy[]) { 
    if (button == MOUSE_LEFT && state == MOUSE_CLICKED) {
      char  str[80];  // default behavior: show the clicked position
      sprintf( str, "(%d, %d)", xy[0], xy[1]);
      setCaption( str );
    }
    if (cb_mouse) cb_mouse( this, button, state, xy ); 
  }
};

}	// end of GUIH namespace

// ===================================================================
#endif	// end of linux version
// ===================================================================

#endif	// GUIH_IMAGE_HPP



// ===================================================================
#if 0	// Example code starts
// ===================================================================
// For common functions, refer to GUIH::Window in "guih_common.hpp".
// compile on Linux:  g++ -o guih guih_example.cpp `pkg-config --cflags gtkglext-1.0 --libs gtkglext-1.0` -Wl,-rpath,/usr/local/lib
#include <iostream>
#include "guih_image.hpp"
void draw(GUIH::Window *w) 
{
  w->setColorForDrawing( 1.0, 0.0, 0.0 );
  w->drawRectangle( 50, 50, 200, 100 );
}
void keys(GUIH::Window *w, int key)
{
  printf("Key pressed with ctrl=%d, alt=%d and shift=%d\n", w->ctrl, w->alt, w->shift);
  switch (key) {
  case 'a':  std::cout << 'a' << std::endl;  break;
  case ' ':  std::cout << "Space" << std::endl;  break;
  case GUIH::F1:  std::cout << "F1 function key" << std::endl;  break;
  case GUIH::ARROW_UP:  std::cout << "arrow up" << std::endl;  break;
  case GUIH::ARROW_LEFT:  std::cout << "arrow left" << std::endl;  break;
  }
}
void mouse(GUIH::Window *w, int button, int state, int xy[])
{
  if        (button==GUIH::MOUSE_LEFT && state==GUIH::MOUSE_CLICKED) {
    printf(" clicked  at (%d,%d)\n", xy[0], xy[1] );
  } else if (button==GUIH::MOUSE_LEFT && state==GUIH::MOUSE_DRAGGED) {
    int mpos, mcnt;  for(mpos=0; xy[mpos]>=0; mpos+=2); mcnt = mpos/2;
    printf(" dragged from (%d,%d) to (%d,%d)\n", xy[0], xy[1], xy[mcnt*2-2], xy[mcnt*2-1]);
  }
}
int  main(void)
{
  // You can create windows as many as you want.
  GUIH::ImageWindow w0(300, 300, "Untitled Window", NULL, draw, keys, mouse);
  GUIH::ImageWindow w1(300, 300, "Untitled Window", NULL, draw, keys, mouse);
  w0.showImage("sample.jpg");
  w1.createImage(320, 240);  w1.resize(320, 240);
  // start the main loop -- press ESC to quit
  w0.runMainLoop();
  return EXIT_SUCCESS;
}
// ===================================================================
#endif	// Example code ends
// ===================================================================
