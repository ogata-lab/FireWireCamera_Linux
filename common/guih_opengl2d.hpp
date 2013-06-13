
//
// GUIH                 : GUI in Header files
// GUIH::OpenGL2DWindow : window for 2D OpenGL graphics
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
// GUIH::OpenGL2DWindow
//   This class provides very convenient GUI for 2D and 3D OpenGL graphics.
//   All the details of setting up OpenGL are hidden, and 2D/3D cameras with
//   key bindings for the navigation are also provided by default.
//
//   Member functions:
//	bool showImage( char *fname, bool grayscale );
//	void initCamera( xmin, ymin, xmax, ymax );
//	void initCameraForBBox( xmin, ymin, xmax, ymax );
//      void drawLine(float x0, float y0, float x1, float y1);
//      void drawRectangle(float x, float y, float width, float height);
//
// For API example, read the example code at the bottom.
// For default key bindings, press Shift+F1 during execution.
//


#ifndef GUIH_OPENGL2D_HPP
#define GUIH_OPENGL2D_HPP

#ifdef WIN32					// Microsoft Windows
#include <windows.h>
#include <cmath>
#include <GL/gl.h>
#include <GL/glu.h>
#elif defined(__APPLE__) || defined(MACOSX)	// Apple OS X
#include <AGL/agl.h>
#include <OpenGL/gl.h>
#include <OpenGL/glu.h>
#include <GLUT/glut.h>
#else						// Linux
#include <GL/gl.h>
#include <GL/glu.h>
#endif
#include "guih_common.hpp"
#include "guih_camera2d.hpp"

#ifndef M_PI
#define M_PI (3.14159265358979323846)
#endif


// ===================================================================
#ifdef  WIN32   // beginning of Window version
// ===================================================================

namespace GUIH {
  

class OpenGL2DWindow : public Window
{
public:
  int   txtw, txth, txtid;	// textured image
  unsigned char	*imgdata;
  pixel_t	imgtype;
  HDC       hDC;
  HGLRC		hRC;
  bool		show_axes;
  Camera2D	camera;
  float		bkgrgb[3];
public:
  OpenGL2DWindow () : imgdata(NULL), txtid(0) { 
    show_axes = false;  skey_by_user = false;
    initCamera( 0, 0, 1, 1 );		// set camera
    bkgrgb[0] = bkgrgb[1] = bkgrgb[2] = 1.0;
  }
  OpenGL2DWindow (int w, int h, char *title = NULL, 
		  void (*init) (Window *w) = NULL,
		  void (*draw) (Window *w) = NULL,
		  void (*keys) (Window *w, int k) = NULL, 
		  void (*mouse)(Window *w, int button, int state, int xy[]) = NULL,
		  void (*close)(Window *w) = NULL) : imgdata(NULL), txtid(0) {
    show_axes = false;  skey_by_user = false;
    initCamera( 0.0f, 0.0f, (float)w, (float)h );		// set camera
    bkgrgb[0] = bkgrgb[1] = bkgrgb[2] = 1.0;
    openWindow(w, h, title, init, draw, keys, mouse, close); 
  }
  ~OpenGL2DWindow() { }
public:
  void openWindow(int w, int h, const char *title = NULL, 
		  void (*init) (Window *w) = NULL,
		  void (*draw) (Window *w) = NULL,
		  void (*keys) (Window *w, int k) = NULL, 
		  void (*mouse)(Window *w, int button, int state, int xy[]) = NULL,
		  void (*close)(Window *w) = NULL) {
    closeWindow();
    setCallbacks(draw, keys, mouse, close);
    createBasicObjects( w, h, title );
  }
  void setBackgroundColor(float r, float g, float b) {
    if (!hDC || !hRC) return;
    bkgrgb[0] = r;  bkgrgb[1] = g;  bkgrgb[2] = b;
    wglMakeCurrent( hDC, hRC );
    glClearColor (bkgrgb[0], bkgrgb[1], bkgrgb[2], 0); // clearing color
    redraw();
  }
  void initCamera(float xmin, float ymin, float xmax, float ymax) {
    // Set camera (Note that it can be scaled, depending on window size)
    camera.initPose( xmin, ymin, xmax, ymax );
    redraw();
  }
  void initCameraForBBox(float xmin, float ymin, float xmax, float ymax, float margin=0.0f) {
    // Set camera (taking into account window size and keeping appropriate scale)
    float sx = (xmax-xmin)*(1+margin),   sy = (ymax-ymin)*(1+margin);
    float cx = (xmax+xmin)/2, cy = (ymax+ymin)/2;
    int   w=0, h=0;   getSize(&w, &h, true);
    if (w<=0 || h<=0) fprintf(stderr, "Warning : initCameraForBBox() called before the window is set\n");
    if (sx/w > sy/h)  sy = sx * h / w;  else  sx = sy * w / h;
    camera.initPose( cx-sx/2, cy-sy/2, cx+sx/2, cy+sy/2 );
    redraw();
  }
  void initCameraForImage(int width, int height) {
    // Set camera (fitting entire image to the window)
    camera.initPoseForImage( width, height );
    redraw();
  }
  void setCamera(float xmin, float ymin, float xmax, float ymax) {
    camera.xmin = xmin;  camera.ymin = ymin;
    camera.xmax = xmax;  camera.ymax = ymax;
    redraw();
  }
  void shiftCamera(float dx, float dy) { camera.shift( dx, dy ); redraw(); }
  
  // image -----------------------------------------------------------
  bool showImage( const char *fname ) { 
//     if (!hwnd || !fname) return false;
//     IplImage *img = cvLoadImage( fname, 1 );
//     if (!img) return false;
//     bool ret = showImage(img); 
//     cvReleaseImage( &img );
//     return ret;
  }
  bool showImage( int w, int h, pixel_t type, void *data, bool use_texture=false ) { 
//     if (!hwnd || !arr) return false;
//     IplImage *src, stub;  src = cvGetImage(arr, &stub);
//     if (src->depth != IPL_DEPTH_8U || src->nChannels != 3) return false;
//     resize( src->width, src->height, true );
//     initCamera( 0, 0, src->width, src->height );
//     // create an OpenGL texture for the image
//     if (txtid == 0) txtid = serial_number+1;
//     txtw = src->width;  txth = src->height;
//     wglMakeCurrent( hDC, hRC );
//     glEnable(GL_TEXTURE_2D);
//     glBindTexture(GL_TEXTURE_2D, txtid);
//     glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);  // GL_LINEAR
//     glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);  // GL_LINEAR
//     glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_DECAL); // GL_MODULATE);
//     if (strcmp(src->channelSeq, "RGB")==0)
//       glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, src->width, src->height, 0, GL_RGB, GL_UNSIGNED_BYTE, src->imageData);
//     else
//       glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, src->width, src->height, 0, GL_BGR_EXT, GL_UNSIGNED_BYTE, src->imageData);
//     glDisable(GL_TEXTURE_2D);
//     redraw();
    return true;
  }
  void removeImage(void) { imgdata = NULL; txtid = 0; }

  // -----------------------------------------------------------------
  // callback functions
  // -----------------------------------------------------------------

  bool runOnInit(void) { 
    // enable OpenGL on the window
    PIXELFORMATDESCRIPTOR pfd;
    ZeroMemory( &pfd, sizeof( pfd ) );	// set the pixel format for the DC
    pfd.nSize = sizeof( pfd );
    pfd.nVersion = 1;
    pfd.dwFlags = PFD_DRAW_TO_WINDOW | PFD_SUPPORT_OPENGL | PFD_DOUBLEBUFFER;
    pfd.iPixelType = PFD_TYPE_RGBA;
    pfd.cColorBits = 24;
    pfd.cDepthBits = 16;
    //pfd.cStencilBits = 0;
    //pfd.cAuxBuffers = 0;
    pfd.iLayerType = PFD_MAIN_PLANE;
    //pfd.bReserved = 0;
    //pfd.dwLayerMask = 0;
    //pfd.dwVisibleMask = 0;
    //pfd.dwDamageMask = 0;
    hDC = GetDC(hwnd);
    int format = ChoosePixelFormat( hDC, &pfd );
    SetPixelFormat( hDC, format, &pfd );
    hRC = wglCreateContext( hDC );	// create and enable the render context (RC)
    wglMakeCurrent( hDC, hRC );    // setup OpenGL
    glClearColor (0.0, 0.0, 0.0, 0.0);	// clearing color
    glColor3ub(255,255,255);
    if (cb_init) cb_init(this);
    return true;
  }
  void runOnDestroy(void) { 
    // This function is called automatically by a event, and cannot be cancelled
    if (cb_close) cb_close(this); 
    // disable OpenGL
    wglMakeCurrent( hDC, NULL );
    wglDeleteContext( hRC );
    ReleaseDC( hwnd, hDC );
  }
  bool runOnDraw(void) {
    wglMakeCurrent( hDC, hRC );
    Camera2D *cam = &(camera);
    glClear(GL_COLOR_BUFFER_BIT);			// clear
    glMatrixMode(GL_PROJECTION);   glLoadIdentity();
    gluOrtho2D(cam->xmin, cam->xmax, cam->ymin, cam->ymax);
    glMatrixMode(GL_MODELVIEW);
    if (show_axes)  cam->renderXYAxes();		// render axes
    glColor3ub(100, 100, 100);
    if (cb_draw) cb_draw( this );			// user callback
	else {				// example
	  RECT rect;
      GetClientRect( hwnd, &rect );
      glBegin(GL_LINE_LOOP);
      glVertex2f(rect.right*0.8f, rect.bottom*0.6f);
      glVertex2f(rect.right*0.2f, rect.bottom*0.6f);
      glVertex2f(rect.right*0.7f, rect.bottom*0.2f);
      glVertex2f(rect.right*0.5f, rect.bottom*0.82f);
      glVertex2f(rect.right*0.3f, rect.bottom*0.2f);
      glEnd();
    }
    glFinish();   SwapBuffers(hDC);
	return true;
  }
  
  void runOnKeyboard(int key);
  void runOnMouse(int button, int state, int xy[]);
  
  void getPositionFromMouseClick(int mxy[2], float wxy[2]) {
    RECT rect; GetWindowRect( hwnd, &rect );
    wxy[0] = camera.xmin + (camera.xmax - camera.xmin) * mxy[0] / abs(rect.right - rect.left);
    wxy[1] = camera.ymax - (camera.ymax - camera.ymin) * mxy[1] / abs(rect.top - rect.bottom);
  }
  
  void drawLine(float x0, float y0, float x1, float y1);
  void drawCircle(double x, double y, double radius, int slices=0, bool filled=false);
  void drawRectangle(double x, double y, double width, double height);
};

}	// end of GUIH namespace

// ===================================================================
#elif defined(__APPLE__) || defined(MACOSX)	// beginning of OS X version
// ===================================================================

namespace GUIH {

class OpenGL2DWindow : public Window
{
public:  // for OpenGL window only
  AGLContext  aglContext;
  int   txtw, txth, txtid;	// textured image
  unsigned char	*imgdata;
  pixel_t	imgtype;
  Camera2D	camera;
  bool		show_axes;
  bool		use_texture_for_image;
  float		bkgrgb[3];
public:
  OpenGL2DWindow () : aglContext(NULL), txtid(0), imgdata(NULL) { 
    show_axes = false;  skey_by_user = false;
    initCamera( 0, 0, 1, 1 );  	// set the camera
    bkgrgb[0] = bkgrgb[1] = bkgrgb[2] = 1.0;
  }
  OpenGL2DWindow (int w, int h, char *title = NULL, 
		  void (*init) (Window *w) = NULL,
		  void (*draw) (Window *w) = NULL,
		  void (*keys) (Window *w, int k) = NULL, 
		  void (*mouse)(Window *w, int button, int state, int xy[]) = NULL,
		  void (*close)(Window *w) = NULL) 
    : aglContext(NULL), txtid(0), imgdata(NULL) {
    show_axes = false;  skey_by_user = false;
    initCamera( 0, 0, w, h );  	// set the camera
    bkgrgb[0] = bkgrgb[1] = bkgrgb[2] = 1.0;
    openWindow(w, h, title, init, draw, keys, mouse, close); 
  }
  ~OpenGL2DWindow() {}
  void initOpenGL2D(void) {
    // This function must be called after createBasicWindow()
    // setup OpenGL context using AGL
    GLint attributes[] = { AGL_RGBA, AGL_DOUBLEBUFFER, AGL_DEPTH_SIZE, 24, AGL_NONE };
    // obtain a pixel format object
    AGLPixelFormat myAGLPixelFormat = aglChoosePixelFormat (NULL, 0, attributes);
    // bind the pixel format object to a rendering contex 
    aglContext = aglCreateContext (myAGLPixelFormat, NULL);
    aglDestroyPixelFormat( myAGLPixelFormat );
    // bind the window to the rendering context (attach the rendering context to the Carbon window)
    aglSetWindowRef (aglContext, window);
    // make the rendering context the current context
    aglSetCurrentContext (aglContext);
    //if (aglGetError()!=AGL_NO_ERROR) return (OSStatus)aglGetError();
    glDisable(GL_DEPTH_TEST);
    glDisable(GL_LIGHTING);
    glMatrixMode (GL_PROJECTION);   glLoadIdentity();
    gluOrtho2D(camera.xmin, camera.xmax, camera.ymin, camera.ymax);
    glMatrixMode (GL_MODELVIEW);    glLoadIdentity();
    glClearColor (bkgrgb[0], bkgrgb[1], bkgrgb[2], 0); // clearing color
  }
public:  
  void openWindow(int w, int h, const char *title = NULL,
		  void (*init) (Window *w) = NULL,
		  void (*draw) (Window *w) = NULL,
		  void (*keys) (Window *w, int k) = NULL, 
		  void (*mouse)(Window *w, int button, int state, int xy[]) = NULL,
		  void (*close)(Window *w) = NULL) {
    cb_init = init;
    setCallbacks(draw, keys, mouse, close);
    createBasicWindow( w, h, title );	// create the window first
    initOpenGL2D();			// initialize OpengGL 2D
    ShowWindow( window );		// show the window
    SelectWindow( window );
  }
  void setBackgroundColor(float r, float g, float b) {
    if (!aglContext) return;
    aglSetCurrentContext (aglContext);
    bkgrgb[0] = r;  bkgrgb[1] = g;  bkgrgb[2] = b;
    glClearColor (bkgrgb[0], bkgrgb[1], bkgrgb[2], 0); // clearing color
  }
  void initCamera(float xmin, float ymin, float xmax, float ymax) {
    // Set camera (Note that it can be scaled, depending on window size)
    camera.initPose( xmin, ymin, xmax, ymax );
    redraw();
  }
  void initCameraForBBox(float xmin, float ymin, float xmax, float ymax, float margin=0.0f) {
    // Set camera (taking into account window size and maintaining world scale)
    float sx = (xmax-xmin)*(1+margin),   sy = (ymax-ymin)*(1+margin);
    float cx = (xmax+xmin)/2, cy = (ymax+ymin)/2;
    int   w=0, h=0;   getSize(&w, &h, true);
    if (w<=0 || h<=0) fprintf(stderr, "Warning : initCameraForBBox() called before the window is set\n");
    if (sx/w > sy/h)  sy = sx * h / w;  else  sx = sy * w / h;
    camera.initPose( cx-sx/2, cy-sy/2, cx+sx/2, cy+sy/2 );
    redraw();
  }
  void initCameraForImage(int width, int height) {
    // Set camera (fitting entire image to the window)
    camera.initPoseForImage( width, height );
    redraw();
  }
  void setCamera(float xmin, float ymin, float xmax, float ymax) {
    camera.xmin = xmin;  camera.ymin = ymin;
    camera.xmax = xmax;  camera.ymax = ymax;
    redraw();
  }
  void shiftCamera(float dx, float dy) { camera.shift( dx, dy ); }
  
  // image -----------------------------------------------------------
  bool showImage( const char *fname, bool use_texture=false ) { return false; }
  bool showImage( int w, int h, pixel_t type, void *data, bool use_texture=false ) {
    return false;
  }
  void removeImage(void) { imgdata = NULL; txtid = 0; }

  bool captureScreenInRGBBuffer(unsigned char *buffer) { 
    // Capture the screen in RGB pixel format
    int w, h;  getSize(&w, &h, true);
    if (!aglContext) return false;
    aglSetCurrentContext (aglContext);
    glReadBuffer(GL_FRONT);
    glReadPixels( 0, 0, w, h, GL_RGB, GL_UNSIGNED_BYTE, buffer );
    if (glGetError() != 0) return false;
    // invert the image vertically (it's upside down)
    unsigned char rgb[3], *p0, *p1;
    for (int row = 0; row < h/2; row++) {
      p0 = buffer + row * 3 * w;
      p1 = buffer + (h-1-row) * 3 * w;
      for (int col = 0; col < w; col++, p0+=3, p1+=3) {
	memcpy( rgb, p0, 3*sizeof(unsigned char) );
	memcpy( p0,  p1, 3*sizeof(unsigned char) );
	memcpy( p1, rgb, 3*sizeof(unsigned char) );
      }
    }
    return true;
  }

  void draw_caption(void) {
    // Draw caption text using GLUT.
    //   Note that Window::draw_caption() does not work for OpenGL windows on OS X.
    //   This function shouldn't be called by user.
    if (!window || !aglContext || !caption[0]) return;
    int w, h;  getSize( &w, &h, true );
    glMatrixMode (GL_PROJECTION);
    glPushMatrix();  glLoadIdentity();  gluOrtho2D(0, w, 0, h);
    glColor4f(1.0, 0.0, 0.0, 0.0);
    glMatrixMode(GL_MODELVIEW);
    glPushMatrix();  glLoadIdentity();
#if true
    // fonts:  GLUT_BITMAP_9_BY_15,  GLUT_BITMAP_8_BY_13, 
    //         GLUT_BITMAP_TIMES_ROMAN_10,  GLUT_BITMAP_TIMES_ROMAN_24,  
    //         GLUT_BITMAP_HELVETICA_10,  GLUT_BITMAP_HELVETICA_12,  GLUT_BITMAP_HELVETICA_18 
    glRasterPos2f(10, 10);
    for (int i=0; caption[i]; i++)  glutBitmapCharacter( GLUT_BITMAP_HELVETICA_12, caption[i] );
#else
    // fonts:  GLUT_STROKE_ROMAN, GLUT_STROKE_MONO_ROMAN
    glTranslatef(10, 10, 0.0);
    glScalef( 0.1f, 0.1f, 0.1f);  // stroke scale
    for (int i=0; caption[i]; i++)  glutStrokeCharacter( GLUT_STROKE_ROMAN, caption[i] );
#endif
    glPopMatrix();
    glMatrixMode (GL_PROJECTION); 
    glPopMatrix();
    glMatrixMode (GL_MODELVIEW);
  }
  
  void draw_image(void) {
    if (use_texture_for_image && txtid > 0) {	// render textured image
      glEnable(GL_TEXTURE_2D);
      glBindTexture(GL_TEXTURE_2D, txtid);
      glBegin(GL_QUADS);
      glTexCoord2f(0.0, 0.0);  glVertex2f(  0 , txth);
      glTexCoord2f(1.0, 0.0);  glVertex2f(txtw, txth);
      glTexCoord2f(1.0, 1.0);  glVertex2f(txtw,   0 );
      glTexCoord2f(0.0, 1.0);  glVertex2f(  0 ,   0 );
      glEnd();
      glDisable(GL_TEXTURE_2D);
    } else if (!use_texture_for_image && imgdata) {
      int x, y, y2;  unsigned char *src = imgdata;
      glBegin(GL_QUADS);
      switch (imgtype) {
      case PIXEL_RGB:  
	for (y = 0; y < txth; y++) {
	  y2 = txth - 1 - y;
	  for (x = 0; x < txtw; x++, src+=3) {
	    glColor3ub(src[0], src[1], src[2]);
	    glVertex2i(x  , y2+1);	glVertex2i(x+1, y2+1);
	    glVertex2i(x+1, y2);	glVertex2i(x  , y2);
	  }
	}
	break;
      case PIXEL_GRAY: 
	for (y = 0; y < txth; y++) {
	  y2 = txth - 1 - y;
	  for (x = 0; x < txtw; x++, src+=1) {
	    glColor3ub(src[0], src[0], src[0]);
	    glVertex2i(x  , y2+1);	glVertex2i(x+1, y2+1);
	    glVertex2i(x+1, y2);	glVertex2i(x  , y2);
	  }
	}
	break;
      default: break;
      }
      glEnd();
    }
  }
  
  // -----------------------------------------------------------------
  // callback functions
  // -----------------------------------------------------------------

  bool runOnInit(WindowRef w) { 
    if (!aglContext) return false;
    aglSetCurrentContext (aglContext);
    if (cb_init)  cb_init( this );
    return true;
  }
  void runOnDestroy(void) { 
    // This function is called automatically by a event, and cannot be cancelled
    if (cb_close) cb_close(this); 
    // clean up OpenGL
    if (aglContext) { aglDestroyContext( aglContext ); aglContext = NULL; }
  }
  bool runOnDraw(void) { 
    if (!aglContext) return false;
    aglSetCurrentContext (aglContext);
    glClear(GL_COLOR_BUFFER_BIT);			// clear
    glMatrixMode (GL_PROJECTION);   glLoadIdentity();
    Camera2D *cam = &(camera);
    gluOrtho2D(cam->xmin, cam->xmax, cam->ymin, cam->ymax);
    glMatrixMode (GL_MODELVIEW);
    draw_image();
    if (show_axes)  cam->renderXYAxes();		// render axes
    glColor3ub(100, 100, 100);
    if (cb_draw) cb_draw( this );			// user callback
    draw_caption();
    aglSwapBuffers( aglGetCurrentContext() );
    return true;
  }
  
  void runOnKeyboard(int key);
  void runOnMouse(int button, int state, int xy[]);
  
  void getPositionFromMouseClick(int mxy[2], float wxy[2]) {
    Rect rect; GetWindowBounds(window, kWindowContentRgn, &rect); 
    if (txtid > 0) {
      wxy[0] = camera.xmin + (camera.xmax - camera.xmin) * mxy[0] / (rect.right-rect.left);
      wxy[1] = (camera.ymaxb - camera.ymax) + (camera.ymax - camera.ymin) * mxy[1] / (rect.bottom-rect.top);
    } else {
      wxy[0] = camera.xmin + (camera.xmax - camera.xmin) * mxy[0] / (rect.right-rect.left);
      wxy[1] = camera.ymax - (camera.ymax - camera.ymin) * mxy[1] / (rect.bottom-rect.top);
    }
  }
  
  void drawLine(float x0, float y0, float x1, float y1);
  void drawCircle(double x, double y, double radius, int slices=0, bool filled=false);
  void drawRectangle(double x, double y, double width, double height);
};

}	// end of GUIH namespace

// ===================================================================
#else		// beginning of Linux/Unix version
// ===================================================================

namespace GUIH {

class OpenGL2DWindow : public Window
{
public:  // for OpenGL window only
  int   txtw, txth, txtid;	// textured image
  unsigned char	*imgdata;
  pixel_t	imgtype;
  GdkPixbuf	*pbuf;
  Camera2D	camera;
  bool		show_axes;
  bool		use_texture_for_image;
  float		bkgrgb[3];
public:
  OpenGL2DWindow () : txtid(0), imgdata(NULL), pbuf(NULL) { 
    show_axes = false;  skey_by_user = false;
    initCamera( 0, 0, 1, 1 );  	// set the camera
    bkgrgb[0] = bkgrgb[1] = bkgrgb[2] = 1.0;
  }
  OpenGL2DWindow (int w, int h, char *title = NULL, 
		  void (*init) (Window *w) = NULL,
		  void (*draw) (Window *w) = NULL,
		  void (*keys) (Window *w, int k) = NULL, 
		  void (*mouse)(Window *w, int button, int state, int xy[]) = NULL,
		  void (*close)(Window *w) = NULL) 
    : txtid(0), imgdata(NULL), pbuf(NULL) {
    show_axes = false;  skey_by_user = false;
    initCamera( 0, 0, w, h );  	// set the camera
    bkgrgb[0] = bkgrgb[1] = bkgrgb[2] = 1.0;
    openWindow(w, h, title, init, draw, keys, mouse, close); 
  }
  ~OpenGL2DWindow() { if (pbuf) gdk_pixbuf_unref(pbuf); }
public:  
  void openWindow(int w, int h, const char *title = NULL,
		  void (*init) (Window *w) = NULL,
		  void (*draw) (Window *w) = NULL,
		  void (*keys) (Window *w, int k) = NULL, 
		  void (*mouse)(Window *w, int button, int state, int xy[]) = NULL,
		  void (*close)(Window *w) = NULL) {
    closeWindow();
    createBasicWidgets( w, h, title );	// create the window and other widgets
    
    // setup OpenGL (Note gtk_gl_init() was omitted, but it seems OK.)
    GdkGLConfigMode mode = (GdkGLConfigMode)(GDK_GL_MODE_RGB|GDK_GL_MODE_DOUBLE);
    GdkGLConfig *glconfig = gdk_gl_config_new_by_mode (mode);
    gtk_widget_set_gl_capability (area, glconfig, NULL, TRUE, GDK_GL_RGBA_TYPE);
    
    RegisterBasicEvents( init );	// setup default event handlers
    setCallbacks(draw, keys, mouse, close);
    gtk_widget_show( window );  	// show the window
  }
  void setBackgroundColor(float r, float g, float b) {
    GdkGLContext  *glcontext  = gtk_widget_get_gl_context (area);
    GdkGLDrawable *gldrawable = gtk_widget_get_gl_drawable(area);
    if (!gdk_gl_drawable_gl_begin (gldrawable, glcontext)) return;  // OpenGL begin
    bkgrgb[0] = r;  bkgrgb[1] = g;  bkgrgb[2] = b;
    glClearColor (bkgrgb[0], bkgrgb[1], bkgrgb[2], 0); // clearing color
    gdk_gl_drawable_gl_end (gldrawable);	// OpenGL END
  }
  void initCamera(float xmin, float ymin, float xmax, float ymax) {
    // Set camera (Note that it can be scaled, depending on window size)
    camera.initPose( xmin, ymin, xmax, ymax );
  }
  void initCameraForBBox(float xmin, float ymin, float xmax, float ymax, float margin=0.0f) {
    // Set camera (taking into account window size and maintaining world scale)
    float sx = (xmax-xmin)*(1+margin),   sy = (ymax-ymin)*(1+margin);
    float cx = (xmax+xmin)/2, cy = (ymax+ymin)/2;
    int   w=0, h=0;   getSize(&w, &h, true);
    if (w<=0 || h<=0) fprintf(stderr, "Warning : initCameraForBBox() called before the window is set\n");
    if (sx/w > sy/h)  sy = sx * h / w;  else  sx = sy * w / h;
    camera.initPose( cx-sx/2, cy-sy/2, cx+sx/2, cy+sy/2 );
  }
  void initCameraForImage(int width, int height) {
    // Set camera (fitting entire image to the window)
    camera.initPoseForImage( width, height );
  }
  void setCamera(float xmin, float ymin, float xmax, float ymax) {
    camera.xmin = xmin;  camera.ymin = ymin;
    camera.xmax = xmax;  camera.ymax = ymax;
  }
  void shiftCamera(float dx, float dy) { camera.shift( dx, dy ); }
  
  // image -----------------------------------------------------------
  bool showImage( char *fname, bool use_texture=false ) { 
    if (!window || !fname) return false;
    if (this->pbuf) gdk_pixbuf_unref( pbuf );
    this->pbuf = gdk_pixbuf_new_from_file( fname, NULL );  // RGB
    if (pbuf == NULL) return false;
    return showImage( gdk_pixbuf_get_width(pbuf), gdk_pixbuf_get_height(pbuf), PIXEL_RGB,
		      (unsigned char*)gdk_pixbuf_get_pixels(pbuf), use_texture );
  }
  bool showImage( int w, int h, pixel_t type, void *data, bool use_texture=false ) { 
    if (!window || !data) return false;
    txtw = w;  txth = h;
    resize( w, h, true );  initCamera( 0, 0, w, h );
    if (use_texture) {	// create an OpenGL texture for the image
      if (txtid == 0) txtid = serial_number+1;
      glEnable(GL_TEXTURE_2D);
      glBindTexture(GL_TEXTURE_2D, txtid);
      glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);  // GL_LINEAR
      glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);  // GL_LINEAR
      glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_DECAL); // GL_MODULATE);
      switch (type) {
      case PIXEL_GRAY:
	glTexImage2D(GL_TEXTURE_2D, 0, 1, txtw, txth, 0, GL_LUMINANCE, GL_UNSIGNED_BYTE, (unsigned char*)data);
	break;
      case PIXEL_RGB:
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, txtw, txth, 0, GL_RGB, GL_UNSIGNED_BYTE, (unsigned char*)data);
	break;
      default: std::cerr << "Error (GUIH::OpenGL2DWindow::showImage): invalid pixel type" << std::endl;  break;
      }
      glDisable(GL_TEXTURE_2D);
      use_texture_for_image = true;
    } else {		// copy the image for later rendering
      imgdata = (unsigned char*) data;
      use_texture_for_image = false;
    }
    this->imgtype = type;
    redraw();
    return true;
  }
  void removeImage(void) { imgdata = NULL; txtid = 0; }

  // -----------------------------------------------------------------
  // callback functions
  // -----------------------------------------------------------------

  bool runOnInit(GtkWidget *widget) { 
    GdkGLContext  *glcontext  = gtk_widget_get_gl_context (widget);
    GdkGLDrawable *gldrawable = gtk_widget_get_gl_drawable(widget);
    if (!gdk_gl_drawable_gl_begin (gldrawable, glcontext)) return FALSE;  // OpenGL begin
    glClearColor (bkgrgb[0], bkgrgb[1], bkgrgb[2], 0); // clearing color
    glColor3ub(100, 100, 100);
    if (cb_init)  cb_init( this );
    gdk_gl_drawable_gl_end (gldrawable);	// OpenGL END
    return true;
  }
  void runOnDestroy(void) { 
    // This function is called automatically by a event, and cannot be cancelled
    if (cb_close) cb_close(this); 
    // disable OpenGL
  }
  bool runOnDraw(GtkWidget *widget, GdkEventExpose *event) {
    GdkGLContext  *glcontext  = gtk_widget_get_gl_context (widget);
    GdkGLDrawable *gldrawable = gtk_widget_get_gl_drawable(widget);
    if (!gdk_gl_drawable_gl_begin (gldrawable, glcontext)) return FALSE;  // OpenGL begin
    Camera2D *cam = &(camera);
    glClear(GL_COLOR_BUFFER_BIT);			// clear
    glMatrixMode (GL_PROJECTION);   glLoadIdentity();
    gluOrtho2D(cam->xmin, cam->xmax, cam->ymin, cam->ymax);
    glMatrixMode (GL_MODELVIEW);
    if (use_texture_for_image && txtid > 0) {	// render textured image
      glEnable(GL_TEXTURE_2D);
      glBindTexture(GL_TEXTURE_2D, txtid);
      glBegin(GL_QUADS);
      glTexCoord2f(0.0, 0.0);  glVertex2f(  0 , txth);
      glTexCoord2f(1.0, 0.0);  glVertex2f(txtw, txth);
      glTexCoord2f(1.0, 1.0);  glVertex2f(txtw,   0 );
      glTexCoord2f(0.0, 1.0);  glVertex2f(  0 ,   0 );
      glEnd();
      glDisable(GL_TEXTURE_2D);
    } else if (!use_texture_for_image && imgdata) {
      int x, y, y2;  unsigned char *src = imgdata;
      glBegin(GL_QUADS);
      switch (imgtype) {
      case PIXEL_RGB:  
	for (y = 0; y < txth; y++) {
	  y2 = txth - 1 - y;
	  for (x = 0; x < txtw; x++, src+=3) {
	    glColor3ub(src[0], src[1], src[2]);
	    glVertex2i(x  , y2+1);	glVertex2i(x+1, y2+1);
	    glVertex2i(x+1, y2);	glVertex2i(x  , y2);
	  }
	}
	break;
      case PIXEL_GRAY: 
	for (y = 0; y < txth; y++) {
	  y2 = txth - 1 - y;
	  for (x = 0; x < txtw; x++, src+=1) {
	    glColor3ub(src[0], src[0], src[0]);
	    glVertex2i(x  , y2+1);	glVertex2i(x+1, y2+1);
	    glVertex2i(x+1, y2);	glVertex2i(x  , y2);
	  }
	}
	break;
      default: break;
      }
      glEnd();
    }
    if (show_axes)  cam->renderXYAxes();		// render axes
    glColor3ub(100, 100, 100);
    if (cb_draw) cb_draw( this );			// user callback
    if (gdk_gl_drawable_is_double_buffered (gldrawable))
      gdk_gl_drawable_swap_buffers (gldrawable);	// flush
    else glFlush ();
    gdk_gl_drawable_gl_end (gldrawable);		// OpenGL end
    return true;
  }
  
  void runOnKeyboard(int key);
  void runOnMouse(int button, int state, int xy[]);
  
  void getPositionFromMouseClick(int mxy[2], float wxy[2]) {
    if (txtid > 0) {
      wxy[0] = camera.xmin + (camera.xmax - camera.xmin) * mxy[0] / area->allocation.width;
      wxy[1] = (camera.ymaxb - camera.ymax) + (camera.ymax - camera.ymin) * mxy[1] / area->allocation.height;
    } else {
      wxy[0] = camera.xmin + (camera.xmax - camera.xmin) * mxy[0] / area->allocation.width;
      wxy[1] = camera.ymax - (camera.ymax - camera.ymin) * mxy[1] / area->allocation.height;
    }
  }
  
  void drawLine(float x0, float y0, float x1, float y1);
  void drawCircle(double x, double y, double radius, int slices=0, bool filled=false);
  void drawRectangle(double x, double y, double width, double height);
};

}	// end of GUIH namespace

// ===================================================================
#endif	// end of linux version
// ===================================================================


namespace GUIH {
  
  void OpenGL2DWindow::runOnKeyboard(int key) 
  {
    bool skey = (key==ESCAPE_KEY || key>200);
    if (skey && !skey_by_user) {  // Special Keys with default behaviors
      switch (key) {
      case ARROW_LEFT:  camera.shift( -0.1, 0 );  break;
      case ARROW_RIGHT: camera.shift( +0.1, 0 );  break;
      case ARROW_UP:    if (ctrl) camera.zoomIn();  else camera.shift(0, +0.1);  break;
      case ARROW_DOWN:  if (ctrl) camera.zoomOut(); else camera.shift(0, -0.1);  break;
      case HOME_KEY: 
	if (ctrl&&shift) { }                            // 
	else if (ctrl)   { camera.printInfo(); }        // print camera info
	else if (shift)  {                              // Show camera view range
	  char buffer[160];
	  sprintf(buffer, "Camera (%.2f %.2f) - (%.2f %.2f)", camera.xmin, camera.ymin, camera.xmax, camera.ymax);
	  setCaption(buffer);
	}
	else { camera.restore(); }                      // initial viewpoint
	break;
      case END_KEY:
	if (ctrl&&shift) {              // change background color
	  float v=(bkgrgb[0]==0 ? 1 : 0);  setBackgroundColor( v, v, v );
	} else if (ctrl) {              // increase the window size
	  int width=0, height=0;  getSize( &width, &height, true );
	  if (width<1600) resize( 2*width, 2*height );
	} else if (shift) {             // decrease the window size
	  int width=0, height=0;  getSize( &width, &height, true );
	  resize( width/2, height/2 );
	} else show_axes = !show_axes;  // turn on/off axes
	break;
      case PAGE_UP:   break;
      case PAGE_DOWN: break;
      default: break;
      }
    }
    if (key == F1 && shift) {
      std::cerr << "---------------------------------------------------------------------" << std::endl;
      std::cerr << "  GUIH::OpenGL2DWindow control keys                                  " << std::endl;
      std::cerr << "    Home       : Go back to the initial camera pose                  " << std::endl;
      std::cerr << "               : ^ PrintInfo $ ShowInfo                              " << std::endl;
      std::cerr << "    End        : Turn on/off axes                                    " << std::endl;
      std::cerr << "               : ^ IncWindow $ DecWindow ^$ ChageBkgColor            " << std::endl;
      std::cerr << "    Ins        : Capture the screen in 'guih_capture.png'            " << std::endl;
      std::cerr << "    Del        : Clear caption                                       " << std::endl;
      std::cerr << "    Arrow keys : Shift the camera (panning)                          " << std::endl;
      std::cerr << "               : ^ Zoom In/Out                                       " << std::endl;
      std::cerr << "---------------------------------------------------------------------" << std::endl;
    }
    // run user key press callback function
    if (cb_keys && (!skey || skey_by_user)) cb_keys( this, key );
  }
  
  void OpenGL2DWindow::runOnMouse(int button, int state, int xy[]) 
  { 
    if (button == MOUSE_LEFT && state == MOUSE_CLICKED) {
      float wxy[2];  char str[80];  // default behavior: show the clicked position
      getPositionFromMouseClick( xy, wxy );
      sprintf( str, "MouseClick: uv(%d, %d) xy(%.2f, %.2f)", xy[0], xy[1], wxy[0], wxy[1]);
      setCaption( str );
    }
    if (cb_mouse) cb_mouse( this, button, state, xy ); 
  }
  
  void OpenGL2DWindow::drawLine(float x0, float y0, float x1, float y1) 
  {
    ////if (linewidth > 0) glSetWidth(line_width);
    glBegin(GL_LINES);
    glVertex2f( x0, y0 );  glVertex2f( x1, y1 );  
    glEnd();
  }
  
  void OpenGL2DWindow::drawCircle(double x, double y, double radius, int slices, bool filled) 
  {
    double tx, ty, theta;
    if (slices < 4) slices = 32;
    glBegin( filled ? GL_POLYGON : GL_LINE_LOOP );
    for (int i = 0; i < slices; i++) {
      theta = i * (2 * M_PI) / slices;
      tx = x + (double)(radius * cos( theta ));
      ty = y + (double)(radius * sin( theta ));
      glVertex2d( tx, ty );
    }
    glEnd();
  }
  
  void OpenGL2DWindow::drawRectangle(double x, double y, double width, double height) {
    double hw = width/2, hh = height/2;
    glBegin(GL_LINE_STRIP);
    glVertex2d( x - hw, y - hh );
    glVertex2d( x + hw, y - hh );
    glVertex2d( x + hw, y + hh );
    glVertex2d( x - hw, y + hh );
    glVertex2d( x - hw, y - hh );
    glEnd();
  }
  
}	// end of GUIH namespace


#endif	// GUIH_OPENGL2D_HPP



// ===================================================================
#if 0	// Example code starts
// ===================================================================
// For common functions, refer to GUIH::Window in "guih_common.hpp".
// compile on Linux:  g++ -o guih guih_example.cpp `pkg-config --cflags gtkglext-1.0 --libs gtkglext-1.0` -Wl,-rpath,/usr/local/lib
#include <iostream>
#include "guih_opengl2d.hpp"
void draw(GUIH::Window *w) {
  // draw a star using standard OpenGL commands
  glBegin(GL_LINE_LOOP);  
  glVertex2f(+3, +1);  glVertex2f(-3, +1);  glVertex2f(+2, -3);
  glVertex2f( 0, +3.2);  glVertex2f(-2, -3);
  glEnd();
  // draw a rectangle using a preset function of GUIH
  GUIH::OpenGL2DWindow *w2d = (GUIH::OpenGL2DWindow *) w;
  w2d->drawRectangle( 5, 5, 8.0, 3.0 );  // cx, cy, width, height
}
void mouse(GUIH::Window *w, int button, int state, int xy[])
{
  if        (button==GUIH::MOUSE_LEFT && state==GUIH::MOUSE_CLICKED) {
    float wxy[2];  w->getPositionFromMouseClick( xy, wxy );
    printf(" clicked  at (%d,%d),  which is (%d,%d) in world,\n", xy[0], xy[1], wxy[0], wxy[1]);
  } else if (button==GUIH::MOUSE_LEFT && state==GUIH::MOUSE_DRAGGED) {
    int mpos, mcnt;  for(mpos=0; xy[mpos]>=0; mpos+=2); mcnt = mpos/2;
    printf(" dragged from (%d,%d) to (%d,%d)\n", xy[0], xy[1], xy[mcnt*2-2], xy[mcnt*2-1]);
  }
}
int  main(void)
{
  // You can create windows as many as you want.
  GUIH::OpenGL2DWindow w0(300, 300, "GUIH OpenGL2D 0", NULL, draw, NULL, mouse);
  w0.showImage("../../sample.jpg", true);  // camera is set to image width/height
  GUIH::OpenGL2DWindow w1(300, 300, "GUIH OpenGL2D 1", NULL, draw, NULL, mouse);
  w1.initCamera( -10.0, -10.0, +20.0, +20.0 );  // xmin, ymin, xmax, ymax
  w1.show_axes = true;
  std::cout << "Press Shift+F1 to see how to use Arrow keys" << std::endl;
  // start the main loop -- press ESC to quit
  w0.runMainLoop();
  return EXIT_SUCCESS;
}
// ===================================================================
#endif	// Example code ends
// ===================================================================
