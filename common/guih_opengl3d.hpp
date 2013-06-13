
//
// GUIH                 : GUI in Header files
// GUIH::OpenGL3DWindow : window for 3D OpenGL graphics
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
// GUIH::OpenGL3DWindow
//   This class provides very convenient GUI for 2D and 3D OpenGL graphics.
//   All the details of setting up OpenGL are hidden, and 2D/3D cameras with
//   key bindings for the navigation are also provided by default.
//
//   Member functions:
//	void setFloorXY(float size, int ngrid=20, float x=0, float y=0, float z=0);
//	void setFloor  (float size, int ngrid=20, float *center=NULL, int axis=2);
//      void initCamera( fx,fy,fz, ax,ay,az, ux,uy,uz, lx,ly,lz );
//      void initCamera( center[], dx,dy,dz, ux,uy,uz, lx,ly,lz );
//      void drawLine (float p0[3], float p1[3], int line_width=0);
//      void drawBox  (float minx, float miny, float minz, float maxx, float maxy, float maxz);
//      void drawDisk (double inner_radius, double outer_radius);
//      void drawCone (double radius, double height, int slices=20, bool closed=false);
//      void drawCube (double radius);
//      void drawSphere  (double radius, int slices=20, int stacks=20);
//      void drawCylinder(double radius, double height, int slices=20, bool closed=false);
//      void drawTorus   (float o_radius, float i_radius, int o_steps, int i_steps);
//      void drawTetrahedron (float radius);
//      void drawOctahedron  (float radius);
//      void drawDodecahedron(float radius);
//      void drawIcosahedron (float radius);
//      void drawArrow (float from[3], float dir[3], float len=0, float radius=0);
//
// For API example, read the example code at the bottom.
// For default key bindings, press Shift+F1 during execution.
//


#ifndef GUIH_OPENGL3D_HPP
#define GUIH_OPENGL3D_HPP

#include <iostream>
#ifdef WIN32					// Microsoft Windows
#include <windows.h>
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
#include "guih_camera3d.hpp"

  
// ===================================================================
#ifdef  WIN32   // beginning of Window version
// ===================================================================

namespace GUIH {
  
  
class OpenGL3DWindow : public Window
{
public:
  HDC              hDC;
  HGLRC            hRC;
  bool	           show_axes;
  bool             show_focal_point;
  bool             show_floor, navigation_mode;
  int              timg, twidth, theight;	// textured image
  float floor_size, floor_center[3];
  int   floor_grid, floor_axis;
  Camera3D	camera;
  float		fov_near;	// near clipping plane
  float		fov_far;	// far  clipping plane
  float		fov_angle;	// viewing angle
  GLUquadric	*quadric;
  float		bkgrgb[3];
public:
  OpenGL3DWindow () { 
    show_axes = false;  show_focal_point = true;  skey_by_user = false;
    show_floor = false;  navigation_mode = false;  floor_size = 0;
    initCamera( 0,0,20, 0,0,0, 0,1,0, -0.2f,+0.4f,52.0f );
    bkgrgb[0] = bkgrgb[1] = bkgrgb[2] = 1.0;
    // for predefined objects to render
    quadric = gluNewQuadric();
    gluQuadricDrawStyle(quadric, GLU_FILL); // GLU_FILL, GLU_LINE, GLU_SILHOUETTE, GLU_POINT
    gluQuadricNormals(quadric, GLU_FLAT);   // GLU_NONE, GLU_FLAT, GLU_SMOOTH
  }
  OpenGL3DWindow (int w, int h, char *title = NULL, 
		  void (*init) (Window *w) = NULL,
		  void (*draw) (Window *w) = NULL,
		  void (*keys) (Window *w, int k) = NULL, 
		  void (*mouse)(Window *w, int button, int state, int xy[]) = NULL,
		  void (*close)(Window *w) = NULL) {
    show_axes = false;  show_focal_point = true;  skey_by_user = false;
    show_floor = false;  navigation_mode = false;  floor_size = 0;
    initCamera( 0,0,20, 0,0,0, 0,1,0, -0.2f,+0.4f,52.0f );
    bkgrgb[0] = bkgrgb[1] = bkgrgb[2] = 1.0;
    openWindow(w, h, title, init, draw, keys, mouse, close); 
    // for predefined objects to render
    quadric = gluNewQuadric();
    gluQuadricDrawStyle(quadric, GLU_FILL); // GLU_FILL, GLU_LINE, GLU_SILHOUETTE, GLU_POINT
    gluQuadricNormals(quadric, GLU_FLAT);   // GLU_NONE, GLU_FLAT, GLU_SMOOTH
  }
  ~OpenGL3DWindow() { }
  inline bool isReady(void) { return (hDC && hRC); }
  
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
  void initCamera(float fx, float fy, float fz,		// from
		  float ax, float ay, float az,		// at
		  float ux, float uy, float uz,		// up
		  float lx=0, float ly=0, float lz=0) {	// light
    camera.initPose( fx, fy, fz, ax, ay, az, ux, uy, uz, lx, ly, lz );
    navigation_mode = false;
    setPerspective();
    redraw();
  }
  void initCameraOverFloor(float fx, float fy, float fz,  float ax, float ay, float az,
			   float floor_equation[4]=NULL, float lx=0, float ly=0, float lz=0) {	// light
    camera.initPoseOverFloor( fx, fy, fz, ax, ay, az, floor_equation, lx, ly, lz );
    navigation_mode = true;
    setPerspective();
    redraw();
  }
  
  void setPerspective(float nearv=0.01, float farv=0, float angle=60);
  void setFloorXY(float size, int ngrid=20, float x=0, float y=0, float z=0);
  void setFloor(float size, int ngrid=20, float *center=NULL, int axis=2);

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
    wglMakeCurrent( hDC, hRC );
    // setup OpenGL
    glEnable (GL_LIGHTING);			// lighting
    GLfloat  ambientLight[]  = { 0.0f, 0.0f, 0.0f, 0.0f };
    GLfloat  diffuseLight[]  = { 0.6f, 0.6f, 0.6f, 0.0f };
    GLfloat  specularLight[] = { 0.0f, 0.0f, 0.0f, 0.0f};
    glLightfv(GL_LIGHT0, GL_AMBIENT,  ambientLight);
    glLightfv(GL_LIGHT0, GL_DIFFUSE,  diffuseLight);
    glLightfv(GL_LIGHT0, GL_SPECULAR, specularLight);
    glEnable(GL_LIGHT0);			// light 0
    glFrontFace(GL_CCW); 			// other global settings
    //glShadeModel (GL_SMOOTH);
    glColorMaterial (GL_FRONT, GL_AMBIENT_AND_DIFFUSE); 
    glEnable (GL_COLOR_MATERIAL);		// coloring
    glClearColor (0.0, 0.0, 0.0, 0.0);	// clearing color
    glEnable(GL_DEPTH_TEST);   		// depth buffering
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
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);	// clear
    glLoadIdentity();
    gluLookAt( camera.From[0], camera.From[1], camera.From[2],	// set 3D camera
	       camera.At[0],   camera.At[1],   camera.At[2],
	       camera.Up[0],   camera.Up[1],   camera.Up[2]);
    glLightfv( GL_LIGHT0, GL_POSITION, camera.Light );
    if (show_axes)   camera.renderXYZAxes();	// render axes
    if (show_focal_point) camera.renderLookAtPosition();	// render the center (LookAt)
    if (show_floor)  camera.renderFloor( floor_size, floor_grid, floor_center, floor_axis );
    glColor3ub(100, 100, 100);
    if (cb_draw) cb_draw( this );			// user callback
    glFinish();
    SwapBuffers(hDC);
    return true;
  }
  
  void runOnKeyboard(int key);
  void runOnMouse(int button, int state, int xy[]);
  
  void getRayFromMouseClick(int mxy[2], float from[3], float dir[3]) {
    // calculate the ray from the viewpoint 
    //   that goes through the point on the image plane
    RECT rect;  GetClientRect( hwnd, &rect );
    camera.ScreenToWorld( abs(rect.right - rect.left), abs(rect.top - rect.bottom), mxy[0], mxy[1], from, dir );
  }
  
  inline void setQuadricDrawStyle(int style) { gluQuadricDrawStyle(quadric, style); }
  inline void setQuadricNormals(int normals) { gluQuadricNormals(quadric, normals); }
  
  void drawLine(float p0[3], float p1[3], int line_width=0);
  void drawBox(float minx, float miny, float minz, float maxx, float maxy, float maxz);
  void drawDisk(double inner_radius, double outer_radius);
  void drawCone(double radius, double height, int slices=20, bool closed=false);
  void drawSphere(double radius, int slices=20, int stacks=20);
  void drawSphereAt(double x, double y, double z, double radius, int slices=20, int stacks=20);
  void drawCylinder(double radius, double height, int slices=20, bool closed=false);
  void drawCube(double radius);
  void drawCubeAt(double x, double y, double z, double radius);
  void drawTorus(float o_radius, float i_radius, int o_steps, int i_steps);
  void drawTetrahedron(float radius);
  void drawOctahedron(float radius);
  void drawDodecahedron(float radius);
  void drawIcosahedron(float radius);
  void drawArrow(float from[3], float dir[3], float len=0, float radius=0);

};

}	// end of GUIH namespace

// ===================================================================
#elif defined(__APPLE__) || defined(MACOSX)	// beginning of OS X version
// ===================================================================

namespace GUIH {

class OpenGL3DWindow : public Window
{
public:  // for OpenGL window only
  AGLContext  aglContext;
  bool	 show_axes, show_focal_point, show_floor, navigation_mode;
  float  floor_size, floor_center[3];
  int    floor_grid, floor_axis;
  Camera3D	camera;		//
  float		fov_near;	// near clipping plane
  float		fov_far;	// far  clipping plane
  float		fov_angle;	// viewing angle
  GLUquadric	*quadric;	// surface for pre-defined objects
  float		bkgrgb[3];
  
public:
  OpenGL3DWindow () : aglContext(NULL) { 
    show_axes = false;  show_focal_point = true;  skey_by_user = false;
    show_floor = false;  navigation_mode = false;  floor_size = 0;
    initCamera( 0.0,0,20, 0,0,0, 0,1,0, -0.2,+0.4,52 );	// set the camera
    bkgrgb[0] = bkgrgb[1] = bkgrgb[2] = 1.0;
  }
  OpenGL3DWindow (int w, int h, char *title = NULL, 
		  void (*init) (Window *w) = NULL,
		  void (*draw) (Window *w) = NULL,
		  void (*keys) (Window *w, int k) = NULL, 
		  void (*mouse)(Window *w, int button, int state, int xy[]) = NULL,
		  void (*close)(Window *w) = NULL) : aglContext(NULL) {
    show_axes = false;  show_focal_point = true;  skey_by_user = false;
    show_floor = false;  navigation_mode = false;  floor_size = 0;
    initCamera( 0.0,0,20, 0,0,0, 0,1,0, -0.2,+0.4,52 );	// set the camera
    bkgrgb[0] = bkgrgb[1] = bkgrgb[2] = 1.0;
    setCallbacks(draw, keys, mouse, close);
    openWindow(w, h, title, init, draw, keys, mouse, close); 
  }
  ~OpenGL3DWindow() { if (quadric) gluDeleteQuadric(quadric); }
  inline bool isReady(void) { return (aglContext != NULL); }
  
public:  
  void initOpenGL3D(void) {
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
    glEnable (GL_DEPTH_TEST);		// depth buffering
    glEnable (GL_LIGHTING);		// lighting
    GLfloat  ambientLight[]  = { 0.0f, 0.0f, 0.0f, 0.0f };
    GLfloat  diffuseLight[]  = { 0.6f, 0.6f, 0.6f, 0.0f };
    GLfloat  specularLight[] = { 0.0f, 0.0f, 0.0f, 0.0f};
    glLightfv(GL_LIGHT0, GL_AMBIENT,  ambientLight);
    glLightfv(GL_LIGHT0, GL_DIFFUSE,  diffuseLight);
    glLightfv(GL_LIGHT0, GL_SPECULAR, specularLight);
    glEnable(GL_LIGHT0);			// light 0
    glFrontFace(GL_CCW); 			// other global settings
    //glShadeModel (GL_SMOOTH);
    glColorMaterial (GL_FRONT, GL_AMBIENT_AND_DIFFUSE); 
    glEnable (GL_COLOR_MATERIAL);		// coloring
    glClearColor (bkgrgb[0], bkgrgb[1], bkgrgb[2], 0); // clearing color
    glColor3ub(100, 100, 100);
    setPerspective();
    quadric = gluNewQuadric();              // for predefined objects to render
    gluQuadricDrawStyle(quadric, GLU_FILL); // GLU_FILL, GLU_LINE, GLU_SILHOUETTE, GLU_POINT
    gluQuadricNormals(quadric, GLU_FLAT);   // GLU_NONE, GLU_FLAT, GLU_SMOOTH
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
    createBasicWindow( w, h, title );	// create the window first (hidden)
    initOpenGL3D();			// initialize OpenGL 3D
    ShowWindow( window );		// show the window
    SelectWindow( window );
  }
  void setBackgroundColor(float r, float g, float b) {
    bkgrgb[0] = r;  bkgrgb[1] = g;  bkgrgb[2] = b;
    if (!aglContext) return;
    aglSetCurrentContext (aglContext);
    glClearColor (bkgrgb[0], bkgrgb[1], bkgrgb[2], 0); // clearing color
    redraw();
  }
  void initCamera(float fx, float fy, float fz,		// from
		  float ax, float ay, float az,		// at
		  float ux, float uy, float uz,		// up
		  float lx=0, float ly=0, float lz=0) {	// light
    camera.initPose( fx, fy, fz, ax, ay, az, ux, uy, uz, lx, ly, lz );
    navigation_mode = false;
    setPerspective();
    redraw();
  }
  void initCameraOverFloor(float fx, float fy, float fz,  float ax, float ay, float az,
			   float floor_equation[4]=NULL, float lx=0, float ly=0, float lz=0) {	// light
    camera.initPoseOverFloor( fx, fy, fz, ax, ay, az, floor_equation, lx, ly, lz );
    navigation_mode = true;
    setPerspective();
    redraw();
  }
  
  void setPerspective(float nearv=0.01, float farv=0, float angle=60);
  void setFloorXY(float size, int ngrid=20, float x=0, float y=0, float z=0);
  void setFloor(float size, int ngrid=20, float *center=NULL, int axis=2);

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
    glDisable(GL_DEPTH_TEST);  glDisable(GL_LIGHTING);
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
    glEnable(GL_DEPTH_TEST);  glEnable(GL_LIGHTING);
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
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);	// clear
    glLoadIdentity();
    gluLookAt( camera.From[0], camera.From[1], camera.From[2],// set 3D camera
	       camera.At[0],   camera.At[1],   camera.At[2],
	       camera.Up[0],   camera.Up[1],   camera.Up[2]);
    glLightfv( GL_LIGHT0, GL_POSITION, camera.Light );
    if (show_axes)   camera.renderXYZAxes();		// render axes
    if (show_focal_point) camera.renderLookAtPosition();	// render the center (LookAt)
    if (show_floor)  camera.renderFloor( floor_size, floor_grid, floor_center, floor_axis );
    glColor3ub(100, 100, 100);
    if (cb_draw) cb_draw( this );		// user callback
    draw_caption();
    aglSwapBuffers( aglGetCurrentContext() );
    return true;
  }
  
  void runOnKeyboard(int key);
  void runOnMouse(int button, int state, int xy[]);
  
  void getRayFromMouseClick(int mxy[2], float from[3], float dir[3]) {
    // calculate the ray from the viewpoint 
    //   that goes through the point on the image plane
    Rect rect; GetWindowBounds(window, kWindowContentRgn, &rect); 
    camera.ScreenToWorld( (rect.right-rect.left), (rect.bottom-rect.top), mxy[0], mxy[1], from, dir );
  }
  
  inline void setQuadricDrawStyle(int style) { gluQuadricDrawStyle(quadric, style); }
  inline void setQuadricNormals(int normals) { gluQuadricNormals(quadric, normals); }
  
  void drawLine(float p0[3], float p1[3], int line_width=0);
  void drawBox(float minx, float miny, float minz, float maxx, float maxy, float maxz);
  void drawDisk(double inner_radius, double outer_radius);
  void drawCone(double radius, double height, int slices=20, bool closed=false);
  void drawSphere(double radius, int slices=20, int stacks=20);
  void drawSphereAt(double x, double y, double z, double radius, int slices=20, int stacks=20);
  void drawCylinder(double radius, double height, int slices=20, bool closed=false);
  void drawCube(double radius);
  void drawCubeAt(double x, double y, double z, double radius);
  void drawTorus(float o_radius, float i_radius, int o_steps, int i_steps);
  void drawTetrahedron(float radius);
  void drawOctahedron(float radius);
  void drawDodecahedron(float radius);
  void drawIcosahedron(float radius);
  void drawArrow(float from[3], float dir[3], float len=0, float radius=0);

};
  
}	// end of GUIH namespace

// ===================================================================
#else		// beginning of Linux/Unix version
// ===================================================================

namespace GUIH {
  

class OpenGL3DWindow : public Window
{
public:  // for OpenGL window only
  bool	show_axes, show_focal_point, show_floor, navigation_mode;
  float floor_size, floor_center[3];
  int   floor_grid, floor_axis;
  Camera3D	camera;		//
  float		fov_near;	// near clipping plane
  float		fov_far;	// far  clipping plane
  float		fov_angle;	// viewing angle
  GLUquadric	*quadric;	// surface for pre-defined objects
  float		bkgrgb[3];
  
public:
  OpenGL3DWindow () { 
    show_axes = false;  show_focal_point = true;  skey_by_user = false;
    show_floor = false;  navigation_mode = false;  floor_size = 0;
    initCamera( 0,0,20, 0,0,0, 0,1,0, -0.2,+0.4,52 );	// set the camera
    bkgrgb[0] = bkgrgb[1] = bkgrgb[2] = 1.0;
    setPerspective();
    // for predefined objects to render
    quadric = gluNewQuadric();
    gluQuadricDrawStyle(quadric, GLU_FILL); // GLU_FILL, GLU_LINE, GLU_SILHOUETTE, GLU_POINT
    gluQuadricNormals(quadric, GLU_FLAT);   // GLU_NONE, GLU_FLAT, GLU_SMOOTH
  }
  OpenGL3DWindow (int w, int h, char *title = NULL, 
		  void (*init) (Window *w) = NULL,
		  void (*draw) (Window *w) = NULL,
		  void (*keys) (Window *w, int k) = NULL, 
		  void (*mouse)(Window *w, int button, int state, int xy[]) = NULL,
		  void (*close)(Window *w) = NULL) {
    show_axes = false;  show_focal_point = true;  skey_by_user = false;
    show_floor = false;  navigation_mode = false;  floor_size = 0;
    openWindow(w, h, title, init, draw, keys, mouse, close); 
    initCamera( 0,0,20, 0,0,0, 0,1,0, -0.2,+0.4,52 );	// set the camera
    bkgrgb[0] = bkgrgb[1] = bkgrgb[2] = 1.0;
    setPerspective();
    // for predefined objects to render
    quadric = gluNewQuadric();
    gluQuadricDrawStyle(quadric, GLU_FILL); // GLU_FILL, GLU_LINE, GLU_SILHOUETTE, GLU_POINT
    gluQuadricNormals(quadric, GLU_FLAT);   // GLU_NONE, GLU_FLAT, GLU_SMOOTH
  }
  ~OpenGL3DWindow() { if (quadric) gluDeleteQuadric(quadric); }
  inline bool isReady(void) { return true; }
  
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
    gtk_gl_init(0, NULL);
    GdkGLConfigMode mode;
    mode = (GdkGLConfigMode)(GDK_GL_MODE_RGB|GDK_GL_MODE_DEPTH|GDK_GL_MODE_DOUBLE);
    GdkGLConfig *glconfig = gdk_gl_config_new_by_mode (mode);
    gtk_widget_set_gl_capability (area, glconfig, NULL, TRUE, GDK_GL_RGBA_TYPE);
    
    RegisterBasicEvents( init );	// setup default event handlers
    setCallbacks(draw, keys, mouse, close);
    gtk_widget_show( window );		// show the window
  }
  void setBackgroundColor(float r, float g, float b) {
    GdkGLContext  *glcontext  = gtk_widget_get_gl_context (area);
    GdkGLDrawable *gldrawable = gtk_widget_get_gl_drawable(area);
    if (!gdk_gl_drawable_gl_begin (gldrawable, glcontext)) return;  // OpenGL begin
    bkgrgb[0] = r;  bkgrgb[1] = g;  bkgrgb[2] = b;
    glClearColor (bkgrgb[0], bkgrgb[1], bkgrgb[2], 0); // clearing color
    gdk_gl_drawable_gl_end (gldrawable);	// OpenGL END
  }
  void initCamera(float fx, float fy, float fz,		// from
		  float ax, float ay, float az,		// at
		  float ux, float uy, float uz,		// up
		  float lx=0, float ly=0, float lz=0) {	// light
    camera.initPose( fx, fy, fz, ax, ay, az, ux, uy, uz, lx, ly, lz );
    navigation_mode = false;
    setPerspective();
    redraw();
  }
  void initCameraOverFloor(float fx, float fy, float fz,  float ax, float ay, float az,
			   float floor_equation[4]=NULL, float lx=0, float ly=0, float lz=0) {	// light
    camera.initPoseOverFloor( fx, fy, fz, ax, ay, az, floor_equation, lx, ly, lz );
    navigation_mode = true;
    setPerspective();
    redraw();
  }
  
  void setPerspective(float near=0.01, float far=0, float angle=60);
  void setFloorXY(float size, int ngrid=20, float x=0, float y=0, float z=0);
  void setFloor(float size, int ngrid=20, float *center=NULL, int axis=2);

  // -----------------------------------------------------------------
  // callback functions
  // -----------------------------------------------------------------

  bool runOnInit(GtkWidget *widget) { 
    GdkGLContext  *glcontext  = gtk_widget_get_gl_context (widget);
    GdkGLDrawable *gldrawable = gtk_widget_get_gl_drawable(widget);
    if (!gdk_gl_drawable_gl_begin (gldrawable, glcontext)) return FALSE;  // OpenGL begin
    glEnable (GL_LIGHTING);			// lighting
    GLfloat  ambientLight[]  = { 0.0f, 0.0f, 0.0f, 0.0f };
    GLfloat  diffuseLight[]  = { 0.6f, 0.6f, 0.6f, 0.0f };
    GLfloat  specularLight[] = { 0.0f, 0.0f, 0.0f, 0.0f};
    glLightfv(GL_LIGHT0, GL_AMBIENT,  ambientLight);
    glLightfv(GL_LIGHT0, GL_DIFFUSE,  diffuseLight);
    glLightfv(GL_LIGHT0, GL_SPECULAR, specularLight);
    glEnable(GL_LIGHT0);			// light 0
    glFrontFace(GL_CCW); 			// other global settings
    //glShadeModel (GL_SMOOTH);
    glColorMaterial (GL_FRONT, GL_AMBIENT_AND_DIFFUSE); 
    glEnable (GL_COLOR_MATERIAL);		// coloring
    glClearColor (bkgrgb[0], bkgrgb[1], bkgrgb[2], 0); // clearing color
    glEnable(GL_DEPTH_TEST);   		// depth buffering
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
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);	// clear
    glLoadIdentity();
    gluLookAt( camera.From[0], camera.From[1], camera.From[2],// set 3D camera
	       camera.At[0],   camera.At[1],   camera.At[2],
	       camera.Up[0],   camera.Up[1],   camera.Up[2]);
    glLightfv( GL_LIGHT0, GL_POSITION, camera.Light );
    if (show_axes)   camera.renderXYZAxes();		// render axes
    if (show_focal_point) camera.renderLookAtPosition();	// render the center (LookAt)
    if (show_floor)  camera.renderFloor( floor_size, floor_grid, floor_center, floor_axis );
    glColor3ub(100, 100, 100);
    if (cb_draw) cb_draw( this );		// user callback
    if (gdk_gl_drawable_is_double_buffered (gldrawable))
      gdk_gl_drawable_swap_buffers (gldrawable);	// flush
    else glFlush ();
    gdk_gl_drawable_gl_end (gldrawable);		// OpenGL end
    return true;
  }
  
  void runOnKeyboard(int key);
  void runOnMouse(int button, int state, int xy[]);
  
  void getRayFromMouseClick(int mxy[2], float from[3], float dir[3]) {
    // calculate the ray from the viewpoint 
    //   that goes through the point on the image plane
    camera.ScreenToWorld( area->allocation.width, area->allocation.height, mxy[0], mxy[1], from, dir );
  }
  
  inline void setQuadricDrawStyle(int style) { gluQuadricDrawStyle(quadric, style); }
  inline void setQuadricNormals(int normals) { gluQuadricNormals(quadric, normals); }
  
  void drawLine(float p0[3], float p1[3], int line_width=0);
  void drawBox(float minx, float miny, float minz, float maxx, float maxy, float maxz);
  void drawDisk(double inner_radius, double outer_radius);
  void drawCone(double radius, double height, int slices=20, bool closed=false);
  void drawSphere(double radius, int slices=20, int stacks=20);
  void drawSphereAt(double x, double y, double z, double radius, int slices=20, int stacks=20);
  void drawCylinder(double radius, double height, int slices=20, bool closed=false);
  void drawCube(double radius);
  void drawCubeAt(double x, double y, double z, double radius);
  void drawTorus(float o_radius, float i_radius, int o_steps, int i_steps);
  void drawTetrahedron(float radius);
  void drawOctahedron(float radius);
  void drawDodecahedron(float radius);
  void drawIcosahedron(float radius);
  void drawArrow(float from[3], float dir[3], float len=0, float radius=0);

};

}	// end of GUIH namespace


// ===================================================================
#endif	// end of linux version
// ===================================================================


namespace GUIH {
  
  void OpenGL3DWindow::setPerspective(float nearv, float farv, float angle) 
  {
    // This function sets near/far clipping plane and viewing angle.
    // Call this function (possibly with default arguments) whenever the depth 
    // of camera (distance from camera to target) is changed significantly.
    // Note that you should not set 'near' to 0. make it as big as possible.
    if (!isReady()) return;
    int   w=0, h=0;  getSize( &w, &h, true );
    if (nearv >= farv) {
      float dist = camera.getDistance();
      nearv = dist * 0.01f;   farv = dist * 50.0f;
    }
    glMatrixMode (GL_PROJECTION);
    glLoadIdentity ();
    gluPerspective(angle, (float)w/(float)h, nearv, farv );
    //glFrustum (-1.0, 1.0, -h/w, h/w, near, far);
    glMatrixMode (GL_MODELVIEW);
    fov_near  = nearv;
    fov_far   = farv;
    fov_angle = angle;
    redraw();
  }
  void OpenGL3DWindow::setFloorXY(float size, int ngrid, float x, float y, float z) 
  {
    float floor_center[3] = { x, y, z };
    setFloor( size, ngrid, floor_center, 2 );
  }
  void OpenGL3DWindow::setFloor(float size, int ngrid, float *center, int axis) 
  {
    // set 'size' to zero, to remove the floor
    if (center) memcpy(floor_center, center, 3 * sizeof(float));
    else        memset(floor_center, 0, 3 * sizeof(float));
    floor_size = size;
    floor_grid = ngrid;
    floor_axis = axis;
    show_floor = true;
  }
  
  void OpenGL3DWindow::runOnKeyboard(int key) 
  {
    bool skey = (key==ESCAPE_KEY || key>200);
    float move = camera.getDistance() / 10;
    if (skey && !skey_by_user) {  // Special Keys with default behaviors
      switch (key) {
      case ARROW_LEFT:      
	if (alt && !shift) camera.turn(0, 0, +5);		// roll to the left
	else if  (!ctrl) {
	  if (shift) camera.move(0, +move, 0);			// move left
	  else       camera.circle( +5, 0 );			// circling left
	} else if (ctrl) {
	  if (shift) camera.moveOverFloor( 0, +move, 0 );	// move left w.r.t. floor
	  else       camera.turnOverFloor( +5, 0 );		// turn left w.r.t. floor
	}
	break;
      case ARROW_RIGHT:     
	if (alt && !shift) camera.turn(0, 0, -5);		// roll to the right
	else if  (!ctrl) {
	  if (shift) camera.move(0, -move, 0);			// move right
	  else       camera.circle( -5, 0 );			// circling right
	} else if (ctrl) {
	  if (shift) camera.moveOverFloor( 0, -move, 0 );	// move right w.r.t. floor
	  else       camera.turnOverFloor( -5, 0 );		// turn right w.r.t. floor
	}
	break;
      case ARROW_UP:        
	if (alt && !shift)  camera.zoomIn();			// zoom in
	else if  (!ctrl) {	// inspection
	  if (alt&&shift) camera.move( +move, 0, 0 );		// move forward
	  else if (shift) camera.move( 0, 0, +move );		// move up
	  else            camera.circle( 0, +5 );		// circling up
	} else if (ctrl) {	// navigation
	  if (alt&&shift) camera.circleOverFloor( 0, +5 );	// circling up w.r.t. floor
	  else if (shift) camera.moveOverFloor( 0, 0, +move );	// move up w.r.t. floor
	  else            camera.moveOverFloor( +move, 0, 0 );	// move forward w.r.t. floor
	}
	break;
      case ARROW_DOWN:      
	if (alt && !shift)  camera.zoomOut();			// zoom out
	else if  (!ctrl) {	// inspection
	  if (alt&&shift) camera.move( -move, 0, 0 );		// move backward
	  else if (shift) camera.move( 0, 0, -move );		// move down
	  else            camera.circle( 0, -5 );		// circling up
	} else if (ctrl) {	// navigation
	  if (alt&&shift) camera.circleOverFloor( 0, -5 );	// circling up w.r.t. floor
	  else if (shift) camera.moveOverFloor( 0, 0, -move );	// move down w.r.t. floor
	  else            camera.moveOverFloor( -move, 0, 0 );	// move backward w.r.t. floor
	}
	break;
      case HOME_KEY:
	if (ctrl&&shift) { }
	else if (ctrl)   { camera.printInfo(); }        // print camera info
	else if (shift)  { show_floor = !show_floor; }  // show floor
	else camera.restore();                          // initial viewpoint
	break;
      case END_KEY:
	if (ctrl&&shift) { 	// change background color
	  float v=(bkgrgb[0]==0 ? 1 : 0);  setBackgroundColor( v, v, v );
	} else if (ctrl) { 	// increase the window size
	  int width=0, height=0;  getSize( &width, &height, true );
	  if (width<1600) resize( 2*width, 2*height );
	} else if (shift) {	// decrease the window size
	  int width=0, height=0;  getSize( &width, &height, true );
	  resize( width/2, height/2 );
	} else {		// turn on/off axes and target point
	  show_axes = !show_axes;  show_focal_point = !show_focal_point; 
	}
	break;
      case PAGE_UP: break;
      case PAGE_DOWN: break;
      }
    }
    if (key == F1 && shift) {
      std::cerr << "---------------------------------------------------------------------" << std::endl;
      std::cerr << "  GUIH::OpenGL3DWindow control keys                                  " << std::endl;
      std::cerr << "    Home (% -) : Go back to the initial viewpoint                    " << std::endl;
      std::cerr << "               : ^ PrintInfo $ ShowFloor                             " << std::endl;
      std::cerr << "    End  (% =) : Turn on/off axes and target point                   " << std::endl;
      std::cerr << "               : ^ IncWindow $ DecWindow ^$ ChageBkgColor            " << std::endl;
      std::cerr << "    Ins  (% i) : Capture the screen in 'guih_capture.png'            " << std::endl;
      std::cerr << "    Del        : Clear caption                                       " << std::endl;
      std::cerr << "    Arrow keys : Move or turn camera, zoom in/out                    " << std::endl;
      std::cerr << "      -Ctrl    : [ Inspection mode ]                                 " << std::endl;
      std::cerr << "        -Shift : move camrea (circling) around the target            " << std::endl;
      std::cerr << "        +Shift : shift left/right, up/down, forward/backward(+alt)   " << std::endl;
      std::cerr << "      +Ctrl    : [ navigation mode ]                                 " << std::endl;
      std::cerr << "        -Shift : turn left/right or go forward/backward              " << std::endl;
      std::cerr << "        +Shift : shift left/right, up/down, rotate up/down(+alt)     " << std::endl;
      std::cerr << "      +Alt     : zoom in/out or roll left/right                      " << std::endl;
      std::cerr << "---------------------------------------------------------------------" << std::endl;
    }
    // run user key press callback function
    if (key > 0 && cb_keys && (!skey || skey_by_user)) cb_keys( this, key );
  }
  
  void OpenGL3DWindow::runOnMouse(int button, int state, int xy[]) 
  { 
    if (button == MOUSE_LEFT && state == MOUSE_CLICKED) {
      // float from[3], dir[3];  char str[80];  // default behavior: show the clicked position
      // getRayFromMouseClick( xy, from, dir );
      // sprintf( str, "MouseClick: uv(%d %d)  ray:(%.2f %.2f %.2f)->(%.2f %.2f %.2f)", 
      // 	       xy[0], xy[1], from[0], from[1], from[2], dir[0], dir[1], dir[2] );
      // setCaption( str );
    }
    if (cb_mouse) cb_mouse( this, button, state, xy ); 
  }
  void OpenGL3DWindow::drawLine(float p0[3], float p1[3], int line_width) 
  {
    ////if (linewidth > 0) glSetWidth(line_width);
    glDisable(GL_LIGHTING);
    glBegin(GL_LINES);
    glVertex3fv( p0 );  glVertex3fv( p1 );  
    glEnd();
    glEnable(GL_LIGHTING);
  }
  
  void OpenGL3DWindow::drawBox(float minx, float miny, float minz, float maxx, float maxy, float maxz) 
  {
    float v[24] = {
      minx, miny, minz,  maxx, miny, minz,  maxx, maxy, minz,  minx, maxy, minz,
      minx, miny, maxz,  maxx, miny, maxz,  maxx, maxy, maxz,  minx, maxy, maxz };
    float *n, n0[18] = {	// normal vectors of the faces
      0, -1, 0,   1, 0, 0,   		// lower/right side
      0, 1, 0,    -1, 0, 0,		// upper/left side
      0, 0, -1,   0, 0, +1 };		// bottom/top
    int *f, f0[24] = {		// vertex indices of the faces
      0, 3, 15, 12,    3, 6, 18, 15,	// lower/right side
      6, 9, 21, 18,    9, 0, 12, 21,	// upper/left side
      9, 6, 3, 0,     12, 15, 18, 21 };	// bottom/top
    n = n0;  f = f0;
    glBegin(GL_QUADS);
    for (int i = 0; i < 6; i++, n+=3, f+=4) {
      glNormal3fv( n );
      glVertex3fv( v + f[0] );   glVertex3fv( v + f[1] );
      glVertex3fv( v + f[2] );   glVertex3fv( v + f[3] );
    }
    glEnd();
  }
  
  void OpenGL3DWindow::drawDisk(double inner_radius, double outer_radius) 
  {
    gluDisk(quadric, inner_radius, outer_radius, 20, 2);
  }
  
  void OpenGL3DWindow::drawCone(double radius, double height, int slices, bool closed) 
  {
    gluCylinder( quadric, radius, 0, height, slices, 2 );
    if (closed) gluDisk(quadric, 0, radius, 20, 2);
  }
  
  void OpenGL3DWindow::drawSphere(double radius, int slices, int stacks) 
  {
    gluSphere( quadric, radius, slices, stacks );
  }
  
  void OpenGL3DWindow::drawSphereAt(double x, double y, double z, double radius, int slices, int stacks) 
  {
    glPushMatrix();
    glTranslated( x, y, z );
    gluSphere( quadric, radius, slices, stacks );
    glPopMatrix();
  }
  
  void OpenGL3DWindow::drawCylinder(double radius, double height, int slices, bool closed) 
  {
    gluCylinder( quadric, radius, radius, height, slices, 2 );
    if (closed) {
      glPushMatrix();
      gluDisk(quadric, 0, radius, 20, 2);
      glTranslated(0,0,height);
      glRotatef(180,1,0,0); // So that normals point right way!
      gluDisk(quadric, 0, radius, 20, 2);
      glPopMatrix();
    }
  }
  
  void OpenGL3DWindow::drawCube(double radius) 
  {
    static GLuint gllist = 0;
    static float last_radius = 0;
    if (radius != last_radius) {
      static float v[24], v0[24] = {
	-1.0, -1.0, -1.0,  +1.0, -1.0, -1.0,  +1.0, +1.0, -1.0,  -1.0, +1.0, -1.0,
	-1.0, -1.0, +1.0,  +1.0, -1.0, +1.0,  +1.0, +1.0, +1.0,  -1.0, +1.0, +1.0 };
      static float *n, n0[18] = {	// normal vectors of the faces
	0, -1, 0,   1, 0, 0,   		// lower/right side
	0, 1, 0,    -1, 0, 0,		// upper/left side
	0, 0, -1,   0, 0, +1 };		// bottom/top
      static int *f, f0[24] = {	// vertex indices of the faces
	0, 3, 15, 12,  3, 6, 18, 15,	// lower/right side
	6, 9, 21, 18,  9, 0, 12, 21,	// upper/left side
	9, 6, 3, 0,   12, 15, 18, 21 };	// bottom/top
      if (radius != last_radius) {
	for (int i = 0; i < 24; i++) v[i] = v0[i] * (float)radius;
	last_radius = (float)radius;
      }
      if (gllist == 0) gllist = glGenLists(1);
      glNewList( gllist, GL_COMPILE_AND_EXECUTE );
      n = n0;  f = f0;
      glBegin(GL_QUADS);
      for (int i = 0; i < 6; i++, n+=3, f+=4) {
	glNormal3fv( n );
	glVertex3fv( v + f[0] );   glVertex3fv( v + f[1] );
	glVertex3fv( v + f[2] );   glVertex3fv( v + f[3] );
      }
      glEnd();
      glEndList();
    } else glCallList( gllist );
  }
  
  void OpenGL3DWindow::drawCubeAt(double x, double y, double z, double radius)
  {
    glPushMatrix();
    glTranslated( x, y, z );
    drawCube( radius );
    glPopMatrix();
  }
  
  void OpenGL3DWindow::drawTorus(float o_radius, float i_radius, int o_steps, int i_steps) 
  {
    static GLuint gllist = 0;
    static float last_o_radius = 0, last_i_radius = 0;
    static int   last_o_steps  = 0, last_i_steps  = 0;
    if (o_radius != last_o_radius || o_steps != last_o_steps || 
	i_radius != last_i_radius || i_steps != last_i_steps ) {
      if (gllist == 0) gllist = glGenLists(1);
      last_o_radius = o_radius;  last_o_steps = o_steps;
      last_i_radius = i_radius;  last_i_steps = i_steps;
      glNewList( gllist, GL_COMPILE_AND_EXECUTE );
      glBegin(GL_QUADS);
      float a, b, n[3], v[3];
      float astep = 3.141592f*2 / o_steps,  aoff = astep * 0.5f;
      float bstep = 3.141592f*2 / i_steps,  boff = bstep * 0.5f;
      for (int i = 0; i < o_steps; i++) {
	a = i * astep;
	for (int j = 0; j < i_steps; j++) {
	  b = j * bstep;
	  n[0] = cos(a) * cos(b);	// nx
	  n[1] = sin(a) * cos(b);	// ny
	  n[2] = sin(b);		// nz
	  glNormal3fv( n );
	  float a_nb[4] = { a-aoff, a+aoff, a+aoff, a-aoff };
	  float b_nb[4] = { b-boff, b-boff, b+boff, b+boff };
	  for (int k = 0; k < 4; k++) {
	    v[0] = cos(a_nb[k]) * o_radius + cos(a_nb[k]) * cos(b_nb[k]) * i_radius; // x
	    v[1] = sin(a_nb[k]) * o_radius + sin(a_nb[k]) * cos(b_nb[k]) * i_radius; // y
	    v[2] = sin(b_nb[k]) * i_radius;  // z
	    glVertex3fv( v );
	  }
	}
      }
      glEnd();
      glEndList();
    } else glCallList( gllist );
  }

  void OpenGL3DWindow::drawTetrahedron(float radius) 
  {
    static GLuint gllist = 0;
    static float last_radius = 0;
    if (radius != last_radius) {
      static float v[12], v0[12] = {
	-1.0000f, -1.0000f, -1.0000f,     +1.0000f, +1.0000f, -1.0000f, 
	+1.0000f, -1.0000f, +1.0000f,     -1.0000f, +1.0000f, +1.0000f
      };
      static float n[12] = {
	0.58f, 0.58f, 0.58f,     0.58f, -0.58f, -0.58f, 
	-0.58f, -0.58f, 0.58f,     -0.58f, 0.58f, -0.58f
      };
      static int f[12] = { 3, 6, 9,  3, 0, 6,  9, 6, 0,  0, 3, 9 };
      if (radius != last_radius) {
	for (int i = 0; i < 12; i++) v[i] = v0[i] * radius;
	last_radius = radius;
      }
      if (gllist == 0) gllist = glGenLists(1);
      glNewList( gllist, GL_COMPILE_AND_EXECUTE );
      glBegin(GL_TRIANGLES);
      for (int i = 0; i < 12; i+=3) {
	glNormal3fv(n+i);
	glVertex3fv(v+f[i+0]);  glVertex3fv(v+f[i+1]);  glVertex3fv(v+f[i+2]);
      }
      glEnd();
      glEndList();
    } else glCallList( gllist );
  }

  void OpenGL3DWindow::drawOctahedron(float radius) 
  {
    static GLuint gllist = 0;
    static float last_radius = 0;
    if (radius != last_radius) {
      static float v[18], v0[18] = {
	+1.0000f, +0.0000f, +0.0000f,     +0.0000f, -1.0000f, +0.0000f, 
	-1.0000f, +0.0000f, +0.0000f,     +0.0000f, +1.0000f, +0.0000f, 
	+0.0000f, +0.0000f, +1.0000f,     +0.0000f, +0.0000f, -1.0000f
      };
      static float n[24] = {
	0.58f, -0.58f, 0.58f,     -0.58f, -0.58f, 0.58f, 
	-0.58f, 0.58f, 0.58f,     0.58f, 0.58f, 0.58f, 
	0.58f, -0.58f, -0.58f,     -0.58f, -0.58f, -0.58f, 
	-0.58f, 0.58f, -0.58f,     0.58f, 0.58f, -0.58f
      };
      static int f[24] = { 
	12, 0, 3,     12, 3, 6,     12, 6, 9,     12, 9, 0, 
	15, 3, 0,     15, 6, 3,     15, 9, 6,     15, 0, 9
      };
      if (radius != last_radius) {
	for (int i = 0; i < 18; i++) v[i] = v0[i] * radius;
	last_radius = radius;
      }
      if (gllist == 0) gllist = glGenLists(1);
      glNewList( gllist, GL_COMPILE_AND_EXECUTE );
      glBegin(GL_TRIANGLES);
      for (int i = 0; i < 24; i+=3) {
	glNormal3fv(n+i);
	glVertex3fv(v+f[i+0]);  glVertex3fv(v+f[i+1]);  glVertex3fv(v+f[i+2]);
      }
      glEnd();
      glEndList();
    } else glCallList( gllist );
  }

  void OpenGL3DWindow::drawDodecahedron(float radius) 
  {
    static GLuint gllist = 0;
    static float last_radius = 0;
    if (radius != last_radius) {
      static float v[60], v0[60] = {
	-0.5774f, -0.5774f, +0.5774f,     +0.9342f, +0.3568f, +0.0000f, 
	+0.9342f, -0.3568f, +0.0000f,     -0.9342f, +0.3568f, +0.0000f, 
	-0.9342f, -0.3568f, +0.0000f,     +0.0000f, +0.9342f, +0.3568f, 
	+0.0000f, +0.9342f, -0.3568f,     +0.3568f, +0.0000f, -0.9342f, 
	-0.3568f, +0.0000f, -0.9342f,     +0.0000f, -0.9342f, -0.3568f, 
	+0.0000f, -0.9342f, +0.3568f,     +0.3568f, +0.0000f, +0.9342f, 
	-0.3568f, +0.0000f, +0.9342f,     +0.5774f, +0.5774f, -0.5774f, 
	+0.5774f, +0.5774f, +0.5774f,     -0.5774f, +0.5774f, -0.5774f, 
	-0.5774f, +0.5774f, +0.5774f,     +0.5774f, -0.5774f, -0.5774f, 
	+0.5774f, -0.5774f, +0.5774f,     -0.5774f, -0.5774f, -0.5774f
      };
      static float *n, n0[36] = {
	0.85f, -0.00f, 0.53f,     0.85f, -0.00f, -0.53f, 
	-0.85f, 0.00f, -0.53f,     -0.85f, -0.00f, 0.53f, 
	-0.53f, 0.85f, 0.00f,     0.53f, 0.85f, -0.00f, 
	0.53f, -0.85f, 0.00f,     -0.53f, -0.85f, -0.00f, 
	0.00f, -0.53f, -0.85f,     0.00f, 0.53f, -0.85f, 
	-0.00f, 0.53f, 0.85f,     0.00f, -0.53f, 0.85f
      };
      static int *f, f0[60] = { 
	3, 6, 54, 33, 42,     3, 39, 21, 51, 6,     9, 12, 57, 24, 45, 
	9, 48, 36, 0, 12,     9, 45, 18, 15, 48,     3, 42, 15, 18, 39, 
	6, 51, 27, 30, 54,     12, 0, 30, 27, 57,     21, 24, 57, 27, 51, 
	18, 45, 24, 21, 39,     15, 42, 33, 36, 48,     30, 0, 36, 33, 54
      };
      if (radius != last_radius) {
	for (int i = 0; i < 60; i++) v[i] = v0[i] * radius;
	last_radius = radius;
      }
      if (gllist == 0) gllist = glGenLists(1);
      glNewList( gllist, GL_COMPILE_AND_EXECUTE );
      n = n0;  f = f0;
      for (int i = 0; i < 12; i++, n+=3, f+=5) {
	glNormal3fv( n );
	glBegin(GL_POLYGON);
	glVertex3fv(v+f[0]);  glVertex3fv(v+f[1]);  glVertex3fv(v+f[2]);  
	glVertex3fv(v+f[3]);  glVertex3fv(v+f[4]);  
	glEnd();
      }
      glEndList();
    } else glCallList( gllist );
  }

  void OpenGL3DWindow::drawIcosahedron(float radius) 
  {
    static GLuint gllist = 0;
    static float last_radius = 0;
    if (radius != last_radius) {
      static float v[36], v0[36] = {	// 12 vertices
	+0.0000f, -0.5257f, +0.8507f,     +0.8507f, +0.0000f, +0.5257f, 
	+0.8507f, +0.0000f, -0.5257f,     -0.8507f, +0.0000f, -0.5257f, 
	-0.8507f, +0.0000f, +0.5257f,     -0.5257f, +0.8507f, +0.0000f, 
	+0.5257f, +0.8507f, +0.0000f,     +0.5257f, -0.8507f, +0.0000f, 
	-0.5257f, -0.8507f, +0.0000f,     +0.0000f, -0.5257f, -0.8507f, 
	+0.0000f, +0.5257f, -0.8507f,     +0.0000f, +0.5257f, +0.8507f
      };
      static float n[60] = {		// 20 faces
	+0.36f, +0.00f, +0.93f,    -0.36f, +0.00f, +0.93f,    -0.58f, -0.58f, +0.58f,
	+0.00f, -0.93f, +0.36f,    +0.58f, -0.58f, +0.58f,    +0.58f, +0.58f, +0.58f,
	+0.00f, +0.93f, +0.36f,    -0.58f, +0.58f, +0.58f,    -0.93f, +0.36f, +0.00f,
	-0.93f, -0.36f, +0.00f,    -0.58f, -0.58f, -0.58f,    +0.00f, -0.93f, -0.36f,
	+0.58f, -0.58f, -0.58f,    +0.93f, -0.36f, +0.00f,    +0.93f, +0.36f, +0.00f,
	+0.00f, +0.93f, -0.36f,    +0.58f, +0.58f, -0.58f,    +0.36f, +0.00f, -0.93f,
	-0.36f, +0.00f, -0.93f,    -0.58f, +0.58f, -0.58f
      };
      static int f[60] = {		// 20 faces
	0, 3, 33,  0, 33, 12,  0, 12, 24,  0, 24, 21,  0, 21, 3,
	3,18,33, 33,18,15, 33,15,12,  12,15,9,  12,9,24,  24,9,27,  24,27,21,  21,27,6,  21,6,3,  3,6,18,
	30,15,18,  30,18,6,  30,6,27,  30,27,9,  30,9,15
      };
      if (radius != last_radius) {
	for (int i = 0; i < 36; i++) v[i] = v0[i] * radius;
	last_radius = radius;
      }
      if (gllist == 0) gllist = glGenLists(1);
      glNewList( gllist, GL_COMPILE_AND_EXECUTE );
      glBegin(GL_TRIANGLES);
      for (int i = 0; i < 60; i+=3) {
	glNormal3fv(n+i);
	glVertex3fv(v+f[i+0]);  glVertex3fv(v+f[i+1]);  glVertex3fv(v+f[i+2]);
      }
      glEnd();
      glEndList();
    } else glCallList( gllist );
  }
  
  void OpenGL3DWindow::drawArrow(float from[3], float dir[3], float len, float radius) 
  {
    if (len <= 0)    len    = GUIH_V3D_LENGTH( dir );
    if (len == 0) return;
    if (radius <= 0) radius = len * 0.05f;
    float bodylen = len * 0.6f;
    float headlen = len * 0.4f;
    float radius2 = radius * 1.5f;
    int   slices = 20;
    
    // calculate local coordinate system where dir[3] becomes z-axis
    float lz[3], lx[3], ly[3], tmp[16];
    GUIH_V3D_COPY( lz, dir );
    GUIH_V3D_NORMALIZE( lz );
    if (lz[0] == 1.0) GUIH_V3D_SET( tmp, 0, 1, 0 );
    else              GUIH_V3D_SET( tmp, 1, 0, 0 );
    GUIH_V3D_CROSS( lx, lz, tmp );  GUIH_V3D_NORMALIZE( lx );
    GUIH_V3D_CROSS( ly, lz, lx  );  GUIH_V3D_NORMALIZE( ly );
    GUIH_M4D_SET( tmp,
		  lx[0], lx[1], lx[2], 0.0, 
		  ly[0], ly[1], ly[2], 0.0, 
		  lz[0], lz[1], lz[2], 0.0, 
		  from[0], from[1], from[2], 1.0 );
    // draw the arrow pointing in the z-axis direction
    glPushMatrix();
    glMultMatrixf( tmp );
    gluCylinder( quadric, radius, radius, bodylen, slices, 2 );	// body
    gluDisk( quadric, 0, radius, 20, 2 );				// bottom
    glPushMatrix();
    glTranslated( 0, 0, bodylen );
    gluCylinder( quadric, radius2, 0, headlen, slices, 2 );	// head
    gluDisk( quadric, 0, radius2, 20, 2 );			// neck
    glPopMatrix();
    glPopMatrix();
  }


}	// end of GUIH namespace


#endif	// GUIH_OPENGL3D_HPP



// ===================================================================
#if 0	// Example code starts
// ===================================================================
// For common functions, refer to GUIH::Window in "guih_common.hpp".
// compile on Linux:  g++ -o guih guih_example.cpp `pkg-config --cflags gtkglext-1.0 --libs gtkglext-1.0` -Wl,-rpath,/usr/local/lib
#include <iostream>
#include "guih_opengl3d.hpp"
void draw(GUIH::Window *w) {
  // draw a star using standard OpenGL commands
  glDisable(GL_LIGHTING);
  glBegin(GL_LINE_LOOP);  
  glVertex2f(+3, +1);  glVertex2f(-3, +1);  glVertex2f(+2, -3);
  glVertex2f( 0, +3.2);  glVertex2f(-2, -3);
  glEnd();
  glEnable(GL_LIGHTING);
  // draw a torus using a preset function of GUIH
  GUIH::OpenGL3DWindow *w3d = (GUIH::OpenGL3DWindow *) w;
  w3d->drawTorus( 5, 2, 20, 12 );  // out_radius, in_radius, slices, steps
}
void mouse(GUIH::Window *w, int button, int state, int xy[])
{
  if        (button==GUIH::MOUSE_LEFT && state==GUIH::MOUSE_CLICKED) {
    float from[3], dir[3];  w->getRayFromMouseClick( xy, from, dir );
    printf(" clicked  at (%d,%d)  ray: (%.1f %.1f %.1f)->(%.2f %.2f %.2f) in world,\n", 
	   xy[0], xy[1], from[0], from[1], from[2], dir[0], dir[1], dir[2] );
  } else if (button==GUIH::MOUSE_LEFT && state==GUIH::MOUSE_DRAGGED) {
    int mpos, mcnt;  for(mpos=0; xy[mpos]>=0; mpos+=2); mcnt = mpos/2;
    printf(" dragged from (%d,%d) to (%d,%d)\n", xy[0], xy[1], xy[mcnt*2-2], xy[mcnt*2-1]);
  }
}
int  main(void)
{
  // You can create windows as many as you want.
  GUIH::OpenGL3DWindow w(300, 300, "GUIH OpenGL3D", NULL, draw, NULL, mouse);
  w.initCamera( 0, 0, 20,  0,0,0, 0,1,0 );  // from[3], at[3], up[3]
  w.show_axes = true;
  std::cout << "Press Shift+F1 to see how to use Arrow keys" << std::endl;
  // start the main loop -- press ESC to quit
  w.runMainLoop();
  return EXIT_SUCCESS;
}
// ===================================================================
#endif	// Example code ends
// ===================================================================
