
//
// GUIH           : GUI in Header files
// GUIH::Camera2D : 2D camera for GUIH::OpenGL2DWindow
//
// (c) 2006  Jaeil Choi
// last modified in May, 2007
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
// GUIH::Camera2D  Features:
//   User-space translation
//   Zoom in/out
// GUIH::Camera2D  Usage:
//   GUIH::Camera3D cam;
//   cam.initPose(0,0,10,  0,0,0,  0,1,0,  50,50,150);
//   //            from     at      up      light
//   cam.zoomIn();		// zoom in
//   cam.shift(0,1,0);		// translation of camera
//   cam.restore();		// return to initialization values
//   cam.printInfo();		// show the current view information
//   //
//   glLoadIdentity();
//   gluLookAt( cam.From[0], cam.From[1], cam.From[2],
//	        cam.At[0],   cam.At[1],   cam.At[2],
//	        cam.Up[0],   cam.Up[1],   cam.Up[2]);
//
// If you use GUIH::OpenGL2DWindow, the camera is already set.
// You only need to specify the initial pose of the camera at the beginning.
//

#ifndef GUIH_CAMERA2D_HPP
#define GUIH_CAMERA2D_HPP


#ifdef WIN32					// Microsoft Windows
#include <windows.h>
#include <GL/gl.h>
#include <GL/glu.h>
#elif defined(__APPLE__) || defined(MACOSX)	// Apple OS X
#include <AGL/agl.h>
#include <OpenGL/gl.h>
#include <OpenGL/glu.h>
#else						// Linux
#include <GL/gl.h>
#include <GL/glu.h>
#endif
#include <iostream>

namespace GUIH {
  
class Camera2D {

public:
  float xmin, xmax, ymin, ymax;	// current viewing area
  float xminb, xmaxb, yminb, ymaxb;  // backup of the initial viewing area

public:
  Camera2D() {
    xmin = xminb = ymin = yminb = 0.0;
    xmax = xmaxb = ymax = ymaxb = 0.0;
  }
  ~Camera2D() {}

  void initPose(float xmin, float ymin, float xmax, float ymax) {
    this->xmin = xmin;  this->ymin = ymin;
    this->xmax = xmax;  this->ymax = ymax;
    backupMinMax();
  }
  void initPoseForBB(float xmin, float ymin, float xmax, float ymax, float margin=0.0f) {
    // Obsolete: Note that it only works only for windows with same width and height.
    float sx, sy, size;
    sx = xmax - xmin;  sy = ymax - ymin;  
    size = (sx > sy ? sx : sy) * (1.0f + margin) / 2.0f;
    sx = (xmax + xmin)/2;  sy = (ymax + ymin)/2;  
    this->xmin = sx - size;  this->ymin = sy - size;
    this->xmax = sx + size;  this->ymax = sy + size;
    backupMinMax();
  }
  void initPoseForImage(int width, int height) {
    xmin = 0;  xmax = width-0.5f;   ymin = 0;  ymax = height-0.5f;
    backupMinMax();
  }
  void backupMinMax(void) { xminb = xmin; xmaxb = xmax; yminb = ymin; ymaxb = ymax; }
  void restore(void)      { xmin = xminb; xmax = xmaxb; ymin = yminb; ymax = ymaxb; }
  
  void zoomIn(void)  { zoom(0.80f); }
  void zoomOut(void) { zoom(1.25f); }
  void zoom(float rate) {
    float xdiff = xmax - xmin;
    float ydiff = ymax - ymin;
    xmin += xdiff * (1.0f-rate)/2.0f;    xmax -= xdiff * (1.0f-rate)/2.0f;
    ymin += ydiff * (1.0f-rate)/2.0f;    ymax -= ydiff * (1.0f-rate)/2.0f;
  }
  
  void shift(float x_ratio, float y_ratio) {
    float xdiff = xmax - xmin;
    float ydiff = ymax - ymin;
    xmin += xdiff * x_ratio;  xmax += xdiff * x_ratio;
    ymin += ydiff * y_ratio;  ymax += ydiff * y_ratio;
  }

  void printInfo(void) {
    printf( "GUIH::Camera2D Information \n");
    printf( "  x : %7.2f - %7.2f     (Initial: %7.2f %7.2f)\n", xmin, xmax, xminb, xmaxb);
    printf( "  y : %7.2f - %7.2f     (Initial: %7.2f %7.2f)\n", ymin, ymax, yminb, ymaxb);
  }
  
  void ScreenToWorld( int w, int h, int sx, int sy, float xy[] ) {
    // in screen space upper-left corner is (0,0).
    xy[0] = xmin + (float)sx / (float)w * (xmax - xmin);
    xy[1] = ymax - (float)sy / (float)h * (ymax - ymin);
  }
  
  // ===================================================================
  // Axes
  // ===================================================================

  void renderXYAxes(bool bUseColors = true) {
    int i;
    float step, pos;
    
    glLineWidth(1.0);
    glColor3f(0.5, 0.5, 0.5);
    glBegin(GL_LINE_STRIP);		// X axis
    if (bUseColors) glColor3f(1, 0.2f, 0.2f);
    for (i = -10; i <= 10; i++) glVertex2f((float)(i*i*i), 0.0f);
    glEnd();
    glBegin(GL_LINE_STRIP);		// Y axis
    if (bUseColors) glColor3f(0.2f, 1, 0.2f);
    for (i = -10; i <= 10; i++) glVertex2f(0.0f, (float)(i*i*i));
    glEnd();
    
    glBegin(GL_LINES);			// grids on the axes
    step = 0.1f;
    while (step <= 10) {
      if      (step < 1) glColor3f(0.2f, 0.2f, 0.2f);
      else if (step > 1) glColor3f(0.2f, 0.2f, 0.5f);
      else               glColor3f(0.2f, 0.4f, 0.4f);
      for (i = -10; i <= +10; i++) {
	pos = i * step;
	glVertex2f( pos, +step/20 );  glVertex2f( pos, -step/20 ); // X
	glVertex2f( +step/20, pos );  glVertex2f( -step/20, pos ); // Y
      }
      step *= 10;
    }
    glEnd();
    glColor3f(0.5, 0.5, 0.5);
  }

};
  
  
}	// end of GUIH namespace


#endif	// GUIH_CAMERA2D_HPP 
  
