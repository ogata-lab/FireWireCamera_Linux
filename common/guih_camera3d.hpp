
//
// GUIH           : GUI in Header files
// GUIH::Camera3D : 3D camera for GUIH::OpenGL2DWindow
//
// (c) 2006  Jaeil Choi
// last modified in Aug, 2006
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
// GUIH::Camera3D  Features:
//   User-space rotation
//   User-space translation
//   Zoom in/out
//   view-dependent light rotation
// GUIH::Camera3D  Usage:
//   GUIH::Camera3D cam;
//   cam.initPose(0,0,10,  0,0,0,  0,1,0,  50,50,150);   // from(3), at(3), up(3), light(3)
//   cam.initPoseOverFloor(0,0,10,  0,0,0,  floor_eq,  50,50,150);  // from(3), at(3), floor_eq(*), light(3)
//   // controlling viewpoint
//   cam.zoom();			// zoom in/out
//   cam.move(forwad,left,up);		// translation of camera
//   cam.turn(yaw,pitch,roll);		// turning camera direction
//   cam.circle(left,up);		// rotation of camera around target point
//   cam.moveOverFloor(forwad,left,up);	// translation of camera w.r.t. floor
//   cam.turnOverFloor(yaw,pitch,roll);	// turning camera direction w.r.t. floor
//   cam.circleOverFloor(left,up);	// rotation of camera around target point w.r.t. floor
//   cam.restore();		// return to the initial pose
//   cam.printInfo();		// show the current view information
//   // setting OpenGL using the current viewpont
//   glLoadIdentity();
//   gluLookAt( cam.From[0], cam.From[1], cam.From[2],
//	        cam.At[0],   cam.At[1],   cam.At[2],
//	        cam.Up[0],   cam.Up[1],   cam.Up[2]);
//   glLightfv( GL_LIGHT0, GL_POSITION, cam.Light );
// 
// If you use GUIH::OpenGL3DWindow, the camera is already set.
// You only need to specify the initial pose of the camera at the beginning.
//


#ifndef GUIH_CAMERA3D_HPP
#define GUIH_CAMERA3D_HPP


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
#include <cstdio>
#include <cmath>

namespace GUIH {
  
#define GUIH_V3D_SET(r,x,y,z)        do { (r)[0] = (x);  (r)[1] = (y);  (r)[2] = (z); } while(0)
#define GUIH_V3D_COPY(d,s)           do { (d)[0] = (s)[0];  (d)[1] = (s)[1];  (d)[2] = (s)[2]; } while(0)
#define GUIH_V3D_ADD(r,a,b)          do { (r)[0] = (a)[0] + (b)[0];  (r)[1] = (a)[1] + (b)[1];  (r)[2] = (a)[2] + (b)[2]; } while(0)
#define GUIH_V3D_SUB(r,a,b)          do { (r)[0] = (a)[0] - (b)[0];  (r)[1] = (a)[1] - (b)[1];  (r)[2] = (a)[2] - (b)[2]; } while(0)
#define GUIH_V3D_CROSS(r,a,b)        do { (r)[0] = (a)[1]*(b)[2]-(a)[2]*(b)[1]; (r)[1] = (a)[2]*(b)[0]-(a)[0]*(b)[2]; (r)[2] = (a)[0]*(b)[1]-(a)[1]*(b)[0]; } while(0)
#define GUIH_V3D_LENGTH(a)           (sqrt((a)[0]*(a)[0] + (a)[1]*(a)[1] + (a)[2]*(a)[2]))
#define GUIH_V3D_NORMALIZE(a)        do { float len = GUIH_V3D_LENGTH(a);  if (len > 0) GUIH_V3D_DIV_VALUE( a, len ); } while(0)
#define GUIH_V3D_DIV_VALUE(a,v)      do { (a)[0] /= (v);  (a)[1] /= (v);  (a)[2] /= (v); } while(0)
#define GUIH_M3D_SET(m, a00, a01, a02,  a10, a11, a12,  a20, a21, a22) \
                       do { (m)[0] = a00;  (m)[1] = a01;  (m)[2] = a02; \
                            (m)[3] = a10;  (m)[4] = a11;  (m)[5] = a12; \
                            (m)[6] = a20;  (m)[7] = a21;  (m)[8] = a22; } while(0)
#define GUIH_M3D_NORM(M)   (sqrt((M)[0]*(M)[0] + (M)[1]*(M)[1] + (M)[2]*(M)[2] + (M)[3]*(M)[3] + (M)[4]*(M)[4] + (M)[5]*(M)[5] + (M)[6]*(M)[6] + (M)[7]*(M)[7] + (M)[8]*(M)[8]))
#define GUIH_M3D_TRANS(R, A) \
                       do { (R)[0] = (A)[0];  (R)[1] = (A)[3];  (R)[2] = (A)[6]; \
                            (R)[3] = (A)[1];  (R)[4] = (A)[4];  (R)[5] = (A)[7]; \
                            (R)[6] = (A)[2];  (R)[7] = (A)[5];  (R)[8] = (A)[8]; } while(0)
#define GUIH_M3D_MUL_MV(r, M, a) \
                       do { (r)[0] = (M)[0] * (a)[0] + (M)[1] * (a)[1] + (M)[2] * (a)[2]; \
                            (r)[1] = (M)[3] * (a)[0] + (M)[4] * (a)[1] + (M)[5] * (a)[2]; \
                            (r)[2] = (M)[6] * (a)[0] + (M)[7] * (a)[1] + (M)[8] * (a)[2]; } while(0)
#define GUIH_M3D_MUL_MM(R, A, B) \
                       do { (R)[0] = (A)[0] * (B)[0] + (A)[1] * (B)[3] + (A)[2] * (B)[6]; \
                            (R)[1] = (A)[0] * (B)[1] + (A)[1] * (B)[4] + (A)[2] * (B)[7]; \
                            (R)[2] = (A)[0] * (B)[2] + (A)[1] * (B)[5] + (A)[2] * (B)[8]; \
                            (R)[3] = (A)[3] * (B)[0] + (A)[4] * (B)[3] + (A)[5] * (B)[6]; \
                            (R)[4] = (A)[3] * (B)[1] + (A)[4] * (B)[4] + (A)[5] * (B)[7]; \
                            (R)[5] = (A)[3] * (B)[2] + (A)[4] * (B)[5] + (A)[5] * (B)[8]; \
                            (R)[6] = (A)[6] * (B)[0] + (A)[7] * (B)[3] + (A)[8] * (B)[6]; \
                            (R)[7] = (A)[6] * (B)[1] + (A)[7] * (B)[4] + (A)[8] * (B)[7]; \
                            (R)[8] = (A)[6] * (B)[2] + (A)[7] * (B)[5] + (A)[8] * (B)[8]; } while(0)
#define GUIH_V4D_SET(r,x,y,z,w)      do { (r)[0] = (x);  (r)[1] = (y);  (r)[2] = (z);  (r)[3] = (w); } while(0)
#define GUIH_V4D_COPY(d,s)           do { (d)[0] = (s)[0];  (d)[1] = (s)[1];  (d)[2] = (s)[2];  (d)[3] = (s)[3]; } while(0)
#define GUIH_M4D_SET(m, a00, a01, a02, a03,  a10, a11, a12, a13,  a20, a21, a22, a23,  a30, a31, a32, a33) \
                       do { (m)[0] = a00;  (m)[1] = a01;  (m)[2] = a02;  (m)[3] = a03;  \
                            (m)[4] = a10;  (m)[5] = a11;  (m)[6] = a12;  (m)[7] = a13;  \
                            (m)[8] = a20;  (m)[9] = a21;  (m)[10] = a22; (m)[11] = a23; \
                            (m)[12] = a30; (m)[13] = a31; (m)[14] = a32; (m)[15] = a33; } while(0)

class CVec4D {		// Vector in 4D
 public:
  float *A;
 private:
  float A2[4];
 public:
  CVec4D (void) { A = A2; memset(A, 0, 4*sizeof(float)); }
  CVec4D (float x, float y, float z, float w) {
    A = A2;
    A[0] = x;  A[1] = y;  A[2] = z;  A[3] = w;
  }
  inline void assign( float *array ) { A = array; }
  inline void assign( float x, float y, float z, float w) { A[0] = x; A[1] = y; A[2] = z; A[3] = w; }
  inline float getDistance (const CVec4D &v) { 
    float dx = A[0]-v.A[0], dy = A[1]-v.A[1], dz = A[2]-v.A[2], dw = A[3]-v.A[3];
    return sqrt( dx * dx + dy * dy + dz * dz + dw * dw); 
  }
  // +, -, *, /
  void add( CVec4D a, CVec4D b ) { for (int i=0; i<4; i++) A[i] = a.A[i] + b.A[i]; }
  void sub( CVec4D a, CVec4D b ) { for (int i=0; i<4; i++) A[i] = a.A[i] - b.A[i]; }
  void mult( float t ) { for (int i=0; i<4; i++) A[i] *= t; }
  // +=, -=, *=, /=, []
  inline void operator+= (const CVec4D &v) { for (int i=0; i<4; i++) A[i] += v.A[i]; }
  inline float &operator[] (const int i) { return A[i]; }
  inline void operator= (const CVec4D &v)  { memcpy (A, v.A, 4*sizeof(float)); }
  inline void operator= (const float i)    { A[0] = i; A[1] = i; A[2] = i; A[3] = i; }
  void printInfo(char *cmmt=NULL) { printf("%s(%.4f %.4f %.4f  %.4f)\n", (cmmt ? cmmt : ""), A[0], A[1], A[2], A[3]); }
};

/* ostream &operator<< (ostream &os, const CVec4D &v) { */
/*   cout << "( " << setw(10) << v.A[0] << " " << setw(10) << v.A[1] << " " << setw(10) << v.A[2]<< " " << setw(10) << v.A[3] << " )"; */
/*   return os; */
/* } */
  
  
class CMat4D {		// Matrix in 4D
 public:
  float *A;
 private:
  float A2[16];
 public:
  CMat4D (void) { A = A2; memset(A, 0, 16*sizeof(float)); }
  CMat4D (float xx, float xy, float xz, float xw,  float yx, float yy, float yz, float yw,
	  float zx, float zy, float zz, float zw,  float wx, float wy, float wz, float ww) {
    A = A2;
    A[0] = xx; A[1] = xy; A[2] = xz; A[3] = xw;
    A[4] = yx; A[5] = yy; A[6] = yz; A[7] = yw;
    A[8] = zx; A[9] = zy; A[10] = zz; A[11] = zw;
    A[12] = wx; A[13] = wy; A[14] = wz; A[15] = ww;
  }
  
  inline void assign( float xx, float xy, float xz, float xw,
		      float yx, float yy, float yz, float yw,
		      float zx, float zy, float zz, float zw,
		      float wx, float wy, float wz, float ww ) { 
    A[0] = xx; A[1] = xy; A[2] = xz; A[3] = xw;
    A[4] = yx; A[5] = yy; A[6] = yz; A[7] = yw;
    A[8] = zx; A[9] = zy; A[10] = zz; A[11] = zw;
    A[12] = wx; A[13] = wy; A[14] = wz; A[15] = ww; 
  }
  inline void assignLookAt (float fx, float fy, float fz, float ax, float ay, float az, float ux, float uy, float uz) {
    float cx[3], cy[3], cz[3], u[3];  // camera frame axes in world coordinate system
    // Note that camera frame is defined as :  cx:right, cy:up, cz:backward
    // The transformation is for : world -> camera
    GUIH_V3D_SET(u,  ux, uy, uz);
    GUIH_V3D_SET(cz,  fx-ax, fy-ay, fz-az);	GUIH_V3D_NORMALIZE( cz );
    GUIH_V3D_CROSS( cx,  u, cz );		GUIH_V3D_NORMALIZE( cx );
    GUIH_V3D_CROSS( cy,  cz, cx );		GUIH_V3D_NORMALIZE( cy );
    assign(cx[0], cx[1], cx[2], - fx*cx[0] - fy*cx[1] - fz*cx[2],
	   cy[0], cy[1], cy[2], - fx*cy[0] - fy*cy[1] - fz*cy[2],
	   cz[0], cz[1], cz[2], - fx*cz[0] - fy*cz[1] - fz*cz[2],
	   0,    0,    0,    1);
  }
  inline float* operator[] (const int i) { return &A[i*4]; }
  inline void operator= (const CMat4D &m) { memcpy (A, m.A, 16*sizeof(float)); }
  // Matrix = Matrix * Matrix
  void mult( CMat4D &a, CMat4D &b ) {
    float t[16];  memcpy( t, a.A, 16*sizeof(float) );
    A[0] = t[0] * b.A[0] + t[1] * b.A[4] + t[2] * b.A[8] + t[3] * b.A[12];
    A[1] = t[0] * b.A[1] + t[1] * b.A[5] + t[2] * b.A[9] + t[3] * b.A[13];
    A[2] = t[0] * b.A[2] + t[1] * b.A[6] + t[2] * b.A[10] + t[3] * b.A[14];
    A[3] = t[0] * b.A[3] + t[1] * b.A[7] + t[2] * b.A[11] + t[3] * b.A[15];
    A[4] = t[4] * b.A[0] + t[5] * b.A[4] + t[6] * b.A[8] + t[7] * b.A[12];
    A[5] = t[4] * b.A[1] + t[5] * b.A[5] + t[6] * b.A[9] + t[7] * b.A[13];
    A[6] = t[4] * b.A[2] + t[5] * b.A[6] + t[6] * b.A[10] + t[7] * b.A[14];
    A[7] = t[4] * b.A[3] + t[5] * b.A[7] + t[6] * b.A[11] + t[7] * b.A[15];
    A[8] = t[8] * b.A[0] + t[9] * b.A[4] + t[10] * b.A[8] + t[11] * b.A[12];
    A[9] = t[8] * b.A[1] + t[9] * b.A[5] + t[10] * b.A[9] + t[11] * b.A[13];
    A[10] = t[8] * b.A[2] + t[9] * b.A[6] + t[10] * b.A[10] + t[11] * b.A[14];
    A[11] = t[8] * b.A[3] + t[9] * b.A[7] + t[10] * b.A[11] + t[11] * b.A[15];
    A[12] = t[12] * b.A[0] + t[13] * b.A[4] + t[14] * b.A[8] + t[15] * b.A[12];
    A[13] = t[12] * b.A[1] + t[13] * b.A[5] + t[14] * b.A[9] + t[15] * b.A[13];
    A[14] = t[12] * b.A[2] + t[13] * b.A[6] + t[14] * b.A[10] + t[15] * b.A[14];
    A[15] = t[12] * b.A[3] + t[13] * b.A[7] + t[14] * b.A[11] + t[15] * b.A[15];
  }
  void transpose(CMat4D &m) {
    A[0] = m.A[0];  A[1] = m.A[4];  A[2] = m.A[8];  A[3] = m.A[12];  
    A[4] = m.A[1];  A[5] = m.A[5];  A[6] = m.A[9];  A[7] = m.A[13];  
    A[8] = m.A[2];  A[9] = m.A[6];  A[10] = m.A[10];  A[11] = m.A[14];  
    A[12] = m.A[3];  A[13] = m.A[7];  A[14] = m.A[11];  A[15] = m.A[15];  
  }
  void transpose(void) {
    float temp;
    temp = A[1];  A[1] = A[4];  A[4] = temp;
    temp = A[2];  A[2] = A[8];  A[8] = temp;
    temp = A[3];  A[3] = A[12];  A[12] = temp;
    temp = A[6];  A[6] = A[9];  A[9] = temp;
    temp = A[7];  A[7] = A[13];  A[13] = temp;
    temp = A[11];  A[11] = A[14];  A[14] = temp;
  }
  void setRotation(char axis, float degree) {
    float rad = degree * 3.14159265358979323846f / 180;
    switch (axis) {
    case 'X': assign( 1,         0,         0,         0,
		      0,         cos(rad),  -sin(rad), 0,
		      0,         sin(rad),  cos(rad),  0,
		      0,         0,         0,         1);  break;
    case 'Y': assign( cos(rad),  0,         sin(rad),  0,
		      0,         1,         0,         0,
		      -sin(rad), 0,         cos(rad),  0,
		      0,         0,         0,         1);  break;
    case 'Z': assign( cos(rad),  -sin(rad), 0,         0,
		      sin(rad),  cos(rad),  0,         0,
		      0,         0,         1,         0,
		      0,         0,         0,         1);  break;
    }
  }
  void printInfo(char *cmmt=NULL) {
    printf("[ %7.2f %7.2f %7.2f %7.2f ]  %s\n", A[0], A[1], A[2], A[3], (cmmt ? cmmt : "") );
    printf("[ %7.2f %7.2f %7.2f %7.2f ]\n", A[4], A[5], A[6], A[7] );
    printf("[ %7.2f %7.2f %7.2f %7.2f ]\n", A[8], A[9], A[10], A[11] );
    printf("[ %7.2f %7.2f %7.2f %7.2f ]\n", A[12], A[13], A[14], A[15] );
  }
			   
};


}	// end of GUIH namespace

/* ================================================================ */
/* Camera in 3D                                                     */
/* ================================================================ */

namespace GUIH {
  
class Camera3D {

public:
  float From[4], At[4], Up[4];	// current view vectors   (read only)
  float Light[4];		// rotated light position (read only)
  float floor[4];		// plane equation of the floor 
  
private:
  CVec4D vcFrom, vcAt, vcUp;	// current view vectors
  CVec4D viFrom, viAt, viUp;	// initial view vectors
  CMat4D mRotat, mRotatB;	// rotation & its backup    (camera frame -> world frame)
  CVec4D vTPos,  vTPosB;	// target position (focal point) & its backup (in world frame)
  float  fScale;		// scaling
  CVec4D vcLight, viLight;	// light position
  
public:
  Camera3D() {
    vcFrom.assign(From);  vcAt.assign(At);  
    vcUp.assign(Up);      vcLight.assign(Light);
    GUIH_V4D_SET( floor, 0, 0, 1, 0 );
  }
  Camera3D(float fx, float fy, float fz,		// from
	   float ax, float ay, float az,		// at
	   float ux, float uy, float uz,		// up
	   float lx = 0, float ly = 0, float lz = 0) {	// light
    vcFrom.assign(From);  vcAt.assign(At);  
    vcUp.assign(Up);      vcLight.assign(Light);
    initPose( fx, fy, fz,  ax, ay, az,  ux, uy, uz,  lx, ly, lz );
    GUIH_V4D_SET( floor, 0, 0, 1, 0 );
  }
  ~Camera3D() { }
  
public:
  // get the distance from the camera to the target position (vTPos)
  float getDistance(void) { return vcFrom.getDistance(vcAt); }
  // restore initial viewpoints
  void  restore(void) { mRotat = mRotatB;   vTPos = vTPosB;  fScale = 1.0f;  updateViewPoints(); }
  
  // -----------------------------------------------------------------
  // control functions (basic)
  // -----------------------------------------------------------------
public:
  void initPose(float fx, float fy, float fz,			// from
		float ax, float ay, float az,			// at
		float ux, float uy, float uz,			// up
		float lx = 0, float ly = 0, float lz = 0) {	// light
    // set current view vectors
    if (fx == ax && fy == ay && fz == az) fz += 1.0;
    vcFrom.assign(fx, fy, fz, 1);
    vcAt.assign  (ax, ay, az, 1);
    vcUp.assign  (ux, uy, uz, 1);
    // calculate initial trasformations & view vectors
    viFrom.sub( vcFrom, vcAt );
    CMat4D mTmp;  // world -> camera transformaion
    mTmp.assignLookAt(viFrom[0], viFrom[1], viFrom[2],
		      0.0f,      0.0f,      0.0f,
		      vcUp[0],   vcUp[1],   vcUp[2]);
    // camera -> world rotation (camera frame axes in 3 column vectors)
    //   This is the rotation matrix all the additon rotation will be accumulated.
    mRotat.assign(mTmp[0][0], mTmp[1][0], mTmp[2][0], 0.0f,
		  mTmp[0][1], mTmp[1][1], mTmp[2][1], 0.0f,
		  mTmp[0][2], mTmp[1][2], mTmp[2][2], 0.0f,
		  0.0f,       0.0f,       0.0f,       1.0f);
    vTPos = vcAt;
    multVM( viFrom, viFrom, mRotat );	// set back to camera space initial position
    multVM( viUp, vcUp, mRotat );	// set back to camera space initial direction
    viAt   = 0;
    fScale = 1.0f;
    // set light position
    if (lx == 0 && ly == 0 && lz == 0)  vcLight.assign( fx*0.995f, fy*1.01f, fz, 1 );
    else                                vcLight.assign( lx, ly, lz, 1 );
    multVM( viLight, vcLight, mRotat );
    mRotatB = mRotat;  vTPosB = vTPos;  // make backup copy
  }
  
  // Change camera position in forward/backward direction,
  //   without changing the target point and camera pose.
  inline void zoomIn(void)   { zoom(0.80f);  updateViewPoints(); }
  inline void zoomOut(void)  { zoom(1.25f);  updateViewPoints(); }
  inline void zoom(float s)  { fScale *= s;  updateViewPoints(); }
  
  void move(float f, float l, float u) {  // forward, left, up 
    // Change the camera position w.r.t. the current camera pose.
    //   Note that this function only need to move the target point.
    CVec4D vec;
    // Note the 'mRotat' is rotation from camera frame to world frame, 
    //   i.e., its column vectors are camera axes in world frame.
    vec.assign( -f*mRotat[0][2], -f*mRotat[1][2], -f*mRotat[2][2], 1.0 );  vTPos += vec;
    vec.assign( -l*mRotat[0][0], -l*mRotat[1][0], -l*mRotat[2][0], 1.0 );  vTPos += vec;
    vec.assign(  u*mRotat[0][1],  u*mRotat[1][1],  u*mRotat[2][1], 1.0 );  vTPos += vec;
    updateViewPoints();
  }
  
  void turn(float yaw, float pitch, float roll) {
    // Change the camera orientation w.r.t. the current camera pose.
    //   Note that this function moves the target point, too.
    CMat4D Rot, Ry, Rp, Rr;
    Ry.setRotation( 'Y', yaw );
    Rp.setRotation( 'X', pitch );  Rot.mult( Ry, Rp );
    Rr.setRotation( 'Z', roll );   Rot.mult( Rot, Rr );
    mRotat.mult( mRotat, Rot );		// update the camera pose
    makeRotatRotation();
    float dist = getDistance();
    CVec4D vec( 0.0, 0.0, -dist, 1.0 );	// target direction in new camera frame
    multMV( vec, mRotat, vec );		// target direction in world frame
    vTPos.add( vcFrom, vec );		// update the target position
    updateViewPoints();
  }
  
  void circle(float left, float up) {
    // Change camera position and orientation by rotating the camera 
    //   around the target point 'vTPos'.
    // Note that 'mRotat' is the rotation from camera frame to world frame,
    //   and all the new rotations should be added from the right hand side of it.
    CMat4D tmp;
    tmp.setRotation( 'Y', -left );  mRotat.mult( mRotat, tmp );
    tmp.setRotation( 'X', -up   );  mRotat.mult( mRotat, tmp );
    makeRotatRotation();
    updateViewPoints();
  }
  
  // -----------------------------------------------------------------
  // control functions for navigation over a plane (floor)
  // -----------------------------------------------------------------
public:
  void initPoseOverFloor(float fx, float fy, float fz,			// from
			 float ax, float ay, float az,			// at
			 float floor_equation[4]=NULL,			// 
			 float lx = 0, float ly = 0, float lz = 0) {	// light
    if (floor_equation) GUIH_V4D_COPY(floor, floor_equation);	// save floor equation
    else                GUIH_V4D_SET (floor, 0, 0, 1, 0);
    GUIH_V3D_NORMALIZE( floor );
    initPose( fx, fy, fz, ax, ay, az, floor[0], floor[1], floor[2], lx, ly, lz );
  }
  
  void moveOverFloor(float f, float l, float u) {
    // Change the camera position    with respect to the plane equation 'floor[4]'.
    //   Note that this function only need to move the target point.
    CMat4D F;   getTemporaryFrameOnFloor( F );
    CVec4D vec;
    vec.assign( f*F[0][0], f*F[1][0], f*F[2][0], 1.0 );  vTPos += vec;  //vec.printInfo("fv");
    vec.assign( l*F[0][1], l*F[1][1], l*F[2][1], 1.0 );  vTPos += vec;  //vec.printInfo("lv");
    vec.assign( u*F[0][2], u*F[1][2], u*F[2][2], 1.0 );  vTPos += vec;  //vec.printInfo("uv");
    updateViewPoints();
  }
  
  void turnOverFloor(float yaw, float pitch) {
    // Change the camera orientation w.r.t. the plane equation 'floor[4]'.
    //   Note that this function moves the target point, too.
    CMat4D Rwf, Rfw, Rcw, Rwc, Rot, Ry, Rp, tmp;
    // a frame on the floor (based on camera orientation)
    getTemporaryFrameOnFloor( Rwf );  Rfw.transpose( Rwf );
    // mRotat' = mRotat * Rcf * Rot * Rfc
    //         = mRotat * (Rcw*Rwf) * Rot * (Rfw*Rwc)     Note 'mRotat == Rwc'
    Ry.setRotation( 'Z', yaw );
    Rp.setRotation( 'Y', pitch );
    Rot.mult( Ry, Rp );			// rotation in the temporary frame on the floor
    Rwc = mRotat;  Rcw.transpose( Rwc );
    Ry.mult( Rcw, Rwf );  Rp.mult( Rfw, Rwc );
    tmp.mult( Ry, Rot );  tmp.mult( tmp, Rp );
    mRotat.mult( mRotat, tmp );    	// update the camera pose
    makeRotatRotation();
    float dist = getDistance();
    CVec4D vec( 0.0, 0.0, -dist, 1.0 );	// target direction in new camera frame
    multMV( vec, mRotat, vec );		// target direction in world frame
    vTPos.add( vcFrom, vec );		// update the target position
    updateViewPoints();
  }
  
  void circleOverFloor (float left, float up) {
    // Change camera position and orientation by rotating the camera 
    //   around the target point, w.r.t. the plane equation 'floor[4]'.
    CMat4D Rwf, Rfw, Rcw, Rwc, Rot, Ry, Rp, tmp;
    // a frame on the floor (based on camera orientation)
    getTemporaryFrameOnFloor( Rwf );  Rfw.transpose( Rwf );
    // mRotat' = mRotat * Rcf * Rot * Rfc
    //         = mRotat * (Rcw*Rwf) * Rot * (Rfw*Rwc)     Note 'mRotat == Rwc'
    Ry.setRotation( 'Z', -left );
    Rp.setRotation( 'Y', +up   );
    Rot.mult( Ry, Rp );			// rotation in the temporary frame on the floor
    Rwc = mRotat;  Rcw.transpose( Rwc );
    Ry.mult( Rcw, Rwf );  Rp.mult( Rfw, Rwc );
    tmp.mult( Ry, Rot );  tmp.mult( tmp, Rp );
    mRotat.mult( mRotat, tmp );    	// update the camera pose
    makeRotatRotation();
    updateViewPoints();
  }
  
  // -----------------------------------------------------------------
  
private:
  void updateViewPoints(void) {
    // Update current view vectors 'vcFrom', 'vcAt', 'vcUp', 'vcLight',
    //   from initial view vectors and additional transformations 'mRotat' and 'vTPos'.
    CVec4D tmp, tmp2;
    //vcFrom =  + vTPos;
    tmp.sub( viFrom, viAt );  tmp.mult( fScale );
    tmp.add( viAt, tmp );
    multMV( tmp2, mRotat, tmp );
    vcFrom.add( tmp2, vTPos );          // vcFrom  = mRotat*(viAt+(viFrom-viAt)*fScale) + vTPos
    vcAt.add( viAt, vTPos );            // vcAt    = viAt + vTPos
    multMV( vcUp, mRotat, viUp );       // vcUp    = mRotat * viUp
    multMV( vcLight, mRotat, viLight ); // vcLight = mRotat * viLight   update positional light
  }
  
  void makeRotatRotation(void) {
    // Make 'mRotat' rotation by enforcing othonormality,
    //   in order to avoid the accumulation of numerical error.
    float  R[9], Rt[9], RtR[9], Tmp[9], error;
    GUIH_M3D_SET( R, mRotat[0][0], mRotat[0][1], mRotat[0][2], mRotat[1][0], mRotat[1][1], mRotat[1][2], mRotat[2][0], mRotat[2][1], mRotat[2][2] );
    for (int i = 0; i < 3; i++) {
      GUIH_M3D_TRANS( Rt, R );
      GUIH_M3D_MUL_MM( RtR, Rt, R );
      RtR[0] -= 1;  RtR[4] -= 1;  RtR[8] -= 1;
      error = GUIH_M3D_NORM( RtR );
      GUIH_M3D_MUL_MM( Tmp, R, RtR );
      for (int j = 0; j < 9; j++) R[j] = R[j] - 0.5f * Tmp[j];
      if (error < 1.0e-6) break;
    }
    mRotat.assign( R[0], R[1], R[2], 0.0,  R[3], R[4], R[5], 0.0,  R[6], R[7], R[8], 0.0,  0.0, 0.0, 0.0, 1.0 );
  }
  
  void getTemporaryFrameOnFloor(CMat4D &F) {
    // Calculate a temporary coordinate system (rotation only) on the floor,
    //   based on current camera orientation. (X:CameraForward, Y:CameraLeft, Z:FloorUp)
    // The resulting transformation is from world to temporary frame.
    float fv[3], lv[3], uv[3];
    // Note the 'mRotat' is rotation from camera frame to world frame, 
    //   i.e., its column vectors are camera axes in world frame.
    GUIH_V3D_SET  ( fv, -mRotat[0][2], -mRotat[1][2], -mRotat[2][2] );  // camera -Z direction
    GUIH_V3D_COPY ( uv, floor );	// up      vector in floor(world) frame
    GUIH_V3D_CROSS( lv, uv, fv );	// left    vector in floor(world) frame
    GUIH_V3D_NORMALIZE( lv );
    GUIH_V3D_CROSS( fv, lv, uv );	// forward vector in floor(world) frame
    F.assign( fv[0], lv[0], uv[0], 0.0,  fv[1], lv[1], uv[1], 0.0,  fv[2], lv[2], uv[2], 0.0,  0.0, 0.0, 0.0, 1.0 );
  }
  
  // -----------------------------------------------------------------

public:
  void printInfo(void) {
    std::cout << "GUIH::Camera3D Information" << std::endl;
    printf( "  From : %7.2f %7.2f %7.2f     (Initial: %7.2f %7.2f %7.2f)\n", vcFrom[0], vcFrom[1], vcFrom[2], viFrom[0], viFrom[1], viFrom[2]);
    printf( "  At   : %7.2f %7.2f %7.2f     (Initial: %7.2f %7.2f %7.2f)\n", vcAt[0], vcAt[1], vcAt[2], viAt[0], viAt[1], viAt[2]);
    printf( "  Up   : %7.2f %7.2f %7.2f     (Initial: %7.2f %7.2f %7.2f)\n", vcUp[0], vcUp[1], vcUp[2], viUp[0], viUp[1], viUp[2]);
    printf( "  Light: %7.2f %7.2f %7.2f     (Initial: %7.2f %7.2f %7.2f)\n", vcLight[0], vcLight[1], vcLight[2], viLight[0], viLight[1], viLight[2]);
    printf( "  Rot  : [ %5.2f %5.2f %5.2f ]    Trs : [ %6.2f ]    Scl : [ %6.2f ]\n", mRotat[0][0], mRotat[0][1], mRotat[0][2], vTPos[0], fScale );
    printf( "         [ %5.2f %5.2f %5.2f ]          [ %6.2f ]              \n", mRotat[1][0], mRotat[1][1], mRotat[1][2], vTPos[1] );
    printf( "         [ %5.2f %5.2f %5.2f ]          [ %6.2f ]              \n", mRotat[2][0], mRotat[2][1], mRotat[2][2], vTPos[2] );
  }
  
  void printForPOVRay(void) {
    printf("camera {\n");
    printf("  location <%f, %f, %f>\n", vcFrom[0], vcFrom[1], -vcFrom[2]);
    printf("  look_at  <%f, %f, %f>\n", vcAt[0],   vcAt[1],   -vcAt[2]);
    printf("  sky      <%f, %f, %f>\n", vcUp[0],   vcUp[1],   -vcUp[2]);
    printf("}\n");
    printf("light_source {\n");
    printf("  <%f, %f, %f> color White\n", vcLight[0], vcLight[1], -vcLight[2]);
    printf("}\n");
  }
  
  void ScreenToWorld( int w, int h, int sx, int sy, float from[3], float dir[3], float fov=60) {
    // interpret a point on screen as a ray from view point in world space
    // save the ray in 'from[]' and 'dir[]'.
    // Note that, in screen space, upper-left corner is (0,0).
    float x_max, y_max, dir2[3], xa[3], ya[3], za[3], M[9];
    // get view point
    GUIH_V3D_COPY( from, From );
    // get screen space direction
    x_max = sin( fov / 2  * 3.14159265358979323846f / 180 );
    y_max = x_max * ( (float)h / (float)w );
    dir2[0] = ((sx - (w/2.0f)) / (w/2.0f)) * x_max;
    dir2[1] = (((h/2.0f) - sy) / (h/2.0f)) * y_max;
    dir2[2] = -cos( fov / 2 * 3.14159265358979323846f / 180 );
    // convert the direction into world space
    GUIH_V3D_SUB( za, From, At );   GUIH_V3D_NORMALIZE( za );
    GUIH_V3D_CROSS( xa, Up, za );   GUIH_V3D_NORMALIZE( xa );
    GUIH_V3D_CROSS( ya, za, xa );   GUIH_V3D_NORMALIZE( ya );
    GUIH_M3D_SET( M,  xa[0], ya[0], za[0],  xa[1], ya[1], za[1],  xa[2], ya[2], za[2] );
    GUIH_M3D_MUL_MV( dir, M, dir2 );
    GUIH_V3D_NORMALIZE( dir );
    //cout << "dir2 " << GUIH_V3D_COUT(dir2) << endl;
  }
  
  // ===================================================================
  // render Axes
  // ===================================================================

  void renderXYZAxes(void) {
    int i;
    float step, pos;
  
    glDisable (GL_LIGHTING);
    glLineWidth(1.0);
  
    glBegin(GL_LINE_STRIP);		// X axis
    glColor3f(1, 0.2f, 0.2f);
    for (i = -10; i <= 10; i++) glVertex3f((float)i*i*i, 0, 0);
    glEnd();
    glBegin(GL_LINE_STRIP);		// Y axis
    glColor3f(0.2f, 1, 0.2f);
    for (i = -10; i <= 10; i++) glVertex3f(0, (float)i*i*i, 0);
    glEnd();
    glBegin(GL_LINE_STRIP);		// Z axis
    glColor3f(0.2f, 0.2f, 1);
    for (i = -10; i <= 10; i++) glVertex3f(0, 0, (float)i*i*i);
    glEnd();
  
    glBegin(GL_LINES);			// grids on the axes
    glColor3f(0.5f, 0.5f, 0.5f);
    step = 0.1f;
    while (step <= 10) {
      if      (step < 1) glColor3f(0.2f, 0.2f, 0.2f);
      else if (step > 1) glColor3f(0.2f, 0.2f, 0.5f);
      else               glColor3f(0.2f, 0.4f, 0.4f);
      for (i = -10; i <= +10; i++) {
	if (i == 0) continue;
	pos = i * step;
	glVertex3f( pos, +step/20, 0 );  glVertex3f( pos, -step/20, 0 ); // X
	glVertex3f( +step/20, pos, 0 );  glVertex3f( -step/20, pos, 0 ); // Y
	//glVertex3f( 0, +step/20, pos );  glVertex3f( 0, -step/20, pos ); // Z
	//glVertex3f( +step/20, 0, pos );  glVertex3f( -step/20, 0, pos ); // Z
      }
      step *= 10;
    }
    glEnd();
  
    glEnable (GL_LIGHTING);
  }

  // ===================================================================
  // render a tiny sphere at 'LookAt' position
  // ===================================================================

  void renderLookAtPosition(void) {
    glPushMatrix();
    glColor3f(0.5, 0, 0.5);
    glDisable (GL_LIGHTING);
    glPointSize(3);
    glBegin(GL_POINTS);
    glVertex3f(At[0], At[1], At[2]);
    glEnd();
    glPointSize(1);
    glEnable (GL_LIGHTING);
    glPopMatrix();
  }
  
  // ===================================================================
  // render the floor
  // ===================================================================
  
  int createFloorTexture(int texture_id=255) {
    int i, j, grid = 32;
    glBindTexture(GL_TEXTURE_2D, texture_id);
    GLubyte *tex = (GLubyte*) malloc(grid * grid * sizeof(GLubyte));
    for (i = 0; i < grid; i++)
      for (j = 0; j < grid; j++)
	tex[i*grid+j] = (GLubyte) ((((i>>4)&0x1)^((j>>4)&0x1)) * 255);
    glTexImage2D(GL_TEXTURE_2D, 0, 1, grid, grid, 0, GL_RED, GL_UNSIGNED_BYTE, tex);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    free(tex);
    return texture_id;
  }
  
  void renderFloor(float size, int grid, float *center=NULL, int axis=2) {
    float origin[3] = {0,0,0};
    if (center == NULL) center = origin;
    if (size <= 0) return;
    int  n = grid/2;
    static int checker_tex = 0;
    if (checker_tex == 0) { checker_tex = createFloorTexture(255); }
    glEnable(GL_TEXTURE_2D); // you can't turn on/off between begin and end calls
    glBindTexture(GL_TEXTURE_2D, checker_tex);
    glBegin(GL_QUADS);
    glColor3ub(100,100,100);
    switch (axis) {
    case 0:			// floor on X axis (YZ plane)
      glNormal3f(1, 0, 0);
      glTexCoord2i(0, 0);   glVertex3f(center[0], center[1]-size, center[2]+size);
      glTexCoord2i(n, 0);   glVertex3f(center[0], center[1]-size, center[2]-size);
      glTexCoord2i(n, n);   glVertex3f(center[0], center[1]+size, center[2]-size);
      glTexCoord2i(0, n);   glVertex3f(center[0], center[1]+size, center[2]+size);
      break;
    case 1:			// floor on Y axis (XZ plane)
      glNormal3f(0, 1, 0);
      glTexCoord2i(0, 0);   glVertex3f(center[0]-size, center[1], center[2]+size);
      glTexCoord2i(n, 0);   glVertex3f(center[0]+size, center[1], center[2]+size);
      glTexCoord2i(n, n);   glVertex3f(center[0]+size, center[1], center[2]-size);
      glTexCoord2i(0, n);   glVertex3f(center[0]-size, center[1], center[2]-size);
      break;
    case 2:			// floor on Z axis (XY plane)
      glNormal3f(0, 0, 1);
      glTexCoord2i(0, 0);   glVertex3f(center[0]-size, center[1]-size, center[2]);
      glTexCoord2i(n, 0);   glVertex3f(center[0]+size, center[1]-size, center[2]);
      glTexCoord2i(n, n);   glVertex3f(center[0]+size, center[1]+size, center[2]);
      glTexCoord2i(0, n);   glVertex3f(center[0]-size, center[1]+size, center[2]);
      break;
    }
    glEnd();
    glDisable(GL_TEXTURE_2D);
  }
  
 private:
  void multVM( CVec4D r, CVec4D v, CMat4D m ) {  // row_vector * matrix
    float  t[4];  memcpy(t, v.A, 4*sizeof(float));  // It's OK to use same 'r' and 'v'
    r.A[0] = t[0] * m.A[0] + t[1] * m.A[4] + t[2] * m.A[8] + t[3] * m.A[12];
    r.A[1] = t[0] * m.A[1] + t[1] * m.A[5] + t[2] * m.A[9] + t[3] * m.A[13];
    r.A[2] = t[0] * m.A[2] + t[1] * m.A[6] + t[2] * m.A[10] + t[3] * m.A[14];
    r.A[3] = t[0] * m.A[3] + t[1] * m.A[7] + t[2] * m.A[11] + t[3] * m.A[15];
  }
  void multMV( CVec4D r, CMat4D m, CVec4D v ) {  // matrix * col_vector
    float  t[4];  memcpy(t, v.A, 4*sizeof(float));  // It's OK to use same 'r' and 'v'
    r.A[0] = m.A[0] * t[0] +m.A[1] * t[1] +m.A[2] * t[2] +m.A[3] * t[3];
    r.A[1] = m.A[4] * t[0] +m.A[5] * t[1] +m.A[6] * t[2] +m.A[7] * t[3];
    r.A[2] = m.A[8] * t[0] +m.A[9] * t[1] +m.A[10] * t[2] +m.A[11] * t[3];
    r.A[3] = m.A[12] * t[0] +m.A[13] * t[1] +m.A[14] * t[2] +m.A[15] * t[3];
  }
  
};

  
}	// end of GUIH namespace


#endif	// GUIH_CAMERA3D_HPP
  
