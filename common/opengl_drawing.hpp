
//
// OpenGL::Drawing class 
//
// Jaeil Choi
// last modified in Mar, 2007
//
// general
//   void setQuadricDrawStyle(int style);
//   void setQuadricNormals(int normals);
// 2D
//   void drawRectangle(T_t x, T_t y, T_t width, T_t height);
//   void drawCircle(T_t x, T_t y, T_t radius, int slices=0, bool filled=true);
//   void drawArrow(T_t xy[2], T_t dir[2], bool outgoing=true, T_t length=0, T_t width=0);
// 3D
//   void renderLine(T_t p0[3], T_t p1[3], int line_width=0);
//   void renderCube(T_t radius);
//   void renderBox(T_t minx, T_t miny, T_t minz, T_t maxx, T_t maxy, T_t maxz);
//   void renderDisk(T_t in_radius, T_t out_radius);
//   void renderCone(T_t radius, T_t height, int slices=20, bool closed=false);
//   void renderSphere(T_t radius, int slices=20, int stacks=20);
//   void renderCylinder(T_t radius, T_t height, int slices=20, bool closed=false);
//   void renderTorus(bool solid, T_t o_radius, T_t i_radius, int o_steps, int i_steps);
//   void renderTetrahedron(bool solid, T_t radius);
//   void renderOctahedron(bool solid, T_t radius);
//   void renderDodecahedron(bool solid, T_t radius);
//   void renderIcosahedron(bool solid, T_t radius);
//    

#ifndef OPENGL_DRAWING_HPP
#define OPENGL_DRAWING_HPP

#include <iostream>
#include <cmath>
#ifdef WIN32					// Microsoft Windows
#include <windows.h>
#include <GL/gl.h>
#include <GL/glu.h>
#elif defined(__APPLE__) || defined(MACOSX)	// Apple OS X
// #include <AGL/agl.h>
#include <OpenGL/gl.h>
#include <OpenGL/glu.h>
// #include <GLUT/glut.h>
#else						// Linux
#include <GL/gl.h>
#include <GL/glu.h>
#endif

namespace OpenGL {
  
#define OGL2V_SET(r,a,b)     do { (r)[0] = (a);  (r)[1] = (b); } while(0)
#define OGL3V_SET(r,x,y,z)   do { (r)[0] = (x);  (r)[1] = (y);  (r)[2] = (z); } while(0)
#define OGL3V_COPY(d,s)      do { (d)[0] = (s)[0];  (d)[1] = (s)[1];  (d)[2] = (s)[2]; } while(0)
#define OGL3V_ADD(r,a,b)     do { (r)[0] = (a)[0] + (b)[0];  (r)[1] = (a)[1] + (b)[1];  (r)[2] = (a)[2] + (b)[2]; } while(0)
#define OGL3V_SUB(r,a,b)     do { (r)[0] = (a)[0] - (b)[0];  (r)[1] = (a)[1] - (b)[1];  (r)[2] = (a)[2] - (b)[2]; } while(0)
#define OGL3V_DIV_VALUE(a,v) do { (a)[0] /= (v);  (a)[1] /= (v);  (a)[2] /= (v); } while(0)
#define OGL3V_CROSS(r,a,b)   do { (r)[0] = (a)[1]*(b)[2]-(a)[2]*(b)[1]; (r)[1] = (a)[2]*(b)[0]-(a)[0]*(b)[2]; (r)[2] = (a)[0]*(b)[1]-(a)[1]*(b)[0]; } while(0)
#define OGL3V_LENGTH(a)      (sqrt((a)[0]*(a)[0] + (a)[1]*(a)[1] + (a)[2]*(a)[2]))
#define OGL3V_NORMALIZE(a)   do { float len = OGL3V_LENGTH(a);  if (len > 0) OGL3V_DIV_VALUE( a, len ); } while(0)
#define OGL4V_SET(r,a,b,c,d) do { (r)[0] = (a);  (r)[1] = (b);  (r)[2] = (c);  (r)[3] = (d); } while(0)
#define OGL4M_SET(m, a00,a01,a02,a03, a10,a11,a12,a13, a20,a21,a22,a23, a30,a31,a32,a33) \
                       do { (m)[0] = a00;  (m)[1] = a01;  (m)[2] = a02;  (m)[3] = a03;  \
                            (m)[4] = a10;  (m)[5] = a11;  (m)[6] = a12;  (m)[7] = a13;  \
                            (m)[8] = a20;  (m)[9] = a21;  (m)[10] = a22; (m)[11] = a23; \
                            (m)[12] = a30; (m)[13] = a31; (m)[14] = a32; (m)[15] = a33; } while(0)
  
template<class T_t>
class Drawing {
public:
  GLUquadric *quadric;
  
  Drawing() : quadric(NULL) { }
  ~Drawing() { if (quadric) gluDeleteQuadric(quadric); }
  
  void init(void) {
    if (!quadric) quadric = gluNewQuadric();
    gluQuadricDrawStyle(quadric, GLU_FILL); // GLU_FILL, GLU_LINE, GLU_SILHOUETTE, GLU_POINT
    gluQuadricNormals(quadric, GLU_FLAT);   // GLU_NONE, GLU_FLAT, GLU_SMOOTH
  }
  
public:
  void setColor(int r, int g, int b)  { glColor3ub( r, g, b ); }
  void setQuadricDrawStyle(int style) { if (!quadric) init(); gluQuadricDrawStyle(quadric, style); }
  void setQuadricNormals(int normals) { if (!quadric) init(); gluQuadricNormals(quadric, normals); }
  
  // =================================================================
  // 2D
  // =================================================================
public:
  void drawPoint(T_t x, T_t y) {
    //glPointSize(1.0f);
    glBegin(GL_POINTS);
    glVertex2d( x, y );
    glEnd();
  }
  void drawLine(T_t sx, T_t sy, T_t ex, T_t ey, bool solid=true) {
    //glLineWidth(1.0);
    if (!solid) { glEnable(GL_LINE_STIPPLE); glLineStipple(3, 0x6666); }
    glBegin(GL_LINES);
    glVertex2d( sx, sy );
    glVertex2d( ex, ey );
    glEnd();
    if (!solid) glDisable(GL_LINE_STIPPLE);
  }
  void drawRect(T_t x, T_t y, T_t width, T_t height, bool solid=true) {
    if (!solid) { glEnable(GL_LINE_STIPPLE); glLineStipple(3, 0x6666); }
    glBegin(GL_LINE_LOOP);
    glVertex2d( x, y );
    glVertex2d( x + width, y );
    glVertex2d( x + width, y + height );
    glVertex2d( x, y + height );
    glEnd();
    if (!solid) glDisable(GL_LINE_STIPPLE);
  }
  void drawRectCenteredAt(T_t x, T_t y, T_t width, T_t height, bool solid=true) {
    T_t hw = width/2, hh = height/2;
    if (!solid) { glEnable(GL_LINE_STIPPLE); glLineStipple(3, 0x6666); }
    glBegin(GL_LINE_LOOP);
    glVertex2d( x - hw, y - hh );
    glVertex2d( x + hw, y - hh );
    glVertex2d( x + hw, y + hh );
    glVertex2d( x - hw, y + hh );
    glVertex2d( x - hw, y - hh );
    glEnd();
    if (!solid) glDisable(GL_LINE_STIPPLE);
  }
  void drawQuad(T_t x0, T_t y0, T_t x1, T_t y1, T_t x2, T_t y2, T_t x3, T_t y3, bool solid=true) {
    if (!solid) { glEnable(GL_LINE_STIPPLE); glLineStipple(3, 0x6666); }
    glBegin(GL_LINE_LOOP);
    glVertex2d( x0, y0 );
    glVertex2d( x1, y1 );
    glVertex2d( x2, y2 );
    glVertex2d( x3, y3 );
    glEnd();
    if (!solid) glDisable(GL_LINE_STIPPLE);
  }
  void drawCircle(T_t x, T_t y, T_t radius, int slices=0, bool filled=false) {
    // render a circle
    T_t tx, ty, theta;
    if (slices < 4) slices = 32;
    glBegin( filled ? GL_POLYGON : GL_LINE_LOOP );
    for (int i = 0; i < slices; i++) {
      theta = i * (2 * M_PI) / slices;
      tx = x + (T_t)(radius * cos( theta ));
      ty = y + (T_t)(radius * sin( theta ));
      glVertex2d( tx, ty );
    }
    glEnd();
  }
  void drawArc(T_t x, T_t y, T_t radius, T_t angs, T_t ange, int slices=0, bool filled=false) {
    glBegin( filled ? GL_POLYGON : GL_LINE_STRIP );
    T_t tx, ty, theta, ang=(ange-angs);
    if (slices < 4) slices = 32;
    for (int i = 0; i <= slices; i++) {
      theta = (angs + i * ang / slices) * M_PI/180;
      tx = x + (T_t)(radius * cos( theta ));
      ty = y + (T_t)(radius * sin( theta ));
      glVertex2d( tx, ty );
    }
    glEnd();
  }
  void drawEllipse(T_t xy[2], T_t cov[2*2], int nsigma=1, bool filled=false) {
    // faster eigenvalue decomposition of 2x2 symmetric matrix A 
    // Jacobi method -- finding a rotation matrix that makes A diagonal
    //   [ c -s ] [ a b ] [ c +s ] = [  cca-2csb+ssd   ccb+cs(a-d)-ssb ]
    //   [ +s c ] [ b d ] [ -s c ]   [ ccb+cs(a-d)-ssb  ccd+2csb+ssa   ]
    //   ccb+cs(a-d)-ssb = 0  =>  1 + (a-d)/b t - tt = 0, where t = s/c
    //   Using smaller t, we get c = 1/sqrt(1+tt), s = ct
    //   Finally, evec = [ c +s ]  and  eval = [ cca-2csb+ssd ]
    //                   [ -s c ]              [ ccd+2csb+ssa ]
    double evec[2*2], eval[2];
    if (cov[1] == 0) {
      OGL4V_SET( evec, 1.0, 0.0, 0.0, 1.0 );
      OGL2V_SET( eval, cov[0], cov[3] );
    } else {
      double bc   = (cov[0] - cov[3]) / cov[1];
      double tanv = 0.5 * ( bc - sqrt( bc*bc + 4 ) );
      double cosv = 1.0 / sqrt( 1.0 + tanv * tanv );
      double sinv = cosv * tanv;
      double cc = cosv * cosv,  ss = sinv * sinv, cs = cosv * sinv;
      OGL4V_SET( evec, cosv, +sinv, -sinv, cosv );
      OGL2V_SET( eval, 
		 cc*cov[0] - 2*cs*cov[1] + ss*cov[3], 
		 cc*cov[3] + 2*cs*cov[1] + ss*cov[0] );
    }
    // render the ellipse
    int    i, n = 20;
    double wh[2], angle = 0, inc = 2*M_PI/n;
    double xl, yl, xw, yw;
    OGL2V_SET( wh, sqrt(eval[0])*nsigma, sqrt(eval[1])*nsigma );
    glBegin( filled ? GL_POLYGON : GL_LINE_LOOP );
    for (i = 0; i < n; i++, angle += inc) {
      xl = wh[0] * cos(angle);  yl = wh[1] * sin(angle);
      xw = xy[0] + evec[0] * xl + evec[1] * yl;	// [ x ] + [ c +s ] [ xl ]
      yw = xy[1] + evec[2] * xl + evec[3] * yl;	// [ y ]   [ -s c ] [ yl ]
      glVertex2d( xw, yw );
    }
    glEnd();
    // // draw the bounding box of the ellipse
    // double Si[2*2];
    // MTX::MatrixSolver<double> solver;
    // solver.inverseByDeterminant( 2, cov, Si );
    // int hw = (int)(nsigma / sqrt(Si[0] - Si[1] * Si[1] / Si[3]));
    // int hh = (int)(nsigma / sqrt(Si[3] - Si[1] * Si[1] / Si[0]));
    // win->drawRectangle( (int)uv[0]-hw, (int)uv[1]-hh, hw*2, hh*2 );
  }
  
  void drawCross(T_t x, T_t y, T_t size) {
    glBegin(GL_LINES);
    glVertex2d( x-size, y ); glVertex2d( x+size, y ); 
    glVertex2d( x, y-size ); glVertex2d( x, y+size ); 
    glEnd();
  }
  void drawArrow(T_t x0, T_t y0, T_t x1, T_t y1, T_t width=0) {
    T_t pa[2]={ x0, y0 }, dir[2];
    dir[0] = x1 - x0;  dir[1] = y1 - y0;
    drawArrow( pa, dir, true, 0, width );
  }
  void drawArrow(T_t xy[2], T_t dir[2], bool outgoing=true, T_t length=0, T_t width=0) {
    // render an arrow
    T_t ndir[2], len = sqrt(dir[0]*dir[0]+dir[1]*dir[1]);
    length = ( (length > 0) ? length : len );
    width = ( (width > 0) ? width  : length/4.0 );
    if (len == 0) return; else { ndir[0] = dir[0]/len; ndir[1] = dir[1]/len; }
    glPushMatrix();
    glTranslated(xy[0], xy[1], 0.0);
    double theta = acos(ndir[0] * 0 + ndir[1] * 1) * 180 / 3.14159265358979323846f;
    if (outgoing) {
      glRotated(theta, 0, 0, (ndir[0] > 0 ? -1 : +1));
    } else {
      glRotated(theta + 180, 0, 0, (ndir[0] > 0 ? -1 : +1));
      glTranslated(0, -length, 0);
    }
    glBegin(GL_LINES); glVertex2d(0, 0); glVertex2d(0, length); glEnd();
    glBegin(GL_POLYGON);	// arrow body (size: width/3)
    glVertex2d(-width/6, 0);  
    glVertex2d(+width/6, 0);
    glVertex2d(+width/6, length*0.66);
    glVertex2d(-width/6, length*0.66);  
    glEnd();
    glBegin(GL_TRIANGLES);	// arrow head (size: width)
    glVertex2d(0, length);
    glVertex2d(-width/2, length*0.66);
    glVertex2d(+width/2, length*0.66);
    glEnd();
    glPopMatrix();
  }
  
  void drawImage(int w, int h, void *data, char *type=NULL, bool with_points=false) {
    if (type==NULL) type = (char*)"RGB";
    int m=0, step=0;
    enum { PIX_GRAY=1, PIX_GRAYA, PIX_RGB, PIX_RGBA, PIX_FLOAT };
    if      (strcmp(type, "GRAY")==0)  { m = PIX_GRAY;  step = 1; }
    else if (strcmp(type, "GRAYA")==0) { m = PIX_GRAYA; step = 2; }
    else if (strcmp(type, "RGB")==0)   { m = PIX_RGB;   step = 3; }
    else if (strcmp(type, "RGBA")==0)  { m = PIX_RGBA;  step = 4; }
    else if (strcmp(type, "FLOAT")==0) { m = PIX_FLOAT; step = 4; }
    else return;
    if (with_points) glBegin(GL_POINTS);
    else             glBegin(GL_QUADS);
    unsigned char *pp = (unsigned char*)data;
    float fv=0;
    for (int y = h-1; y >= 0; y--) 
      for (int x = 0; x < w; x++, pp+=1) { 
	switch (m) {
	case PIX_GRAY : glColor3ub(pp[0], pp[0], pp[0]); break;
	case PIX_GRAYA: glColor4ub(pp[0], pp[0], pp[0], pp[1]); break;
	case PIX_RGB  : glColor3ub(pp[0], pp[1], pp[2]); break;
	case PIX_RGBA : glColor4ub(pp[0], pp[1], pp[2], pp[3]); break;
	case PIX_FLOAT: fv = ((float*)pp)[0]; glColor3f(fv, fv, fv); break;
	}
	if (with_points) glVertex2i(x, y);
	else {
	  glVertex2i(x  , y+1);  glVertex2i(x+1, y+1);
	  glVertex2i(x+1, y);    glVertex2i(x  , y);
	}
      }
    glEnd();
  }
  
  // =================================================================
  // 3D
  // =================================================================
public:
  void renderLine(T_t p0[3], T_t p1[3], int line_width=0) {
    ////if (linewidth > 0) glSetWidth(line_width);
    glDisable(GL_LIGHTING);
    glBegin(GL_LINES);
    glVertex3d( p0[0], p0[1], p0[2] );  
    glVertex3d( p1[0], p1[1], p1[2] );  
    glEnd();
    glEnable(GL_LIGHTING);
  }
  void renderCube(T_t radius, bool use_gllist=false) {
    static T_t last_radius = 0;
    if (use_gllist == false) {
      static double v[24], v0[24] = {
	-1.0, -1.0, -1.0,  +1.0, -1.0, -1.0,  +1.0, +1.0, -1.0,  -1.0, +1.0, -1.0,
	-1.0, -1.0, +1.0,  +1.0, -1.0, +1.0,  +1.0, +1.0, +1.0,  -1.0, +1.0, +1.0 };
      static double *n, n0[18] = {	// normal vectors of the faces
	0, -1, 0,   1, 0, 0,   		// lower/right side
	0, 1, 0,    -1, 0, 0,		// upper/left side
	0, 0, -1,   0, 0, +1 };		// bottom/top
      static int *f, f0[24] = {		// vertex indices of the faces
	0, 3, 15, 12,  3, 6, 18, 15,	// lower/right side
	6, 9, 21, 18,  9, 0, 12, 21,	// upper/left side
	9, 6, 3, 0,   12, 15, 18, 21 };	// bottom/top
      n = n0;  f = f0;
      if (radius != last_radius) {
	for (int i = 0; i < 24; i++) v[i] = v0[i] * radius;
	last_radius = radius;
      }
      glBegin(GL_QUADS);
      for (int i = 0; i < 6; i++, n+=3, f+=4) {
	glNormal3dv( n );
	glVertex3dv( v + f[0] );   glVertex3dv( v + f[1] );
	glVertex3dv( v + f[2] );   glVertex3dv( v + f[3] );
      }
      glEnd();
    } else {
      static GLuint gllist = 0;
      if (gllist && radius==last_radius) glCallList( gllist );
      else {
	if (gllist == 0) gllist = glGenLists(1);
	glNewList( gllist, GL_COMPILE_AND_EXECUTE );
	renderCube( radius, false );
	glEndList();
      }
    }
  }
  void renderCubeAt(double x, double y, double z, T_t radius, bool use_gllist=false) {
    glPushMatrix();
    glTranslated( x, y, z );
    renderCube( radius, use_gllist );
    glPopMatrix();
  }
  
  void renderBox(T_t minx, T_t miny, T_t minz, T_t maxx, T_t maxy, T_t maxz) {
    double v[24] = {
      minx, miny, minz,  maxx, miny, minz,  maxx, maxy, minz,  minx, maxy, minz,
      minx, miny, maxz,  maxx, miny, maxz,  maxx, maxy, maxz,  minx, maxy, maxz };
    double *n, n0[18] = {	// normal vectors of the faces
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
      glNormal3dv( n );
      glVertex3dv( v + f[0] );   glVertex3dv( v + f[1] );
      glVertex3dv( v + f[2] );   glVertex3dv( v + f[3] );
    }
    glEnd();
  }
  inline void renderDisk(T_t inner_radius, T_t outer_radius) {
    if (!quadric) init(); 
    gluDisk(quadric, inner_radius, outer_radius, 20, 2);
  }
  inline void renderCone(T_t radius, T_t height, int slices=20, bool closed=false) {
    if (!quadric) init(); 
    gluCylinder( quadric, radius, 0, height, slices, 2 );
    if (closed) gluDisk(quadric, 0, radius, 20, 2);
  }
  inline void renderSphere(T_t radius, int slices=20, int stacks=20) {
    if (!quadric) init(); 
    gluSphere( quadric, radius, slices, stacks );
  }
  void renderSphereAt(double x, double y, double z, T_t radius, int slices=20, int stacks=20) {
    glPushMatrix();
    glTranslated( x, y, z );
    renderSphere( radius, slices, stacks );
    glPopMatrix();
  }
  void renderCylinder(T_t radius, T_t height, int slices=20, bool closed=false) {
    if (!quadric) init(); 
    gluCylinder( quadric, radius, radius, height, slices, 2 );
    if (closed) {
      glPushMatrix();
      glRotatef(180,1,0,0); // So that normals point right way!
      gluDisk(quadric, 0, radius, 20, 2);
      glRotatef(180,1,0,0); // So that normals point right way!
      glTranslated(0,0, (double)height);
      gluDisk(quadric, 0, radius, 20, 2);
      glPopMatrix();
    }
  }
  
  void renderTorus(bool solid, T_t o_radius, T_t i_radius, int o_steps, int i_steps) {
    static GLuint gllist = 0;
    static T_t last_o_radius = 0, last_i_radius = 0;
    static int   last_o_steps  = 0, last_i_steps  = 0;
    if (o_radius != last_o_radius || o_steps != last_o_steps || 
	i_radius != last_i_radius || i_steps != last_i_steps ) {
      if (gllist == 0) gllist = glGenLists(1);
      last_o_radius = o_radius;  last_o_steps = o_steps;
      last_i_radius = i_radius;  last_i_steps = i_steps;
      glNewList( gllist, GL_COMPILE_AND_EXECUTE );
      glBegin(GL_QUADS);
      double a, b, n[3], v[3];
      double astep = M_PI*2 / o_steps,  aoff = astep * 0.5;
      double bstep = M_PI*2 / i_steps,  boff = bstep * 0.5;
      for (int i = 0; i < o_steps; i++) {
	a = i * astep;
	for (int j = 0; j < i_steps; j++) {
	  b = j * bstep;
	  n[0] = cos(a) * cos(b);	// nx
	  n[1] = sin(a) * cos(b);	// ny
	  n[2] = sin(b);		// nz
	  glNormal3dv( n );
	  double a_nb[4] = { a-aoff, a+aoff, a+aoff, a-aoff };
	  double b_nb[4] = { b-boff, b-boff, b+boff, b+boff };
	  for (int k = 0; k < 4; k++) {
	    v[0] = cos(a_nb[k]) * o_radius + cos(a_nb[k]) * cos(b_nb[k]) * i_radius; // x
	    v[1] = sin(a_nb[k]) * o_radius + sin(a_nb[k]) * cos(b_nb[k]) * i_radius; // y
	    v[2] = sin(b_nb[k]) * i_radius;  // z
	    glVertex3dv( v );
	  }
	}
      }
      glEnd();
      glEndList();
    } else glCallList( gllist );
  }

  void renderTetrahedron(bool solid, T_t radius) {
    static GLuint gllist = 0;
    static T_t last_radius = 0;
    if (radius != last_radius) {
      static double v[12], v0[12] = {
	-1.0000f, -1.0000f, -1.0000f,     +1.0000f, +1.0000f, -1.0000f, 
	+1.0000f, -1.0000f, +1.0000f,     -1.0000f, +1.0000f, +1.0000f
      };
      static double n[12] = {
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
	glNormal3dv(n+i);
	glVertex3dv(v+f[i+0]);  glVertex3dv(v+f[i+1]);  glVertex3dv(v+f[i+2]);
      }
      glEnd();
      glEndList();
    } else glCallList( gllist );
  }

  void renderOctahedron(bool solid, T_t radius) {
    static GLuint gllist = 0;
    static T_t last_radius = 0;
    if (radius != last_radius) {
      static double v[18], v0[18] = {
	+1.0000f, +0.0000f, +0.0000f,     +0.0000f, -1.0000f, +0.0000f, 
	-1.0000f, +0.0000f, +0.0000f,     +0.0000f, +1.0000f, +0.0000f, 
	+0.0000f, +0.0000f, +1.0000f,     +0.0000f, +0.0000f, -1.0000f
      };
      static double n[24] = {
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
	glNormal3dv(n+i);
	glVertex3dv(v+f[i+0]);  glVertex3dv(v+f[i+1]);  glVertex3dv(v+f[i+2]);
      }
      glEnd();
      glEndList();
    } else glCallList( gllist );
  }

  void renderDodecahedron(bool solid, T_t radius) {
    static GLuint gllist = 0;
    static T_t last_radius = 0;
    if (radius != last_radius) {
      static double v[60], v0[60] = {
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
      static double *n, n0[36] = {
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
	glNormal3dv( n );
	glBegin(GL_POLYGON);
	glVertex3dv(v+f[0]);  glVertex3dv(v+f[1]);  glVertex3dv(v+f[2]);  
	glVertex3dv(v+f[3]);  glVertex3dv(v+f[4]);  
	glEnd();
      }
      glEndList();
    } else glCallList( gllist );
  }

  void renderIcosahedron(bool solid, T_t radius) {
    static GLuint gllist = 0;
    static T_t last_radius = 0;
    if (radius != last_radius) {
      static double v[36], v0[36] = {	// 12 vertices
	+0.0000f, -0.5257f, +0.8507f,     +0.8507f, +0.0000f, +0.5257f, 
	+0.8507f, +0.0000f, -0.5257f,     -0.8507f, +0.0000f, -0.5257f, 
	-0.8507f, +0.0000f, +0.5257f,     -0.5257f, +0.8507f, +0.0000f, 
	+0.5257f, +0.8507f, +0.0000f,     +0.5257f, -0.8507f, +0.0000f, 
	-0.5257f, -0.8507f, +0.0000f,     +0.0000f, -0.5257f, -0.8507f, 
	+0.0000f, +0.5257f, -0.8507f,     +0.0000f, +0.5257f, +0.8507f
      };
      static double n[60] = {		// 20 faces
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
	glNormal3dv(n+i);
	glVertex3dv(v+f[i+0]);  glVertex3dv(v+f[i+1]);  glVertex3dv(v+f[i+2]);
      }
      glEnd();
      glEndList();
    } else glCallList( gllist );
  }

  void renderArrow(T_t from[3], T_t dir[3], bool outgoing=true, 
		   T_t len=0, T_t radius=0, int r=-1, int g=-1, int b=-1) {
    if (!quadric) init(); 
    if (len <= 0)    len    = OGL3V_LENGTH( dir );
    if (len == 0) return;
    if (radius <= 0) radius = len * 0.05;
    double bodylen = len * 0.6;
    double headlen = len * 0.4;
    double radius2 = radius * 1.5;
    int    slices = 20;
    if (r >= 0 && g >= 0 && b >= 0) glColor3ub( r, g, b );
    
    // calculate local coordinate system where dir[3] becomes z-axis
    double lz[3], lx[3], ly[3], tmp[16];
    OGL3V_COPY( lz, dir );
    OGL3V_NORMALIZE( lz );
    if (lz[0] == 1.0) OGL3V_SET( tmp, 0, 1, 0 );
    else              OGL3V_SET( tmp, 1, 0, 0 );
    OGL3V_CROSS( lx, lz, tmp );  OGL3V_NORMALIZE( lx );
    OGL3V_CROSS( ly, lz, lx  );  OGL3V_NORMALIZE( ly );
    OGL4M_SET( tmp,
	     lx[0], lx[1], lx[2], 0.0, 
	     ly[0], ly[1], ly[2], 0.0, 
	     lz[0], lz[1], lz[2], 0.0, 
	     from[0], from[1], from[2], 1.0 );
    // draw the arrow pointing in the z-axis direction
    glPushMatrix();
    glMultMatrixd( tmp );
    gluCylinder( quadric, radius, radius, bodylen, slices, 2 );	// body
    gluDisk( quadric, 0, radius, 20, 2 );			// bottom
    glPushMatrix();
    glTranslated( 0, 0, bodylen );
    gluCylinder( quadric, radius2, 0, headlen, slices, 2 );	// head
    gluDisk( quadric, 0, radius2, 20, 2 );			// neck
    glPopMatrix();
    glPopMatrix();
  }
  
  void renderAxes(T_t len) {
    T_t o[3]={0,0,0}, x[3]={1,0,0}, y[3]={0,1,0}, z[3]={0,0,1};
    renderArrow( o, x, true, len, 0, 255, 0, 0 );
    renderArrow( o, y, true, len, 0, 0, 255, 0 );
    renderArrow( o, z, true, len, 0, 0, 0, 255 );
  }
};
  
}	// namespace OpenGL

#endif  // OPENGL_DRAWING_HPP
