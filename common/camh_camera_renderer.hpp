
// 
// CAMH::CameraRenderer<> class template
//
// Jaeil Choi
// last modified in Sep, 2009
//
// -------------------------------------------------------------------
// This file is a part of CAMH (Camera classes in Header) library
//   - CAMH::Camera<>		in camh_camera.hpp
//   - CAMH::CameraFileIO<>	in camh_camera_fileio.hpp
//   - CAMH::CameraRenderer<>	in camh_camera_renderer.hpp
//   - CAMH::Stereo<>		in camh_stereo.hpp
// 
// This class is designed for:
//   - visualization of camera pose (or world pose in camera frame)
//   - visualization of an arbitrary pose
//   - visualization of 2D/3D point matrix in camera frame
//   - visualization of a pose in camera frame
// 

#ifndef CAMH_CAMERA_RENDERER_HPP
#define CAMH_CAMERA_RENDERER_HPP

#include <iostream>
#include <cfloat>
#include "opengl_drawing.hpp"
#include "camh_camera.hpp"


namespace CAMH {

  
template<class T_t>
class CameraRenderer 
{
public:
  T_t			p_size;
private:
  OpenGL::Drawing<T_t>  oglr;
  
public:
  CameraRenderer() : p_size(0) { }
  ~CameraRenderer() { }
  void clear() { p_size = 0; }
  
  // -----------------------------------------------------------------
  // render camera
  // -----------------------------------------------------------------
public:
  void  renderCamera(CAMH::Camera<T_t> *c, T_t width) {
    if (c == NULL) return;
    glMatrixMode(GL_MODELVIEW);
    glPushMatrix();
    T_t Twc[3];  double M[4*4];
    c->getCameraPosition( Twc );
    G4M_SET( M,
	     c->Rcw[0], c->Rcw[1], c->Rcw[2], 0.0,
	     c->Rcw[3], c->Rcw[4], c->Rcw[5], 0.0,
	     c->Rcw[6], c->Rcw[7], c->Rcw[8], 0.0,
	     Twc[0],    Twc[1],    Twc[2],    1.0 );
    glMultMatrixd( M );
    // draw the camera (as in the camera coordinate system)
    //   located at (0,0,0), pointing to (0,0,1), upward(0,-1,0)
    renderCameraModel( width );
    glPopMatrix();
  }
  
  void  renderCameraAxes(CAMH::Camera<T_t> *c, T_t length) {
    // visualize camera coordinates system with 3 XYZ axis arrows
    if (c == NULL) return;
    T_t pos[3];
    c->getCameraPosition(pos);
    oglr.renderArrow(pos, c->Rcw+0, true, length, length/12, 255, 0, 0);
    oglr.renderArrow(pos, c->Rcw+3, true, length, length/12, 0, 255, 0);
    oglr.renderArrow(pos, c->Rcw+6, true, length, length/12, 0, 0, 255);
  }
  
  void  renderWorldAxesInCameraSpace(CAMH::Camera<T_t> *c, T_t length) {
    if (c == NULL) return;
    T_t R[9] = { c->Rcw[0], c->Rcw[3], c->Rcw[6], 
		 c->Rcw[1], c->Rcw[4], c->Rcw[7], 
		 c->Rcw[2], c->Rcw[5], c->Rcw[8] };
    oglr.renderArrow(c->Tcw, R+0, true, length, length/12, 255, 0, 0);
    oglr.renderArrow(c->Tcw, R+3, true, length, length/12, 0, 255, 0);
    oglr.renderArrow(c->Tcw, R+6, true, length, length/12, 0, 0, 255);
  }
  
  void renderCameraModel(T_t width) {
    T_t bw=width*0.5,   bh=width*0.75,  bl=width*1.5;
    T_t tw=width*0.25,  th=bh*0.1,      tl=width*0.5;
    glColor3ub(255,   0,   0);
    oglr.renderBox( -tw, -bh-th, -bl+tl, +tw, -bh+th, 0-tl ); // button on top
    glColor3ub(255, 255, 0);
    oglr.renderBox( -bw, -bh, -bl, +bw, +bh, 0.0 );           // body
    if (!oglr.quadric) oglr.init();
    gluCylinder( oglr.quadric, 0.00, bh, bh, 16, 1 );         // lens
    //oglr.renderBox( -0.01, -0.034, -0.04, +0.01, -0.03, -0.02 ); // button on top
    //oglr.renderBox( -0.02, -0.03, -0.06, +0.02, +0.03, +0.00 );  // body
    //gluCylinder( oglr.quadric, 0.00, 0.03, 0.03, 16, 1 );        // lens
  }
  
  // -----------------------------------------------------------------
  // render image (or target object)
  // -----------------------------------------------------------------
public:  
  void  renderCameraImage(CAMH::Camera<T_t> *c, T_t z_dist, int text_id) {
    // render the perspective image using texture in 3D space.
    if (c == NULL || text_id < 0) return;
    // calculate the 4 corners of the screen in 3D
    float  pos[3], ct[3], d, v0[3], v1[3], v2[3], v3[3];  T_t pos2[3];
    c->getCameraPosition(pos2);  
    G3V_SET( pos, (float)pos2[0], (float)pos2[1], (float)pos2[2] );
    G3V_SCALED_ADD( ct, pos, z_dist, c->Rcw + 6 );
    d = fabs(z_dist) / 1.732050;
    G3V_SCALED_ADD(v0, ct, -d, c->Rcw+0);  G3V_SCALED_ADD(v0, v0, +d, c->Rcw+3);
    G3V_SCALED_ADD(v1, ct, +d, c->Rcw+0);  G3V_SCALED_ADD(v1, v1, +d, c->Rcw+3);
    G3V_SCALED_ADD(v2, ct, +d, c->Rcw+0);  G3V_SCALED_ADD(v2, v2, -d, c->Rcw+3);
    G3V_SCALED_ADD(v3, ct, -d, c->Rcw+0);  G3V_SCALED_ADD(v3, v3, -d, c->Rcw+3);
    // bind texture
    glEnable(GL_TEXTURE_2D);
    glBindTexture(GL_TEXTURE_2D, text_id);
    // render the image on screen
    glBegin(GL_QUADS);
    glNormal3f(-c->Rcw[6], -c->Rcw[7], -c->Rcw[8]);
    glTexCoord2f(0.0, 0.0);   glVertex3fv(v0);
    glTexCoord2f(1.0, 0.0);   glVertex3fv(v1);
    glTexCoord2f(1.0, 1.0);   glVertex3fv(v2);
    glTexCoord2f(0.0, 1.0);   glVertex3fv(v3);
    glEnd();
    glDisable(GL_TEXTURE_2D);
    // render the border of the image
    if (false) {
      float values[4];
      glGetFloatv(GL_COLOR_CLEAR_VALUE, values);
      if (values[0] == 0) glColor3f(1.0, 1.0, 1.0);
      else                glColor3f(0.0, 0.0, 0.0);
      glDisable(GL_LIGHTING);
      glEnable(GL_POLYGON_OFFSET_LINE);
      glPolygonOffset (1.0, 1.0);
      glBegin(GL_LINE_LOOP);
      glVertex3fv(v0);  glVertex3fv(v1);  glVertex3fv(v2);  glVertex3fv(v3);
      glEnd();
      glDisable(GL_POLYGON_OFFSET_LINE);
      glEnable(GL_LIGHTING);
    }
  }
  
  int   prepareImageTexture(void *image, int w, int h, int row_stride, int text_id) {
    // 'img' must be a unsigned char array, with 3 RGB channels.
    // Prepare the texture of the image, and bind it to given texture id.
    // Note that 'w' and 'h' must be power of 2 (128, 256, 512, 1024, etc).
    if (image == NULL) return -1;
    unsigned char *img = (unsigned char*)image;
    // create RGB data in 'w'x'h' texture format (It's crude scaling)
    int   i, j, nw, nh;
    int   size[] = { 64, 128, 256, 512, 1024, 2048, 4096, 8192 };
    T_t scalex, scaley;
    unsigned char *data, *dst, *src;
    for (i = 0; i < 7; i++) if (size[i] >= w) { nw = size[i]; break; }
    for (i = 0; i < 7; i++) if (size[i] >= h) { nh = size[i]; break; }
    data = (unsigned char*) calloc(nw*nh, 3*sizeof(unsigned char));
    scalex = w / (T_t)nw;
    scaley = h / (T_t)nh;
    for (dst = data, i = 0; i < nh; i++) {
      for (j = 0; j < nw; j++, dst += 3) {
	if (i >= h || j >= w) {
	  dst[0] = dst[1] = dst[2] = 0;
	} else {
	  // for OpenGL textures, it should be from bottom to top
	  src = img + ( (h - 1 - int(i*scaley)) * row_stride + int(j*scalex) );
	  dst[0] = src[0];  dst[1] = src[1];  dst[2] = src[2];
	}
      }
    }
    // create texture of the perspective image
    glEnable(GL_TEXTURE_2D);
    glBindTexture(GL_TEXTURE_2D, text_id);
    //glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
    //glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
    //glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
    //glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR); // GL_NEAREST);
    glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_DECAL); // GL_MODULATE);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, w, h, 0, GL_RGB, GL_UNSIGNED_BYTE, data);
    glDisable(GL_TEXTURE_2D);
  
    free(data);
    return text_id;
  }
  
  // -----------------------------------------------------------------
  // render an arbitrary pose
  // -----------------------------------------------------------------
public:
  void  renderPose(T_t Rco[9], T_t Tco[3], T_t viz_size) {
    // Given an object pose in camera frame, 
    //   visualize the axes of the object frame in camera space.
    T_t R[9] = { Rco[0], Rco[3], Rco[6], 
		 Rco[1], Rco[4], Rco[7], 
		 Rco[2], Rco[5], Rco[8] };
    oglr.renderArrow(Tco, R+0, true, viz_size, viz_size/12, 255, 0, 0);
    oglr.renderArrow(Tco, R+3, true, viz_size, viz_size/12, 0, 255, 0);
    oglr.renderArrow(Tco, R+6, true, viz_size, viz_size/12, 0, 0, 255);
  }
  
#ifdef USE_MTX_MATRIX
  // -----------------------------------------------------------------
  // render 3D points
  // -----------------------------------------------------------------
public:
  void  render3DPointMatrix(MTX::Matrix<T_t> &m, int r=200, int g=200, int b=200) {
    if (m.nRow < 2 || m.nCol < 2) return;
    int   pidx, pos;
    T_t x, y, z;
  
    // decide point size
    if (p_size <= 0) p_size = DecidePointSize(m);
  
    // render point size
    glColor3ub( (unsigned char)r, (unsigned char)g, (unsigned char)b );
    for (pidx = 0; pidx < m.nCol; pidx++) {
      pos  = pidx;	x = m.data[pos];
      pos += m.nCol;	y = m.data[pos];
      if (m.nRow > 2) {
	pos += m.nCol;	z = m.data[pos];
	if (m.nRow == 4 && m.data[pos + m.nCol] != 1.0) continue;
      } else		z = 0.0;
      glPushMatrix();
      glTranslatef(x, y, z);
      oglr.renderCube(p_size);
      glPopMatrix();
    }
    glFlush();
  }
 private:
  T_t DecidePointSize(MTX::Matrix<T_t> &m) {
    T_t  xmin, xmax, ymin, ymax, xsize, ysize;
    xmin=+FLT_MAX;  xmax=-FLT_MAX;  ymin=+FLT_MAX;  ymax=-FLT_MAX;
    for (int pidx = 0; pidx < m.nCol; pidx++) {
      T_t x = m.data[pidx];
      T_t y = m.data[pidx + m.nCol];
      if (x < xmin) xmin = x;   if (x > xmax) xmax = x;
      if (y < ymin) ymin = y;   if (y > ymax) ymax = y;
    }
    xsize = xmax - xmin;  ysize = ymax - ymin;
    return ((xsize + ysize) / (8.0 * (int)sqrt(m.nCol)));
  }
  
  // -----------------------------------------------------------------
  // render 2D points
  // -----------------------------------------------------------------
 public:
  void  render2DPointMatrix(MTX::Matrix<T_t> &m, int height, int r=200, int g=200, int b=200, char type, int size) {
    if (m.nRow < 2 || m.nCol < 2) return;
    if (size <= 0) size = (int)DecidePointSize(m);	// decide point size
    glColor3ub( (unsigned char)r, (unsigned char)g, (unsigned char)b );
    for (int pidx = 0; pidx < m.nCol; pidx++) {
      T_t u = m.data[pidx];
      T_t v = m.data[pidx + m.nCol];
      render2DPoint( u, v, height, type, size );
      if (pidx == 0 || pidx == 1) oglr.drawCircle( u, (height-1)-v, 4 );
    }
    glFlush();
  }
  void  renderReprojections(MTX::Matrix<T_t> &m, CAMH::Camera<T_t> *cam, bool use_distortion, int r=200, int g=200, int b=200, char type, int size) {
    if (m.nRow < 2 || m.nCol < 2) return;
    if (size <= 0) size = (int)DecidePointSize(m);	// decide point size
    glColor3ub( (unsigned char)r, (unsigned char)g, (unsigned char)b );
    T_t xyz[3], uv[2];
    for (int pidx = 0; pidx < m.nCol; pidx++) {
      xyz[0] = m.data[pidx];
      xyz[1] = m.data[pidx + m.nCol];
      xyz[2] = m.data[pidx + 2 * m.nCol];
      cam->World2Pixel( xyz, uv, use_distortion);
      render2DPoint( uv[0], uv[1], cam->wh[1], type, size );
    }
    glFlush();
  }
#endif USE_MTX_MATRIX
  
 private:
  void  render2DPoint(T_t u, T_t v, int height, char type, int size) {
    // Note the different coordinate system:
    // camera  (0,0) +------> X   OpenGL:   Y ^        
    // &image:       |                        |
    //               |                        |
    //             Y V                  (0,0) +------> X 
    double du = u, dv = (height - 1) - v;
    double dsize = (type=='x' ? (size*1.4*0.5) : size*0.5);
    switch (type) {
    case '.':	// pixel
      glBegin(GL_POINTS);
      glVertex2d(du, dv);
      glEnd();
      break;
    case '+':	// cross
      glBegin(GL_LINES);
      glVertex2d(du-dsize, dv);  glVertex2d(du+dsize, dv);  
      glVertex2d(du, dv-dsize);  glVertex2d(du, dv+dsize);  
      glEnd();
      break;
    case 'x':	// 'x'
      glBegin(GL_LINES);
      glVertex2d(du-dsize, dv-dsize);  glVertex2d(du+dsize, dv+dsize);  
      glVertex2d(du-dsize, dv+dsize);  glVertex2d(du+dsize, dv-dsize);  
      glEnd();
      break;
    case 'r':	// rectangle
      glBegin(GL_QUADS);
      glVertex2d(du - dsize, dv - dsize);
      glVertex2d(du - dsize, dv + dsize);
      glVertex2d(du + dsize, dv + dsize);
      glVertex2d(du + dsize, dv - dsize);
      glEnd();
      break;
    default: break;
    }
  }
  
};


}	// end of CAMH namespace

#endif // CAMH_CAMERA_RENDERER_HPP

