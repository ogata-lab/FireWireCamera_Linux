
// 
// CAMH::Stereo<> class template
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
//   - setting up and using a stereo camera
//   - conversion betwen left and right camera
// 
// StereoCamera coordinate system is defined as:
//   X-axis : direction from Left camera center to Right camera center
//   Y-axis : perpendicular to X and Left Camera's Z (cross-product of them)
//   Z-axis : perpendicular to X, on the plane defind by X and Left camera's Z
//   origin : at the center of Left and Right camera centers
// 

#ifndef CAMH_STEREO_HPP
#define CAMH_STEREO_HPP

#include <iostream>
#include "vm_macros.h"
#include "camh_camera.hpp"

#ifdef USE_GUIH_CONFIG
#include "guih_config.hpp"
#endif
#ifdef USE_GEOMETRY_3D
#include "geometry3d.hpp"
#endif

namespace CAMH {

template <class T_t>
class Stereo 
{
public:
  // extrinsic parameters of the stereo camera
  int		wh[2];		// image resolution
  T_t		Rcw[9], Tcw[3];
  
  // stereo 
  Camera<T_t>	*camL;		// left  camera
  Camera<T_t>	*camR;		// right camera
  
private:
  T_t		baseline;       // 
  Camera<T_t>	camL_dummy;	// left  camera dummy place holder
  Camera<T_t>	camR_dummy;	// right camera dummy place holder
  T_t		Rqc[9], Tqc[3];	// xform from camera frame to left  cam
  T_t		Rpc[9], Tpc[3];	// xform from camera frame to right cam
public:
  T_t		stmRT[6];	// affine xform (R[4] and T[2]) for stereo matching
  
public:
  Stereo() : camL(&camL_dummy), camR(&camR_dummy), baseline(0) { 
    G3M_SET( Rcw, 1, 0, 0, 0, 1, 0, 0, 0, 1 );
    G3V_SET( Tcw, 0, 0, 0 );
  }
  Stereo(Camera<T_t> *cL, Camera<T_t> *cR, T_t baseline=0) { 
    setCamera( cL, cR );
    if (baseline > 0) setStereoTransformWithBaseline(baseline);
    G3M_SET( Rcw, 1, 0, 0, 0, 1, 0, 0, 0, 1 );
    G3V_SET( Tcw, 0, 0, 0 );
  }
  ~Stereo() {}
  inline bool isReady(void) { return (baseline>0 && camL->fc[0]>0 && camR->fc[0]>0); }
  inline void setCamera(Camera<T_t> *cL, Camera<T_t> *cR) {
    camL = (cL ? cL : &camL_dummy);
    camR = (cR ? cR : &camR_dummy);
  }
  
  // -----------------------------------------------------------------
  // 
  // -----------------------------------------------------------------
public:
  bool setupKnownCamera(const char *name) {
    double Rqc[9], Tqc[3], Rpc[9], Tpc[3];
    if (!name || !name[0]) return false;
    else if (strcmp(name,"BBOU0")==0 || strcmp(name,"Bumblebee_OsakaUniv_0")==0) {
      setStereoIntrinsic( 640, 480,
			  552.8220, 564.9148, 0.0, 322.1500, 242.6670, 9.42823e-07,
			  549.6514, 560.0063, 0.0, 319.9728, 239.9296, 9.4221e-07 );
      G3M_SET( Rqc, 0.9996, 0.0045, -0.0291, -0.0045, 1.0000, 0.0001, 0.0291, -0.0000, 0.9996);
      G3V_SET( Tqc, 0.0608, -0.0003, 0.0018 );
      G3M_SET( Rpc, 0.9998, 0.0105, -0.0161, -0.0104, 0.9999, 0.0064, 0.0162, -0.0062, 0.9998);
      G3V_SET( Tpc, -0.0608, 0.0006, -0.0010 );
      setStereoTransformation( Rqc, Tqc, Rpc, Tpc );
      G6V_SET( stmRT, 0.9921, 0.0091, -0.0091, 0.9921,  0.00, 1.54 );
    } else if (strcmp(name,"BBOU1")==0 || strcmp(name,"Bumblebee_OsakaUniv_1")==0) {
      setStereoIntrinsic( 640, 480,	// Model: Bumblebee2 (BB2-03s2c)  S/N = 6345059
			  539.9279, 547.3207, 0.0, 308.7435, 232.7912, 1.23766e-06,
			  536.6165, 547.6781, 0.0, 320.2642, 240.2368, 1.13871e-06 );
      G3M_SET( Rqc, 0.9991, -0.0080, -0.0417, 0.0080, 1.0000, -0.0003, 0.0417, 0.0000, 0.9991);
      G3V_SET( Tqc, 0.0603, 0.0005, 0.0025 );
      G3M_SET( Rpc, 0.9986, -0.0133, -0.0512, 0.0124, 0.9998, -0.0181, 0.0514, 0.0174, 0.9985);
      G3V_SET( Tpc, -0.0602, -0.0007, -0.0031 );
      setStereoTransformation( Rqc, Tqc, Rpc, Tpc );
      G6V_SET( stmRT, 1.0115, -0.0025, 0.0025, 1.0115,  0.00, -2.69 );
    } else if (strcmp(name,"BBOU7")==0 || strcmp(name,"Bumblebee_OsakaUniv_7")==0) {
      setStereoIntrinsic( 640, 480,	// Model: Bumblebee2 (BB2-03s2c)  S/N = 8220020
			  536.2552, 542.3329, 0.0, 320.0111, 240.0174, 1.17436e-06,
			  537.9141, 543.9839, 0.0, 319.9916, 239.9882, 1.02629e-06 );
      G3M_SET( Rqc, 0.9984, 0.0031, 0.0556, -0.0031, 1.0000, -0.0002, -0.0556, 0.0000, 0.9985);
      G3V_SET( Tqc, 0.0602, -0.0002, -0.0034 );
      G3M_SET( Rpc, 0.9983, 0.0004, 0.0588, -0.0012, 0.9999, 0.0134, -0.0588, -0.0134, 0.9982);
      G3V_SET( Tpc, -0.0602, 0.0001, 0.0035 );
      setStereoTransformation( Rqc, Tqc, Rpc, Tpc );
      G6V_SET( stmRT, 1.0026, -0.0026, 0.0026, 1.0026, 0.00, 7.20 );
    } else return false;
    return true;
  }
  
  // -----------------------------------------------------------------
  // set the intrinsic parameters of the left/right cameras
  // -----------------------------------------------------------------
public:
  void setStereoIntrinsic(int imgw, int imgh,
			  T_t fxL, T_t fyL, T_t afL, T_t cxL, T_t cyL, T_t dtL,
			  T_t fxR, T_t fyR, T_t afR, T_t cxR, T_t cyR, T_t dtR) {
    G2V_SET( wh, imgw, imgh );
    camL->setIntrinsic( imgw, imgh, fxL, fyL, afL, cxL, cyL, dtL );
    camR->setIntrinsic( imgw, imgh, fxR, fyR, afR, cxR, cyR, dtR );
  }
  void changeStereoResolution(double ratio) {
    G2V_SET( wh, (int)(wh[0]*ratio), (int)(wh[1]*ratio) );
    camL->changeCameraResolution( ratio );
    camR->changeCameraResolution( ratio );
  }
  
  // -----------------------------------------------------------------
  // set the extrinsic pose of the left/right cameras
  // -----------------------------------------------------------------
public:
  void setStereoExtrinsic(T_t Rcw[9], T_t Tcw[3]) {
    if (Rcw != this->Rcw) G3M_COPY( this->Rcw, Rcw );
    if (Tcw != this->Tcw) G3V_COPY( this->Tcw, Tcw );
    G3M_XFORM_MERGE( camL->Rcw, camL->Tcw,  Rqc, Tqc,  Rcw, Tcw );
    G3M_XFORM_MERGE( camR->Rcw, camR->Tcw,  Rpc, Tpc,  Rcw, Tcw );
  }
  void setStereoExtrinsicWithCameraPose(T_t Rwc[9], T_t Twc[3]) {
    G3M_XFORM_INVERSE( Rcw, Tcw, Rwc, Twc );
    setStereoExtrinsic( Rcw, Tcw );
  }
  
  // -----------------------------------------------------------------
  // set relative pose of left/right camera
  // -----------------------------------------------------------------
public:
  void setStereoTransformWithBaseline(T_t baseline) {
    // Set the relative pose of left/right camera, 
    //   using baseline size and assuming the cameras are aligned perfectly.
    this->baseline = baseline;
    G3M_SET_ID( Rqc );  G3V_SET( Tqc, +baseline/2, 0, 0 );  // camera -> left
    G3M_SET_ID( Rpc );  G3V_SET( Tpc, -baseline/2, 0, 0 );  // camera -> right
    T_t Rcw[9]={1,0,0, 0,1,0, 0,0,1}, Tcw[3]={0,0,0};
    setStereoExtrinsic( Rcw, Tcw );
  }
  void setStereoTransformWithCameraXForm(T_t Rqw[9], T_t Tqw[3], T_t Rpw[9], T_t Tpw[3]) {
    // Set the relative pose of left/right camera, using absolute pose of each camera.
    double Rwq[9], Twq[3], Rwp[9], Twp[3];
    G3M_XFORM_INVERSE( Rwq, Twq, Rqw, Tqw );
    G3M_XFORM_INVERSE( Rwp, Twp, Rpw, Tpw );
    double zL[3]; G3V_SET( zL, Rwq[2], Rwq[5], Rwq[8] );
    this->baseline = G3V_DISTANCE( Twq, Twp );
    double x[3], y[3], z[3], Rwc[9], Twc[3];      // Stereo camera coordinate system
    G3V_SUB( x, Twp, Twq );  G3V_NORMALIZE( x );  // X-axis: direction from Left camera center to Right camera center
    G3V_CROSS( y, zL, x );   G3V_NORMALIZE( y );  // Y-axis: perpendicular to X and Left Camera's Z (cross-product of them)
    G3V_CROSS( z,  x, y );   G3V_NORMALIZE( z );  // Z-axis: perpendicular to X, on the plane defind by X and Left camera's Z
    G3V_SCALED_ADD( Twc, Twq, baseline/2, x );    // origin: at the center of Left and Right camera centers
    G3M_SET( Rwc,  x[0], y[0], z[0],  x[1], y[1], z[1],  x[2], y[2], z[2] );
    G3M_XFORM_INVERSE( this->Rcw, this->Tcw, Rwc, Twc );
    G3M_XFORM_MERGE( this->Rqc, this->Tqc, Rqw, Tqw, Rwc, Twc );
    G3M_XFORM_MERGE( this->Rpc, this->Tpc, Rpw, Tpw, Rwc, Twc );
    makeRotationRotation( this->Rqc );
    makeRotationRotation( this->Rpc );
    setStereoExtrinsic( Rcw, Tcw );
  }
  void setStereoTransformWithCameraPoses(T_t Rwq[9], T_t Twq[3], T_t Rwp[9], T_t Twp[3]) {
    // Set the relative pose of left/right camera, using absolute pose of each camera.
    double zL[3];  G3V_SET( zL, Rwq[2], Rwq[5], Rwq[8] );
    this->baseline = G3V_DISTANCE( Twq, Twp );
    double x[3], y[3], z[3], Rwc[9], Twc[3];      // Stereo camera coordinate system
    G3V_SUB( x, Twp, Twq );  G3V_NORMALIZE( x );  // X-axis: direction from Left camera center to Right camera center
    G3V_CROSS( y, zL, x );   G3V_NORMALIZE( y );  // Y-axis: perpendicular to X and Left Camera's Z (cross-product of them)
    G3V_CROSS( z,  x, y );   G3V_NORMALIZE( z );  // Z-axis: perpendicular to X, on the plane defind by X and Left camera's Z
    G3V_SCALED_ADD( Twc, Twq, baseline/2, x );    // origin: at the center of Left and Right camera centers
    G3M_SET( Rwc,  x[0], y[0], z[0],  x[1], y[1], z[1],  x[2], y[2], z[2] );
    double Rqw[9], Tqw[3], Rpw[9], Tpw[3];
    G3M_XFORM_INVERSE( Rqw, Tqw, Rwq, Twq );
    G3M_XFORM_INVERSE( Rpw, Tpw, Rwp, Twp );
    G3M_XFORM_INVERSE( this->Rcw, this->Tcw, Rwc, Twc );
    G3M_XFORM_MERGE( this->Rqc, this->Tqc, Rqw, Tqw, Rwc, Twc );
    G3M_XFORM_MERGE( this->Rpc, this->Tpc, Rpw, Tpw, Rwc, Twc );
    makeRotationRotation( this->Rqc );
    makeRotationRotation( this->Rpc );
    setStereoExtrinsic( Rcw, Tcw );
  }
  void setStereoTransformation(T_t Rqc[9], T_t Tqc[3], T_t Rpc[9], T_t Tpc[3]) {
    // Set the relative pose of left/right camera, using relative pose of each camera.
    memcpy( this->Rqc, Rqc, 9*sizeof(double) );
    memcpy( this->Tqc, Tqc, 3*sizeof(double) );
    makeRotationRotation( this->Rqc );
    memcpy( this->Rpc, Rpc, 9*sizeof(double) );
    memcpy( this->Tpc, Tpc, 3*sizeof(double) );
    makeRotationRotation( this->Rpc );
    this->baseline = G3V_LENGTH( Tqc ) + G3V_LENGTH( Tpc );
    setStereoExtrinsic( Rcw, Tcw );
  }
  
  void getCameraPose (T_t Rwc[9], T_t Twc[3]) { G3M_XFORM_INVERSE( Rwc, Twc, Rcw, Tcw ); }
  void getCameraPoseL(T_t Rcq[9], T_t Tcq[3]) { G3M_XFORM_INVERSE( Rcq, Tcq, Rqc, Tqc ); }
  void getCameraPoseR(T_t Rcp[9], T_t Tcp[3]) { G3M_XFORM_INVERSE( Rcp, Tcp, Rpc, Tpc ); }
  
  void getPitchYawRoll(T_t pyr[3], T_t Rrc[9], T_t Rwr[9]=NULL) {
    // Calculate pitch/yaw/roll angles (in degrees) of camera.
    //   Camera pose is given by 'Rrc[9]' in a reference frame 'Rwr[9]'.
    //   Pitch : rotation along X axis. positive when it looks upward (opposite to Enon)
    //   Yaw   : rotation along Y axis. positive when it turns to the right
    //   Roll  : rotation along Z axis. positive when it tilts to the right (clock-wise)
    T_t cX[3]={Rrc[0],Rrc[3],Rrc[6]};  G3V_NORMALIZE( cX );
    T_t cY[3]={Rrc[1],Rrc[4],Rrc[7]};  G3V_NORMALIZE( cY );
    T_t cZ[3]={Rrc[2],Rrc[5],Rrc[8]};  G3V_NORMALIZE( cZ );
    T_t X[3]={1,0,0}, Y[3]={0,1,0}, Z[3]={0,0,1}, nv[3], s, dot;
    if (Rwr) {
      G3V_SET( X, Rwr[0], Rwr[3], Rwr[6] );
      G3V_SET( Y, Rwr[1], Rwr[4], Rwr[7] );
      G3V_SET( Z, Rwr[2], Rwr[5], Rwr[8] );
    }
    // for pitch: project cZ on ZX, dot (cZ.cZ') with sign (cZ'.-Y)
    // for yaw  : project cX on XY, dot (cX.cX') with sign (cX'.-Z)
    // for roll : project cY on YZ, dot (cY.cY') with sign (cY'.-X)
    s   = G3V_DOT(cZ, Y);  G3V_SCALED_ADD(nv, cZ,-s,Y);  G3V_NORMALIZE(nv);
    dot = G3V_DOT(nv,cZ);  if (dot>1) dot = 1.0;
    pyr[0] = (G3V_DOT(nv,Y)<0 ? +1:-1) * acos( dot ) * 180/M_PI;
    s   = G3V_DOT(cX, Z);  G3V_SCALED_ADD(nv, cX,-s,Z);  G3V_NORMALIZE(nv);
    dot = G3V_DOT(nv,cX);  if (dot>1) dot = 1.0;
    pyr[1] = (G3V_DOT(nv,Z)<0 ? +1:-1) * acos( dot ) * 180/M_PI;
    s   = G3V_DOT(cY, X);  G3V_SCALED_ADD(nv, cY,-s,X);  G3V_NORMALIZE(nv);
    dot = G3V_DOT(nv,cY);  if (dot>1) dot = 1.0;
    pyr[2] = (G3V_DOT(nv,X)<0 ? +1:-1) * acos( dot ) * 180/M_PI;
  }
  void makeRotationRotation(T_t R[9]) {
    T_t  RtR[9], Re[9];
    for (int i = 0; i < 10; i++) {
      G3M_MUL_MtM( RtR, R, R );
      RtR[0] -= 1;  RtR[4] -= 1;  RtR[8] -= 1;
      G3M_MUL_MM( Re, R, RtR );
      for (int j = 0; j < 9; j++) R[j] = R[j] - 0.5 * Re[j];
      if (G3M_NORM( RtR ) < 1.0e-6) break;
    }
  }
  
  // -----------------------------------------------------------------
  // Stereo Geometry
  // -----------------------------------------------------------------
  void convertLeft2Right(T_t uvL[2], T_t uvR[2], T_t distance) {
    // Convert left image position 'uvL[2]' to right image 'uvR[2]', with given distance.
    T_t xyz[3], dir[3];
    camL->Pixel2WorldRay( uvL, dir, xyz );
    G3V_SCALED_ADD( xyz, xyz, distance, dir );
    camR->World2Pixel( xyz, uvR );
  }
  
  void convertRight2Left(T_t uvR[2], T_t uvL[2], T_t distance) {
    // Convert right image position 'uvR[2]' to left image 'uvL[2]', with given distance.
    T_t xyz[3], dir[3];
    camR->Pixel2WorldRay( uvR, dir, xyz );
    G3V_SCALED_ADD( xyz, xyz, distance, dir );
    camL->World2Pixel( xyz, uvL );
  }
  
#ifdef USE_GEOMETRY_3D
  bool triangulatePosition(T_t uvL[2], T_t uvR[2], T_t xyz[3]) {
    T_t xyzL[3], xyzR[3], dirL[3], dirR[3];
    camL->Pixel2WorldRay( uvL, dirL, xyzL );
    camR->Pixel2WorldRay( uvR, dirR, xyzR );
//     printf("RayL (%.2f %.2f %.2f)->(%.2f %.2f %.2f)  RayR (%.2f %.2f %.2f)->(%.2f %.2f %.2f)\n", 
// 	   xyzL[0], xyzL[1], xyzL[2], dirL[0], dirL[1], dirL[2], 
// 	   xyzR[0], xyzR[1], xyzR[2], dirR[0], dirR[1], dirR[2]);
    Geometry3D<double> geom;
    int mode = -1;   // -1: OnLeftLine,  +1: OnRightLine,  0: InBetween
    int ret = geom.getPointByIntersectingTwoLines( xyzL, dirL, xyzR, dirR, xyz, mode );  //// new function
    return (ret > 0);
  }
  bool triangulateLine(T_t lineL[3], T_t lineR[3], T_t xyz[3], T_t dir[3]) {
    T_t planeL[4], planeR[4];
    convertImageLine2WorldPlane( &camL, lineL, planeL );
    convertImageLine2WorldPlane( &camR, lineR, planeR );
    Geometry3D<double> geom;
    return geom.getLineByIntersectingTwoPlanes(planeL, planeR, xyz, dir );  //// new function
  }
  void convertImageLine2WorldPlane(Camera<T_t> *cp, T_t line[3], T_t plane[4]) {
    // 
    T_t  tmp[2]={1,0}, uv1[2], uv2[2], cxyz[3], dir1[3], dir2[3], p1[3], p2[3];
    if (fabs(G2V_DOT(tmp,line)) < 0.7) { // the line is horizontal
      G2V_SET( uv1, 0, -line[2]/line[1] );
      G2V_SET( uv2, cp->wh[0]-1, -(line[2] + line[0]*(cp->wh[0]-1))/line[1] );
    } else {				// the line is vertical
      G2V_SET( uv1, -line[2]/line[0], 0 );
      G2V_SET( uv2, -(line[2] + line[1]*(cp->wh[1]-1))/line[0], cp->wh[1]-1 );
    }
    cp->Pixel2WorldRay( uv1, dir1, cxyz );   G3V_ADD( p1, cxyz, dir1 );
    cp->Pixel2WorldRay( uv2, dir2, cxyz );   G3V_ADD( p2, cxyz, dir2 );
    Geometry3D<double> geom;
    geom.getPlaneFromTwoLines( plane, cxyz, p1, p2 );	//// new function
  }
#endif
  
  void convertFloorL2R(T_t uvL[2], T_t uvR[2], T_t plane[4], T_t dist[1]=NULL) {
    T_t xyz[3];
    if (!camL->Pixel2PointOnPlane( uvL, xyz, plane, dist ))  G2V_COPY( uvR, uvL );
    else camR->World2Pixel( xyz, uvR );
  }
  void convertFloorR2L(T_t uvR[2], T_t uvL[2], T_t plane[4], T_t dist[1]=NULL) {
    T_t xyz[3];
    if (!camR->Pixel2PointOnPlane( uvR, xyz, plane, dist ))  G2V_COPY( uvL, uvR );
    else camL->World2Pixel( xyz, uvL );
  }
  
  // -----------------------------------------------------------------
  // Stereo Matching
  // -----------------------------------------------------------------
public:
  T_t convertXYZ2Disp(T_t xyz[3], T_t uvL[2]=NULL, T_t uvR[2]=NULL) {
    T_t uvl[2];  if (!uvL) uvL = uvl;
    T_t uvr[2];  if (!uvR) uvR = uvr;
    camL->World2Pixel( xyz, uvL );
    camR->World2Pixel( xyz, uvR );
    return (uvL[0] - uvR[0]);
  }
  T_t convertDist2Disp(T_t u, T_t v, T_t distance, T_t duv[2]=NULL) {
    if (!isReady()) return 0;
    T_t uvL[2]={u,v}, uvR[2]={0,0}, dir[3], from[3], xyz[3];
    camL->Pixel2WorldRay( uvL, dir, from );
    G3V_SCALED_ADD( xyz, from, distance, dir );
    camR->World2Pixel( xyz, uvR );
    if (duv) G2V_SUB( duv, uvL, uvR );
    return (T_t)(uvL[0] - uvR[0]);
  }
#ifdef USE_GEOMETRY_3D
  T_t convertDisp2Dist(T_t u, T_t v, T_t disp, T_t xyz[3]=NULL) {
    T_t uvL[2]={u,v}, uvR[2]={u-disp,v};
    T_t dirL[3], xyzL[3], dirR[3], xyzR[3], wxyz[3];
    camL->Pixel2WorldRay( uvL, dirL, xyzL );
    camR->Pixel2WorldRay( uvR, dirR, xyzR );
    Geometry3D<double> geom;
    int mode = -1;   // -1: OnLeftLine,  +1: OnRightLine,  0: InBetween
    if (!xyz) xyz = wxyz;
    int ret = geom.getPointByIntersectingTwoLines( xyzL, dirL, xyzR, dirR, xyz, mode );
    return (ret>0 ? G3V_DISTANCE( xyzL, xyz ) : -1);
  }
#endif
  
  // -----------------------------------------------------------------
  // File I/O
  // -----------------------------------------------------------------
public:
  bool writeStereoConfigFile(const char *filename, const char *section) {
    // Write stereo camera information in the specified section of a configuration file.
    // Note that this function will overwrite the existing file, if any.
    FILE *fp = fopen( filename, "w+" );
    if (!fp) { printf("Error (Stereo::writeStereoConfigFile): failed to open file '%s'\n", filename); return false; }
    T_t Rwc[9], Twc[3], Rcq[9], Tcq[3], Rcp[9], Tcp[3], Ang[3], AngL[3], AngR[3];
    T_t Rwr[9] = { 0, -1, 0,  0,0,-1,  1,0,0 };
    getCameraPose ( Rwc, Twc );   getPitchYawRoll( Ang,  Rwc, Rwr );
    getCameraPoseL( Rcq, Tcq );   getPitchYawRoll( AngL, Rcq );
    getCameraPoseR( Rcp, Tcp );   getPitchYawRoll( AngR, Rcp );
    fprintf( fp, "\n" );
    fprintf( fp, "[%s]\n", section );
    fprintf( fp, "StereoResolution  = %d x %d\n", wh[0], wh[1] );
    fprintf( fp, "StereoBaseline    = %g\n", baseline );
    fprintf( fp, "# Located at: Twc[3]=(%7.3f,%7.3f,%7.3f)  PYR[3]=(%4.1f %4.1f %4.1f)\n", Twc[0], Twc[1], Twc[2], Ang[0], Ang[1], Ang[2]);
    fprintf( fp, "StereoExtrinsic0  = %9.4f %9.4f %9.4f  %9.4f\n", Rcw[0], Rcw[1], Rcw[2], Tcw[0] );
    fprintf( fp, "StereoExtrinsic1  = %9.4f %9.4f %9.4f  %9.4f\n", Rcw[3], Rcw[4], Rcw[5], Tcw[1] );
    fprintf( fp, "StereoExtrinsic2  = %9.4f %9.4f %9.4f  %9.4f\n", Rcw[6], Rcw[7], Rcw[8], Tcw[2] );
    fprintf( fp, "# Left  camera at Tcq[3]=(%6.3f,%6.3f,%6.3f)  PYR[3]=(%4.1f %4.1f %4.1f)\n", Tcq[0], Tcq[1], Tcq[2], AngL[0], AngL[1], AngL[2]);
    fprintf( fp, "StereoLIntrinsic0 = %9.4f %9.4f %9.4f \n", camL->fc[0], camL->afc  , camL->cc[0] );
    fprintf( fp, "StereoLIntrinsic1 = %9.4f %9.4f %9.4f \n",         0.0, camL->fc[1], camL->cc[1] );
    fprintf( fp, "StereoLDistortion =  %g\n", camL->dt[0] );
    fprintf( fp, "StereoLTransform0 = %9.4f %9.4f %9.4f  %9.4f\n", Rqc[0], Rqc[1], Rqc[2], Tqc[0] );
    fprintf( fp, "StereoLTransform1 = %9.4f %9.4f %9.4f  %9.4f\n", Rqc[3], Rqc[4], Rqc[5], Tqc[1] );
    fprintf( fp, "StereoLTransform2 = %9.4f %9.4f %9.4f  %9.4f\n", Rqc[6], Rqc[7], Rqc[8], Tqc[2] );
    fprintf( fp, "# Right camera at Tcp[3]=(%6.3f,%6.3f,%6.3f)  PYR[3]=(%4.1f %4.1f %4.1f)\n", Tcp[0], Tcp[1], Tcp[2], AngR[0], AngR[1], AngR[2]);
    fprintf( fp, "StereoRIntrinsic0 = %9.4f %9.4f %9.4f \n", camR->fc[0], camR->afc  , camR->cc[0] );
    fprintf( fp, "StereoRIntrinsic1 = %9.4f %9.4f %9.4f \n",         0.0, camR->fc[1], camR->cc[1] );
    fprintf( fp, "StereoRDistortion =  %g\n", camR->dt[0] );
    fprintf( fp, "StereoRTransform0 = %9.4f %9.4f %9.4f  %9.4f\n", Rpc[0], Rpc[1], Rpc[2], Tpc[0] );
    fprintf( fp, "StereoRTransform1 = %9.4f %9.4f %9.4f  %9.4f\n", Rpc[3], Rpc[4], Rpc[5], Tpc[1] );
    fprintf( fp, "StereoRTransform2 = %9.4f %9.4f %9.4f  %9.4f\n", Rpc[6], Rpc[7], Rpc[8], Tpc[2] );
    fprintf( fp, "StereoMatchAXform = R[4]={%7.4f %7.4f %7.4f %7.4f}  T[2]={%5.2f %5.2f}\n", stmRT[0], stmRT[1], stmRT[2], stmRT[3], stmRT[4], stmRT[5] );
    fclose( fp );
    return true;
  }
#ifdef USE_GUIH_CONFIG
  bool readStereoConfigFile(const char *filename, const char *section, 
		      bool read_int=true, bool read_ext=true) {
    // Read stereo camera information from the specified section in a configuration file.
    char  res[80], ex0[200], ex1[200], ex2[200], stmRTs[200];
    char  cLin0[200], cLin1[200], cLtr0[200], cLtr1[200], cLtr2[200];
    char  cRin0[200], cRin1[200], cRtr0[200], cRtr1[200], cRtr2[200];
    double base, cLdst, cRdst;
    GUIH::Config cfg;
    cfg.set( section, "StereoResolution", CFG_STRING, res );
    cfg.set( section, "StereoBaseline",   CFG_DOUBLE, &base );
    cfg.set( section, "StereoExtrinsic0", CFG_STRING, ex0 );
    cfg.set( section, "StereoExtrinsic1", CFG_STRING, ex1 );
    cfg.set( section, "StereoExtrinsic2", CFG_STRING, ex2 );
    cfg.set( section, "StereoLIntrinsic0", CFG_STRING, cLin0 );
    cfg.set( section, "StereoLIntrinsic1", CFG_STRING, cLin1 );
    cfg.set( section, "StereoLDistortion", CFG_DOUBLE, &cLdst );
    cfg.set( section, "StereoLTransform0", CFG_STRING, cLtr0 );
    cfg.set( section, "StereoLTransform1", CFG_STRING, cLtr1 );
    cfg.set( section, "StereoLTransform2", CFG_STRING, cLtr2 );
    cfg.set( section, "StereoRIntrinsic0", CFG_STRING, cRin0 );
    cfg.set( section, "StereoRIntrinsic1", CFG_STRING, cRin1 );
    cfg.set( section, "StereoRDistortion", CFG_DOUBLE, &cRdst );
    cfg.set( section, "StereoRTransform0", CFG_STRING, cRtr0 );
    cfg.set( section, "StereoRTransform1", CFG_STRING, cRtr1 );
    cfg.set( section, "StereoRTransform2", CFG_STRING, cRtr2 );
    cfg.set( section, "StereoMatchAXform", CFG_STRING, stmRTs );
    if (!cfg.process( filename, section ) ||
	!cfg.processed(section, "StereoResolution") ||
	!cfg.processed(section, "StereoLTransform2") ||
	!cfg.processed(section, "StereoRTransform2")) {
      printf("Error (Camera::readStereoConfigFile): cannot read camera parameters '%s' in '%s'\n", section, filename);
      return false;
    }
    double RL[9], TL[3], RR[9], TR[3], R[9], T[3];
    sscanf( res,   "%d x %d", wh+0, wh+1 );
    baseline = base;
    sscanf( cLtr0, "%lf %lf %lf %lf", RL+0, RL+1, RL+2, TL+0 );
    sscanf( cLtr1, "%lf %lf %lf %lf", RL+3, RL+4, RL+5, TL+1 );
    sscanf( cLtr2, "%lf %lf %lf %lf", RL+6, RL+7, RL+8, TL+2 );
    sscanf( cRtr0, "%lf %lf %lf %lf", RR+0, RR+1, RR+2, TR+0 );
    sscanf( cRtr1, "%lf %lf %lf %lf", RR+3, RR+4, RR+5, TR+1 );
    sscanf( cRtr2, "%lf %lf %lf %lf", RR+6, RR+7, RR+8, TR+2 );
    setStereoTransformation( RL, TL, RR, TR );
    // read intrinsic parameters
    if (read_int) {
      double tmp, fL[2], cL[2], aL, fR[2], cR[2], aR;
      sscanf( cLin0, "%lf %lf %lf", fL+0, &aL,  cL+0 );
      sscanf( cLin1, "%lf %lf %lf", &tmp, fL+1, cL+1 );
      sscanf( cRin0, "%lf %lf %lf", fR+0, &aR,  cR+0 );
      sscanf( cRin1, "%lf %lf %lf", &tmp, fR+1, cR+1 );
      setStereoIntrinsic( wh[0], wh[1],
			  (T_t)fL[0], (T_t)fL[1], (T_t)aL, (T_t)cL[0], (T_t)cL[1], (T_t)cLdst,
			  (T_t)fR[0], (T_t)fR[1], (T_t)aR, (T_t)cR[0], (T_t)cR[1], (T_t)cRdst );
    }
    // read extrinsic parameters
    if (read_ext) {
      sscanf( ex0, "%lf %lf %lf %lf", R+0, R+1, R+2, T+0 );
      sscanf( ex1, "%lf %lf %lf %lf", R+3, R+4, R+5, T+1 );
      sscanf( ex2, "%lf %lf %lf %lf", R+6, R+7, R+8, T+2 );
      setStereoExtrinsic( R, T );
    }
    // affine transformation for rectification in stereo matching
    int read = sscanf( stmRTs, "R[4]={%lf %lf %lf %lf}  T[2]={%lf %lf}", R+0, R+1, R+2, R+3, T+0, T+1 );
    if (read==6) G6V_SET( stmRT, (T_t)R[0], (T_t)R[1], (T_t)R[2], (T_t)R[3], (T_t)T[0], (T_t)T[1] );
    else         G6V_SET( stmRT, 1, 0, 0, 1, 0, 0 );
    return true;
  }
#endif
  
  // -----------------------------------------------------------------
  // 
  // -----------------------------------------------------------------
public:
  void printInfo(const char *cmmt=NULL) {
    T_t Rwc[9], Twc[3], Rcq[9], Tcq[3], Rcp[9], Tcp[3], Ang[3], AngL[3], AngR[3];
    getCameraPose ( Rwc, Twc );   getPitchYawRoll( Ang,  Rwc );
    getCameraPoseL( Rcq, Tcq );   getPitchYawRoll( AngL, Rcq );
    getCameraPoseR( Rcp, Tcp );   getPitchYawRoll( AngR, Rcp );
    printf("Stereo Camera  %s  %s\n", (cmmt ? cmmt : ""), (isReady() ? "":"NOT_READY"));
    printf("  Located at: Twc[3]=(%7.3f,%7.3f,%7.3f)  PYR[3]=(%4.1f %4.1f %4.1f)\n", Twc[0], Twc[1], Twc[2], Ang[0], Ang[1], Ang[2]);
    printf("    [ %7.4f %7.4f %7.4f %9.4f ] [Rwc Twc]\n", Rwc[0], Rwc[1], Rwc[2], Twc[0]);
    printf("    [ %7.4f %7.4f %7.4f %9.4f ]\n", Rwc[3], Rwc[4], Rwc[5], Twc[1]);
    printf("    [ %7.4f %7.4f %7.4f %9.4f ]\n", Rwc[6], Rwc[7], Rwc[8], Twc[2]);
    printf("  Resolution: (%d x %d)   Baseline: %.3f \n", wh[0], wh[1], baseline);
    if (fabs(baseline-(Tcp[0]-Tcq[0])) > 0.001) printf("  Warning: calculated baseline is %.3f\n", (Tcp[0]-Tcq[0]));
    // Left camera ---------------------------------------------------
    printf("  Left  camera  at  Tcq[3]=(%6.3f,%6.3f,%6.3f)  PYR[3]=(%4.1f %4.1f %4.1f)\n", Tcq[0], Tcq[1], Tcq[2], AngL[0], AngL[1], AngL[2]);
    printf("    [ %7.2f %7.2f %7.2f ] [ %7.4f %7.4f %7.4f %9.4f ]  [K] [Rqw Tqw]\n", camL->fc[0], camL->afc, camL->cc[0],  camL->Rcw[0], camL->Rcw[1], camL->Rcw[2], camL->Tcw[0]);
    printf("    [ %7.2f %7.2f %7.2f ] [ %7.4f %7.4f %7.4f %9.4f ]\n", 0.0, camL->fc[1], camL->cc[1], camL->Rcw[3], camL->Rcw[4], camL->Rcw[5], camL->Tcw[1]);
    printf("    [ %7.2f %7.2f %7.2f ] [ %7.4f %7.4f %7.4f %9.4f ]\n", 0.0, 0.0, 1.0, camL->Rcw[6], camL->Rcw[7], camL->Rcw[8], camL->Tcw[2]);
    if (camL->dt_type == 'I') {	// invertible distortion model (as in Andrew Davison's vSLAM)
      printf("    Distortion parameter: [ %g ]  (rd = r / sqrt(1 + 2*k*r^2)) \n", camL->dt[0]);
    } else {		// standard distortion model (as in Matlab toolbox)
      printf("    Distortion parameters: [ %8.4f %8.4f (%8.4f %8.4f) %8.4f ]\n", camL->dt[0], camL->dt[1], camL->dt[2], camL->dt[3], camL->dt[4]);
    }
//     printf("    [ %7.4f %7.4f %7.4f %9.4f ]  [Rcq Tcq]\n", Rcq[0], Rcq[1], Rcq[2], Tcq[0]);
//     printf("    [ %7.4f %7.4f %7.4f %9.4f ]\n", Rcq[3], Rcq[4], Rcq[5], Tcq[1]);
//     printf("    [ %7.4f %7.4f %7.4f %9.4f ]\n", Rcq[6], Rcq[7], Rcq[8], Tcq[2]);
    // Right camera --------------------------------------------------
    printf("  Right camera  at  Tcp[3]=(%6.3f,%6.3f,%6.3f)  PYR[3]=(%4.1f %4.1f %4.1f)\n", Tcp[0], Tcp[1], Tcp[2], AngR[0], AngR[1], AngR[2]);
    printf("    [ %7.2f %7.2f %7.2f ] [ %7.4f %7.4f %7.4f %9.4f ]  [K] [Rpw Tpw]\n", camR->fc[0], camR->afc, camR->cc[0],  camR->Rcw[0], camR->Rcw[1], camR->Rcw[2], camR->Tcw[0]);
    printf("    [ %7.2f %7.2f %7.2f ] [ %7.4f %7.4f %7.4f %9.4f ]\n", 0.0, camR->fc[1], camR->cc[1], camR->Rcw[3], camR->Rcw[4], camR->Rcw[5], camR->Tcw[1]);
    printf("    [ %7.2f %7.2f %7.2f ] [ %7.4f %7.4f %7.4f %9.4f ]\n", 0.0, 0.0, 1.0, camR->Rcw[6], camR->Rcw[7], camR->Rcw[8], camR->Tcw[2]);
    if (camR->dt_type == 'I') {	// invertible distortion model (as in Andrew Davison's vSLAM)
      printf("    Distortion parameter: [ %g ]  (rd = r / sqrt(1 + 2*k*r^2)) \n", camR->dt[0]);
    } else {		// standard distortion model (as in Matlab toolbox)
      printf("    Distortion parameters: [ %8.4f %8.4f (%8.4f %8.4f) %8.4f ]\n", camR->dt[0], camR->dt[1], camR->dt[2], camR->dt[3], camR->dt[4]);
    }
//     printf("    [ %7.4f %7.4f %7.4f %9.4f ]  [Rcp Tcp]\n", Rcp[0], Rcp[1], Rcp[2], Tcp[0]);
//     printf("    [ %7.4f %7.4f %7.4f %9.4f ]\n", Rcp[3], Rcp[4], Rcp[5], Tcp[1]);
//     printf("    [ %7.4f %7.4f %7.4f %9.4f ]\n", Rcp[6], Rcp[7], Rcp[8], Tcp[2]);
    printf("  StereoMatching Xform: R=[ %5.2f %5.2f %5.2f %5.2f ] T=[ %4.1f %4.1f ]\n", stmRT[0], stmRT[1], stmRT[2], stmRT[3], stmRT[4], stmRT[5]);
  }
};


}	// end of CAMH namespace

#endif // CAMH_STEREO_HPP
