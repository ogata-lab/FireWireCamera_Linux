
//
// MTH::Rotation<> class template
//
// Jaeil Choi
// Last modified in Dec, 2005
//

#ifndef MTH_ROTATION_HPP
#define MTH_ROTATION_HPP
#define USE_MTH_ROTATION

#include <iostream>
#include <cmath>
#include "vm_macros.h"

namespace MTH {
  
template <class T>
class Rotation {
 public:
  char  info[256];
 public:
  Rotation()  { }
  ~Rotation() { }
  
public:
  void setRotation(T R[9], T x, T y, T z, T angle) {
    T r[3], q[4];
    G3V_SET( r, x, y, z );  G3V_NORMALIZE( r );
    Rv2R( r, angle, R );
  }
  
  // ================================================================
  // iterative fitting to rotation
  // ================================================================
public:  
  void MakeRotation(T R[9]) {
    // Converge R to the nearest orthogonal matrix
    //   Let error F(R) = Rt R - I, then
    //   R(i+1) = R(i) - (F')^{-1} F
    //          = R(i) - (2 Rt)^{-1} (Rt R - I)
    //          = R(i) - 0.5 Rt^{-1} (Rt R - I)
    //         ~= R(i) - 0.5 R (Rt R - I)
    T  Rt[9], RtR[9], Tmp[9], error;
    for (int i = 0; i < 10; i++) {
      G3M_TRANS( Rt, R );
      G3M_MUL_MM( RtR, Rt, R );
      RtR[0] -= 1;  RtR[4] -= 1;  RtR[8] -= 1;
      error = G3M_NORM( RtR );
      G3M_MUL_MM( Tmp, R, RtR );
      for (int j = 0; j < 9; j++) R[j] = R[j] - 0.5 * Tmp[j];
      //cout << "MakeRotation (" << i << ") error = " << error << endl;
      if (error < 1.0e-6) break;
    }
  }
  
  // ================================================================
  // Rotation Angles
  // ================================================================
public:  
  void ZYX2R(T x, T y, T z, T R[9], bool degree=true) {
    // R = M(Z) * M(Y) * M(X)
    if (degree) { x = x*M_PI/180;  x = x*M_PI/180;  x = x*M_PI/180; }
    T cx=cos(x), sx=sin(x), cy=cos(y), sy=sin(y), cz=cos(z), sz=sin(z);
    G3M_SET( R,
	     cy*cz, cx*sz+sx*sy*cz,  -cx*sy*cz+sx*sz,
	     -cy*sz, cx*cz-sx*sy*sz, sx*cz+cx*sy*sz,
	     sy, -sx*cy, cx*cy );
  }
  void R2ZYX(T R[9], T *x, T *y, T *z) {
    *y = asin( R[6] );
    *x = asin( -R[7] / cos(*y) );
    *z = asin( -R[3] / cos(*y) );
  }
  
  void XYZ2R(T x, T y, T z, T R[9], bool degree=true) {
    // R = M(X) * M(Y) * M(Z)
    if (degree) { x = x*M_PI/180;  x = x*M_PI/180;  x = x*M_PI/180; }
    T cx=cos(x), sx=sin(x), cy=cos(y), sy=sin(y), cz=cos(z), sz=sin(z);
    G3M_SET( R,
	     cy*cz,  cy*sz,  -sy,
	     sx*sy*cz-cx*sz,  sx*sy*sz+cx*cz,  sx*cy,
	     cx*sy*cz+sx*sz,  cx*sy*sz-sx*cz,  cx*cy );
  }
  void R2XYZ(T R[9], T *x, T *y, T *z) {
    *y = asin( -R[2] );
    *x = asin( R[5] / cos(*y) );
    *z = asin( R[1] / cos(*y) );
  }
  
  // ================================================================
  // Rotation vector and Rodrigues rotation formula
  // ================================================================
public:  
  // O. Faugeras. "Three-dimensional Computer Vision: A Geometric Viewpoint", MIT 1993
  
  void R2Rv(T R[9], T r[3], T *theta=NULL) {
    // convert rotation matrix to rotation axis and angle (radian)
    T q[4], angle;
    R2Q( R, q );           // rotation matrix -> normalized quaternion
    Q2Rv( q, r, &angle, false );  // quaternion -> rotation axis and angle
    if (theta) *theta = angle;
    else G3V_MUL_VALUE( r, angle );
  }
  
  void Rv2R(T r[3], T R[9]) {
    // convert rotation axis (|r|=angle) to rotation matrix
    T    len,  ra[3];
    if ((len = G3V_LENGTH(r)) == 0) {
      G3M_SET( R,  1, 0, 0,  0, 1, 0,  0, 0, 1 );
    } else {
      G3V_SET( ra, r[0]/len, r[1]/len, r[2]/len );
      Rv2R( ra, len, R, false );
    }
  }
  
  void Rv2R(T r[3], T theta, T R[9], bool normalize=true) {
    // convert rotation axis and angle to rotation matrix
    if (normalize) G3V_NORMALIZE_WITH_TYPE( r, T );
    T  c, s, u;
    c = cos(theta);  s = sin(theta);  u = 1 - c;
    G3M_SET( R,
	     r[0]*r[0]*u + c,      r[0]*r[1]*u - r[2]*s,  r[0]*r[2]*u + r[1]*s,
	     r[0]*r[1]*u + r[2]*s, r[1]*r[1]*u + c,       r[1]*r[2]*u - r[0]*s,
	     r[0]*r[2]*u - r[1]*s, r[1]*r[2]*u + r[0]*s,  r[2]*r[2]*u + c );
  }
  
  // ================================================================
  // Quaternion   q[] = (w, x, y, z)
  // ================================================================
public:  
  inline T QNorm(T q[4]) {
    return sqrt( q[0]*q[0] + q[1]*q[1] + q[2]*q[2] + q[3]*q[3] );
  }
  void QNormalize(T q[4]) {
    T len = QNorm( q );
    if (len < 0.000001) return;
    q[0] /= len;  q[1] /= len;  q[2] /= len;  q[3] /= len;  
  }
  
  void R2Q(T R[9], T q[4]) {
    // convert rotation matrix to normalized quaternion
    T  s, tr = 1 + R[0] + R[4] + R[8];
    if (tr > 0.00000001) {
      s = sqrt(tr) * 2;
      q[0] = 0.25 * s;
      q[1] = ( R[7] - R[5] ) / s;  
      q[2] = ( R[2] - R[6] ) / s;
      q[3] = ( R[3] - R[1] ) / s;
    } else {
      if (R[0] > R[4] && R[0] > R[8]) {
	s = sqrt( 1.0 + R[0] - R[4] - R[8] ) * 2;
	q[0] = ( R[7] - R[5] ) / s;
	q[1] = 0.25 * s;
	q[2] = ( R[3] + R[1] ) / s;
	q[3] = ( R[2] + R[6] ) / s;
      } else if (R[4] > R[8]) {
	s = sqrt( 1.0 + R[4] - R[0] - R[8] ) * 2;
	q[0] = ( R[2] - R[6] ) / s;
	q[1] = ( R[3] + R[1] ) / s;
	q[2] = 0.25 * s;
	q[3] = ( R[7] + R[5] ) / s;
      } else {
	s = sqrt( 1.0 + R[8] - R[0] - R[4] ) * 2;
	q[0] = ( R[3] - R[1] ) / s;
	q[1] = ( R[2] + R[6] ) / s;
	q[2] = ( R[7] + R[5] ) / s;
	q[3] = 0.25 * s;
      }
    }
    QNormalize(q);
  }
  void Q2R(T q[4], T R[9], bool normalize=true) {
    // convert quaternion to rotation matrix
    if (normalize) QNormalize( q );
    T xx = q[1]*q[1];  T xy = q[1]*q[2];  T xz = q[1]*q[3];
    T xw = q[1]*q[0];  T yy = q[2]*q[2];  T yz = q[2]*q[3];
    T yw = q[2]*q[0];  T zz = q[3]*q[3];  T zw = q[3]*q[0];
    R[0] = 1 - 2 * ( yy + zz );
    R[1] =     2 * ( xy - zw );
    R[2] =     2 * ( xz + yw );
    R[3] =     2 * ( xy + zw );
    R[4] = 1 - 2 * ( xx + zz );
    R[5] =     2 * ( yz - xw );
    R[6] =     2 * ( xz - yw );
    R[7] =     2 * ( yz + xw );
    R[8] = 1 - 2 * ( xx + yy );
  }
  void Rv2Q(T axis[3], T angle, T q[4], bool normalize=true) {
    // convert rotation axis and angle to quaternion
    //   q = [ cos(angle/2)            ]
    //       [ sin(angle/2) * x/norm ]
    //       [ sin(angle/2) * y/norm ]
    //       [ sin(angle/2) * z/norm ]
    if (normalize) G3V_NORMALIZE( axis );
    T sinv = sin( angle / 2 );
    T cosv = cos( angle / 2 );
    q[0] = cosv;
    q[1] = sinv * axis[0];
    q[2] = sinv * axis[1];
    q[3] = sinv * axis[2];
  }
  void Wt2Q(T w[3], T dt, T q[4]) {
    // convert angular velocity and delte_t to quaternion
    //   q(w*t) = [ cos(angle/2)            ]
    //            [ sin(angle/2) * wx/wnorm ]
    //            [ sin(angle/2) * wy/wnorm ]
    //            [ sin(angle/2) * wz/wnorm ], where angle = sqrt(wxwx+wywy+wzwz)*dt
    // Note this is just a linear approximation, assuming ||w|| and dt are small.
    T norm = sqrt( w[0]*w[0] + w[1]*w[1] + w[2]*w[2] );
    if (norm > 1.0e-8) {
      T angle = norm * dt;
      T sinv  = sin( angle / 2.0 ) / norm;
      T cosv  = cos( angle / 2.0 );
      q[0] = cosv;
      q[1] = sinv * w[0];
      q[2] = sinv * w[1];
      q[3] = sinv * w[2];
    } else {
      q[0] = 1.0;  q[1] = q[2] = q[3] = 0.0;
    }
  }
  void Q2Rv(T q[4], T axis[3], T *angle, bool normalize=true) {
    // convert quaternion to rotation axis and angle(radian)
    if (normalize) QNormalize( q );
    T cos = q[0];
    T sin = sqrt( 1 - cos * cos );
    if ( fabs(sin) < 0.0005 ) sin = 1;
    *angle = acos( cos ) * 2;
    axis[0] = q[1] / sin;
    axis[1] = q[2] / sin;
    axis[2] = q[3] / sin;
  }
  void QMult(T qr[4], T qa[4], T qb[4]) {
    // multiply two quaternions 'qa' and 'qb'
    // w = w1w2 - x1x2 - y1y2 - z1z2 = w1w2 - ra.rb
    // x = w1x2 + x1w2 + y1z2 - z1y2
    // y = w1y2 + y1w2 + z1x2 - x1z2
    // z = w1z2 + z1w2 + x1y2 - y1x2
    qr[0] = qa[0]*qb[0] - qa[1]*qb[1] - qa[2]*qb[2] - qa[3]*qb[3];
    qr[1] = qa[0]*qb[1] + qa[1]*qb[0] + qa[2]*qb[3] - qa[3]*qb[2];
    qr[2] = qa[0]*qb[2] + qa[2]*qb[0] + qa[3]*qb[1] - qa[1]*qb[3];
    qr[3] = qa[0]*qb[3] + qa[3]*qb[0] + qa[1]*qb[2] - qa[2]*qb[1];
  }
  void QInverse(T qr[4], T q[4]) {
    // inverse of a quaternion 'q'
    // qr = conj(qa) / ||q||
    //    = [w -x -y -z] / ||q||
    T len = QNorm(q);
    qr[0] = +q[0] / len;
    qr[1] = -q[1] / len;
    qr[2] = -q[2] / len;
    qr[3] = -q[3] / len;
  }
  void QRotateVector(T r[3], T q[4], T v[3]) {
    // rotate a vector 'v' using a quaternion 'q'
    // [1 v'] = q * [1 v] * q^{-1}
    //        = [w x y z] * [1 v] * [w -x -y -z]
    T qw, qx, qy, qz, q_v, tmx, tmy, tmz;
    qw = q[0];  qx = q[1];  qy = q[2];  qz = q[3];
    q_v  = qx * v[0] + qy * v[1] + qz * v[2];
    tmx  = qw * v[0] + qy * v[2] - qz * v[1];
    tmy  = qw * v[1] + qz * v[0] - qx * v[2];
    tmz  = qw * v[2] + qx * v[1] - qy * v[0];
    // Note w' = w*w + x*x + y*y + z*z = ||q|| = 1  for rotation
    r[0] = (q_v * qx + tmx * qw - tmy * qz + tmz * qy);
    r[1] = (q_v * qy + tmy * qw - tmz * qx + tmx * qz);
    r[2] = (q_v * qz + tmz * qw - tmx * qy + tmy * qx);
  }
  void QJacobian_dq3_dq1(T q2[4], T dq3_dq1[4*4]) {
    // calculate the Jacobian dq3/dq1 (4x4 matrix), where q3 = q2 x q1
    G4M_SET( dq3_dq1,
	     +q2[0], -q2[1], -q2[2], -q2[3],
	     +q2[1], +q2[0], -q2[3], +q2[2],
	     +q2[2], +q2[3], +q2[0], -q2[1],
	     +q2[3], -q2[2], +q2[1], +q2[0] );
  }
  void QJacobian_dq3_dq2(T q1[4], T dq3_dq2[4*4]) {
    // calculate the Jacobian dq3/dq2 (4x4 matrix), where q3 = q2 x q1
    G4M_SET( dq3_dq2,
	     +q1[0], -q1[1], -q1[2], -q1[3],
	     +q1[1], +q1[0], +q1[3], -q1[2],
	     +q1[2], -q1[3], +q1[0], +q1[1],
	     +q1[3], +q1[2], -q1[1], +q1[0] );
  }
  void QJacobian_dqwt_dw(T w[3], T dt, T dqwt_dw[4*3]) {
    // calculate the Jacobian dq(w*dt)/dw (4x3 matrix), where w is angular velocity
    //   qwt = [ cos(angle/2)            ]
    //         [ sin(angle/2) * wx/wnorm ]
    //         [ sin(angle/2) * wy/wnorm ]
    //         [ sin(angle/2) * wz/wnorm ], where angle = sqrt(wxwx+wywy+wzwz)*dt
    T wnorm = sqrt(w[0] * w[0] + w[1] * w[1] + w[2] * w[2]);
    if (wnorm < 1.0e-8) { memset(dqwt_dw, 0, 4*3*sizeof(T)); return; }
    T angle = wnorm * dt;
    T sinv  = sin( angle / 2 );
    T cosv  = cos( angle / 2 );
    T comm1 = -0.5 * sinv * dt / wnorm;
    T comm2 = cosv * (dt/2) / (wnorm*wnorm) - sinv / (wnorm*wnorm*wnorm);
    T added = sinv / wnorm;
    dqwt_dw[0] = comm1 * w[0];			// 0th row
    dqwt_dw[1] = comm1 * w[1];
    dqwt_dw[2] = comm1 * w[2];
    dqwt_dw[3] = comm2 * w[0] * w[0] + added;	// 1st row
    dqwt_dw[4] = comm2 * w[0] * w[1];
    dqwt_dw[5] = comm2 * w[0] * w[2];
    dqwt_dw[6] = comm2 * w[1] * w[0];		// 2nd row
    dqwt_dw[7] = comm2 * w[1] * w[1] + added;
    dqwt_dw[8] = comm2 * w[1] * w[2];
    dqwt_dw[9] = comm2 * w[2] * w[0];		// 3rd row
    dqwt_dw[10]= comm2 * w[2] * w[1];
    dqwt_dw[11]= comm2 * w[2] * w[2] + added;
  }
  void QJacobian_dqinv_dq(T dqinv_dq[4*4]) {
    // calculate the Jacobian		// FROM : dqbar_by_dq()
    G4M_SET( dqinv_dq,
	     1.0,  0.0,  0.0,  0.0,
	     0.0, -1.0,  0.0,  0.0,
	     0.0,  0.0, -1.0,  0.0,
	     0.0,  0.0,  0.0, -1.0 );
  }
  void QJacobian_dRv_dq(T q[4], T v[3], T dqa_dq[3*4]) {
    // calculate the Jacobian d(q*v)/dq  (3x4 matrix)  // FROM : dRq_times_a_by_dq()
    T M[9];
    G3M_SET( M,	// dR_dq0,
	     +2*q[0], -2*q[3], +2*q[2],
	     +2*q[3], +2*q[0], -2*q[1],
	     -2*q[2], +2*q[1], +2*q[0] );
    // the 0th column
    dqa_dq[0*4+0] = M[0*3+0] * v[0] + M[0*3+1] * v[1] + M[0*3+2] * v[2];
    dqa_dq[1*4+0] = M[1*3+0] * v[0] + M[1*3+1] * v[1] + M[1*3+2] * v[2];
    dqa_dq[2*4+0] = M[2*3+0] * v[0] + M[2*3+1] * v[1] + M[2*3+2] * v[2];
    G3M_SET( M,	// dR_dqx,
	     +2*q[1], +2*q[2], +2*q[3],
	     +2*q[2], -2*q[1], -2*q[0],
	     +2*q[3], +2*q[0], -2*q[1] );
    // the 1st column
    dqa_dq[0*4+1] = M[0*3+0] * v[0] + M[0*3+1] * v[1] + M[0*3+2] * v[2];
    dqa_dq[1*4+1] = M[1*3+0] * v[0] + M[1*3+1] * v[1] + M[1*3+2] * v[2];
    dqa_dq[2*4+1] = M[2*3+0] * v[0] + M[2*3+1] * v[1] + M[2*3+2] * v[2];
    G3M_SET( M,	// dR_dqy,
	     -2*q[2], +2*q[1], +2*q[0],
	     +2*q[1], +2*q[2], +2*q[3],
	     -2*q[0], +2*q[3], -2*q[2] );
    // the 2nd column
    dqa_dq[0*4+2] = M[0*3+0] * v[0] + M[0*3+1] * v[1] + M[0*3+2] * v[2];
    dqa_dq[1*4+2] = M[1*3+0] * v[0] + M[1*3+1] * v[1] + M[1*3+2] * v[2];
    dqa_dq[2*4+2] = M[2*3+0] * v[0] + M[2*3+1] * v[1] + M[2*3+2] * v[2];
    G3M_SET( M,	// dR_dqz,
	     -2*q[3], -2*q[0], +2*q[1],
	     +2*q[0], -2*q[3], +2*q[2],
	     +2*q[1], +2*q[2], +2*q[3] );
    // the 3rd column
    dqa_dq[0*4+3] = M[0*3+0] * v[0] + M[0*3+1] * v[1] + M[0*3+2] * v[2];
    dqa_dq[1*4+3] = M[1*3+0] * v[0] + M[1*3+1] * v[1] + M[1*3+2] * v[2];
    dqa_dq[2*4+3] = M[2*3+0] * v[0] + M[2*3+1] * v[1] + M[2*3+2] * v[2];
  }
  void QJacobian_dqnorm_dq(T q[4], T dqnorm_dq[16]) {
    T qq = q[0] * q[0] + q[1] * q[1] + q[2] * q[2] + q[3] * q[3];
    T qq2 = qq * qq,  qq3 = qq * qq * qq;
    G4M_SET ( dqnorm_dq,
	      +(1-q[0]*q[0]/qq2)/qq, -(q[0] * q[1]) / qq3,  -(q[0] * q[2]) / qq3,  -(q[0] * q[3]) / qq3,
	      -(q[1] * q[0]) / qq3,  +(1-q[1]*q[1]/qq2)/qq, -(q[1] * q[2]) / qq3,  -(q[1] * q[3]) / qq3,
	      -(q[2] * q[0]) / qq3,  -(q[2] * q[1]) / qq3,  +(1-q[2]*q[2]/qq2)/qq, -(q[2] * q[3]) / qq3,
	      -(q[3] * q[0]) / qq3,  -(q[3] * q[1]) / qq3,  -(q[3] * q[2]) / qq3,  +(1-q[3]*q[3]/qq2)/qq );
  }
  
  // =================================================================
  // 
  // =================================================================
public:
  char* getInfoR(T R[9]) {
    T Rv[3], angle;
    R2Rv( R, Rv, &angle );
    sprintf(info, "angle=%.1f(d) on axis=(%.2f %.2f %.2f)", angle*180/M_PI, Rv[0], Rv[1], Rv[2]);
    return info;
  }
  char* getInfoQ(T q[4]) {
    T Rv[3], angle;
    Q2Rv( q, Rv, &angle );
    sprintf(info, "angle=%.1f(d) on axis=(%.2f %.2f %.2f)", angle*180/M_PI, Rv[0], Rv[1], Rv[2]);
    return info;
  }
};


}	// end of namespace MTH

#endif // CROTATION_HPP


//
// test code
//
// #define RA_COUT(r, a)  "(" << setprecision(4) << r[0] << " " << setprecision(4) << r[1] << " " << setprecision(4) << r[2] << " ; " << angle*180/M_PI << ")"
// #define Q_COUT(q)  "(" << setprecision(4) << q[0] << " " << setprecision(4) << q[1] << " " << setprecision(4) << q[2] << " " << setprecision(4) << q[3] << ")"
// void test_code(void)
// {
//   MTH::Rotation<float> rot;
//   float r[3], angle, R[9], q[4], v[3];
//   G3V_SET( r, -1.0, 0.2, 0.4 );  G3V_NORMALIZE(r);  angle = 18 * M_PI/180;
//   cout << "Rv " << RA_COUT( r, angle ) << endl;
//   rot.Rv2R( r, angle, R );
//   rot.R2Rv( R, r, &angle );
//   cout << "Rv " << RA_COUT( r, angle ) << endl;
  
//   cout << endl;
//   rot.setRotation(R, 1, 0, 0, 90 * M_PI/180);
//   cout << "R  " << G3M_COUT(R) << endl;
//   rot.R2Q( R, q );
//   cout << "Qt " << Q_COUT( q ) << endl;
//   rot.Q2R( q, R );
//   cout << "R  " << G3M_COUT(R) << endl;
//   rot.R2Q( R, q );
//   cout << "Qt " << Q_COUT( q ) << endl;
//   rot.R2Rv( R, r, &angle );
//   cout << "Rv " << RA_COUT( r, angle ) << endl;
// }

