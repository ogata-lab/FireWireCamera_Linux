
//
// Geometry3D<T> class template
//
// Jaeil Choi
// last modified in June, 2003
//
// Note that functions of template classes are included 
// in the compiled codes only if they are called.
//

#ifndef CGEOMETRY3D_HPP
#define CGEOMETRY3D_HPP
#ifndef USE_GEOMETRY_3D
#define USE_GEOMETRY_3D
#endif

#include <iostream>
#include <iomanip>
#include <cfloat>
#include <cmath>
#include "vm_macros.h"
using namespace std;


template <class T>
class Geometry3D {
  
 public:
  
  // ================================================================ 
  // Line / Line Segment
  //
  // T     getRandomSampleOnLineSegment(T v0[3], T v1[3], T p[3]);
  // T     getDistanceLinesegPoint( T v0[3], T v1[3], T p[3], T *projection = NULL);
  // void  getGradientLineLength(T v0[3], T v1[3], T result[3]);
  // ================================================================ 
  
  int getPointByIntersectingTwoLines(T ap[3], T adir[3], T bp[3], T bdir[3], T xp[3], int mode=0) {
    // Find the intersection between two lines, 'ap[]+s*adir[]' and 'bp[]+s*bdir[]'.
    //   Note that 'adir[]' and 'bdir[]' are assumed to be normalized already.
    T nx[3], ny[3], nz[3];
    // find the plane of 'ap[3]', 'bp[3]', and 'xyz[3]', and its rotation matrix
    G3V_SUB  ( nx, bp, ap );  G3V_NORMALIZE( nx );
    if      (mode < 0)  G3V_CROSS( nz, nx, adir );	// intersection on line 'a'
    else if (mode > 0)  G3V_CROSS( nz, nx, bdir );	// intersection on line 'b'
    else {						// intersection in between
      G3V_AVERAGE2( ny, adir, bdir );  G3V_NORMALIZE( ny );
      G3V_CROSS( nz, nx, ny );
    }
    G3V_NORMALIZE( nz );
    G3V_CROSS( ny, nz, nx );
    T Rno[9], Tno[3], Ron[9], Ton[3];
    G3M_SET( Ron, nx[0], ny[0], nz[0], nx[1], ny[1], nz[1], nx[2], ny[2], nz[2] );
    G3V_SET( Ton, ap[0], ap[1], ap[2] );
    G3M_XFORM_INVERSE( Rno, Tno, Ron, Ton );
    // convert the two 3D lines to the new coordinate system of the plane
    T apn[3], adirn[3], bpn[3], bdirn[3], lineA[3], lineB[3];
    G3M_XFORM_MVV( apn, Rno, ap, Tno );  G3M_MUL_MV( adirn, Rno, adir );
    G3M_XFORM_MVV( bpn, Rno, bp, Tno );  G3M_MUL_MV( bdirn, Rno, bdir );
    // get 2D lines on the new plane by removing Z components
    G2V_LINE_EQ_NP( lineA, -adirn[1], +adirn[0], apn[0], apn[1] );
    G2V_LINE_EQ_NP( lineB, -bdirn[1], +bdirn[0], bpn[0], bpn[1] );
    // intersect the two lines
    T xyz[3] = {0,0,0};
    T det = ( lineA[0] * lineB[1] - lineA[1] * lineB[0] );
    if (fabs(det) < 0.0001) return 0;	//  0 : DO NOT INTERSECT
    xyz[0] = - (+ lineB[1] * lineA[2] - lineA[1] * lineB[2]) / det;
    xyz[1] = - (- lineB[0] * lineA[2] + lineA[0] * lineB[2]) / det;
    // get the intersection point back to 3D using the transformation 'Ron[9]' and 'Ton[3]'
    G3M_XFORM_MVV( xp, Ron, xyz, Ton );	// +1 : INTERSECTS IN THE FRONT
    return (xyz[1] >= 0 ? +1 : -1);	// -1 : INTERSECTS IN THE BACK
  }
  
  T getRandomSampleOnLineSegment(T v0[3], T v1[3], T p[3]) {
    // get a uniformly random sample p[] on the line segment,
    T w = (T) random() / (T) RAND_MAX;
    G3V_WEIGHTED_ADD( p, (w), v0, (1-w), v1 );
    // return the weight of v0 for the sample.
    return w;
  }
  
  T getDistanceLinePoint( T v[3], T dir[3], T p[3] ) {
    // Calculate the distance from a line 'v[]+s*dir[]' to a point 'p[3]'.
    T vp[3], pp[3];
    G3V_SUB( vp, p, v );		// vector to 'p', from 'v'
    T s = G3V_DOT( vp, dir );		// project  'vp'  onto 'dir'
    G3V_SCALED_ADD( pp, vp, -s, dir );	// vector to 'p', perpendicular to 'dir'
    return G3V_LENGTH( pp );
  }
  
  T getDistanceLinesegPoint( T v0[3], T v1[3], T p[3], T *projection = NULL) {
    // calculate distance from vertex p to edge v0v1
    T v0_p[3], v0_v1[3], dist2_v0_p, line_length, projected_length, dist;
    // to avoid numerical error, use lower address as reference point
    if ((unsigned int)v0 > (unsigned int)v1) { T *temp = v0;  v0 = v1;  v1 = temp; }
    //if (debug) printf("  distance of (%.8f %.8f %.8f) from edge between (%.2f %.2f %.2f) and (%.2f %.2f %.2f)\n", p[0], p[1], p[2], v1[0], v1[1], v1[2], v2[0], v2[1], v2[2]);
    // calculate vector from v0 to p, and its squared distance
    G3V_SUB( v0_p, p, v0 );
    dist2_v0_p = (T) G3V_LENGTH2( v0_p );
    // calculate normalized vector from v0 to v1
    G3V_SUB( v0_v1, v1, v0 );
    line_length = (T) G3V_LENGTH( v0_v1 );
    if (line_length > 0) G3V_DIV_VALUE( v0_v1, line_length );
    // projected lenth of vector v0_p onto v0_v1 (dot product)
    projected_length = G3V_DOT( v0_v1, v0_p );
    // calculate final distance from vertex p to line v0_v1
    if (projected_length > line_length) {
      dist = G3V_DISTANCE( p, v1 );
      if (projection) G3V_COPY( projection, v1 );
    } else if (projected_length < 0) {
      dist = G3V_DISTANCE( p, v0 );
      if (projection) G3V_COPY( projection, v0 );
    } else {
      dist = (float)sqrt(dist2_v0_p - projected_length * projected_length);
      if (projection) {
	G3V_MUL_VALUE( v0_v1, projected_length );
	G3V_ADD( projection, v0, v0_v1 );
      }
    }
    //if (debug) printf("  projected dist = %.8f   distance from v1 = %.8f  =>  dist = %.8f\n", projected_length, dist2_v0_p, dist);
    return dist;
  }
  
  void getGradientLineLength(T v0[3], T v1[3], T result[3]) {
    // L(v0)  = || v0 - v1 ||
    G3V_SUB( result, v0, v1 );
    G3V_NORMALIZE( result );
  }
  
  // ================================================================ 
  // Plane
  //
  // bool getPlane(T plane[4], T p0[3], T p1[3], T p2[3]);
  // int  intersectPlaneAndLineSegment( T plane[4], T p0[3], T p1[3], T result[3] = NULL);
  // int  intersectPlaneAndLineDir( T plane[4], T from[3], T dir[3], T result[3] = NULL);
  // ================================================================ 
  
  bool getPlane(T plane[4], T p0[3], T p1[3], T p2[3]) {
    // calculate the equation of the plane containing given three points
    // assuming the points are CCW order
    T v1[3], v2[3];
    G3V_SUB( v1, p1, p0 );  G3V_NORMALIZE( v1 );
    G3V_SUB( v2, p2, p0 );  G3V_NORMALIZE( v2 );
    G3V_CROSS( plane, v1, v2 );  G3V_NORMALIZE( plane );
    plane[3] = - G3V_DOT( plane, p0 );
    return true;
  }
  bool getPlaneFromTwoLines(T plane[4], T p[3], T dir1[3], T dir2[3]) {
    T p1[3];  G3V_ADD( p1, p, dir1 );
    T p2[3];  G3V_ADD( p2, p, dir2 );
    return getPlane( p, p1, p2 );
  }
  T DistanceFromPointToPlane(T point[3], T p0[3], T p1[3], T p2[3]) {
    T plane[4];
    getPlane( plane, p0, p1, p2 );
    return ( G3V_DOT( plane, point ) + plane[3] );
  }
  
  int intersectPlaneAndLineSegment(T plane[4], T p0[3], T p1[3], T result[3] = NULL) {
    int intersect;
    float dir[3];    // the direction of the line
    G3V_SUB( dir, p1, p0 );
    
    // check whether they intersect or not
    T p0side = G3V_PLANE_POINT( plane, p0 );
    T p1side = G3V_PLANE_POINT( plane, p1 );
    if      (p0side == 0 && p1side == 0) intersect = 0;   // does not intersect
    else if (p0side <= 0 && p1side >= 0) intersect = +1;  // from front to back
    else if (p0side >= 0 && p1side <= 0) intersect = -1;  // from back to front
    else intersect = 0;
    //printf(" p0 (%.2f %.2f %.2f : %.2f)  p1 (%.2f %.2f %.2f : %.2f)  plane (%.2f %.2f %.2f %.2f)  \n", p0[0], p0[1], p0[2], p0side, p1[0], p1[1], p1[2], p1side, plane[0], plane[1], plane[2], plane[3]);
    // find the coordinates of the intersection point
    if (intersect != 0 && result != NULL) {
      // N . (p0 + t * dir) + plane[3] = 0  , where N = (plane[0] plane[1] plane[2])
      // t = - ( N . p0 + plane[3] ) / ( N . dir )
      // result = p0 + t * dir
      T t = -( G3V_PLANE_POINT(plane, p0) / G3V_DOT(plane, dir) );
      G3V_SCALED_ADD( result, p0, t, dir );  
    }
    return intersect;
  }
  
  int intersectPlaneAndLineDir( T plane[4], T from[3], T dir[3], T result[3] = NULL) {
    int intersect;
    // check whether they intersect or not
    if      ( G3V_DOT( dir, plane ) < -1.0e-8 ) intersect = -1;  // from front to back
    else if ( G3V_DOT( dir, plane ) > +1.0e-8 ) intersect = +1;  // from back to front
    else intersect = 0;                                     // does not intersect
    // find the coordinates of the intersection point
    if (intersect != 0 && result != NULL) {
      // N . (from + t * dir) + plane[3] = 0  , where N = (plane[0] plane[1] plane[2])
      // t = - ( N . from + plane[3] ) / ( N . dir )
      // result = from + t * dir
      T t = -( G3V_PLANE_POINT(plane, from) / G3V_DOT(plane, dir) );
      G3V_SCALED_ADD( result, from, t, dir );  
    }
    return intersect;
  }
  
  bool getIntersectionOf3Planes(T point[3], T p0[4], T p1[4], T p2[4]) {
    // ax + by + cz + j = 0
    // dx + ey + fz + k = 0
    // gx + hy + iz + l = 0
    T det = ( + p0[0]*(p1[1]*p2[2]-p1[2]*p2[1]) 	// det = a(ei-fh) - b(di-fg) + c(dh-eg)
	      - p0[1]*(p1[0]*p2[2]-p1[2]*p2[0])
	      + p0[2]*(p1[0]*p2[1]-p1[1]*p2[0]) );
    if (det == 0) return false;
    T M[3][3] = { {+(p1[1]*p2[2]-p1[2]*p2[1]), -(p0[1]*p2[2]-p0[2]*p2[1]), +(p0[1]*p1[2]-p0[2]*p1[1])}, 
		  {-(p1[0]*p2[2]-p1[2]*p2[0]), +(p0[0]*p2[2]-p0[2]*p2[0]), -(p0[0]*p1[2]-p0[2]*p1[0])}, 
		  {+(p1[0]*p2[1]-p1[1]*p2[0]), -(p0[0]*p2[1]-p0[1]*p2[0]), +(p0[0]*p1[1]-p0[1]*p1[0])} };
    point[0] = - (M[0][0]*p0[3] + M[0][1]*p1[3] + M[0][2]*p2[3]) / det;
    point[1] = - (M[1][0]*p0[3] + M[1][1]*p1[3] + M[1][2]*p2[3]) / det;
    point[2] = - (M[2][0]*p0[3] + M[2][1]*p1[3] + M[2][2]*p2[3]) / det;
    return true;
  }
  
  bool getLineByIntersectingTwoPlanes(T pa[4], T pb[4], T xyz[3], T dir[3]) {
    // calculate the normal of the intersecting line using cross product
    if (fabs(G3V_DOT(pa, pb)) > 0.999) return false;
    G3V_CROSS( dir, pa, pb );  G3V_NORMALIZE( dir );
    // find a point that lies on both 'pa' and 'pb'
    T  A[4], x[2], b[2], Ai[4], det;
    if (fabs(dir[0]) > fabs(dir[1])) {	// search it on YZ plane (when X=0)
      G2M_SET( A, pa[1], pa[2], pb[1], pb[2] );
      G2V_SET( b, -pa[3], -pb[3] );	// A * x = b
      det = ( A[0] * A[3] - A[1] * A[2] );
      G2M_SET( Ai, +A[3]/det, -A[1]/det, -A[2]/det, +A[0]/det );
      G2M_MUL_MV( x, Ai, b );		// x = Ainv * b
      G3V_SET( xyz, 0, x[0], x[1] );	// the point on the intersecting line
    } else {				// search it on XZ plane (when Y=0)
      // solve A*x = b   =>   x = Ainv * b
      G2M_SET( A, pa[0], pa[2], pb[0], pb[2] );
      G2V_SET( b, -pa[3], -pb[3] );	// A * x = b
      det = ( A[0] * A[3] - A[1] * A[2] );
      G2M_SET( Ai, +A[3]/det, -A[1]/det, -A[2]/det, +A[0]/det );
      G2M_MUL_MV( x, Ai, b );		// x = Ainv * b
      G3V_SET( xyz, x[0], 0, x[1] );	// the point on the intersecting line
    }
    return true;
  }
  
  // ================================================================ 
  // Triangle
  //
  // T     getRandomSampleOnTriangle(T v0[3], T v1[3], T v2[3], T p[3]);
  // void  getTriangleTangentPlane(T v0[3], T v1[3], T v2[3], T plane[4]);
  // void  getEdgeBoundary(T normal[3], T v0[3], T v1[3], T plane[4]);
  // T     DistanceTriangleAndPoint( T v0[3], T v1[3], T v2[3], T p[3], T *projection = NULL);
  // T     intersectTriangleAndRay( T v0[3], T v1[3], T v2[3], T from[3], T dir[3], T *pold = NULL);
  // void  getGradientTriangleArea(T v0[3], T v1[3], T v2[3], T result[3]);
  // ================================================================ 
  
  T getTriangleArea(T v0[3], T v1[3], T v2[3]) {
    T v10[3], v20[3], normal[3];
    G3V_SUB( v10, v1, v0 );
    G3V_SUB( v20, v2, v0 );
    G3V_CROSS( normal, v10, v20 );
    return G3V_LENGTH( normal ) * 0.5;
  }
  
  T getRandomSampleOnTriangle(T v0[3], T v1[3], T v2[3], T p[3]) {
    // get a uniformly random sample p[] on the triangle
    T w1, w2, v0v1[3], v0v2[3], distance;
    w1 = (T) random() / (T) RAND_MAX;
    w2 = (T) random() / (T) RAND_MAX;
    distance = 1 - (w1 + w2);
    if (distance < 0) { w1 = 1 - w1;  w2 = 1 - w2; }
    G3V_SUB( v0v1, v1, v0 );
    G3V_SUB( v0v2, v2, v0 );
    G3V_SCALED_ADD( p, v0, w1, v0v1 );
    G3V_SCALED_ADD( p, p,  w2, v0v2 );
    // return the weight of v0 for the sample.
    return distance * distance;
  }
    
  bool getTriangleNormal(T v0[3], T v1[3], T v2[3], T n[3]) {
    // calculate the plane equation (assuming the vertex order is CCW)
    T v10[3], v20[3];
    G3V_SUB( v10, v1, v0 );
    G3V_SUB( v20, v2, v0 );
    G3V_CROSS( n, v10, v20 );
    T length = G3V_LENGTH( n );
    if (length == 0) return false;
    G3V_DIV_VALUE( n, length );
    return true;
  }
  
  void  getTriangleTangentPlane(T v0[3], T v1[3], T v2[3], T plane[4]) {
    // calculate the plane equation (assuming the vertex order is CCW)
    T va[3], vb[3];
    G3V_SUB( va, v0, v1 );
    G3V_SUB( vb, v2, v1 );
    G3V_CROSS( plane, vb, va );
    G3V_NORMALIZE( plane );
    plane[3] = - G3V_DOT( plane, v0 );
  }
  
  void  getEdgeBoundary(T normal[3], T v0[3], T v1[3], T plane[4]) {
    // assuming the vertex order is CCW
    // plane[0] * x + plane[1] * y + plane[2] * z + plane[3] = 0
    T v0v1[3];
    G3V_SUB( v0v1, v1, v0 );
    // plane[] = cross product of normal[] and v0v1[]
    G3V_CROSS( plane, v0v1, normal );
    G3V_NORMALIZE( plane );
    plane[3] = - G3V_DOT( plane, v0 );
  }
  
  T DistanceTriangleAndPoint( T v0[3], T v1[3], T v2[3], T p[3], 
			      T *xyz_from = NULL) {
    // assuming the vertex order is CCW
    // calculate the plane equation (assuming the vertex order is CCW)
    T plane[4], dist, dist2=0;
    getTriangleTangentPlane( v0, v1, v2, plane );
    // calculate edge boundaries
    T b01[4], b12[4], b20[4];
    getEdgeBoundary( plane, v0, v1, b01 );
    getEdgeBoundary( plane, v1, v2, b12 );
    getEdgeBoundary( plane, v2, v0, b20 );
    // distance between the point and tangent plane
    dist = G3V_PLANE_POINT( plane, p );
    // check if the point lies in the boundaries of edges
    bool in01 = ( G3V_PLANE_POINT( b01, p ) <= 0);
    bool in12 = ( G3V_PLANE_POINT( b12, p ) <= 0);
    bool in20 = ( G3V_PLANE_POINT( b20, p ) <= 0);
    // decide the final distance and nearest point on the triangle
    if ( in01 && in12 && in20) {
      if (xyz_from) {
	G3V_MUL_VALUE( plane, dist );
	G3V_SUB( xyz_from, p, plane );
      }
    } else {
      if      (!in20 && !in01)  { dist2 = G3V_DISTANCE( p, v0 );  if (xyz_from) G3V_COPY( xyz_from, v0 ); }
      else if (!in01 && !in12)  { dist2 = G3V_DISTANCE( p, v1 );  if (xyz_from) G3V_COPY( xyz_from, v1 ); }
      else if (!in12 && !in20)  { dist2 = G3V_DISTANCE( p, v2 );  if (xyz_from) G3V_COPY( xyz_from, v2 ); }
      else if (!in01)  dist2 = getDistanceLinesegPoint(v0, v1, p, xyz_from);
      else if (!in12)  dist2 = getDistanceLinesegPoint(v1, v2, p, xyz_from);
      else if (!in20)  dist2 = getDistanceLinesegPoint(v2, v0, p, xyz_from);
      dist = (dist > 0 ? +dist2 : -dist2);
    }
    return dist;
  }
  
  T intersectTriangleAndRay( T v0[3], T v1[3], T v2[3], T from[3], T dir[3], T *pold = NULL) {
    T plane[4], screen[4], pnew[3], *p;
    if (pold) p = pold; else p = pnew;
    // calculate the intersection point of the plane and the line
    getTriangleTangentPlane( v0, v1, v2, plane );  // assuming the vertex order is CCW
    if ( intersectPlaneAndLineDir( plane, from, dir, p ) == 0 ) return -1;
    //cout << "Tri : v0 " << G3V_COUT(v0) << "  v1 " << G3V_COUT(v1) << "  v2 " << G3V_COUT(v2) << "  plane " << G3V_PCOUT(plane) << endl;
    //cout << "Ray : from " << G3V_COUT(from) << "  dir " << G3V_COUT(dir) << "  intersected at " << G3V_COUT(p) << endl;
    // check if the ray intersects in the front
    G3V_COPY( screen, dir );
    screen[3] = -G3V_DOT( screen, from );
    if ( G3V_PLANE_POINT( screen, p ) < 0 ) return -1;
    // calculate edge boundaries
    T b01[4], b12[4], b20[4];
    getEdgeBoundary( plane, v0, v1, b01 );
    getEdgeBoundary( plane, v1, v2, b12 );
    getEdgeBoundary( plane, v2, v0, b20 );
    // check if the intersection point lies in the boundaries of edges
    if ( G3V_PLANE_POINT( b01, p ) > 0 ) return -1;
    if ( G3V_PLANE_POINT( b12, p ) > 0 ) return -1;
    if ( G3V_PLANE_POINT( b20, p ) > 0 ) return -1;
    // return the distance 
    return G3V_DISTANCE( from, p );
  }
  
  void getGradientTriangleArea(T v0[3], T v1[3], T v2[3], T result[3]) {
    // A(v0)  = 0.5 * || v2v1 x v2v0 ||
    // A'(v0) = 1/|| v2v1 x v2v0 || * [(v2v1 x v2v0) . (v2v1 x v2v0)']
    T v2v1[3], v2v0[3], cross[3], length;
    T e[3][3] = { {1,0,0}, {0,1,0}, {0,0,1} };
    T ecross0[3], ecross1[3], ecross2[3];
    G3V_SUB( v2v1, v1, v2 );
    G3V_SUB( v2v0, v0, v2 );
    G3V_CROSS( cross, v2v1, v2v0 );
    length = G3V_LENGTH( cross );
    if (length == 0) return;
    G3V_CROSS( ecross0, v2v1, e[0] );
    G3V_CROSS( ecross1, v2v1, e[1] );
    G3V_CROSS( ecross2, v2v1, e[2] );
    result[0] = G3V_DOT( cross, ecross0 );
    result[1] = G3V_DOT( cross, ecross1 );
    result[2] = G3V_DOT( cross, ecross2 );
    G3V_DIV_VALUE( result, length );
  }
  
  bool getCircumcenter( T center[3], T v0[3], T v1[3], T v2[3]) {
    T c01[3], c12[3], P1[4], P2[4], P0[4];
    // calculate the centers of edges 
    G3V_AVERAGE2( c01, v0, v1 );
    G3V_AVERAGE2( c12, v1, v2 );
    // calculate the equations of the perpendicular planes of each edge
    getPlane( P0, v0, v1, v2 );					// ax + by + cz + j = 0
    G3V_SUB( P1, v1, v0 );  P1[3] = -G3V_DOT( P1, c01 );	// dx + ey + fz + k = 0
    G3V_SUB( P2, v2, v1 );  P2[3] = -G3V_DOT( P2, c12 );	// gx + hy + iz + l = 0
    // find the center of the circumcircle (intersection)
    return getIntersectionOf3Planes( center, P0, P1, P2 );
  }

  bool getTriangleCircumsphere(T v0[3], T v1[3], T v2[3], T center[3], T *radius) {
    T c01[3], c12[3], P1[4], P2[4], P0[4];
    // calculate the centers of edges 
    G3V_AVERAGE2( c01, v0, v1 );
    G3V_AVERAGE2( c12, v1, v2 );
    // calculate the equations of the perpendicular planes of each edge
    getPlane( P0, v0, v1, v2 );					// ax + by + cz + j = 0
    G3V_SUB( P1, v1, v0 );  P1[3] = -G3V_DOT( P1, c01 );	// dx + ey + fz + k = 0
    G3V_SUB( P2, v2, v1 );  P2[3] = -G3V_DOT( P2, c12 );	// gx + hy + iz + l = 0
    // find the center of the circumcircle (intersection)
    return getIntersectionOf3Planes( center, P0, P1, P2 );
  }

  // ================================================================ 
  // Tetrahedron
  //
  // bool getTetrahedronCircumsphere( T v0[3], T v1[3], T v2[3], T v3[3], T center[3], T &radius ) ;
  // bool getTetrahedronCircumsphere( T v0[3], T v1[3], T v2[3], T v3[3], T center[3], T &radius ) ;
  // bool getTetrahedronIncenter( T v0[3], T v1[3], T v2[3], T v3[3], T center[3], T &radius ) ;
  // T    getTetrahedronAspectRatio(T v0[3], T v1[3], T v2[3], T v3[3]) ;
  // ================================================================ 
  
  T intersectTetrahedronAndLineDir( T v0[3], T v1[3], T v2[3], T v3[3], T from[3], T dir[3], T p[3] ) {
    T plane[4], cp[3], dist, b01[4], b12[4], b20[4];
    T *a=NULL, *b=NULL, *c=NULL, temp, min = FLT_MAX;
    for (int i = 0; i < 4; i++) {
      switch (i) {
      case 0:  a = v0;  b = v1;  c = v2;  break;
      case 1:  a = v1;  b = v3;  c = v2;  break;
      case 2:  a = v2;  b = v3;  c = v0;  break;
      case 3:  a = v3;  b = v1;  c = v0;  break;
      }
      getTriangleTangentPlane( a, b, c, plane );
      if ((temp = G3V_DOT( plane, dir )) == 0) continue;
      dist = - ( G3V_DOT( plane, from ) + plane[3] ) / temp;
      G3V_SCALED_ADD( cp, from, dist, dir );
      // calculate edge boundaries
      getEdgeBoundary( plane, a, b, b01 );
      getEdgeBoundary( plane, b, c, b12 );
      getEdgeBoundary( plane, c, a, b20 );
      // check if the intersection point lies in the boundaries of edges
      if ( G3V_PLANE_POINT( b01, cp ) > 0 ) continue;
      if ( G3V_PLANE_POINT( b12, cp ) > 0 ) continue;
      if ( G3V_PLANE_POINT( b20, cp ) > 0 ) continue;
      if (dist < min) { min = dist;  G3V_COPY( p, cp ); }
    }
    return min;		// if it does not intersect, min = FLT_MAX;
  }
    
  T getTetrahedronVolume(T v0[3], T v1[3], T v2[3], T v3[3] ) {
    T v10[3], v20[3], v30[3], det;
    G3V_SUB( v10, v1, v0 );	//                    ( v10x v10y v10z )
    G3V_SUB( v20, v2, v0 );	//  volume = 1/6 * det( v20x v20y v20z )
    G3V_SUB( v30, v3, v0 );	//                    ( v30x v30y v30z ) 
    det = ( + v10[0] * ( v20[1] * v30[2] - v20[2] * v30[1] )
	    - v10[1] * ( v20[0] * v30[2] - v20[2] * v30[0] )
	    + v10[2] * ( v20[0] * v30[1] - v20[1] * v30[0] ) );
    //if (det < 0) det *= -1;
    return -det / 6.0;
  }
  
  bool testPointInTetrahedron(T p[3], T v0[3], T v1[3], T v2[3], T v3[3]) {
    T vol0 = getTetrahedronVolume( p, v0, v2, v1 );
    T vol1 = getTetrahedronVolume( p, v0, v3, v2 );
    T vol2 = getTetrahedronVolume( p, v0, v1, v3 );
    T vol3 = getTetrahedronVolume( p, v1, v2, v3 );
    return (vol0 >= 0 && vol1 >= 0 && vol2 >= 0 && vol3 >= 0);
  }
  
  void parameterizePointInTetrahedron(T p[3], T v0[3], T v1[3], T v2[3], T v3[3], T param[4]) {
    T A[9], B[3], I[9];
    // a*v0 + b*v1 + c*v2 + d*v3 = p, while a + b + c + d = 1.
    A[0]=v0[0]-v3[0];  A[1]=v1[0]-v3[0];  A[2]=v2[0]-v3[0];  B[0]=p[0]-v3[0];
    A[3]=v0[1]-v3[1];  A[4]=v1[1]-v3[1];  A[5]=v2[1]-v3[1];  B[1]=p[1]-v3[1];
    A[6]=v0[2]-v3[2];  A[7]=v1[2]-v3[2];  A[8]=v2[2]-v3[2];  B[2]=p[2]-v3[2];
    // A^{-1} : inverse of 3x3 matrix
    T det = ( + A[0] * A[4] * A[8] + A[1] * A[5] * A[6] + A[2] * A[3] * A[7]
	      - A[0] * A[5] * A[7] - A[1] * A[3] * A[8] - A[2] * A[4] * A[6] );
    I[0] = +(A[4]*A[8]-A[5]*A[7]);  I[3] = -(A[3]*A[8]-A[5]*A[6]);  I[6] = +(A[3]*A[7]-A[4]*A[6]);  
    I[1] = -(A[1]*A[8]-A[2]*A[7]);  I[4] = +(A[0]*A[8]-A[2]*A[6]);  I[7] = -(A[0]*A[7]-A[1]*A[6]);  
    I[2] = +(A[1]*A[5]-A[2]*A[4]);  I[5] = -(A[0]*A[5]-A[2]*A[3]);  I[8] = +(A[0]*A[4]-A[1]*A[3]);  
    for (int i = 0; i < 9; i++) I[i] /= det;
    // A^{-1} * B
    param[0] = I[0] * B[0] + I[1] * B[1] + I[2] * B[2];
    param[1] = I[3] * B[0] + I[4] * B[1] + I[5] * B[2];
    param[2] = I[6] * B[0] + I[7] * B[1] + I[8] * B[2];
    param[3] = 1 - param[0] - param[1] - param[2];
  }
    
  bool getTetrahedronCircumsphere( T v0[3], T v1[3], T v2[3], T v3[3], T center[3], T *radius ) {
    // | C - Vi | = R     after squaring and subtracting the equation for i = 0,
    // 2 (Vi - Vo).C = Vi.Vi - Vo.Vo	or, equivalently,
    // (Vi - Vo).(C - Vo) = 1/2 (Vi - Vo).(Vi - Vo) = 1/2 |Vi - Vo|^2
    // This is 3 equations for 3 unknowns.
    T v10[3], v20[3], v30[3], det;
    G3V_SUB( v10, v1, v0 );	//                    ( v10x v10y v10z )
    G3V_SUB( v20, v2, v0 );	//  volume = 1/6 * det( v20x v20y v20z )
    G3V_SUB( v30, v3, v0 );	//                    ( v30x v30y v30z ) 
    det = ( + v10[0] * ( v20[1] * v30[2] - v20[2] * v30[1] )
	    - v10[1] * ( v20[0] * v30[2] - v20[2] * v30[0] )
	    + v10[2] * ( v20[0] * v30[1] - v20[1] * v30[0] ) );
    if (det == 0) return false;
    T L10 = G3V_LENGTH2( v10 );
    T L20 = G3V_LENGTH2( v20 );
    T L30 = G3V_LENGTH2( v30 );
    center[0] = v0[0] + ( + (v20[1]*v30[2] - v20[2]*v30[1]) * L10
			  - (v10[1]*v30[2] - v10[2]*v30[1]) * L20
			  + (v10[1]*v20[2] - v10[2]*v20[1]) * L30 ) / (2 * det);
    center[1] = v0[1] + ( - (v20[0]*v30[2] - v20[2]*v30[0]) * L10
			  + (v10[0]*v30[2] - v10[2]*v30[0]) * L20
			  - (v10[0]*v20[2] - v10[2]*v20[0]) * L30 ) / (2 * det);
    center[2] = v0[2] + ( + (v20[0]*v30[1] - v20[1]*v30[0]) * L10
			  - (v10[0]*v30[1] - v10[1]*v30[0]) * L20
			  + (v10[0]*v20[1] - v10[1]*v20[0]) * L30 ) / (2 * det);
    *radius = G3V_DISTANCE( center, v0 );
    return true;
  }

  bool getTetrahedronInsphere( T v0[3], T v1[3], T v2[3], T v3[3], T center[3], T *radius ) {
    // Ni . (C - Vi) = R	where Ni are inward normals.  Equivalently,
    // (Ni,-1).(C,R) = Ni.Vi
    // This is 4 equations for 4 unknowns.
    
    // T M[4][4], b[4];
    // if (getTriangleNormal( v0, v2, v1, M[0] ) == false) return false;
    // if (getTriangleNormal( v1, v2, v3, M[1] ) == false) return false;
    // if (getTriangleNormal( v2, v0, v3, M[2] ) == false) return false;
    // if (getTriangleNormal( v3, v0, v1, M[3] ) == false) return false;
    // M[0][3] = M[1][3] = M[2][3] = M[3][3] = -1;
    // b[0] = G3V_DOT( M[0], v0 );
    // b[1] = G3V_DOT( M[1], v1 );
    // b[2] = G3V_DOT( M[2], v2 );
    // b[3] = G3V_DOT( M[3], v3 );
    T va[3], vb[3], N[3], L0, L1, L2, L3, L;
    G3V_SUB(va,v2,v0);  G3V_SUB(vb,v1,v0);  G3V_CROSS(N,va,vb);  L0 = G3V_LENGTH(N);
    G3V_SUB(va,v2,v1);  G3V_SUB(vb,v3,v1);  G3V_CROSS(N,va,vb);  L1 = G3V_LENGTH(N);
    G3V_SUB(va,v0,v2);  G3V_SUB(vb,v3,v2);  G3V_CROSS(N,va,vb);  L2 = G3V_LENGTH(N);
    G3V_SUB(va,v0,v3);  G3V_SUB(vb,v1,v3);  G3V_CROSS(N,va,vb);  L3 = G3V_LENGTH(N);
    if (L0 == 0 || L1 == 0 || L2 == 0 || L3 == 0) return false;
    L = L0 + L1 + L2 + L3;
    center[0] = ( v0[0]*L3 + v1[0]*L2 + v2[0]*L1 + v3[0]*L0 ) / L ;
    center[1] = ( v0[1]*L3 + v1[1]*L2 + v2[1]*L1 + v3[1]*L0 ) / L ;
    center[2] = ( v0[2]*L3 + v1[2]*L2 + v2[2]*L1 + v3[2]*L0 ) / L ;
    *radius   = abs ( + ( + v1[0] * ( v2[1]*v3[2] - v3[1]*v2[2] )
			  - v2[0] * ( v1[1]*v3[2] - v3[1]*v1[2] )
			  + v3[0] * ( v1[1]*v2[2] - v2[1]*v1[2] ) ) +
		      - ( + v0[0] * ( v2[1]*v3[2] - v3[1]*v2[2] )
			  - v2[0] * ( v0[1]*v3[2] - v3[1]*v0[2] )
			  + v3[0] * ( v0[1]*v2[2] - v2[1]*v0[2] ) ) +
		      + ( + v0[0] * ( v1[1]*v3[2] - v3[1]*v1[2] )
			  - v1[0] * ( v0[1]*v3[2] - v3[1]*v0[2] )
			  + v3[0] * ( v0[1]*v1[2] - v1[1]*v0[2] ) ) +
		      - ( + v0[0] * ( v1[1]*v2[2] - v2[1]*v1[2] )
			  - v1[0] * ( v0[1]*v2[2] - v2[1]*v0[2] )
			  + v2[0] * ( v0[1]*v1[2] - v1[1]*v0[2] ) ) ) / L ;
    return true;
  }

  // =================================================================
  // Measurement of mesh quality
  // =================================================================

  T getTetrahedronEdgeLength(T v0[3], T v1[3], T v2[3], T v3[3], bool max) {
    T d[6], ext=0;
    d[0] = G3V_DISTANCE2( v0, v1 );  d[1] = G3V_DISTANCE2( v0, v2 );
    d[2] = G3V_DISTANCE2( v0, v3 );  d[3] = G3V_DISTANCE2( v1, v2 );
    d[4] = G3V_DISTANCE2( v1, v3 );  d[5] = G3V_DISTANCE2( v2, v3 );
    if (max) {
      for (int i=0; i<6; i++) if (i == 0 || d[i] > ext) ext = d[i];
    } else {
      for (int i=0; i<6; i++) if (i == 0 || d[i] < ext) ext = d[i];
    }
    return (T)(sqrt(ext));
  }
  
  T getTetrahedronAltitude(T v0[3], T v1[3], T v2[3], T v3[3], bool max) {
    T d[4], ext=0;
    d[0] = DistanceFromPointToPlane( v0, v1, v2, v3 );
    d[1] = DistanceFromPointToPlane( v1, v0, v3, v2 );
    d[2] = DistanceFromPointToPlane( v2, v0, v1, v3 );
    d[3] = DistanceFromPointToPlane( v3, v0, v2, v1 );
    if (max) {
      for (int i=0; i<4; i++) if (i == 0 || d[i] > ext) ext = d[i];
    } else {
      for (int i=0; i<4; i++) if (i == 0 || d[i] < ext) ext = d[i];
    }
    return ext;
  }
  
  T getTetrahedronAngle(T v0[3], T v1[3], T v2[3], T v3[3], bool max) {
    T p[4][4], angle, ext;
    getPlane(p[0], v1, v3, v2);
    getPlane(p[1], v0, v2, v3);
    getPlane(p[2], v0, v3, v1);
    getPlane(p[3], v0, v1, v2);
    ext = (max ? M_PI : 0);	// calculate angle between face normals
    for (int i = 0; i < 4; i++)
      for (int j = i+1; j < 4; j++) {
	angle = (T)acos(G3V_DOT(p[i], p[j]));
	if (max) { if (angle < ext) ext = angle; }
	else     { if (angle > ext) ext = angle; }
      }
    return (M_PI - ext);	// return angle between faces
  }
  
  T getTetrahedronAspectRatio(T v0[3], T v1[3], T v2[3], T v3[3], int type) {
    // Note :
    //   negative types : the bigger,  the better
    //   positive types : the smaller, the better
    float numerator, denominator, tmp[3], max_angle, ratio = 0;
    switch (type) {
    case +1:	// min_height / max_edge [ 0 => 0.8165 ]
      numerator   = getTetrahedronAltitude(v0, v1, v2, v3, false);
      denominator = getTetrahedronEdgeLength(v0, v1, v2, v3, true);
      ratio = (denominator==0 ? 0 : (numerator/denominator));
      break;
    case +2:	// min_height / max_edge + 1/4 * max_facenormal_angle
      numerator   = getTetrahedronAltitude(v0, v1, v2, v3, false);
      denominator = getTetrahedronEdgeLength(v0, v1, v2, v3, true);
      max_angle   = getTetrahedronAngle(v0, v1, v2, v3, true);
      if (denominator == 0) return 0;
      ratio = (numerator / denominator) + 0.25 * cos(max_angle);
      break;
    case -1:	// circumradius / min_edge [ 0.61237 <= $$ ]
      if (!getTetrahedronCircumsphere(v0, v1, v2, v3, tmp, &numerator)) return 0;
      denominator = getTetrahedronEdgeLength(v0, v1, v2, v3, false);
      ratio = (denominator==0 ? 0 : (numerator/denominator));
      break;
    case -2:	// circumradius / inradius  [ 3 <= $$ ]
      if (!getTetrahedronCircumsphere(v0, v1, v2, v3, tmp, &numerator)) return 0;
      if (!getTetrahedronInsphere(v0, v1, v2, v3, tmp, &denominator)) return 0;
      ratio = (denominator==0 ? 0 : (numerator/denominator));
      break;
    default:
      cerr << "Warning: Invalid type for getTetrahedronAspectRatio()" << endl;
      ratio = 0;
    }
    return ratio;
  }
  
};

#endif  // CGEOMETRY3D_HPP

