
//
// Geometry2D<T> class template
//
// Jaeil Choi
// last modified in June, 2003
//
// Note that functions of template classes are included 
// in the source file only when they are called.
//

#ifndef CGEOMETRY2D_HPP
#define CGEOMETRY2D_HPP

#include <iostream>
#include <iomanip>
#include <cmath>
#include "vm_macros.h"
using namespace std;

#ifndef M_PI
#define M_PI (3.14159265358979323846)
#endif


template <class T>
class Geometry2D {
  
 public:
  
  // ================================================================
  // Lines & Line Segments
  //
  // T     getRandomSampleOnLineSegment(T v0[2], T v1[2], T p[2]);
  // T     getDistanceLinesegPoint( T v0[2], T v1[2], T p[2], T *projection = NULL);
  // void  getGradientLineLength(T v0[2], T v1[2], T result[2]);
  // ================================================================
  
  bool getLineEquation(T v1[2], T v2[2], T line[3], bool bNormalize=true) {
    G2V_SET( line, v2[1]-v1[1], v1[0]-v2[0] );  
    T len2 = G2V_LENGTH2( line );
    if (len2 == 0) return false;
    if (bNormalize) G2V_DIV_VALUE( line, sqrt(len2) );
    line[2] = -G2V_DOT( line, v1 );
    return true;
  }
    
  bool testPointOntoLineSegment(T p[2], T v1[2], T v2[2], T tol=0.0) {
    // test if a point 'p' can be projected onto the line segment 'v1v2'
    T v12[2], v21[2], v1p[2], v2p[2];
    G2V_SUB( v12, v2, v1 );  G2V_SUB( v1p, p, v1 );
    G2V_SUB( v21, v1, v2 );  G2V_SUB( v2p, p, v2 );
    return (G2V_DOT( v1p, v12 ) >= -tol && G2V_DOT( v2p, v21 ) >= -tol);
  }
  
  bool getPointByIntersection(T p[2], T eqa[3], T eqb[3]) {
    // intersection of infinite lines
    p[0] = p[1] = 0;
    T det = ( eqa[0] * eqb[1] - eqa[1] * eqb[0] );
    if (fabs(det) < 0.0001) return false;
    p[0] = - (+ eqb[1] * eqa[2] - eqa[1] * eqb[2]) / det;
    p[1] = - (- eqb[0] * eqa[2] + eqa[0] * eqb[2]) / det;
    return true;	
  }
  
  void getPointByProjection(T p[2], T eq[3]) {
    // find the closest point (double) on the line, given the line equation
    T  v[2], vo[2], len;
    T  vl[2] = { +eq[1], -eq[0] };		// vl[2] : unit vector of the line direction
    // get a point 'v[2]' on the line
    if (fabs(eq[0]) > fabs(eq[1])) { v[1] = p[1];  v[0] = -(eq[1]*v[1]+eq[2])/eq[0]; }  // ax = - by - c
    else                           { v[0] = p[0];  v[1] = -(eq[0]*v[0]+eq[2])/eq[1]; }  // by = - ax - c
    // project the input point to line using 'v[2]'
    vo[0] = p[0] - v[0];  vo[1] = p[1] - v[1];	// vo[] = p[] - v[]
    len = G2V_DOT( vo, vl );			// len  = (vo . vl)
    G2V_SCALED_ADD( p, v, len, vl );		// p[] = v[] + len * vl[]
  }
  
  bool intersectLines(T v1[2], T v2[2], T v3[2], T v4[2], T p[2]) {
    // calculate the intersection point of two infinite lines, if possible.
#if 0
    T det = (v1[0] - v2[0]) * (v3[1] - v4[1]) - (v1[1] - v2[1]) * (v3[0] - v4[0]);
    if (det == 0) return false;
    T det2 = (v1[0]*v2[1]-v1[1]*v2[0]);
    T det3 = (v3[0]*v4[1]-v3[1]*v4[0]);
    p[0] = ( det2 * (v3[0] - v4[0]) - (v1[0] - v2[0]) * det3 ) / det ;
    p[1] = ( det2 * (v3[1] - v4[1]) - (v1[1] - v2[1]) * det3 ) / det ;
#else
    T den = ((v4[1] - v3[1])*(v2[0] - v1[0]) - (v4[0] - v3[0])*(v2[1] - v1[1]));
    if (den == 0) return false;
    T ua  = ((v4[0] - v3[0])*(v1[1] - v3[1]) - (v4[1] - v3[1])*(v1[0] - v3[0])) / den;
    p[0] = v1[0] + ua * (v2[0] - v1[0]);
    p[1] = v1[1] + ua * (v2[1] - v1[1]);
//     T ub  = ((v2[0] - v1[0])*(v1[1] - v3[1]) - (v2[1] - v1[1])*(v1[0] - v3[0])) / den;
//     p[0] = v3[0] + ub * (v4[0] - v3[0]);
//     p[1] = v3[1] + ub * (v4[1] - v3[1]);
#endif
    return true;
  }
  
  bool intersectLineSegments(T v1[2], T v2[2], T v3[2], T v4[2], T p[2]) {
    // calculate the intersection point of two line segments, if possible.
    return ( intersectLines( v1, v2, v3, v4, p ) &&
	     testPointOntoLineSegment( p, v1, v2 ) &&
	     testPointOntoLineSegment( p, v3, v4 ) );
  }
  
  T getRandomSampleOnLineSegment(T v0[2], T v1[2], T p[2]) {
    // get a uniformly random sample p[] on the line segment,
    T  w = (T) random() / (T) RAND_MAX;
    G2V_WEIGHTED_ADD( p, (w), v0, (1-w), v1 );
    // return the weight of v0 for the sample.
    return w;
  }
  
  T DistanceFromPointToLine( T p[2], T v1[2], T v2[2] ) {
    // calculate the distance from a point 'p' to the line 'v1v2'
    T line[3]={0,0,0};
    getLineEquation( v1, v2, line, true );
    return (G2V_DOT( line, p ) + line[2]);
  }
  
  void ProjectPointOntoLine( T p[2], T v1[2], T v2[2], T result[2] ) {
    T v12[2], v1p[2];
    G2V_SUB( v12, v2, v1 );
    G2V_NORMALIZE( v12 );
    G2V_SUB( v1p, p, v1 );
    G2V_SCALED_ADD( result, v1, G2V_DOT( v1p, v12 ), v12 );
  }
  
  T DistanceFromPointToLineSegment( T p[2], T v1[2], T v2[2], T *closest = NULL) {
    // calculate the distance from a point 'p' to the line segment 'v1v2'
    // returns SIGNED distance (positive if 'p' lies on the right-hand side)
    T v12[2], v21[2], v1p[2], v2p[2], line[3], dist;
    getLineEquation(v1, v2, line);
    dist = DistanceFromPointToLine( p, v1, v2 );
    int sign = (dist >= 0 ? +1 : -1);
    
    G2V_SUB( v12, v2, v1 );  G2V_SUB( v1p, p, v1 );
    G2V_SUB( v21, v1, v2 );  G2V_SUB( v2p, p, v2 );
    if        (G2V_DOT( v1p, v12 ) < 0) {
      if (closest)  G2V_COPY( closest, v1 );
      return sign * G2V_DISTANCE( p, v1 );
    } else if (G2V_DOT( v2p, v21 ) < 0) {
      if (closest)  G2V_COPY( closest, v2 );
      return sign * G2V_DISTANCE( p, v2 );
    } else {
      if (closest)  ProjectPointOntoLine( p, v1, v2, closest );
      return dist;
    }
  }
    
  void getGradientLineLength(T v0[2], T v1[2], T result[2]) {
    // L(v0)  = || v0 - v1 ||
    G2V_SUB( result, v0, v1 );
    G2V_NORMALIZE( result );
  }
  
  bool LeftOfLine(T p[2], T v1[2], T v2[2]) {
    // check if 'p' is on the left-hand side of the vector 'v1'->'v2'
    T v1p[2], v12[2];
    G2V_SUB( v1p, p, v1 );
    G2V_SUB( v12, v2, v1 );
    return ( (v12[0] * v1p[1] - v12[1] * v1p[0]) > 0 );   // rotate v1p[2] to the right
  }
    
  bool RightOfLine(T p[2], T v1[2], T v2[2]) {
    // check if 'p' is on the right-hand side of the vector 'v1'->'v2'
    T v1p[2], v12[2];
    G2V_SUB( v1p, p, v1 );
    G2V_SUB( v12, v2, v1 );
    return ( (v12[0] * v1p[1] - v12[1] * v1p[0]) < 0 );
  }
    
  bool testPointInCircumcircle( T p[2], T v0[2], T v1[2] ) {
    // Test if the point 'p' is inside of the diametral circle of the edge
    // (smallest circumcircle of the line segment)
    T center[2];
    G2V_AVERAGE2( center, v0, v1 );
    return ( G2V_DISTANCE2( center, p ) <= G2V_DISTANCE2( center, v0 ) );
  }
  
  
  // ================================================================ 
  // Triangle
  //
  // T     getRandomSampleOnTriangle(T v0[2], T v1[2], T v2[2], T p[2]);
  // void  getEdgeBoundary(T normal[2], T v0[2], T v1[2], T line[3]);
  // T     DistanceTriangleAndPoint( T v0[2], T v1[2], T v2[2], T p[2], T *projection = NULL);
  // ================================================================ 
  
  T calculateTriangleArea(T v0[2], T v1[2], T v2[2]) {
    T v01[2], v02[2];
    G2V_SUB( v01, v1, v0 );
    G2V_SUB( v02, v2, v0 );
    return (v01[0] * v02[1] - v01[1] * v02[0]) * 0.5;
  }
  
  T getRandomSampleOnTriangle(T v0[2], T v1[2], T v2[2], T p[2]) {
    // get a uniformly random sample p[] on the triangle
    T w1, w2, v0v1[2], v0v2[2], distance;
    w1 = (T) random() / (T) RAND_MAX;
    w2 = (T) random() / (T) RAND_MAX;
    distance = 1 - (w1 + w2);
    if (distance < 0) { w1 = 1 - w1;  w2 = 1 - w2; }
    G2V_SUB( v0v1, v1, v0 );
    G2V_SUB( v0v2, v2, v0 );
    G2V_SCALED_ADD( p, v0, w1, v0v1 );
    G2V_SCALED_ADD( p, p,  w2, v0v2 );
    // return the weight of v0 for the sample.
    return distance * distance;
  }
    
  void  getEdgeBoundary(T v0[2], T v1[2], T line[3]) {
    // assuming the vertex order is CCW
    // line[0] * x + line[1] * y + line[2] * z + line[3] = 0
    T v0v1[2];
    G2V_SUB( v0v1, v1, v0 );
    G2V_CROSS( line, v0v1 );
    G2V_NORMALIZE( line );
    line[2] = - G2V_DOT( line, v0 );
  }
  
  T DistanceTriangleAndPoint( T v0[2], T v1[2], T v2[2], T p[2], 
			      T *projection = NULL) {
    // assuming the vertex order is CCW
    // calculate the plane equation (assuming the vertex order is CCW)
    T dist, b01[3], b12[3], b20[3];
    getEdgeBoundary( v0, v1, b01 );
    getEdgeBoundary( v1, v2, b12 );
    getEdgeBoundary( v2, v0, b20 );
    // check if the point lies in the boundaries of edges
    bool in01 = ( G2V_LINE_POINT( b01, p ) <= 0);
    bool in12 = ( G2V_LINE_POINT( b12, p ) <= 0);
    bool in20 = ( G2V_LINE_POINT( b20, p ) <= 0);
    // decide the final distance and nearest point on the triangle
    if ( in01 && in12 && in20) {
      dist = 0;
      if (projection) G2V_COPY( projection, p );
    } else {
      if      (!in20 && !in01)  { dist = G2V_DISTANCE( p, v0 );  if (projection) G2V_COPY( projection, v0 ); }
      else if (!in01 && !in12)  { dist = G2V_DISTANCE( p, v1 );  if (projection) G2V_COPY( projection, v1 ); }
      else if (!in12 && !in20)  { dist = G2V_DISTANCE( p, v2 );  if (projection) G2V_COPY( projection, v2 ); }
      else if (!in01)  dist = DistanceFromPointToLineSegment(p, v0, v1, projection);
      else if (!in12)  dist = DistanceFromPointToLineSegment(p, v1, v2, projection);
      else if (!in20)  dist = DistanceFromPointToLineSegment(p, v2, v0, projection);
    }
    return dist;
  }
  
  int testPointInTriangle( T p[2], T v0[2], T v1[2], T v2[2] ) {
    // assuming the vertex order is CCW
    // calculate the plane equation (assuming the vertex order is CCW)
    T b01[3], b12[3], b20[3];
    getEdgeBoundary( v0, v1, b01 );
    getEdgeBoundary( v1, v2, b12 );
    getEdgeBoundary( v2, v0, b20 );
    // check if the point lies in the boundaries of edges
    T d01 = G2V_LINE_POINT( b01, p );
    T d12 = G2V_LINE_POINT( b12, p );
    T d20 = G2V_LINE_POINT( b20, p );
    if      (d01 < 0 && d12 < 0 && d20 < 0) return +1;
    else if (d01 > 0 || d12 > 0 || d20 > 0) return 0;
    else return -1;
//     else if ((d01 == 0 || d12 == 0 || d20 == 0) &&
// 	     (d01 <= 0 && d12 <= 0 && d20 <= 0)) return -1;
//     else return 0;
  }
  
  bool getCircumcenter( T center[2], T v0[2], T v1[2], T v2[2] ) {
    T c01[2], c12[2], L01[3], L12[3];
    // calculate the centers of edges 
    G2V_AVERAGE2( c01, v0, v1 );
    G2V_AVERAGE2( c12, v1, v2 );
    // calculate the equations of the perpendicular lines of the edges
    G2V_SUB( L01, v1, v0 );  L01[2] = G2V_DOT( L01, c01 );	// ax + by = e
    G2V_SUB( L12, v2, v1 );  L12[2] = G2V_DOT( L12, c12 ); 	// cx + dy = f
    // find the center of the circumcircle (intersection)
    T det = L01[0] * L12[1] - L01[1] * L12[0];			// det = a * d - b * c
    if (det == 0) return false;
    center[0] = ( L12[1] * L01[2] - L01[1] * L12[2]) / det;	// x  = ( d * e - b * f) / det
    center[1] = (-L12[0] * L01[2] + L01[0] * L12[2]) / det;	// y  = (-c * e + a * f) / det
    return true;
  }
  
  bool testPointInCircumcircle( T p[2], T v0[2], T v1[2], T v2[2], T *ratio=NULL ) {
    // Note that this function is vulnerable to floating point errors.
    T center[2], radius, dist;
    if (!getCircumcenter( center, v0, v1, v2 )) {
      bool v0v1eq = G2V_EQUAL(v0,v1);
      bool v1v2eq = G2V_EQUAL(v1,v2);
      bool v2v0eq = G2V_EQUAL(v2,v0);
      if (v0v1eq && v1v2eq) G2V_COPY( center, v0 );
      else if      (v0v1eq) G2V_AVERAGE2( center, v0, v2 );
      else if      (v1v2eq) G2V_AVERAGE2( center, v0, v1 );
      else if      (v2v0eq) G2V_AVERAGE2( center, v0, v1 );
      else { // three differents points are on a same line
	if (ratio) *ratio = -1; return true;
      }
    }
    // compare the distances from the center of the circumcircle
    radius = (T)G2V_DISTANCE2( center, v0 );  if (radius<=0) return false;
    dist   = (T)G2V_DISTANCE2( center, p  ); 
    if (ratio) *ratio = dist / radius;
    return ( dist <= radius );
  }
  
  float getCircumradiusRatio( T v0[2], T v1[2], T v2[2] ) {
    T center[2], radius2, d01, d12, d02, min2;
    getCircumcenter( center, v0, v1, v2 );
    radius2 = G2V_DISTANCE2( center, v0 );
    d01 = G2V_DISTANCE2( v0, v1 );
    d12 = G2V_DISTANCE2( v1, v2 );
    d02 = G2V_DISTANCE2( v0, v2 );
    min2 = (d01 <= d12) ? ((d01 <= d02) ? d01 : d02) : ((d12 <= d02) ? d12 : d02);
    if (min2 == 0) return 1000; 
    return sqrt( radius2 / min2 );
  }
  
  void getTriangleAngles( T v0[2], T v1[2], T v2[2], float angles[3], bool bRadian = true ) {
    float d[3], a, b, c, angle;
    d[0] = G2V_DISTANCE( v1, v2 );
    d[1] = G2V_DISTANCE( v2, v0 );
    d[2] = G2V_DISTANCE( v0, v1 );
    if (d[0] == 0 || d[1] == 0 || d[2] == 0) return;
    for (int i = 0; i < 3; i++) {
      a = d[i];  b = d[(i+1)%3];  c = d[(i+2)%3];
      angle = (float) acos (( b*b + c*c - a*a ) / ( 2 * b * c ));
      if (angle < 0) angle = (float)M_PI + angle;
      if (!bRadian) angle = angle * 180 / (float)M_PI;
      angles[i] = angle;
    }
  }
  
  // ================================================================ 
  // Parabola
  // ================================================================ 
  
  void getParabolicEquation( T p[2], T v0[2], T v1[2], T eq[6] ) {
    // calculate the equation of parabola 
    // that is equidistant from a point 'p' and a line 'v0v1'
    T line[3]={0,0,0};
    getLineEquation( v0, v1, line, false );
    // equation of parabola : a x^2 + b y^z + c x y + d x + e y + f = 0
    eq[0] = line[0] * line[0] - 1;
    eq[1] = line[1] * line[1] - 1;
    eq[2] = 2 * line[0] * line[1];
    eq[3] = 2 * line[0] * line[2] + 2 * p[0];
    eq[4] = 2 * line[1] * line[2] + 2 * p[1];
    eq[5] = line[2] * line[2] - p[0] * p[0] - p[1] * p[1];
  }
  
  int intersectLineAndParabola( T line[3], T para[6], T p1[2], T p2[2] ) {
    T   x, rt, eq[3];	// a x^2 + b x + c = 0
    int nIntersections;
    if (line[1] == 0) {	// the line is vertical, solve for 'y'
      x = - line[2] / line[0];
      eq[0] = para[1];
      eq[1] = para[2] * x + para[4];
      eq[2] = para[0] * x * x + para[3] * x + para[5];
      if (fabs(eq[0]) < 0.000001) {	// orientation is same (intersect at single point)
	p1[0] = x;
	p1[1] = - eq[2] / eq[1];
	nIntersections = 1;
      } else {
	rt = eq[1] * eq[1] - 4 * eq[0] * eq[2] ;
	if (rt < 0)  nIntersections = 0;		// No intersection
	else if (rt == 0) {				// the line is tangent
	  p1[0] = x;
	  p1[1] = - eq[1] * eq[1] / ( 2 * eq[0] );
	  nIntersections = 1;
	} else {
	  rt = sqrt(rt);				// intersect at two points
	  p1[0] = x;
	  p1[1] = - (eq[1] * eq[1] - rt) / ( 2 * eq[0] );
	  p2[0] = x;
	  p2[1] = - (eq[1] * eq[1] + rt) / ( 2 * eq[0] );
	  nIntersections = 2;
	}
      }
    } else {		// replace 'y', solve for 'x'
      eq[0] = ( + para[0] 
		+ para[1] * line[0]*line[0]/(line[1]*line[1]) 
		- para[2] * line[0]/line[1] );
      eq[1] = ( + para[3]
		+ para[1] * 2*line[0]*line[2]/(line[1]*line[1])
		- para[2] * line[2]/line[1]
		- para[4] * line[0]/line[1] );
      eq[2] = ( + para[5]
		+ para[1] * line[2]*line[2]/(line[1]*line[1])
		- para[4] * line[2]/line[1] );
      // x = ( - b^2 +/- sqrt( b^2 - 4 a c ) ) / (2 a)
      if (fabs(eq[0]) < 0.000001) {	// orientation is same (intersect at single point)
	p1[0] = - eq[2] / eq[1];
	p1[1] = - (line[0] * p1[0] + line[2] ) / line[1];
	nIntersections = 1;
      } else {
	rt = eq[1] * eq[1] - 4 * eq[0] * eq[2];
	if (rt < 0)  nIntersections = 0;		// No intersection
	else if (rt == 0) {				// the line is tangent
	  p1[0] = - eq[1] * eq[1] / ( 2 * eq[0] );
	  p1[1] = - (line[0] * p1[0] + line[2] ) / line[1];
	  nIntersections = 1;
	} else {					// intersect at two points
	  rt = sqrt(rt);
	  p1[0] = - (eq[1] * eq[1] - rt) / ( 2 * eq[0] );
	  p1[1] = - (line[0] * p1[0] + line[2] ) / line[1];
	  p2[0] = - (eq[1] * eq[1] + rt) / ( 2 * eq[0] );
	  p2[1] = - (line[0] * p1[0] + line[2] ) / line[1];
	  nIntersections = 2;
	}
      }
    }
    return nIntersections;
  }
  
};

#endif  // CGEOMETRY2D_HPP

