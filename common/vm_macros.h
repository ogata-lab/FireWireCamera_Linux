
//
// Macros for vectors and matrices
//
// Jaeil Choi
// last modified in July, 2006
//
// After using classes (Vector3D and Matrix3D, for example) for a while,
// I gave up those classes (at least for basic vector/matrix operations),
// and returned to simple macros. Macros are faster and easy to use.
// They may slightly increase the size after compilation, but who really 
// cares the size of the program these days?
//

#ifndef VECTOR_MATRIX_MACROS_H
#define VECTOR_MATRIX_MACROS_H

#include <iostream>
#include <iomanip>
#include <cmath>

#ifndef M_PI
#define M_PI (3.14159265358979323846)
#endif

/* ================================================================ */
/*  2D                                                              */
/* ================================================================ */

#define G2V_SET(r,x,y)          do { (r)[0] = (x);  (r)[1] = (y); } while(0)
#define G2V_COPY(d,s)           do { (d)[0] = (s)[0];  (d)[1] = (s)[1]; } while(0)

#define G2V_ADD(r,a,b)          do { (r)[0] = (a)[0] + (b)[0];  (r)[1] = (a)[1] + (b)[1]; } while(0)
#define G2V_SUB(r,a,b)          do { (r)[0] = (a)[0] - (b)[0];  (r)[1] = (a)[1] - (b)[1]; } while(0)
#define G2V_ADD_VALUES(r,a,b)   do { (r)[0] += (a);  (r)[1] += (b); } while(0)
#define G2V_SUB_VALUES(r,a,b)   do { (r)[0] -= (a);  (r)[1] -= (b); } while(0)
#define G2V_MUL_VALUES(r,a,b)   do { (r)[0] *= (a);  (r)[1] *= (b); } while(0)
#define G2V_DIV_VALUES(r,a,b)   do { (r)[0] /= (a);  (r)[1] /= (b); } while(0)
#define G2V_ADD_VALUE(a,v)      do { (a)[0] += (v);  (a)[1] += (v); } while(0)
#define G2V_SUB_VALUE(a,v)      do { (a)[0] -= (v);  (a)[1] -= (v); } while(0)
#define G2V_MUL_VALUE(a,v)      do { (a)[0] *= (v);  (a)[1] *= (v); } while(0)
#define G2V_DIV_VALUE(a,v)      do { (a)[0] /= (v);  (a)[1] /= (v); } while(0)
#define G2V_SCALED_ADD(r,a,s,b) do { (r)[0] = (a)[0] + (b)[0]*(s);  (r)[1] = (a)[1] + (b)[1]*(s); } while(0)
#define G2V_WEIGHTED_ADD(r,w1,a,w2,b) do { (r)[0] = (a)[0]*(w1) + (b)[0]*(w2);  (r)[1] = (a)[1]*(w1) + (b)[1]*(w2); } while(0)

#define G2V_DOT(a,b)            ((a)[0]*(b)[0] + (a)[1]*(b)[1])
#define G2V_CROSS(r,a)          do { (r)[0] = (a)[1];  (r)[1] = -(a)[0]; } while(0)
#define G2V_LINE_POINT(l,a)     ((l)[0]*(a)[0] + (l)[1]*(a)[1] + (l)[2])
#define G2V_LINE_EQ(eq,a,b)     do { G2V_SET(eq, +((b)[1]-(a)[1]), -((b)[0]-(a)[0])); G2V_NORMALIZE(eq); eq[2] = -G2V_DOT(eq,a); } while(0)
#define G2V_LINE_EQ_NP(eq,nx,ny,px,py)  do { G2V_SET(eq, nx, ny); G2V_NORMALIZE(eq); eq[2] = -(eq[0]*px+eq[1]*py); } while(0)
#define G2V_POINT_ON_LINE(r,p,eq) do { float d=G2V_LINE_POINT(eq,p); G2V_SCALED_ADD(r, p, -d, eq); } while(0)
// normal pointing outward, when (a,b) in given in CCW order

#define G2V_DISTANCE2(a,b)      (((b)[0]-(a)[0])*((b)[0]-(a)[0]) + ((b)[1]-(a)[1])*((b)[1]-(a)[1]))
#define G2V_DISTANCE(a,b)       (sqrt(((b)[0]-(a)[0])*((b)[0]-(a)[0]) + ((b)[1]-(a)[1])*((b)[1]-(a)[1])))
#define G2V_DIST_ABS(a,b)       (fabs((b)[0]-(a)[0]) + fabs((b)[1]-(a)[1]))
#define G2V_LENGTH2(a)          ((a)[0]*(a)[0] + (a)[1]*(a)[1])
#define G2V_LENGTH(a)           (sqrt((a)[0]*(a)[0] + (a)[1]*(a)[1]))
#define G2V_NORMALIZE(a)        do { float len = G2V_LENGTH(a);  if (len > 0) G2V_DIV_VALUE( a, len ); } while(0)
#define G2V_NORMALIZED(a)       do { double len = G2V_LENGTH(a);  if (len > 0) G2V_DIV_VALUE( a, len ); } while(0)

#define G2V_EQUAL(a,b)          ((a)[0] == (b)[0] && (a)[1] == (b)[1])
#define G2V_AVERAGE2(r,a,b)     do { (r)[0] = ((a)[0] + (b)[0])/2;  (r)[1] = ((a)[1] + (b)[1])/2; } while(0)
#define G2V_AVERAGE3(r,a,b,c)   do { (r)[0] = ((a)[0] + (b)[0] + (c)[0])/3;  (r)[1] = ((a)[1] + (b)[1] + (c)[1])/3; } while(0)
#define G2V_AVERAGE4(r,a,b,c,d) do { (r)[0] = ((a)[0] + (b)[0] + (c)[0] + d[0])/4;  (r)[1] = ((a)[1] + (b)[1] + (c)[1] + d[1])/4; } while(0)
#define G2V_MINV(a)             ( (a)[0]<(a)[1] ? (a)[0] : (a)[1] )
#define G2V_MAXV(a)             ( (a)[0]>(a)[1] ? (a)[0] : (a)[1] )
#define G2V_LOWER_LIMIT(a,l)    do { if ((a)[0]<(l)[0]) (a)[0]=(l)[0]; if ((a)[1]<(l)[1]) (a)[1]=(l)[1]; } while(0)
#define G2V_UPPER_LIMIT(a,u)    do { if ((a)[0]>(u)[0]) (a)[0]=(u)[0]; if ((a)[1]>(u)[1]) (a)[1]=(u)[1]; } while(0)
#define G2V_ABS_LIMIT(a,u)      do { if ((a)[0]>(u)[0]) (a)[0]=(u)[0]; if ((a)[1]>(u)[1]) (a)[1]=(u)[1]; if ((a)[0]<-(u)[0]) (a)[0]=-(u)[0]; if ((a)[1]<-(u)[1]) (a)[1]=-(u)[1]; } while(0)
#define G2V_LIMIT(a,l,u)        do { G2V_LLIMIT(a,l); G2V_ULIMIT(a,u); } while(0)
#define G2V_SWAP(a,b,t)         do { G2V_COPY(t, a); G2V_COPY(a, b); G2V_COPY(b, t); } while(0)
#define G2V_BBOX_INITIALIZE(a, bmin, bmax) do { G2V_COPY(bmin,a); G2V_COPY(bmax,a); } while(0)
#define G2V_BBOX_UPDATE(a, bmin, bmax) do { if ((a)[0]<(bmin)[0]) (bmin)[0] = (a)[0];  if ((a)[1]<(bmin)[1]) (bmin)[1] = (a)[1];  if ((a)[0]>(bmax)[0]) (bmax)[0] = (a)[0];  if ((a)[1]>(bmax)[1]) (bmax)[1] = (a)[1]; } while(0)
#define G2V_BBOX_LENGTH(bmin, bmax) ((bmax)[0]-(bmin)[0] > (bmax)[1]-(bmin)[1] ? (bmax)[0]-(bmin)[0] : (bmax)[1]-(bmin)[1]);
#define G2V_BBOX_CENTER(c, bmin, bmax) do { (c)[0] = ((bmin)[0]+(bmax)[0])/2; (c)[1] = ((bmin)[1]+(bmax)[1])/2; } while(0)
#define G2V_BBOX_INSIDE(x, y, bmin, bmax, tol) ( (x) >= (bmin)[0]-tol && (x) <= (bmax)[0]+tol && (y) >= (bmin)[1]-tol && (y) <= (bmax)[1]+tol )
#define G2V_BBOX_OVERLAP(amin, amax, bmin, bmax, tol)	   \
  ( G2V_BBOX_INSIDE( (amin)[0], (amin)[1], bmin, bmax, tol ) || \
    G2V_BBOX_INSIDE( (amax)[0], (amax)[1], bmin, bmax, tol ) || \
    G2V_BBOX_INSIDE( (amin)[0], (amax)[1], bmin, bmax, tol ) || \
    G2V_BBOX_INSIDE( (amax)[0], (amin)[1], bmin, bmax, tol ) || \
    G2V_BBOX_INSIDE( (bmin)[0], (bmin)[1], amin, amax, tol ) )

#define G2V_ROTATION_ANGLE(v) \
  (((v)[0] >= 0) ? asin((v)[1]) : \
   (((v)[1] >= 0) ? M_PI - asin((v)[1]) : -M_PI - asin((v)[1])))	// -PI ~ +PI

#define G2V_COUT(a)    "(" << setprecision(4) << (a)[0] << " " << (a)[1] << ")"
#define G2V_LCOUT(a)   "(" << setprecision(4) << (a)[0] << " " << (a)[1] << " " << (a)[2] << ")"



#define G2M_SET(M, a00, a01,  a10, a11) \
                       do { (M)[0] = a00;  (M)[1] = a01; \
                            (M)[2] = a10;  (M)[3] = a11; } while(0)
#define G2M_SET_ID(M)  do { (M)[0] = 1;  (M)[1] = 0; \
                            (M)[2] = 0;  (M)[3] = 1; } while(0)
#define G2M_NORM(M)   (sqrt((M)[0]*(M)[0] + (M)[1]*(M)[1] + (M)[2]*(M)[2] + (M)[3]*(M)[3]))
#define G2M_COPY(R,A)  \
                       do { (R)[0] = (A)[0];  (R)[1] = (A)[1]; \
                            (R)[2] = (A)[2];  (R)[3] = (A)[3]; } while(0)
#define G2M_ADD(R, A, B) \
                       do { (R)[0] = (A)[0] + (B)[0];  (R)[1] = (A)[1] + (B)[1]; \
                            (R)[2] = (A)[2] + (B)[2];  (R)[3] = (A)[3] + (B)[3]; } while(0)
#define G2M_SUB(R, A, B) \
                       do { (R)[0] = (A)[0] - (B)[0];  (R)[1] = (A)[1] - (B)[1]; \
                            (R)[2] = (A)[2] - (B)[2];  (R)[3] = (A)[3] - (B)[3]; } while(0)
#define G2M_SCALED_ADD(R, A, s, B) \
                       do { (R)[0] = (A)[0] + (s)*(B)[0];  (R)[1] = (A)[1] + (s)*(B)[1]; \
                            (R)[2] = (A)[2] + (s)*(B)[2];  (R)[3] = (A)[3] + (s)*(B)[3]; } while(0)
#define G2M_TRANS(R,A)  \
                       do { (R)[0] = (A)[0];  (R)[1] = (A)[2]; \
                            (R)[2] = (A)[1];  (R)[3] = (A)[3]; } while(0)
#define G2M_MUL_VALUE(A,v)  \
                       do { (A)[0] *= (v);  (A)[1] *= (v); \
                            (A)[2] *= (v);  (A)[3] *= (v); } while(0)
#define G2M_MUL_MV(r, M, a) \
                       do { (r)[0] = (M)[0] * (a)[0] + (M)[1] * (a)[1]; \
                            (r)[1] = (M)[2] * (a)[0] + (M)[3] * (a)[1]; } while(0)
#define G2M_MUL_VM(r, M, a) \
                       do { (r)[0] = (M)[0] * (a)[0] + (M)[2] * (a)[1]; \
                            (r)[1] = (M)[1] * (a)[0] + (M)[3] * (a)[1]; } while(0)
#define G2M_MUL_VMV(a, M, b) \
                       ( ((M)[0] * (b)[0] + (M)[1] * (b)[1]) * (a)[0] + \
                         ((M)[2] * (b)[0] + (M)[3] * (b)[1]) * (a)[1] )
#define G2M_MUL_MM(R, A, B) \
                       do { (R)[0] = (A)[0] * (B)[0] + (A)[1] * (B)[2]; \
                            (R)[1] = (A)[0] * (B)[1] + (A)[1] * (B)[3]; \
                            (R)[2] = (A)[2] * (B)[0] + (A)[3] * (B)[2]; \
                            (R)[3] = (A)[2] * (B)[1] + (A)[3] * (B)[3]; } while(0)
#define G2M_MUL_MtM(R, A, B) \
                       do { (R)[0] = (A)[0] * (B)[0] + (A)[2] * (B)[2]; \
                            (R)[1] = (A)[0] * (B)[1] + (A)[2] * (B)[3]; \
                            (R)[2] = (A)[1] * (B)[0] + (A)[3] * (B)[2]; \
                            (R)[3] = (A)[1] * (B)[1] + (A)[3] * (B)[3]; } while(0)
#define G2M_MUL_MMt(R, A, B) \
                       do { (R)[0] = (A)[0] * (B)[0] + (A)[1] * (B)[1]; \
                            (R)[1] = (A)[0] * (B)[2] + (A)[1] * (B)[3]; \
                            (R)[2] = (A)[2] * (B)[0] + (A)[3] * (B)[1]; \
                            (R)[3] = (A)[2] * (B)[2] + (A)[3] * (B)[3]; } while(0)
#define G2M_DETERMINANT(A) (+ A[0] * A[3] - A[1] * A[2])
#define G2M_COUT(a)    "(" << setprecision(4) << (a)[0] << " " << (a)[1] << "  " << (a)[2] << " " << (a)[3] << ")"

#define G2M_ROTATION_ANGLE(R) \
  (((R)[0] >= 0) ? asin((R)[2]) : \
   (((R)[2] >= 0) ? M_PI - asin((R)[2]) : -M_PI - asin((R)[2])))	// -PI ~ +PI

#define G2M_XFORM_INVERSE(Rdst, Tdst,   Rsrc, Tsrc) \
                      do { G2M_TRANS( Rdst, Rsrc ); G2M_MUL_MV(Tdst, Rdst, Tsrc); G2V_MUL_VALUE(Tdst, -1); } while(0)
#define G2M_XFORM_MVV(c, R, a, T) \
                      do { G2M_MUL_MV(c, R, a); G2V_ADD( c, c, T ); } while(0)
#define G2M_XFORM_MERGE(Rac, Tac,   Rab, Tab,   Rbc, Tbc) \
                      do { G2M_MUL_MM(Rac, Rab, Rbc); G2M_XFORM_MVV(Tac, Rab, Tbc, Tab); } while(0)
#define G2M_XFORM_PRINTF(R, T) \
                      do { printf("  [ %5.2f %5.2f  %6.2f ] \n",(R)[0], (R)[1], (T)[0]); printf("  [ %5.2f %5.2f  %6.2f ] \n",(R)[2], (R)[3], (T)[1]); } while(0)

/* ================================================================ */
/*  3D                                                              */
/* ================================================================ */

#define G3V_SET(r,x,y,z)        do { (r)[0] = (x);  (r)[1] = (y);  (r)[2] = (z); } while(0)
#define G3V_COPY(d,s)           do { (d)[0] = (s)[0];  (d)[1] = (s)[1];  (d)[2] = (s)[2]; } while(0)

#define G3V_ADD(r,a,b)          do { (r)[0] = (a)[0] + (b)[0];  (r)[1] = (a)[1] + (b)[1];  (r)[2] = (a)[2] + (b)[2]; } while(0)
#define G3V_SUB(r,a,b)          do { (r)[0] = (a)[0] - (b)[0];  (r)[1] = (a)[1] - (b)[1];  (r)[2] = (a)[2] - (b)[2]; } while(0)
#define G3V_ADD_VALUES(r,a,b,c) do { (r)[0] += (a);  (r)[1] += (b);  (r)[2] += (c); } while(0)
#define G3V_SUB_VALUES(r,a,b,c) do { (r)[0] -= (a);  (r)[1] -= (b);  (r)[2] -= (c); } while(0)
#define G3V_MUL_VALUES(r,a,b,c) do { (r)[0] *= (a);  (r)[1] *= (b);  (r)[2] *= (c); } while(0)
#define G3V_DIV_VALUES(r,a,b,c) do { (r)[0] /= (a);  (r)[1] /= (b);  (r)[2] /= (c); } while(0)
#define G3V_ADD_VALUE(a,v)      do { (a)[0] += (v);  (a)[1] += (v);  (a)[2] += (v); } while(0)
#define G3V_SUB_VALUE(a,v)      do { (a)[0] -= (v);  (a)[1] -= (v);  (a)[2] -= (v); } while(0)
#define G3V_MUL_VALUE(a,v)      do { (a)[0] *= (v);  (a)[1] *= (v);  (a)[2] *= (v); } while(0)
#define G3V_DIV_VALUE(a,v)      do { (a)[0] /= (v);  (a)[1] /= (v);  (a)[2] /= (v); } while(0)
#define G3V_SCALED_ADD(r,a,s,b) do { (r)[0] = (a)[0] + (b)[0]*(s);  (r)[1] = (a)[1] + (b)[1]*(s);  (r)[2] = (a)[2] + (b)[2]*(s); } while(0)
#define G3V_WEIGHTED_ADD(r,w1,a,w2,b) do { (r)[0] = (a)[0]*(w1) + (b)[0]*(w2);  (r)[1] = (a)[1]*(w1) + (b)[1]*(w2);  (r)[2] = (a)[2]*(w1) + (b)[2]*(w2); } while(0)

#define G3V_DOT(a,b)            ((a)[0]*(b)[0] + (a)[1]*(b)[1] + (a)[2]*(b)[2])
#define G3V_CROSS(r,a,b)        do { (r)[0] = (a)[1]*(b)[2]-(a)[2]*(b)[1]; (r)[1] = (a)[2]*(b)[0]-(a)[0]*(b)[2]; (r)[2] = (a)[0]*(b)[1]-(a)[1]*(b)[0]; } while(0)
#define G3V_PLANE_POINT(p,a)    (p[0]*(a)[0] + p[1]*(a)[1] + p[2]*(a)[2] + p[3])

#define G3V_DISTANCE2(a,b)      (((b)[0]-(a)[0])*((b)[0]-(a)[0]) + ((b)[1]-(a)[1])*((b)[1]-(a)[1]) + ((b)[2]-(a)[2])*((b)[2]-(a)[2]))
#define G3V_DISTANCE(a,b)       (sqrt(((b)[0]-(a)[0])*((b)[0]-(a)[0]) + ((b)[1]-(a)[1])*((b)[1]-(a)[1]) + ((b)[2]-(a)[2])*((b)[2]-(a)[2])))
#define G3V_DIST_FABS(a,b)      (fabs((b)[0]-(a)[0]) + fabs((b)[1]-(a)[1]) + fabs((b)[2]-(a)[2]))
#define G3V_DIST_ABS(a,b)       (abs((b)[0]-(a)[0]) + abs((b)[1]-(a)[1]) + abs((b)[2]-(a)[2]))
#define G3V_LENGTH2(a)          ((a)[0]*(a)[0] + (a)[1]*(a)[1] + (a)[2]*(a)[2])
#define G3V_LENGTH(a)           (sqrt((a)[0]*(a)[0] + (a)[1]*(a)[1] + (a)[2]*(a)[2]))
#define G3V_NORMALIZE(a)        do { float  len = (float)G3V_LENGTH(a);  if (len > 0) G3V_DIV_VALUE( a, len ); } while(0)
#define G3V_NORMALIZE_WITH_TYPE(a,T_t)  do { T_t len = (T_t)G3V_LENGTH(a);  if (len > 0) G3V_DIV_VALUE( a, len ); } while(0)

#define G3V_EQUAL(a,b)          ((a)[0] == (b)[0] && (a)[1] == (b)[1] && (a)[2] == (b)[2])
#define G3V_AVERAGE2(r,a,b)     do { (r)[0] = ((a)[0] + (b)[0])/2;  (r)[1] = ((a)[1] + (b)[1])/2;  (r)[2] = ((a)[2] + (b)[2])/2; } while(0)
#define G3V_AVERAGE3(r,a,b,c)   do { (r)[0] = ((a)[0] + (b)[0] + (c)[0])/3;  (r)[1] = ((a)[1] + (b)[1] + (c)[1])/3;  (r)[2] = ((a)[2] + (b)[2] + (c)[2])/3; } while(0)
#define G3V_AVERAGE4(r,a,b,c,d) do { (r)[0] = ((a)[0] + (b)[0] + (c)[0] + (d)[0])/4;  (r)[1] = ((a)[1] + (b)[1] + (c)[1] + (d)[1])/4;  (r)[2] = ((a)[2] + (b)[2] + (c)[2] + (d)[2])/4; } while(0)
#define G3V_MINV(a)             ( (a)[0]<(a)[1] ? ((a)[0]<(a)[2] ? (a)[0] : (a)[2]) : ((a)[1]<(a)[2] ? (a)[1] : (a)[2]) )
#define G3V_MAXV(a)             ( (a)[0]>(a)[1] ? ((a)[0]>(a)[2] ? (a)[0] : (a)[2]) : ((a)[1]>(a)[2] ? (a)[1] : (a)[2]) )
#define G3V_MAKE_ABOVE(a,l)     do { if ((a)[0]<(l)[0]) (a)[0]=(l)[0]; if ((a)[1]<(l)[1]) (a)[1]=(l)[1]; if ((a)[2]<(l)[2]) (a)[2]=(l)[2]; } while(0)
#define G3V_MAKE_BELOW(a,u)     do { if ((a)[0]>(u)[0]) (a)[0]=(u)[0]; if ((a)[1]>(u)[1]) (a)[1]=(u)[1]; if ((a)[2]>(u)[2]) (a)[2]=(u)[2]; } while(0)
#define G3V_MAKE_BETWEEN(a,l,u) do { G3V_MAKE_ABOVE(a,l); G3V_MAKE_BELOW(a,u); } while(0)
#define G3V_TEST_ABOVE(a,u)     ( (a)[0]>=(u)[0] && (a)[1]>=(u)[1] && (a)[2]>=(u)[2] )
#define G3V_TEST_BELOW(a,l)     ( (a)[0]<=(l)[0] && (a)[1]<=(l)[1] && (a)[2]<=(l)[2] )
#define G3V_SWAP(a,b,t)         do { G3V_COPY(t, a); G3V_COPY(a, b); G3V_COPY(b, t); } while(0)
#define G3V_BBOX_INITIALIZE(a, bmin, bmax) do { G3V_COPY(bmin,a); G3V_COPY(bmax,a); } while(0)
#define G3V_BBOX_UPDATE(a, amin, amax) do { if ((a)[0]<(amin)[0]) (amin)[0] = (a)[0];  if ((a)[1]<(amin)[1]) (amin)[1] = (a)[1];  if ((a)[2]<(amin)[2]) (amin)[2] = (a)[2];  if ((a)[0]>(amax)[0]) (amax)[0] = (a)[0];  if ((a)[1]>(amax)[1]) (amax)[1] = (a)[1];  if ((a)[2]>(amax)[2]) (amax)[2] = (a)[2]; } while(0)

#define G3V_PRINTF(v) do { printf("  %5.2f %5.2f %5.2f \n",(v)[0], (v)[1], (v)[2]); } while(0)
#define G3V_COUT(a)   "(" << setprecision(4) << (a)[0] << " " << (a)[1] << " " << (a)[2] << ")"
#define G3V_PCOUT(a)  "(" << setprecision(4) << (a)[0] << " " << (a)[1] << " " << (a)[2] << " " << (a)[3] << ")"


#define G3M_SET(m, a00,a01,a02,  a10,a11,a12,  a20,a21,a22) \
                       do { (m)[0] = a00;  (m)[1] = a01;  (m)[2] = a02; \
                            (m)[3] = a10;  (m)[4] = a11;  (m)[5] = a12; \
                            (m)[6] = a20;  (m)[7] = a21;  (m)[8] = a22; } while(0)
#define G3M_SET_ID(m)  do { (m)[0] = 1;  (m)[1] = 0;  (m)[2] = 0; \
                            (m)[3] = 0;  (m)[4] = 1;  (m)[5] = 0; \
                            (m)[6] = 0;  (m)[7] = 0;  (m)[8] = 1; } while(0)
#define G3M_NORM(M)   (sqrt((M)[0]*(M)[0] + (M)[1]*(M)[1] + (M)[2]*(M)[2] + (M)[3]*(M)[3] + (M)[4]*(M)[4] + (M)[5]*(M)[5] + (M)[6]*(M)[6] + (M)[7]*(M)[7] + (M)[8]*(M)[8]))
#define G3M_COPY(R, A) \
                       do { (R)[0] = (A)[0];  (R)[1] = (A)[1];  (R)[2] = (A)[2]; \
                            (R)[3] = (A)[3];  (R)[4] = (A)[4];  (R)[5] = (A)[5]; \
                            (R)[6] = (A)[6];  (R)[7] = (A)[7];  (R)[8] = (A)[8]; } while(0)
#define G3M_ADD(R, A, B) \
                       do { (R)[0] = (A)[0] + (B)[0];  (R)[1] = (A)[1] + (B)[1];  (R)[2] = (A)[2] + (B)[2]; \
                            (R)[3] = (A)[3] + (B)[3];  (R)[4] = (A)[4] + (B)[4];  (R)[5] = (A)[5] + (B)[5]; \
                            (R)[6] = (A)[6] + (B)[6];  (R)[7] = (A)[7] + (B)[7];  (R)[8] = (A)[8] + (B)[8]; } while(0)
#define G3M_SUB(R, A, B) \
                       do { (R)[0] = (A)[0] - (B)[0];  (R)[1] = (A)[1] - (B)[1];  (R)[2] = (A)[2] - (B)[2]; \
                            (R)[3] = (A)[3] - (B)[3];  (R)[4] = (A)[4] - (B)[4];  (R)[5] = (A)[5] - (B)[5]; \
                            (R)[6] = (A)[6] - (B)[6];  (R)[7] = (A)[7] - (B)[7];  (R)[8] = (A)[8] - (B)[8]; } while(0)
#define G3M_SCALED_ADD(R, A, s, B) \
                       do { (R)[0] = (A)[0] + (s)*(B)[0];  (R)[1] = (A)[1] + (s)*(B)[1];  (R)[2] = (A)[2] + (s)*(B)[2]; \
                            (R)[3] = (A)[3] + (s)*(B)[3];  (R)[4] = (A)[4] + (s)*(B)[4];  (R)[5] = (A)[5] + (s)*(B)[5]; \
                            (R)[6] = (A)[6] + (s)*(B)[6];  (R)[7] = (A)[7] + (s)*(B)[7];  (R)[8] = (A)[8] + (s)*(B)[8]; } while(0)
#define G3M_TRANS(R, A) \
                       do { (R)[0] = (A)[0];  (R)[1] = (A)[3];  (R)[2] = (A)[6]; \
                            (R)[3] = (A)[1];  (R)[4] = (A)[4];  (R)[5] = (A)[7]; \
                            (R)[6] = (A)[2];  (R)[7] = (A)[5];  (R)[8] = (A)[8]; } while(0)
#define G3M_MUL_VALUE(R, v) \
                       do { (R)[0] *= (v);  (R)[1] *= (v);  (R)[2] *= (v); \
                            (R)[3] *= (v);  (R)[4] *= (v);  (R)[5] *= (v); \
                            (R)[6] *= (v);  (R)[7] *= (v);  (R)[8] *= (v); } while(0)
#define G3M_MUL_MV(r, M, a) \
                       do { (r)[0] = (M)[0] * (a)[0] + (M)[1] * (a)[1] + (M)[2] * (a)[2]; \
                            (r)[1] = (M)[3] * (a)[0] + (M)[4] * (a)[1] + (M)[5] * (a)[2]; \
                            (r)[2] = (M)[6] * (a)[0] + (M)[7] * (a)[1] + (M)[8] * (a)[2]; } while(0)
#define G3M_MUL_VM(r, a, M) \
                       do { (r)[0] = (M)[0] * (a)[0] + (M)[3] * (a)[1] + (M)[6] * (a)[2]; \
                            (r)[1] = (M)[1] * (a)[0] + (M)[4] * (a)[1] + (M)[7] * (a)[2]; \
                            (r)[2] = (M)[2] * (a)[0] + (M)[5] * (a)[1] + (M)[8] * (a)[2]; } while(0)
#define G3M_MUL_VVt(M, a, b) \
                       do { (M)[0] = (a)[0] * (b)[0];  (M)[1] = (a)[0] * (b)[1];  (M)[2] = (a)[0] * (b)[2]; \
                            (M)[3] = (a)[1] * (b)[0];  (M)[4] = (a)[1] * (b)[1];  (M)[5] = (a)[1] * (b)[2]; \
                            (M)[6] = (a)[2] * (b)[0];  (M)[7] = (a)[2] * (b)[1];  (M)[8] = (a)[2] * (b)[2]; } while(0)
#define G3M_MUL_VMV(a, M, b) \
                       ( ((M)[0] * (b)[0] + (M)[1] * (b)[1] + (M)[2] * (b)[2]) * (a)[0] + \
                         ((M)[3] * (b)[0] + (M)[4] * (b)[1] + (M)[5] * (b)[2]) * (a)[1] + \
                         ((M)[6] * (b)[0] + (M)[7] * (b)[1] + (M)[8] * (b)[2]) * (a)[2] )
#define G3M_MUL_MM(R, A, B) \
                       do { (R)[0] = (A)[0] * (B)[0] + (A)[1] * (B)[3] + (A)[2] * (B)[6]; \
                            (R)[1] = (A)[0] * (B)[1] + (A)[1] * (B)[4] + (A)[2] * (B)[7]; \
                            (R)[2] = (A)[0] * (B)[2] + (A)[1] * (B)[5] + (A)[2] * (B)[8]; \
                            (R)[3] = (A)[3] * (B)[0] + (A)[4] * (B)[3] + (A)[5] * (B)[6]; \
                            (R)[4] = (A)[3] * (B)[1] + (A)[4] * (B)[4] + (A)[5] * (B)[7]; \
                            (R)[5] = (A)[3] * (B)[2] + (A)[4] * (B)[5] + (A)[5] * (B)[8]; \
                            (R)[6] = (A)[6] * (B)[0] + (A)[7] * (B)[3] + (A)[8] * (B)[6]; \
                            (R)[7] = (A)[6] * (B)[1] + (A)[7] * (B)[4] + (A)[8] * (B)[7]; \
                            (R)[8] = (A)[6] * (B)[2] + (A)[7] * (B)[5] + (A)[8] * (B)[8]; } while(0)
#define G3M_MUL_MMt(R, A, B) \
                       do { (R)[0] = (A)[0] * (B)[0] + (A)[1] * (B)[1] + (A)[2] * (B)[2]; \
                            (R)[1] = (A)[0] * (B)[3] + (A)[1] * (B)[4] + (A)[2] * (B)[5]; \
                            (R)[2] = (A)[0] * (B)[6] + (A)[1] * (B)[7] + (A)[2] * (B)[8]; \
                            (R)[3] = (A)[3] * (B)[0] + (A)[4] * (B)[1] + (A)[5] * (B)[2]; \
                            (R)[4] = (A)[3] * (B)[3] + (A)[4] * (B)[4] + (A)[5] * (B)[5]; \
                            (R)[5] = (A)[3] * (B)[6] + (A)[4] * (B)[7] + (A)[5] * (B)[8]; \
                            (R)[6] = (A)[6] * (B)[0] + (A)[7] * (B)[1] + (A)[8] * (B)[2]; \
                            (R)[7] = (A)[6] * (B)[3] + (A)[7] * (B)[4] + (A)[8] * (B)[5]; \
                            (R)[8] = (A)[6] * (B)[6] + (A)[7] * (B)[7] + (A)[8] * (B)[8]; } while(0)
#define G3M_MUL_MtM(R, A, B) \
                       do { (R)[0] = (A)[0] * (B)[0] + (A)[3] * (B)[3] + (A)[6] * (B)[6]; \
                            (R)[1] = (A)[0] * (B)[1] + (A)[3] * (B)[4] + (A)[6] * (B)[7]; \
                            (R)[2] = (A)[0] * (B)[2] + (A)[3] * (B)[5] + (A)[6] * (B)[8]; \
                            (R)[3] = (A)[1] * (B)[0] + (A)[4] * (B)[3] + (A)[7] * (B)[6]; \
                            (R)[4] = (A)[1] * (B)[1] + (A)[4] * (B)[4] + (A)[7] * (B)[7]; \
                            (R)[5] = (A)[1] * (B)[2] + (A)[4] * (B)[5] + (A)[7] * (B)[8]; \
                            (R)[6] = (A)[2] * (B)[0] + (A)[5] * (B)[3] + (A)[8] * (B)[6]; \
                            (R)[7] = (A)[2] * (B)[1] + (A)[5] * (B)[4] + (A)[8] * (B)[7]; \
                            (R)[8] = (A)[2] * (B)[2] + (A)[5] * (B)[5] + (A)[8] * (B)[8]; } while(0)
#define G3M_PRINTF(R) do { printf("  %5.2f %5.2f %5.2f \n",(R)[0], (R)[1], (R)[2]); printf("  %5.2f %5.2f %5.2f \n",(R)[3], (R)[4], (R)[5]); printf("  %5.2f %5.2f %5.2f \n",(R)[6], (R)[7], (R)[8]);  } while(0)
#define G3M_MUL_MMM(R, A, B, C) \
                      do { double T[9];  G3M_MUL_MM( T, B, C );  G3M_MUL_MM( R, A, T ); } while(0)
#define G3M_XFORM_INVERSE(Rdst, Tdst,   Rsrc, Tsrc) \
                      do { G3M_TRANS( Rdst, Rsrc ); G3M_MUL_MV(Tdst, Rdst, Tsrc); G3V_MUL_VALUE(Tdst, -1); } while(0)
#define G3M_XFORM_MVV(c, R, a, T) \
                      do { G3M_MUL_MV(c, R, a); G3V_ADD( c, c, T ); } while(0)
#define G3M_XFORM_MERGE(Rac, Tac,   Rab, Tab,   Rbc, Tbc) \
                      do { G3M_MUL_MM(Rac, Rab, Rbc); G3M_XFORM_MVV(Tac, Rab, Tbc, Tab); } while(0)
#define G3M_XFORM_PRINTF(R, T) \
                      do { printf("  [ %7.4f %7.4f %7.4f  %8.4f ] \n",(R)[0], (R)[1], (R)[2], (T)[0]); printf("  [ %7.4f %7.4f %7.4f  %8.4f ] \n",(R)[3], (R)[4], (R)[5], (T)[1]); printf("  [ %7.4f %7.4f %7.4f  %8.4f ] \n",(R)[6], (R)[7], (R)[8], (T)[2]);  } while(0)

#define G3M_DETERMINANT(A) (+ A[0] * A[4] * A[8] + A[1] * A[5] * A[6] + A[2] * A[3] * A[7] - A[0] * A[5] * A[7] - A[1] * A[3] * A[8] - A[2] * A[4] * A[6])
#define G3M_COUT(a)   "(" << setprecision(4) << (a)[0] << " " << (a)[1] << " " << (a)[2] << "  " << (a)[3] << " " << (a)[4] << " " << (a)[5] << "  " << (a)[6] << " " << (a)[7] << " " << (a)[8] << ")"


/* ================================================================ */
/*  4D                                                              */
/* ================================================================ */

#define G4V_SET(r,x,y,z,w)      do { (r)[0] = (x);  (r)[1] = (y);  (r)[2] = (z);  (r)[3] = (w); } while(0)
#define G4V_COPY(d,s)           do { (d)[0] = (s)[0];  (d)[1] = (s)[1];  (d)[2] = (s)[2];  (d)[3] = (s)[3]; } while(0)
#define G4V_ADD(r,a,b)          do { (r)[0] = (a)[0] + (b)[0];  (r)[1] = (a)[1] + (b)[1];  (r)[2] = (a)[2] + (b)[2];  (r)[3] = (a)[3] + (b)[3]; } while(0)
#define G4V_SUB(r,a,b)          do { (r)[0] = (a)[0] - (b)[0];  (r)[1] = (a)[1] - (b)[1];  (r)[2] = (a)[2] - (b)[2];  (r)[3] = (a)[3] - (b)[3]; } while(0)
#define G4V_ADD_VALUE(a,v)      do { (a)[0] += (v);  (a)[1] += (v);  (a)[2] += (v);  (a)[3] += (v); } while(0)
#define G4V_SUB_VALUE(a,v)      do { (a)[0] -= (v);  (a)[1] -= (v);  (a)[2] -= (v);  (a)[3] -= (v); } while(0)
#define G4V_MUL_VALUE(a,v)      do { (a)[0] *= (v);  (a)[1] *= (v);  (a)[2] *= (v);  (a)[3] *= (v); } while(0)
#define G4V_DIV_VALUE(a,v)      do { (a)[0] /= (v);  (a)[1] /= (v);  (a)[2] /= (v);  (a)[3] /= (v); } while(0)
#define G4V_SCALED_ADD(r,a,s,b) do { (r)[0] = (a)[0] + (b)[0]*(s);  (r)[1] = (a)[1] + (b)[1]*(s);  (r)[2] = (a)[2] + (b)[2]*(s);  (r)[3] = (a)[3] + (b)[3]*(s); } while(0)
#define G4V_WEIGHTED_ADD(r,sa,a,sb,b)  do { (r)[0] = (sa)*(a)[0]+(sb)*(b)[0]; (r)[1] = (sa)*(a)[1]+(sb)*(b)[1]; (r)[2] = (sa)*(a)[2]+(sb)*(b)[2]; (r)[3] = (sa)*(a)[3]+(sb)*(b)[3]; } while(0)
#define G4V_SWAP(a,b,t)         do { G4V_COPY(t, a); G4V_COPY(a, b); G4V_COPY(b, t); } while(0)
#define G4V_COUT(a)   "(" << setprecision(4) << (a)[0] << " " << (a)[1] << " " << (a)[2] << " " << (a)[3] << ")"

#define G4M_SET(m, a00,a01,a02,a03, a10,a11,a12,a13, a20,a21,a22,a23, a30,a31,a32,a33) \
                       do { (m)[0] = a00;  (m)[1] = a01;  (m)[2] = a02;  (m)[3] = a03;  \
                            (m)[4] = a10;  (m)[5] = a11;  (m)[6] = a12;  (m)[7] = a13;  \
                            (m)[8] = a20;  (m)[9] = a21;  (m)[10] = a22; (m)[11] = a23; \
                            (m)[12] = a30; (m)[13] = a31; (m)[14] = a32; (m)[15] = a33; } while(0)


/* ================================================================ */
/*                                                                  */
/* ================================================================ */

#define G5V_COPY(p,s)  do { (p)[0] = (s)[0]; (p)[1] = (s)[1]; (p)[2] = (s)[2]; (p)[3] = (s)[3]; (p)[4] = (s)[4]; } while(0)
#define G6V_COPY(p,s)  do { (p)[0] = (s)[0]; (p)[1] = (s)[1]; (p)[2] = (s)[2]; (p)[3] = (s)[3]; (p)[4] = (s)[4]; (p)[5] = (s)[5]; } while(0)
#define G5V_SET(p,v0,v1,v2,v3,v4)    do { (p)[0] = v0; (p)[1] = v1; (p)[2] = v2; (p)[3] = v3; (p)[4] = v4; } while(0)
#define G6V_SET(p,v0,v1,v2,v3,v4,v5) do { (p)[0] = v0; (p)[1] = v1; (p)[2] = v2; (p)[3] = v3; (p)[4] = v4; (p)[5] = v5; } while(0)
#define G5V_COUT(a)   "(" << setprecision(4) << (a)[0] << " " << (a)[1] << " " << (a)[2] << " " << (a)[3] << " " << (a)[4] << ")"
#define G6V_COUT(a)   "(" << setprecision(4) << (a)[0] << " " << (a)[1] << " " << (a)[2] << " " << (a)[3] << " " << (a)[4] << " " << (a)[5] << ")"

#define G7V_SET(p,v0,v1,v2,v3,v4,v5,v6)    do { (p)[0] = v0; (p)[1] = v1; (p)[2] = v2; (p)[3] = v3; (p)[4] = v4; (p)[5] = v5; (p)[6] = v6; } while(0)
#define G8V_SET(p,v0,v1,v2,v3,v4,v5,v6,v7) do { (p)[0] = v0; (p)[1] = v1; (p)[2] = v2; (p)[3] = v3; (p)[4] = v4; (p)[5] = v5; (p)[6] = v6; (p)[7] = v7; } while(0)

#endif  // VECTOR_MATRIX_MACROS_H


/* ================================================================ */
/*                                                                  */
/* ================================================================ */

#define PRINT_ARRAY1(array,size,cmmt,format,Return) \
do { \
  printf("%s ", (cmmt!=NULL ? cmmt : "ARRY")); \
  char *fmt2 = (char*)format; if (fmt2==NULL) fmt2 = (char*)"%12g ";	\
  for (int k=0; k<size; k++) printf(fmt2, (array)[k]); \
  if (Return) printf("\n"); \
} while(0)	// print 1D array in a single line

#define PRINT_ARRAY2(array,w,h,cmmt,format) \
do { \
  int i, nsp=0, len=(cmmt!=NULL?strlen(cmmt):0);  char spaces[41]; \
  if (cmmt!=NULL) { for (i=0; i<40&&i<len; i++) if (cmmt[i]==' ') nsp++; } \
  for (i=0; i<nsp; i++) spaces[i] = ' ';  spaces[i] = '\0'; \
  printf("%s\n", (cmmt!=NULL ? cmmt : "ARRY")); \
  char *fmt = (char*)format; if (fmt==NULL) fmt = (char*)"%12g ";	\
  for (int j=0; j < h; j++) { \
    printf("%s  ", spaces); \
    for (int i=0; i < w; i++) printf(fmt, (array)[j*(w)+i]); \
    printf("\n"); \
  } \
} while(0)	// print 2D array

#define PRINT_ARRAY3(array,w,h,s,cmmt,format) \
do { \
  int i, nsp=0, len=(cmmt!=NULL?strlen(cmmt):0);  char spaces[41]; \
  if (cmmt!=NULL) { for (i=0; i<40&&i<len; i++) if (cmmt[i]==' ') nsp++; } \
  for (i=0; i<nsp; i++) spaces[i] = ' ';  spaces[i] = '\0'; \
  printf("%s\n", (cmmt!=NULL ? cmmt : "ARRY")); \
  char *fmt = (char*)format; if (fmt==NULL) fmt = (char*)"%12g ";	\
  for (int j=0; j < h; j++) { \
    printf("%s  ", spaces); \
    for (int i=0; i < w; i++)  \
      PRINT_ARRAY1( array+(j*(w)+i)*(s), s, "", fmt, false ); \
    printf("\n"); \
  } \
} while(0)	// print 3D array

#define PRINT_ARRAY2_RANGE(array,w,h,cmmt,format,x,y,ww,hh) \
do { \
  int i, nsp=0, len=(cmmt!=NULL?strlen(cmmt):0);  char spaces[41]; \
  if (cmmt!=NULL) { for (i=0; i<40&&i<len; i++) if (cmmt[i]==' ') nsp++; } \
  for (i=0; i<nsp; i++) spaces[i] = ' ';  spaces[i] = '\0'; \
  printf("%s\n", (cmmt!=NULL ? cmmt : "ARRY")); \
  char *fmt = (char*)format; if (fmt==NULL) fmt = (char*)"%12g ";	\
  for (int j=0; j < h; j++) { \
    if (j < y) continue; if (j >= y+hh) break; \
    printf("%s  ", spaces); \
    for (int i=0; i < w; i++) \
      if (i >= x && i < x+ww) printf(fmt, (array)[j*(w)+i]);	\
    printf("\n"); \
  } \
} while(0)

