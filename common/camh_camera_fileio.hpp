
// 
// CAMH::CameraFileIO<> class template
//
// Jaeil Choi
// last modified in Oct, 2004
//
// -------------------------------------------------------------------
// This file is a part of CAMH (Camera classes in Header) library
//   - CAMH::Camera<>		in camh_camera.hpp
//   - CAMH::CameraFileIO<>	in camh_camera_fileio.hpp
//   - CAMH::CameraRenderer<>	in camh_camera_renderer.hpp
//   - CAMH::Stereo<>		in camh_stereo.hpp
// 
// This class is designed for
//   - reading camera calibration parameters from  various file format
//   - writing camera calibration parameters using various file format
//

#ifndef CAMH_FILEIO_HPP
#define CAMH_FILEIO_HPP

#include <iostream>
#include "vm_macros.h"
#include "util_word_breaker.hpp"
#include "guih_config.hpp"
#include "camh_camera.hpp"

namespace CAMH {

template <class T_t>
class CameraFileIO 
{
 public:
  CameraFileIO()  {}
  ~CameraFileIO() {}
  
 public:
  bool writeFile(Camera<T_t> *c, FILE *fp, bool write_motion=true) {
    if (fp == NULL) return false;
    fprintf(fp, "CameraIntrinsic : [K]  (%d x %d)\n", c->wh[0], c->wh[1]);
    fprintf(fp, "\t%-9g %-9g %-9g\n\t%-9g %-9g %-9g\n\t%-9g %-9g %-9g\n",
	    c->fc[0],  c->afc,    c->cc[0], 
	    0.0,      c->fc[1],   c->cc[1],
	    0.0,         0.0,       1.0  );
    if (c->dt_type == 'I') {	// invertible distortion model
      fprintf(fp, "CameraDistortion: %g\n", c->dt[0]);
    } else {		// standard distortion model
      fprintf(fp, "CameraDistortion: %g %g %g %g %g\n", 
	      c->dt[0], c->dt[1], c->dt[2], c->dt[3], c->dt[4]);
    }
    fprintf(fp, "CameraExtrinsic : [Rcw|Tcw]\n");
    fprintf(fp, "\t%-12g %-12g %-12g %-12g \n\t%-12g %-12g %-12g %-12g \n\t%-12g %-12g %-12g %-12g \n",
	    c->Rcw[0], c->Rcw[1], c->Rcw[2], c->Tcw[0], 
	    c->Rcw[3], c->Rcw[4], c->Rcw[5], c->Tcw[1], 
	    c->Rcw[6], c->Rcw[7], c->Rcw[8], c->Tcw[2]);
    if (write_motion && c->nframes > 1) {
      fprintf(fp, "CameraMotion [ %d frames * (Rcw[9]+Tcw[3]) ]\n", c->nframes);
      for (int fidx = 0; fidx < c->nframes; fidx++) {
	T_t *R = c->motion + fidx * (9 + 3);
	T_t *T = R + 9;
	fprintf(fp, "\t%-12g %-12g %-12g  %-12g %-12g %-12g  %-12g %-12g %-12g   %-12g %-12g %-12g\n",
		R[0], R[1], R[2], R[3], R[4], R[5], R[6], R[7], R[8], T[0], T[1], T[2]);
      }
    }
    return true;
  }
  
  bool writeFile(Camera<T_t> *c, char *filename) {
    FILE *fp = fopen(filename, "w+");
    if (fp == NULL) return false;
    bool ret = writeFile(c, fp, true);
    fclose(fp);
    return ret;
  }
  
  Camera<T_t>* readFile(Camera<T_t> *c, FILE *fp, bool read_motion) {
    // read generic calibration file format
    if (fp == NULL) return NULL;
    if (c==NULL) c = new Camera<T_t>;   
    UTIL::WordBreaker wb;
    char   line[1024], word[80], *rest;
    float  fc[2], afc, cc[2], dum, dt[5], R[9], T[3];
    c->wh[0] = c->wh[1] = 0;
    // read intrinsic parameters
    while (fgets(line, 1024, fp)) if (!wb.isEmpty(line, '#')) break;
    sscanf(line, "%s\n", word);
    if (strncmp(word, "CameraIntrinsic : [K]", 14)==0) {
      sscanf(line+21, " (%d x %d)", c->wh+0, c->wh+1);
      fscanf(fp, " %f %f %f\n", fc+0, &afc, cc+0);  
      fscanf(fp, " %f %f %f\n", &dum,  fc+1,  cc+1);
      fscanf(fp, " %f %f %f\n", &dum,  &dum,  &dum);
      G2V_COPY(c->fc, fc);  G2V_COPY(c->cc, cc);
      c->afc = afc;
    } else {
      fprintf(stderr, "Error (CameraFileIO::readFile): invalid intrinsic parameters\n");
      c->clear(); return NULL;
    }
    // read distortion coefficients
    while (fgets(line, 1024, fp)) if (!wb.isEmpty(line, '#')) break;
    rest = wb.parseFirst( line, word );
    if (strncmp(word, "CameraDistortion:", 14)==0) {
      int n = sscanf(rest, " %f %f %f %f %f", dt+0, dt+1, dt+2, dt+3, dt+4);
      G5V_SET( c->dt, dt[0], dt[1], dt[2], dt[3], dt[4] );
      c->dt_type = (n == 1 ? 'I' : 'S');
    } else {
      fprintf(stderr, "Error (CameraFileIO::readFile): invalid distortion coefficients\n");
      c->clear(); return NULL;
    }
    // read extrinsic parameters
    c->createMotionFrames(1);
    while (fgets(line, 1024, fp)) if (!wb.isEmpty(line, '#')) break;
    sscanf(line, "%s\n", word);
    if (strncmp(word, "CameraExtrinsic : [Rcw|Tcw]", 14)==0) {
      fscanf(fp, " %f %f %f %f\n", R+0, R+1, R+2, T+0);
      fscanf(fp, " %f %f %f %f\n", R+3, R+4, R+5, T+1);
      fscanf(fp, " %f %f %f %f\n", R+6, R+7, R+8, T+2);
      G3M_SET( c->Rcw, R[0], R[1], R[2], R[3], R[4], R[5], R[6], R[7], R[8] );
      G3V_SET( c->Tcw, T[0], T[1], T[2] );
    } else {
      fprintf(stderr, "Error (CameraFileIO::readFile): invalid extrinsic parameters\n");
      c->clear(); return NULL;
    }
    c->updatePmatUsingParam();
    // read camera motion, if requested
    if (read_motion) {
      while (fgets(line, 1024, fp)) if (!wb.isEmpty(line, '#')) break;
      int count;
      if (sscanf(line, "CameraMotion [ %d frames * (R[9]+T[3]) ]\n", &count) == 1) {
	// note that the extrinsic parameters read above are lost
	c->createMotionFrames(count);
	for (int fidx = 0; fidx < c->nframes; fidx++) {
	  T_t *cR = c->getR(fidx);
	  T_t *cT = c->getT(fidx);
	  fscanf(fp, "%f %f %f  %f %f %f  %f %f %f  %f %f %f\n",
		 R+0, R+1, R+2, R+3, R+4, R+5, R+6, R+7, R+8, T+0, T+1, T+2);
	  G3M_SET( cR,  R[0], R[1], R[2], R[3], R[4], R[5], R[6], R[7], R[8] );
	  G3V_SET( cT,  T[0], T[1], T[2] );
	}
      }
    }
    return c;
  }
  
  Camera<T_t>* readFile(Camera<T_t> *c, char *filename) {
    // read generic calibration file format
    FILE *fp = fopen(filename, "r");
    if (fp == NULL) {
      //cerr << "Error (CameraFileIO::readFile): failed to open '" << filename << "'" << endl;
      return NULL;
    }
    Camera<T_t>* ret = readFile(c, fp, true);
    fclose(fp);
    return ret;
  }
  
  Camera<T_t>* readCameraConfigFile(Camera<T_t> *c, char *filename, char *section) {
    char  res[80], in0[200], in1[200], ex0[200], ex1[200], ex2[200], dst[200];
    GUIH::Config cfg;
    cfg.set( section, "CameraResolution", CFG_STRING, res );
    cfg.set( section, "CameraIntrinsic0", CFG_STRING, in0 );
    cfg.set( section, "CameraIntrinsic1", CFG_STRING, in1 );
    cfg.set( section, "CameraExtrinsic0", CFG_STRING, ex0 );
    cfg.set( section, "CameraExtrinsic1", CFG_STRING, ex1 );
    cfg.set( section, "CameraExtrinsic2", CFG_STRING, ex2 );
    cfg.set( section, "CameraDistortion", CFG_STRING, dst );
    if (!cfg.process( filename, section )) {
      fprintf(stderr, "Error (CameraFileIO::readCameraConfigFile): cannot read camera parameters\n");
      return NULL;
    }
    if (!cfg.processed(section, "CameraIntrinsic0") && 
	!cfg.processed(section, "CameraExtrinsic0")) return false;
    float fc[2], cc[2], afc, R[9], T[3], dt[5], dummy;
    if (c == NULL) c = new Camera<T_t>;
    sscanf( res, "%d x %d", c->wh+0, c->wh+1 );
    sscanf( in0, "%f %f %f", fc+0,   &afc, cc+0 );
    sscanf( in1, "%f %f %f", &dummy, fc+1,   cc+1 );
    G2V_COPY( c->fc, fc );  G2V_COPY( c->cc, cc );   c->afc = afc;
    sscanf( ex0, "%f %f %f %f", R+0, R+1, R+2, T+0 );
    sscanf( ex1, "%f %f %f %f", R+3, R+4, R+5, T+1 );
    sscanf( ex2, "%f %f %f %f", R+6, R+7, R+8, T+2 );
    G3M_SET( c->Rcw,  R[0], R[1], R[2], R[3], R[4], R[5], R[6], R[7], R[8] );
    G3V_SET( c->Tcw,  T[0], T[1], T[2] );
    int n = sscanf( dst, "%f %f %f %f %f", dt+0, dt+1, dt+2, dt+3, dt+4);
    G5V_SET( c->dt, dt[0], dt[1], dt[2], dt[3], dt[4] );
    c->dt_type = (n > 1 ? 'S' : 'I');
  
    c->updatePmatUsingParam();
    return c;
  }
  
  Camera<T_t>* readGabFile(Camera<T_t> *c, char *filename) {
    // read Gabriel Brostow's calibration file.
    // 'filename' is something like 'calib_gui.calib'
    char   buffer[256];
    int    i_cnt = 0, t_cnt = 0, r_cnt = 0;
    float  K[9], tmp, dt[4], rt[3];
    FILE *fp = fopen(filename, "r");
    if (fp == (FILE*)NULL) {
      fprintf(stderr, "Error (CameraFileIO::readGabFile): cannot open file '%s'\n", filename);
      return NULL;
    }
    if (c == NULL) c = new Camera<T_t>;
    while (fgets(buffer, 256, fp)) {
      if        (strncmp(buffer, "Intrin ", 7) == 0) {
	if (i_cnt >= 3) { fprintf(stderr, "Error (readGabFile): invalid intrinsic\n"); break; }
	sscanf(buffer + 7, "%f %f %f", K+(3*i_cnt+0), K+(3*i_cnt+1), K+(3*i_cnt+2) );
	i_cnt++;
      } else if (strncmp(buffer, "Trans_mm ", 9) == 0) {
	if (t_cnt >= 3) { fprintf(stderr, "Error (readGabFile): invalid translation\n"); break; }
	sscanf( buffer + 9, "%f", &tmp); c->Tcw[t_cnt] = tmp * 0.001f;
	t_cnt++;
      } else if (strncmp(buffer, "Rot ", 4) == 0) {
	if (r_cnt >= 3) { fprintf(stderr, "Error (readGabFile): invalid rotation\n"); break; }
	sscanf( buffer + 4, "%f %f %f", rt+0, rt+1, rt+3); 
	c->Rcw[3*r_cnt+0] = rt[0];  c->Rcw[3*r_cnt+1] = rt[1];  c->Rcw[3*r_cnt+2] = rt[2];
	r_cnt++;
      } else if (strncmp(buffer, "Distort ", 8) == 0) {
	sscanf( buffer + 8, "%f %f %f %f", dt+0, dt+1, dt+2, dt+3); 
	G4V_COPY( c->dt, dt );
      }
      buffer[0] = '\0';
    }
    c->fc[0] = K[0];	// focal length in pixels (X)
    c->fc[1] = K[4];	// focal length in pixels (X)
    c->afc   = K[1];	// skew coefficient
    c->cc[0] = K[2];	// principal point X coordinate
    c->cc[1] = K[5];	// principal point Y coordinate
    fclose(fp);
    if (buffer[0] != '\0' || i_cnt != 3 || t_cnt != 3 || r_cnt != 3)
      return NULL;
    c->updatePmatUsingParam();
    return c;
  }
  
  Camera<T_t>* readMCSCPmatFile(Camera<T_t> *c, char *filename) {
    // read the MCSC(Multi-Camera Self-Calibration) projection matrix file.
    // Usually they are names as "camera%d.Pmat.dat"
    FILE *fp = fopen(filename, "r");
    if (fp == NULL) {
      fprintf(stderr, "Error (CameraFileIO::readMCSCPmatFile): cannot open file '%s'\n", filename);
      return NULL;
    }
    // read the projection matrix 'Pmat' 
    if (c == NULL) c = new Camera<T_t>;
    int n;
    n = fscanf(fp, "%f %f %f %f", c->Pmat+0, c->Pmat+1, c->Pmat+2, c->Pmat+3);
    n = fscanf(fp, "%f %f %f %f", c->Pmat+4, c->Pmat+5, c->Pmat+6, c->Pmat+7);
    n = fscanf(fp, "%f %f %f %f", c->Pmat+8, c->Pmat+9, c->Pmat+10, c->Pmat+11);
    if (n != 4) {
      fprintf(stderr, "Error (CameraFileIO::readMCSCPmatFile): invalid projection matrix in '%s'\n", filename);
      return NULL;
    }
    fclose(fp);
#ifdef MTX_MATRIX_SOLVER
    c->updateParamUsingPmat();
#endif
    return true;
  }
  
  Camera<T_t>* readMCSCRadFile(Camera<T_t> *c, char *filename ) {
    // read the MCSC(Multi-Camera Self-Calibration) radial distortion file.
    // Usually they are names as "*.rad"
    FILE *fp = fopen(filename, "r");
    if (fp == NULL) {
      fprintf(stderr, "Error (CameraFileIO::readMCSCRadFile): cannot open file '%s'\n", filename);
      return NULL;
    }
    if (c == NULL) c = new Camera<T_t>;
    // intrinsic parameters
    float fc2[2], cc2[2], afc, temp, dt[4];
    fscanf(fp, "K11 = %f\n", fc2+0);	// focal length in pixels (X)
    fscanf(fp, "K12 = %f\n", &afc);	// skew (c->alpha * c->fc[1])
    fscanf(fp, "K13 = %f\n", cc2+0);	// principal point X coordinate
    fscanf(fp, "K21 = %f\n", &temp);
    fscanf(fp, "K22 = %f\n", fc2+1);	// focal length in pixels (Y)
    fscanf(fp, "K23 = %f\n", cc2+1);	// principal point Y coordinate
    fscanf(fp, "K31 = %f\n", &temp);
    fscanf(fp, "K32 = %f\n", &temp);
    if (fscanf(fp, "K33 = %f\n\n", &temp) != 1) return false;
    G2V_COPY( c->fc, fc2 );  G2V_COPY( c->cc, cc2 );
    c->afc = afc;
    // check 
    // read non-linear intrinsic parameters (radial distortions)
    fscanf(fp, "kc1 = %f\n", dt+0);	// 2th order coeff. of radial distortion
    fscanf(fp, "kc2 = %f\n", dt+1);	// 4th order coeff. of radial distortion
    fscanf(fp, "kc3 = %f\n", dt+2);	// coefficient of tangential distortion
    fscanf(fp, "kc4 = %f",   dt+3);	// coefficient of tangential distortion
    G4V_COPY( c->dt, dt );
    fclose(fp);
    return c;
  }
  
};

}	// end of CAMH namespace

#endif // CAMH_FILEIO_HPP

