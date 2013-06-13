
// 
// DC1394 class  -- based on libdc1394-2.0.2
//   IEEE1394 (Firewire) digital camera interface 
//
// Jaeil Choi
// last modified in Oct, 2008
//
// This class depends on
//   libdc1394 (ver 2.0.2) and libraw1394 (ver 1.2.0),
//   and should be linked with appropriate libraries, '-lraw1394' and '-ldc1394'.
//


#ifndef DC1394_H
#define DC1394_H

#include <iostream>
#include <dc1394/dc1394.h>
#include "util_color.h"


class DC1394 
{
private:
  dc1394_t		*dc1394;
  dc1394camera_list_t	*dclist;
public:
  dc1394camera_t	*camera;
  int			camera_index;
  int			camera_total;
  unsigned int		w, h;
  dc1394framerate_t	video_fps;
  dc1394video_mode_t	video_mode;
  dc1394color_coding_t	color_coding;
private:
  dc1394error_t		error;
  
public:
  // private members (public only because the capture thread uses them)
  dc1394video_frame_t  *framebuffer;
  pthread_t		thread_tid;
  bool			thread_running;
  void	(*capture_callback) (DC1394*);

public:
  bool			stereo;			// set by 'openCamera()'
  bool			skip_frames;
  bool			grayscale;		// whether the mode is MONO8 or MONO16
  bool			updated;
  bool			verbose;
  char			errmsg[256];
  
public:
  DC1394() : dc1394(NULL), dclist(NULL), camera(NULL), camera_index(0), camera_total(0),
	     framebuffer(NULL), thread_tid(0), thread_running(false), capture_callback(NULL), 
	     skip_frames(true), updated(false), verbose(false) { }
  ~DC1394() { closeCamera(); }
  
  // -----------------------------------------------------------------
  // opening / closing 
  // -----------------------------------------------------------------
public:
  bool openCamera(int cidx=0, int midx=-1, int fps=-1);
  void closeCamera(void);
  void resetCamera(void);
  
  // -----------------------------------------------------------------
  // functions for snapshot
  // -----------------------------------------------------------------
public:
  bool captureFrame(char pixel_type, char *fname,  char *secondary=NULL); // 'R'GB 'G'ray 'Y'UV
  bool captureFrame(char pixel_type, void *buffer, void *secondary=NULL);
  
  // -----------------------------------------------------------------
  // functions for streaming (using capture thread)
  // -----------------------------------------------------------------
public:
  bool startCaptureThread(void(*callback)(DC1394*)=NULL);
  void stopCaptureThread(void);
  bool copyToRGB (void *dest, void *secondary=NULL);	// for user callback function
  bool copyToGray(void *dest, void *secondary=NULL);	// for user callback function
  bool copyToYUV (void *dest, void *secondary=NULL);	// for user callback function
  
  // -----------------------------------------------------------------
  // private functions
  // -----------------------------------------------------------------
public:
  void skipOldFrames(int nLeftOver);
private:
  dc1394video_mode_t getModeByIndex(int midx) { 
    if (midx >= 0) {
      return (dc1394video_mode_t)(midx + (int)DC1394_VIDEO_MODE_160x120_YUV444); 
    } else if (camera) {
      dc1394video_modes_t vm;
      dc1394_video_get_supported_modes( camera, &vm );
      int mid = (int)vm.modes[0] - (int)DC1394_VIDEO_MODE_160x120_YUV444;
      return getModeByIndex( mid );
    } else return DC1394_VIDEO_MODE_160x120_YUV444;
  }
  dc1394framerate_t  getFrameRateByInt(int fps) {
    dc1394framerate_t ret;
    switch (fps) {
    case  1: ret = DC1394_FRAMERATE_1_875; break;
    case  3: ret = DC1394_FRAMERATE_3_75;  break;
    case  7: ret = DC1394_FRAMERATE_7_5;   break;
    case 15: ret = DC1394_FRAMERATE_15;    break;
    case 30: ret = DC1394_FRAMERATE_30;    break;
    case 60: ret = DC1394_FRAMERATE_60;    break;
    default: ret = DC1394_FRAMERATE_15;    break;
    }
    return ret;
  }
  dc1394speed_t getISOSpeedByInt(int speed) {
    dc1394speed_t ret;
    switch (speed) {
    case  100: ret = DC1394_ISO_SPEED_100;  break;
    case  200: ret = DC1394_ISO_SPEED_200;  break;
    case  400: ret = DC1394_ISO_SPEED_400;  break;
    case  800: ret = DC1394_ISO_SPEED_800;  break;
    case 1600: ret = DC1394_ISO_SPEED_1600;  break;
    case 3200: ret = DC1394_ISO_SPEED_3200;  break;
    default:   ret = DC1394_ISO_SPEED_400;  break;
    }
    return ret;
  }
  
  // -----------------------------------------------------------------
  // get/set functions for camera features
  // -----------------------------------------------------------------
public:
  bool getCameraFeature(char feature, int *value, int *min=NULL, int *max=NULL);
  bool setCameraFeature(char feature, int  value);
  bool getCameraWhiteBalance(int *u_b, int *v_r);
  bool setCameraWhiteBalance(int u_v, int v_r);
private:
  bool getFeatureByChar(char feature, dc1394feature_info_t **fip, dc1394feature_t *fidp=NULL);
  
  // -----------------------------------------------------------------
  // miscellaneous
  // -----------------------------------------------------------------
public:
  bool updateVideoInfo(void);
  void printInfo(void);
private:
  char *mode2string(dc1394video_mode_t mode);
  char *color2string(dc1394color_coding_t color);
  char *fps2string(dc1394framerate_t fps);
  
};

#endif // DC1394_H

