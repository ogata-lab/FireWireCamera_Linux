
// 
// DC1394 class  -- based on libdc1394-2.1.0
//   IEEE1394 (Firewire) digital camera interface 
//
// Jaeil Choi
// last modified in Oct, 2008
//
// This class depends on libdc1394 (ver 2.1.0), and 
//   it should be linked with '-ldc1394'.
//

#include <iostream>
#include <pthread.h>
#include "util_pixfmt_converter.hpp"
#include "imgh_common.hpp"
#include "imgh_editor.hpp"
#include "imgh_fileio.hpp"
#include "DC1394.hpp"
using namespace std;


// -------------------------------------------------------------------
// functions for the capture thread
// -------------------------------------------------------------------

static void* capture_thread(void* p)
{
  DC1394 *fcam = (DC1394*) p;
  fcam->thread_running = true;
  
  // capture the frame
  while (true) {
    // get the current image
    fcam->skipOldFrames(1);
    if (dc1394_capture_dequeue( fcam->camera, DC1394_CAPTURE_POLICY_WAIT, 
				&(fcam->framebuffer) ) != DC1394_SUCCESS) {
      cerr << "Error (DC1394 capture_thread): failed to capture" << endl;
      fcam->framebuffer = NULL;  break;
    }
    
    // run user-callback function
    if (fcam->capture_callback) fcam->capture_callback(fcam);
    fcam->updated = true;
    
    // release the frame
    if (fcam->framebuffer) 
      dc1394_capture_enqueue( fcam->camera, fcam->framebuffer );
    fcam->framebuffer = NULL;
    
    if (fcam->thread_running == false) break;	// stopped by other thread
  }
  
  fcam->thread_running = false;
  fcam->thread_tid = 0;
  return NULL;
}

// -------------------------------------------------------------------
// opening / closing 
// -------------------------------------------------------------------

bool DC1394::openCamera(int cidx, int mode, int fps)
{
  closeCamera();
  if (!dc1394) dc1394 = dc1394_new();
  
  // turn off some messages
  //dc1394_log_register_handler(DC1394_LOG_ERROR, NULL, NULL);
  dc1394_log_register_handler(DC1394_LOG_DEBUG, NULL, NULL);
  dc1394_log_register_handler(DC1394_LOG_WARNING, NULL, NULL);
  
  // find cameras
  error = dc1394_camera_enumerate( dc1394, &dclist );
  if (error != DC1394_SUCCESS) {
    if (verbose) {
      cerr << "Error (DC1394::initializeCamera): can't find cameras" << endl;
      cerr << "On Linux, you may want to check that " << endl;
      cerr << "  - the kernel modules `ieee1394',`raw1394' and `ohci1394' are loaded " << endl;
      cerr << "  - you have read/write access to /dev/raw1394" << endl;
      cerr << "  - the camera is connected with power" << endl;
    }
    return false;
  }
  camera_total = dclist->num;
  if (camera_total < 1) { 
    if (verbose) printf("[DC1394] No cameras found\n");
    dc1394_camera_free_list (dclist);
    return false; 
  } else {
    if (verbose) printf("[DC1394] Found %d camera%s\n", dclist->num, (dclist->num>1 ? "s":""));
  }
  if (cidx >= camera_total) {
    if (verbose) printf("[DC1394] camera %d requested, but only %d cameras available\n", cidx, camera_total);
    dc1394_camera_free_list (dclist);
    return false; 
  }
  
  if (verbose) { printf("[DC1394] Connecting to camera %d / %d .. ", cidx, dclist->num); fflush(stdout); }
  camera = dc1394_camera_new (dc1394, dclist->ids[cidx].guid);
  if (!camera) {
    dc1394_log_error("Failed to initialize camera with guid %llx", dclist->ids[cidx].guid);
    dc1394_camera_free_list (dclist);
    return false;
  }
  camera_index = cidx;
  if (verbose) printf("done\n");
  
  /*-----------------------------------------------------------------------
   *  setup capture
   *-----------------------------------------------------------------------*/
  
  video_mode = getModeByIndex( mode );
  video_fps  = getFrameRateByInt( fps );
  dc1394_get_color_coding_from_video_mode( camera, video_mode, &color_coding );
  
  error=dc1394_video_set_iso_speed(camera, DC1394_ISO_SPEED_400);
  if (error != DC1394_SUCCESS) { if (verbose) printf("failed to set iso speed\n"); return false; }

  error=dc1394_video_set_mode(camera, video_mode);
  if (error != DC1394_SUCCESS) { if (verbose) printf("failed to set video mode\n"); return false; }

  error=dc1394_video_set_framerate(camera, video_fps);
  if (error != DC1394_SUCCESS) { if (verbose) printf("failed to set framerate\n"); return false; }

  error=dc1394_capture_setup(camera, 4, DC1394_CAPTURE_FLAGS_DEFAULT);
  if (error != DC1394_SUCCESS) { if (verbose) printf("failed to setup camera\n"); return false; }
  //make sure that the video mode and framerate are\nsupported by your camera\n");

  /*-----------------------------------------------------------------------
   *  have the camera start sending us data
   *-----------------------------------------------------------------------*/
  error=dc1394_video_set_transmission(camera, DC1394_ON);
  //DC1394_ERR_CLN_RTN(error,cleanup_and_exit(camera),"Could not start camera iso transmission\n");
  
  updateVideoInfo();
  
  if (verbose) fprintf(stderr, "[DC1394] camera  opened  (%d-th camera)\n", camera_index);
  dc1394_camera_free_list (dclist);
  return true;
}

void DC1394::closeCamera(void)
{
  if (camera) {
    // stop capture thread, if it's running
    if (thread_tid) stopCaptureThread();
    // close each camera
    dc1394_video_set_transmission(camera, DC1394_OFF);
    dc1394_capture_stop(camera);
    dc1394_camera_free(camera);
    if (verbose) fprintf(stderr, "[DC1394] camera  closed \n");
  }
  camera = NULL;
  if (dc1394) { dc1394_free( dc1394 ); dc1394 = NULL; }
}

void DC1394::resetCamera(void)
{
  dc1394_camera_reset( camera );
  if (verbose) fprintf(stderr, "[DC1394] camera  reset \n");
}

// -------------------------------------------------------------------
// functions for snapshot
// -------------------------------------------------------------------

bool DC1394::captureFrame(char pixel_type, char *fname, char *secondary)
{
  // Capture an image from the camera with a pixel format 'pixel_type' {'G','R','Y'},
  //   and save the image in memory 'dest'. For stereo camera,
  //   left and right images are saved separately in 'fname' and 'secondary'.
  if (!camera) return false;
  IMGH::Image temp( w, h, IMGH::PIXEL_RGB );
  bool ret;
  // capture the frame
  if (!stereo) {
    ret = captureFrame( pixel_type, temp.data );
  } else {
    IMGH::Image tempL( w, h, IMGH::PIXEL_RGB );
    IMGH::Image tempR( w, h, IMGH::PIXEL_RGB );
    ret = captureFrame( pixel_type, tempL.data, tempR.data);
    if (ret) {
      IMGH::ImageEditor imgedt;
      imgedt.mergeh( &tempL, &tempR, &temp );
    }
  }
  if (ret) {
    IMGH::ImageFileIO ifile;
    ifile.writeFile( fname, &temp );
  }
  return ret;
}

bool DC1394::captureFrame(char pixel_type, void *dest, void *secondary)
{
  // Capture an image from the camera with a pixel format 'pixel_type' {'G','R','Y'},
  //   and save the image in memory 'dest'. For stereo camera,
  //   the left image is saved in 'dest', and the right in 'secondary'.
  if (!camera) return false;
  if (dest == NULL) { sprintf(errmsg, "invalid buffer (NULL)");  return false; }
  // Stop the running thread first, to avoid conflict.
  void (*callback_backup)(DC1394*) = NULL;
  if (thread_running) { callback_backup = capture_callback;  stopCaptureThread(); }
  
  // get the current frame
  if (skip_frames) skipOldFrames(0);
  if (dc1394_capture_dequeue( camera, DC1394_CAPTURE_POLICY_WAIT, 
			      &framebuffer ) != DC1394_SUCCESS ) {
    sprintf(errmsg, "failed to dequeue in captureFrame()");
    return false;
  }
  
  video_mode   = framebuffer->video_mode;
  color_coding = framebuffer->color_coding;
  grayscale = (color_coding == DC1394_COLOR_CODING_MONO8 || 
	       color_coding == DC1394_COLOR_CODING_MONO16);
  w = framebuffer->size[0];   h = framebuffer->size[1];
  
  bool ret = true;
  switch (pixel_type) {
  case 'G': case 'g':  ret = copyToGray( dest, secondary );  break;
  case 'R': case 'r':  ret = copyToRGB ( dest, secondary );  break;
  case 'Y': case 'y':  ret = copyToYUV ( dest, secondary );  break;
  }
  
  // release the frame buffer of the camera
  dc1394_capture_enqueue( camera, framebuffer );
  
  // Resume the lastest running thread, if any.
  if (callback_backup) startCaptureThread( callback_backup );
  return ret;
}

void DC1394::skipOldFrames(int nLeftOver)
{
  int frames_behind;
  do {
    if (dc1394_capture_dequeue( camera, DC1394_CAPTURE_POLICY_WAIT, 
				&framebuffer ) != DC1394_SUCCESS ) {
      sprintf(errmsg, "failed to dequeue in captureFrame()");
      return;
    }
    frames_behind = framebuffer->frames_behind;
    dc1394_capture_enqueue( camera, framebuffer );
  } while (frames_behind > nLeftOver);
}

// -------------------------------------------------------------------
// functions for streaming (using capture thread)
// -------------------------------------------------------------------

bool DC1394::startCaptureThread( void(*callback)(DC1394*) )
{
  if (!camera) return false;
  this->capture_callback = callback;
  this->updated = false;
  // create the thread for the capture in a infinite loop
  pthread_create( &thread_tid, NULL, capture_thread, this );
  if (verbose) cerr << "[DC1394] capture thread started (thread_id = " << thread_tid << ")" << endl;
  return true;
}

void DC1394::stopCaptureThread(void)
{
  if (thread_tid == 0) return;
  pthread_t saved_tid = thread_tid;
  // set the capture thread to stop
  thread_running = false;
  // wait until the capture thread is completely stopped
  pthread_join( thread_tid, NULL );
  if (verbose) cerr << "[DC1394]         thread stopped (thread_id = " << saved_tid << ")" << endl;
  thread_tid = 0;
}

bool DC1394::copyToRGB (void *dest, void *secondary)
{
  // color space conversion to 3 x UChar 8 bits
  //   'dest'  : destination RGB buffer 
  //   'dest2' : destination RGB buffer for secondary camera of stereo
  if ( !framebuffer || !dest ) return false;
  UTIL::PixFmtConverter conv;
  static void          *temp=NULL;
  static unsigned long  temp_size=0;
  // NEAREST, SIMPLE, BILINEAR, HQLINEAR, DOWNSAMPLE, EDGESENSE, VNG, AHD
  dc1394bayer_method_t  bayer_decode_method = DC1394_BAYER_METHOD_AHD;  // BILINEAR;
  
  switch (color_coding) {
  case DC1394_COLOR_CODING_MONO8:	// 8 bits per pixel (Gray -> RGB)
    conv.convert_gray8_rgb24( framebuffer->image, dest, w, h );
    break;
  case DC1394_COLOR_CODING_MONO16:	// 16 bits per pixel (Gray -> RGB)
    conv.convert_gray16_rgb24( framebuffer->image, dest, w, h );
    break;
  case DC1394_COLOR_CODING_YUV422:	// (UYVY) 16 bits per pixel (UYVY -> RGBRGB)
    conv.convert_uyvy_rgb24( framebuffer->image, dest, w, h );
    break;
  case DC1394_COLOR_CODING_YUV444:	// (IYU2) 24 bits per pixel (UYV -> RGB)
    conv.convert_uyv_rgb24( framebuffer->image, dest, w, h );
    break;
  case DC1394_COLOR_CODING_YUV411:	// (IYU1) 12 bits per pixel (uyyvyy -> RGBRGBRGBRGB)
    conv.convert_uyyvyy_rgb24( framebuffer->image, dest, w, h );
    break;
  case DC1394_COLOR_CODING_RGB8:
    memcpy ( dest, framebuffer->image, w * h * 3 );
    break;
  case DC1394_COLOR_CODING_RAW8:
    // choose a method from {NEAREST,SIMPLE,BILINEAR,HQLINEAR,DOWNSAMPLE,VNG,AHD,DOWNSAMPLE}
    dc1394_bayer_decoding_8bit((uint8_t*)framebuffer->image, (uint8_t*)dest, w, h, 
			       framebuffer->color_filter, bayer_decode_method);
    break;
  case DC1394_COLOR_CODING_RAW16:
    if (!temp || temp_size != 2*w*h*sizeof(uint8_t)) {
      if (temp) free(temp);
      temp_size = 2 * w * h * sizeof(uint8_t);
      temp = (unsigned char*) malloc( temp_size );
    }
    if (strncmp(camera->model, "Bumblebee2",10)==0) {
      dc1394_deinterlace_stereo ( (uint8_t*)framebuffer->image, (uint8_t*)temp, w, 2*h );
      dc1394_bayer_decoding_8bit( (uint8_t*)temp,       (uint8_t*)dest,         w, h,
				  framebuffer->color_filter, bayer_decode_method );
      if (secondary) 
	dc1394_bayer_decoding_8bit( (uint8_t*)temp + w*h, (uint8_t*)secondary, w, h,
				    framebuffer->color_filter, bayer_decode_method );
    } else {
      // choose a method from {NEAREST,SIMPLE,BILINEAR,HQLINEAR,EDGESENSE,VNG,AHD,DOWNSAMPLE}
      dc1394_bayer_decoding_16bit((uint16_t*)framebuffer->image, (uint16_t*)dest, w, h, 
				  framebuffer->color_filter, bayer_decode_method, 8);
    }
    break;
  default:
    cerr << "Error (DC1394::copyToRGB): invalid pixel format (or not implemented yet)" << endl;
    return false;
  }
  return true;
}

bool DC1394::copyToGray(void *dest, void *secondary)
{
  // color space conversion to 1 x UChar 8 bits
  if ( !framebuffer || !dest ) return false;
  UTIL::PixFmtConverter conv;
  static unsigned char *temp=NULL;
  static unsigned long  temp_size=0;
  // NEAREST, SIMPLE, BILINEAR, HQLINEAR, DOWNSAMPLE, EDGESENSE, VNG, AHD
  dc1394bayer_method_t  bayer_decode_method = DC1394_BAYER_METHOD_AHD;  // BILINEAR;
  
  switch (color_coding) {
  case DC1394_COLOR_CODING_MONO8:	// 8 bits per pixel (Gray8 -> Gray8)
    memcpy ( dest, framebuffer->image, w * h );
    break;
  case DC1394_COLOR_CODING_MONO16:	// 16 bits per pixel (Gray16 -> Gray8)
    conv.convert_gray16_gray8( framebuffer->image, dest, w, h );
    break;
  case DC1394_COLOR_CODING_YUV422:	// (YUY2) 16 bits per pixel (YUYV -> WW)
    conv.convert_uyvy_gray8 ( framebuffer->image, dest, w, h );
    break;
  case DC1394_COLOR_CODING_YUV444:	// (IYU2) 24 bits per pixel (UYV -> Y)
    conv.convert_uyv_gray8 ( framebuffer->image, dest, w, h );
    break;
  case DC1394_COLOR_CODING_YUV411:	// (IYU1) 12 bits per pixel (uyyvyy -> WWWW)
    conv.convert_uyyvyy_gray8 ( framebuffer->image, dest, w, h );
    break;
  case DC1394_COLOR_CODING_RGB8:	// 24 bits per pixel (RGB -> Y)
    conv.convert_rgb24_gray8 ( framebuffer->image, dest, w, h );
    break;
  case DC1394_COLOR_CODING_RAW8:	// 8 bits per pixel (Bayer -> Y)
    if (!temp || temp_size != 3 * w * h * sizeof(unsigned char)) {
      if (temp) free(temp);
      temp_size = 3 * w * h * sizeof(unsigned char);
      temp = (unsigned char*)malloc( temp_size );
    }
    // choose a method from {NEAREST,SIMPLE,BILINEAR,HQLINEAR,DOWNSAMPLE,EDGESENSE,VNG}
    dc1394_bayer_decoding_8bit((uint8_t*)framebuffer->image, (uint8_t*)temp, w, h, 
			       framebuffer->color_filter, bayer_decode_method);
    conv.convert_rgb24_gray8( temp, dest, w, h );
    break;
  case DC1394_COLOR_CODING_RAW16:	// 16 bits per pixel (Bayer -> Y)
    if (!temp || temp_size != 2 * 3 * w * h * sizeof(unsigned char)) {
      if (temp) free(temp);
      temp_size = 2 * 3 * w * h * sizeof(unsigned char);
      temp = (unsigned char*)malloc( temp_size );
    }
    if (strncmp(camera->model, "Bumblebee2",10)==0) {
      dc1394_deinterlace_stereo ( (uint8_t*)framebuffer->image, (uint8_t*)temp, w, 2*h );
      dc1394_bayer_decoding_8bit( (uint8_t*)temp,       (uint8_t*)dest,         w, h,
				  framebuffer->color_filter, bayer_decode_method );
      //conv.convert_rgb24_gray8( temp, dest, w, h );
      if (secondary) {
	dc1394_bayer_decoding_8bit( (uint8_t*)temp + w*h, (uint8_t*)secondary, w, h,
				    framebuffer->color_filter, bayer_decode_method );
	//conv.convert_rgb24_gray8( (uint8_t*)temp+w*h, (uint8_t*)secondary, w, h );
      }
    } else {
      // choose a method from {NEAREST,SIMPLE,BILINEAR,HQLINEAR,DOWNSAMPLE,EDGESENSE,VNG}
      dc1394_bayer_decoding_16bit((uint16_t*)framebuffer->image, (uint16_t*)temp, w, h, 
				  framebuffer->color_filter, bayer_decode_method, 8);
      conv.convert_rgb24_gray8( temp, dest, w, h );
    }
  default:
    cerr << "Error (DC1394::copyToGray): invalid pixel format (or not implemented yet)" << endl;
    return false;
    break;
  }
  return true;
}

bool DC1394::copyToYUV (void *dest, void *secondary)
{
  // color space conversion to 3 x UChar 8 bits
  UTIL::PixFmtConverter conv;
  if ( !framebuffer || !dest ) return false;
  switch (color_coding) {
  case DC1394_COLOR_CODING_MONO8:
  case DC1394_COLOR_CODING_MONO16:
  case DC1394_COLOR_CODING_YUV422:
  case DC1394_COLOR_CODING_YUV444:
  case DC1394_COLOR_CODING_YUV411:
  case DC1394_COLOR_CODING_RGB8:
  default:
    cerr << "Error (DC1394::copyToYUV): invalid pixel format (or not implemented yet)" << endl;
    return false;
  }
  return true;
}
  
// -------------------------------------------------------------------
// get/set functions for camera features
// -------------------------------------------------------------------

bool DC1394::getCameraFeature(char feature, int *value, int *minp, int *maxp)
{
  // get the value of the feature of the camera
  //   feature : 'b'rightness, 'e'xposure, 's'harpness, 'w'hitebalance
  //   value   : non-negative values for manual mode, -1 for auto mode
  dc1394feature_info_t *fip=NULL;
  if (!getFeatureByChar( feature, &fip )) return false;
//   *value = (int)fip->value;
//   if (minp) *minp = (int)fip->min;
//   if (maxp) *maxp = (int)fip->max;
  return true;
}

bool DC1394::setCameraFeature(char feature, int value)
{
  // set the value of the feature of the camera
  //   feature : 'b'rightness, 'e'xposure, 's'harpness, 'w'hitebalance
  //   value   : non-negative values for manual mode, -1 for auto mode
  dc1394feature_t fid;
  dc1394feature_info_t *fip;
  if (!getFeatureByChar( feature, &fip, &fid )) return false;
//   if (value < 0) {	// auto mode
//     if (!fip->auto_capable) return false;
//     error = dc1394_feature_set_mode( camera, fid, DC1394_FEATURE_MODE_AUTO );
//   } else  {		// manual mode. set value
//     if (!fip->manual_capable) return false;
//     error = dc1394_feature_set_mode( camera, fid, DC1394_FEATURE_MODE_MANUAL );
//     if (error != DC1394_SUCCESS) return false;
//     error = dc1394_feature_set_value( camera, fid, value );
//   }
  return (error == DC1394_SUCCESS);
}

bool DC1394::getCameraWhiteBalance(int *u_b, int *v_r)
{
  dc1394feature_info_t *fip;
  if (!getFeatureByChar( 'w', &fip )) return false;
  unsigned int ub, vr;
  error = dc1394_feature_whitebalance_get_value( camera, &ub, &vr );
  *u_b = (int)ub;   *v_r = (int)vr;
  return (error == DC1394_SUCCESS);
}

bool DC1394::setCameraWhiteBalance(int u_b, int v_r)
{
  dc1394feature_t fid;
  dc1394feature_info_t *fip;
  if (!getFeatureByChar( 'w', &fip, &fid)) return false;
//   if (u_b < 0 || v_r < 0) {	// auto mode
//     if (!fip->auto_capable) return false;
//     error = dc1394_feature_set_mode( camera, fid, DC1394_FEATURE_MODE_AUTO );
//   } else {
//     if (!fip->manual_capable) return false;
//     error = dc1394_feature_set_mode( camera, fid, DC1394_FEATURE_MODE_MANUAL );
//     if (error != DC1394_SUCCESS) return false;
//     error = dc1394_feature_whitebalance_set_value( camera, u_b, v_r );
//   }
  return (error == DC1394_SUCCESS);
}

bool DC1394::getFeatureByChar(char feature, dc1394feature_info_t **fip, dc1394feature_t *fidp)
{
  dc1394feature_t  fid=DC1394_FEATURE_BRIGHTNESS;
  switch (feature) {
  case 'b':  case 'B':  fid = DC1394_FEATURE_BRIGHTNESS;  break;
  case 'e':  case 'E':  fid = DC1394_FEATURE_EXPOSURE;  break;
  case 'f':  case 'F':  fid = DC1394_FEATURE_FRAME_RATE;  break;
  case 'g':  case 'G':  fid = DC1394_FEATURE_GAIN;  break;
  case 'h':  case 'H':  fid = DC1394_FEATURE_HUE;  break;
  case 's':  case 'S':  fid = DC1394_FEATURE_SHARPNESS;  break;
  case 'w':  case 'W':  fid = DC1394_FEATURE_WHITE_BALANCE;  break;
  default:  return false;
  }
//   if (fidp) *fidp = fid;
//   // get all the feature set
//   static dc1394featureset_t fset;
//   if (dc1394_get_camera_feature_set(camera,&fset)!=DC1394_SUCCESS) return false;
//   // find the specific feature
//   for (int fidx = 0; true; ) {
//     if (fset.feature[fidx].id == fid) { *fip = fset.feature+fidx; break; }
//     if (++fidx == DC1394_FEATURE_NUM) return false;
//   }
//   // check if the feature is available
//   if (!(*fip)->available) return false;
  return true;
}


// -------------------------------------------------------------------
//
// -------------------------------------------------------------------

bool DC1394::updateVideoInfo(void)
{
  // update the value of video mode by reading a frame
  if (verbose) { printf("[DC1394] updateVideoInfo (dc1394_capture_dequeue) .. "); fflush(stdout); }
  //   (This is where the warning happens: "libdc1394 warning: packet 955 had error.."
  if (dc1394_capture_dequeue( camera, DC1394_CAPTURE_POLICY_WAIT, 
			      &framebuffer ) != DC1394_SUCCESS ) {
    sprintf(errmsg, "failed to dequeue in updateVideoInfo()");
    return false;
  }
  if (verbose) printf("done\n");
  video_mode = framebuffer->video_mode;
  dc1394_capture_enqueue( camera, framebuffer );
  
  // update mode information
  if (video_mode < DC1394_VIDEO_MODE_FORMAT7_0) {	// format 0, 1, 2
    dc1394_get_image_size_from_video_mode  ( camera, video_mode, &w, &h );
    dc1394_get_color_coding_from_video_mode( camera, video_mode, &color_coding );
  } else {					// format 7
    dc1394format7mode_t f7;
    dc1394_format7_get_mode_info( camera, video_mode, &f7 );
    w = f7.size_x;  h = f7.size_y;
    color_coding = f7.color_coding;
  }
  grayscale = (color_coding == DC1394_COLOR_CODING_MONO8 || 
	       color_coding == DC1394_COLOR_CODING_MONO16);
  
  // update frame rate
  error = dc1394_video_get_framerate(camera, &video_fps);
  ////if (error != DC1394_SUCCESS) video_fps = 0;
  ////else dc1394_framerate_as_float( video_fps, &fps );
  
  stereo = (strncmp(camera->model,"Bumblebee",9)==0 && video_mode >= DC1394_VIDEO_MODE_FORMAT7_0);
  
  // decide Bayer pattern
//   if ( strcmp(camera->vendor, "Point Grey Research") == 0) {
//     uint32_t qValue;	// read Bayer Tile Mapping register directly
//     GetCameraControlRegister( camera, 0x1040, &qValue);
//     switch (qValue) {
//     case 0x42474752:  camera->cf = DC1394_COLOR_FILTER_BGGR;  break;
//     case 0x47524247:  camera->cf = DC1394_COLOR_FILTER_GRBG;  break;
//     case 0x52474742:  camera->cf = DC1394_COLOR_FILTER_RGGB;  break;
//     case 0x47425247:  camera->cf = DC1394_COLOR_FILTER_GBRG;  break;
//     case 0x59595959:  // YYYY
//     default:          camera->cf = DC1394_COLOR_FILTER_GRBG; break;
//     }
//   } else camera->cf = DC1394_COLOR_FILTER_GRBG;
  return true;
}

char *DC1394::mode2string(dc1394video_mode_t mode)
{
  char *str;
  switch (mode) {
  case DC1394_VIDEO_MODE_160x120_YUV444: str = "DC1394_VIDEO_MODE_160x120_YUV444"; break;
  case DC1394_VIDEO_MODE_320x240_YUV422: str = "DC1394_VIDEO_MODE_320x240_YUV422"; break;
  case DC1394_VIDEO_MODE_640x480_YUV411: str = "DC1394_VIDEO_MODE_640x480_YUV411"; break;
  case DC1394_VIDEO_MODE_640x480_YUV422: str = "DC1394_VIDEO_MODE_640x480_YUV422"; break;
  case DC1394_VIDEO_MODE_640x480_RGB8:   str = "DC1394_VIDEO_MODE_640x480_RGB8"; break;
  case DC1394_VIDEO_MODE_640x480_MONO8:  str = "DC1394_VIDEO_MODE_640x480_MONO8"; break;
  case DC1394_VIDEO_MODE_640x480_MONO16: str = "DC1394_VIDEO_MODE_640x480_MONO16"; break;
  case DC1394_VIDEO_MODE_800x600_YUV422: str = "DC1394_VIDEO_MODE_800x600_YUV422"; break;
  case DC1394_VIDEO_MODE_800x600_RGB8:   str = "DC1394_VIDEO_MODE_800x600_RGB8"; break;
  case DC1394_VIDEO_MODE_800x600_MONO8:  str = "DC1394_VIDEO_MODE_800x600_MONO8"; break;
  case DC1394_VIDEO_MODE_1024x768_YUV422: str = "DC1394_VIDEO_MODE_1024x768_YUV422"; break;
  case DC1394_VIDEO_MODE_1024x768_RGB8:   str = "DC1394_VIDEO_MODE_1024x768_RGB8"; break;
  case DC1394_VIDEO_MODE_1024x768_MONO8:  str = "DC1394_VIDEO_MODE_1024x768_MONO8"; break;
  case DC1394_VIDEO_MODE_800x600_MONO16:  str = "DC1394_VIDEO_MODE_800x600_MONO16"; break;
  case DC1394_VIDEO_MODE_1024x768_MONO16: str = "DC1394_VIDEO_MODE_1024x768_MONO16"; break;
  case DC1394_VIDEO_MODE_1280x960_YUV422:  str = "DC1394_VIDEO_MODE_1280x960_YUV422"; break;
  case DC1394_VIDEO_MODE_1280x960_RGB8:    str = "DC1394_VIDEO_MODE_1280x960_RGB8"; break;
  case DC1394_VIDEO_MODE_1280x960_MONO8:   str = "DC1394_VIDEO_MODE_1280x960_MONO8"; break;
  case DC1394_VIDEO_MODE_1600x1200_YUV422: str = "DC1394_VIDEO_MODE_1600x1200_YUV422"; break;
  case DC1394_VIDEO_MODE_1600x1200_RGB8:   str = "DC1394_VIDEO_MODE_1600x1200_RGB8"; break;
  case DC1394_VIDEO_MODE_1600x1200_MONO8:  str = "DC1394_VIDEO_MODE_1600x1200_MONO8"; break;
  case DC1394_VIDEO_MODE_1280x960_MONO16:  str = "DC1394_VIDEO_MODE_1280x960_MONO16"; break;
  case DC1394_VIDEO_MODE_1600x1200_MONO16: str = "DC1394_VIDEO_MODE_1600x1200_MONO16"; break;
  case DC1394_VIDEO_MODE_FORMAT7_0: str = "DC1394_VIDEO_MODE_FORMAT7_0"; break;
  case DC1394_VIDEO_MODE_FORMAT7_1: str = "DC1394_VIDEO_MODE_FORMAT7_1"; break;
  case DC1394_VIDEO_MODE_FORMAT7_2: str = "DC1394_VIDEO_MODE_FORMAT7_2"; break;
  case DC1394_VIDEO_MODE_FORMAT7_3: str = "DC1394_VIDEO_MODE_FORMAT7_3"; break;
  case DC1394_VIDEO_MODE_FORMAT7_4: str = "DC1394_VIDEO_MODE_FORMAT7_4"; break;
  case DC1394_VIDEO_MODE_FORMAT7_5: str = "DC1394_VIDEO_MODE_FORMAT7_5"; break;
  case DC1394_VIDEO_MODE_FORMAT7_6: str = "DC1394_VIDEO_MODE_FORMAT7_6"; break;
  case DC1394_VIDEO_MODE_FORMAT7_7: str = "DC1394_VIDEO_MODE_FORMAT7_7"; break;
  default:  str = "DC1394_VIDEO_MODE_UNKNOWN"; break;
  }
  return str;
}

char *DC1394::color2string(dc1394color_coding_t color)
{
  char *str;
  switch (color) {
  case DC1394_COLOR_CODING_MONO8:   str = "MONO8";   break;
  case DC1394_COLOR_CODING_YUV411:  str = "YUV411";  break;
  case DC1394_COLOR_CODING_YUV422:  str = "YUV422";  break;
  case DC1394_COLOR_CODING_YUV444:  str = "YUV444";  break;
  case DC1394_COLOR_CODING_RGB8:    str = "RGB8";    break;
  case DC1394_COLOR_CODING_MONO16:  str = "MONO16";  break;
  case DC1394_COLOR_CODING_RGB16:   str = "RGB16";   break;
  case DC1394_COLOR_CODING_MONO16S: str = "MONO16S"; break;
  case DC1394_COLOR_CODING_RGB16S:  str = "RGB16S";  break;
  case DC1394_COLOR_CODING_RAW8:    str = "RAW8";    break;
  case DC1394_COLOR_CODING_RAW16:   str = "RAW16";   break;
  default: str = "UNKNOWN"; break;
  }
  return str;
}

char *DC1394::fps2string(dc1394framerate_t fps)
{
  char *str;
  switch(fps) {
  case DC1394_FRAMERATE_1_875: str = "1.875 fps"; break;
  case DC1394_FRAMERATE_3_75:  str = "3.75 fps"; break;
  case DC1394_FRAMERATE_7_5:   str = "7.5 fps"; break;
  case DC1394_FRAMERATE_15:    str = "15 fps"; break;
  case DC1394_FRAMERATE_30:    str = "30 fps"; break;
  case DC1394_FRAMERATE_60:    str = "60 fps"; break;
  default: str = "UNKNOWN"; break;
  }
  return str;
}

void DC1394::printInfo(void)
{
  if (!camera) return;
  int  i, j, mid;
  printf("Camera  detail information  (%d/%d  %s)\n", camera_index, camera_total, (thread_tid ? "Running" : "NotRunning"));
  printf("  %-12s : [%s] %s \n", "Make & Model", camera->vendor, camera->model);
  printf("  %-12s : %s  %s \n", "Mode & FPS", mode2string(video_mode), fps2string(video_fps));
  dc1394featureset_t fset;
  if(dc1394_feature_get_all(camera, &fset) !=DC1394_SUCCESS) {
    printf("  unable to get all the features\n");
  } else {
    // dc1394_feature_print_all(&fset, stdout);
//     char *fname, *bar;
//     for (i = 0; i < DC1394_FEATURE_NUM; i++) {
//       dc1394feature_t       fid = fset.feature[i].id;
//       dc1394feature_info_t *fip = fset.feature + i;
//       if (fid == DC1394_FEATURE_BRIGHTNESS ||
// 	  fid == DC1394_FEATURE_EXPOSURE ||
// 	  fid == DC1394_FEATURE_SHARPNESS ) {
// 	switch (fid) {
// 	case DC1394_FEATURE_BRIGHTNESS: fname = "Brightness"; break;
// 	case DC1394_FEATURE_EXPOSURE:   fname = "Exposure";   break;
// 	case DC1394_FEATURE_SHARPNESS:  fname = "Sharpness";  break;
// 	default: continue;
// 	}
// 	if (fip->available) {
// 	  float rate = ((float)fip->value - fip->min) / (fip->max - fip->min);
// 	  if      (rate < 0.1) bar = "+---------";
// 	  else if (rate < 0.2) bar = "-+--------";
// 	  else if (rate < 0.3) bar = "--+-------";
// 	  else if (rate < 0.4) bar = "---+------";
// 	  else if (rate < 0.5) bar = "----+-----";
// 	  else if (rate < 0.6) bar = "-----+----";
// 	  else if (rate < 0.7) bar = "------+---";
// 	  else if (rate < 0.8) bar = "-------+--";
// 	  else if (rate < 0.9) bar = "--------+-";
// 	  else                 bar = "---------+";
// 	  printf("  %-12s : %-6s [%2s %2s]  %4d [%3d %s %4d] \n", fname, 
// 		 (fip->auto_capable && fip->auto_active ? "Auto" : "Manual"),
// 		 (fip->auto_capable ? "AC" : "  "),
// 		 (fip->manual_capable ? "MC" : "  "),
// 		 fip->value, fip->min, bar, fip->max);
// 	} else {
// 	  printf("  %-12s : unavailable \n", fname);
// 	}
//       } else if (fid == DC1394_FEATURE_WHITE_BALANCE) {
// 	if (fip->available) {
// 	  unsigned int u_b, v_r;
// 	  dc1394_feature_whitebalance_get_value( camera, &u_b, &v_r );
// 	  printf("  %-12s : %-6s [%2s %2s]  u_b = %d  v_r = %d \n", "WhiteBalance",
// 		 (fip->auto_capable && fip->auto_active ? "Auto" : "Manual"),
// 		 (fip->auto_capable ? "AC" : "  "),
// 		 (fip->manual_capable ? "MC" : "  "), u_b, v_r);
// 	} else {
// 	  printf("  %-12s : unavailable \n", "WhiteBalance");
// 	}
//       }
//     }
  }
  
  printf("  Supported Modes and Framerates:\n");
  dc1394video_modes_t vm;
  dc1394framerates_t  fr;
  dc1394_video_get_supported_modes( camera, &vm );
  for (i = 0; i < (int)vm.num; i++) {
    mid = (int)vm.modes[i] - (int)DC1394_VIDEO_MODE_160x120_YUV444;
    printf("  %02d %-34s ", mid, mode2string(vm.modes[i]));
    if (vm.modes[i] < DC1394_VIDEO_MODE_FORMAT7_0) {
      // ordinary predetermined modes
      error = dc1394_video_get_supported_framerates( camera, vm.modes[i], &fr );
      if (error != DC1394_SUCCESS) printf(" unknown");
      else {
	for (j = 0; j < (int)fr.num; j++) {
	  float fps;
	  dc1394_framerate_as_float( fr.framerates[j], &fps );
	  printf(" %3g", fps);
	}
      }
    } else {      // format 7
      unsigned int          w, h, bpp=0, psize=0, ppf=0;
      dc1394color_coding_t  cc;
      dc1394_format7_get_max_image_size( camera, vm.modes[i], &w, &h );
      dc1394_format7_get_color_coding( camera, vm.modes[i], &cc );
      dc1394_format7_get_packet_size( camera, vm.modes[i], &psize );
      dc1394_format7_get_packets_per_frame( camera, vm.modes[i], &ppf );
      printf(" %d x %d  %6s  Pkg=%d Pkg/F=%d", w, h, color2string(cc), bpp, ppf);
    }
    printf("\n");
  }
}

