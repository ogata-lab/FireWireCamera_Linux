// -*- C++ -*-
/*!
 * @file  FireWireCamera.cpp
 * @brief Camera components , captrue image form 1394 Firewire  camera 
 * @date $Date$
 *
 * $Id$
 */

#include "FireWireCamera.h"

// Module specification
// <rtc-template block="module_spec">
static const char* firewirecamera_spec[] =
  {
    "implementation_id", "FireWireCamera",
    "type_name",         "FireWireCamera",
    "description",       "Camera components , captrue image form 1394 Firewire  camera ",
    "version",           "1.0.0",
    "vendor",            "Kyohei Iwane  Osaka.Univ",
    "category",          "Category",
    "activity_type",     "PERIODIC",
    "kind",              "DataFlowComponent",
    "max_instance",      "0",
    "language",          "C++",
    "lang_type",         "compile",
    // Configuration variables
    "conf.default.ImageWindow", "off",
    "conf.default.string_output_color_format", "RGB",
    "conf.default.string_camera_param_filename", "camera.yml",
    // Widget
    "conf.__widget__.ImageWindow", "text",
    "conf.__widget__.string_output_color_format", "text",
    "conf.__widget__.string_camera_param_filename", "text",
    // Constraints
    ""
  };
// </rtc-template>

/*!
 * @brief constructor
 * @param manager Maneger Object
 */
FireWireCamera::FireWireCamera(RTC::Manager* manager)
  // <rtc-template block="initializer">
  : RTC::DataFlowComponentBase(manager),
    m_ImageDataOut("ImageData", m_ImageData)

    // </rtc-template>
{
}

/*!
 * @brief destructor
 */
FireWireCamera::~FireWireCamera()
{
}



RTC::ReturnCode_t FireWireCamera::onInitialize()
{
  // Registration: InPort/OutPort/Service
  // <rtc-template block="registration">
  // Set InPort buffers
  
  // Set OutPort buffer
  addOutPort("ImageData", m_ImageDataOut);
  
  // Set service provider to Ports
  
  // Set service consumers to Ports
  
  // Set CORBA Service Ports
  
  // </rtc-template>

  // <rtc-template block="bind_config">
  // Bind variables and configuration variable
  bindParameter("ImageWindow", m_ImageWindow, "off");
  bindParameter("string_output_color_format", m_string_output_color_format, "RGB");
  bindParameter("string_camera_param_filename", m_string_camera_param_filename, "camera.yml");
  // </rtc-template>
  return RTC::RTC_OK;
}

/*
  RTC::ReturnCode_t FireWireCamera::onFinalize()
  {
  return RTC::RTC_OK;
  }
*/

/*
  RTC::ReturnCode_t FireWireCamera::onStartup(RTC::UniqueId ec_id)
  {
  return RTC::RTC_OK;
  }
*/

/*
  RTC::ReturnCode_t FireWireCamera::onShutdown(RTC::UniqueId ec_id)
  {
  return RTC::RTC_OK;
  }
*/


RTC::ReturnCode_t FireWireCamera::onActivated(RTC::UniqueId ec_id)
{
  // camera open 
  // 4 means DC1394_VIDEO_MODE_640x480_RGB8 : Bumblebee -> 27  DC1394_VIDEO_MODE_FORMAT7_2 
  // 正直なぜBumblebee の時がFormat7 なのかはわからない。
  if( !_bbc.openCamera( 0 ,4 )){
    std::cout << "Open Camera Error " << std::endl;
    return RTC::RTC_ERROR;
  }

  //*********************** load camera parameter file**************************
  std::cout << "Load camera parameter from " << m_string_camera_param_filename << std::endl;
  cv::FileStorage fs(m_string_camera_param_filename, cv::FileStorage::READ);
  if (!fs.isOpened()){
    RTC_ERROR(("Unable to open camera parameter file. [%s]",
			m_string_camera_param_filename.c_str()));
    std::cerr << "Could not open camera parameter." << std::endl;
    return RTC::RTC_ERROR;
  }
  //*****************************************************************************

  //**************** intrinsic dostortion  *****************
#if CV_MAJOR_VERSION == 2 && CV_MINOR_VERSION >= 2
  // std::cout << "intrinsic Parameter" << std::endl;
  fs["intrinsic"] >> cam_param.cameraMatrix;
  // std::cout << "distrotion Parameter" << std::endl;
  fs["distortion"] >> cam_param.distCoeffs;
#else  // CV_MAJOR_VERSION == 2 && CV_MINOR_VERSION >= 2
  cv::FileNode node(fs.fs, NULL);
  cv::FileNode fn = node["intrinsic"];
  cv::read(fn, cam_param.cameraMatrix);
  fn = node["distortion"];
  cv::read(fn, cam_param.distCoeffs);
#endif  // CV_MAJOR_VERSION == 2 && CV_MINOR_VERSION >= 2

  // std::cout << "camera" << std::endl;

  //********** camera matrix *****************************//
  std::cout << "camera matrix" << std::endl;
  for(int i=0;i<3;i++){
    for(int j=0;j<3;j++)
      std::cout << cam_param.cameraMatrix.at<double>(i,j) << ", " << std::ends;
    std::cout << std::endl;
  }

  //********** distortion coefficient ********************//
  std::cout << "distrotion coefficient "  << std::endl;
  for(int i=0;i<5;i++)
    std::cout << cam_param.distCoeffs.at<double>(0, i) << "," << std::ends;
  std::cout << std::endl;

  //************ intrinsic matrix ********************
  m_ImageData.data.intrinsic.matrix_element[0] = cam_param.cameraMatrix.at<double>(0, 0);
  m_ImageData.data.intrinsic.matrix_element[1] = cam_param.cameraMatrix.at<double>(0, 1);
  m_ImageData.data.intrinsic.matrix_element[2] = cam_param.cameraMatrix.at<double>(1, 1);
  m_ImageData.data.intrinsic.matrix_element[3] = cam_param.cameraMatrix.at<double>(0, 2);
  m_ImageData.data.intrinsic.matrix_element[4] = cam_param.cameraMatrix.at<double>(1, 2);

  //*********** intrinsic.distortion_coefficient.length ********//
  // std::cout << "distrotion  length " << std::endl;
  m_ImageData.data.intrinsic.distortion_coefficient.length(cam_param.distCoeffs.cols);

  /************************* channels check ***********************/
  if(m_string_output_color_format == "GRAY")
    channels = 1;
  else
    channels = 3;

  //******************
  m_ImageWidth = 640;
  m_ImageHeight = 480;

  //**************** image data length ************************//
  // std::cout << "Image data  length " << std::endl;
  m_ImageData.data.image.raw_data.length(m_ImageWidth*m_ImageHeight * channels);

  //*************** image buffer allocated *******************//
  // std::cout << "image buffer allocated " << std::endl;
  _imgL = (unsigned char*) calloc( _bbc.w * _bbc.h * channels, sizeof(unsigned char));

  _src = cvCreateImage( cvSize( _bbc.w , _bbc.h ) , IPL_DEPTH_8U , channels);

  //********************** width, height, format **********************//
  m_ImageData.data.image.width = _bbc.w;
  m_ImageData.data.image.height = _bbc.h;
  if(channels == 1)
    m_ImageData.data.image.format = Img::CF_GRAY;
  else
    m_ImageData.data.image.format = Img::CF_RGB;

  _size = _bbc.w * _bbc.h;
  cvNamedWindow( "cap");

  //******************* show paramaters *********************//
  std::cout << "camera matrix" << std::endl;
  for(int i=0;i<5;i++)
    std::cout << m_ImageData.data.intrinsic.matrix_element[i] << "," << std::ends;
  std::cout << std::endl;
  std::cout << "distrotion coefficient "  << std::endl;
  for(int i=0;i<5;i++)
    std::cout << m_ImageData.data.intrinsic.distortion_coefficient[i] << "," << std::ends;
  std::cout << std::endl;
  //std::cout << "set distrotion  parameter finish " << std::endl;
  //std::cout << "Activate finish " << std::endl;

  return RTC::RTC_OK;
}


RTC::ReturnCode_t FireWireCamera::onDeactivated(RTC::UniqueId ec_id)
{
  cvDestroyWindow( "cap");
  _bbc.closeCamera();
  
  return RTC::RTC_OK;
}


RTC::ReturnCode_t FireWireCamera::onExecute(RTC::UniqueId ec_id)
{
  //bbc.captureFrame( 'R', imgL );  // this is RGB pixel format
  //memcpy( &m_ImageData.data[0], imgL , _size) ;

  /************************* channels check ***********************/

  if(m_string_output_color_format == "GRAY" && channels == 3 ){

    channels = 1;
    m_ImageData.data.image.raw_data.length(m_ImageWidth*m_ImageHeight );
    _imgL = (unsigned char*) calloc( _bbc.w * _bbc.h , sizeof(unsigned char));
    cvReleaseImage( &_src );
    _src = cvCreateImage( cvSize( _bbc.w , _bbc.h ) , IPL_DEPTH_8U , 1);
    m_ImageData.data.image.format = Img::CF_GRAY;

  }else if(m_string_output_color_format == "RGB" && channels == 1 ){

    channels = 3;
    m_ImageData.data.image.raw_data.length(m_ImageWidth*m_ImageHeight * 3);
    _imgL = (unsigned char*) calloc( _bbc.w * _bbc.h * 3, sizeof(unsigned char));
    cvReleaseImage( &_src );
    _src = cvCreateImage( cvSize( _bbc.w , _bbc.h ) , IPL_DEPTH_8U , 3);
    m_ImageData.data.image.format = Img::CF_RGB;

  }

  //************************imagedata***************************
  if(channels == 1)
    _bbc.captureFrame( 'G', &m_ImageData.data.image.raw_data[0] );
  else
    _bbc.captureFrame( 'R', &m_ImageData.data.image.raw_data[0] );
  //************************************************************
  // std::cout << "get image data finish " << std::endl;

  //**********************intrinsic*****************************
  m_ImageData.data.intrinsic.matrix_element[0] = cam_param.cameraMatrix.at<double>(0, 0);
  m_ImageData.data.intrinsic.matrix_element[1] = cam_param.cameraMatrix.at<double>(0, 1);
  m_ImageData.data.intrinsic.matrix_element[2] = cam_param.cameraMatrix.at<double>(1, 1);
  m_ImageData.data.intrinsic.matrix_element[3] = cam_param.cameraMatrix.at<double>(0, 2);
  m_ImageData.data.intrinsic.matrix_element[4] = cam_param.cameraMatrix.at<double>(1, 2);
  //************************************************************
  // std::cout << "set intrinsic parameter finish " << std::endl;

  //**********************distortion*****************************
  for(int j = 0; j < 5; ++j)
    m_ImageData.data.intrinsic.distortion_coefficient[j] = cam_param.distCoeffs.at<double>(j, 0);
  //*************************************************************

  //*********************length*********************************
  //m_ImageData.data.intrinsic.distortion_coefficient.length(cam_param.distCoeffs.cols);
  //************************************************************


  m_ImageDataOut.write();

  if( m_ImageWindow == "on" ){

    if(channels == 1){
	 for(int i = 0 ; i < _size ; i++)
	   _src->imageData[i]  = m_ImageData.data.image.raw_data[i];
    }else if(channels == 3){
	 for(int i = 0 ; i < _size ; i++){
	   _src->imageData[3*i]   = m_ImageData.data.image.raw_data[3*i+2];
	   _src->imageData[3*i+1] = m_ImageData.data.image.raw_data[3*i+1];
	   _src->imageData[3*i+2] = m_ImageData.data.image.raw_data[3*i];
	 }
    }


    cvShowImage( "cap" , _src );
    char key = cvWaitKey(10);
    if ( key == 's'){
      static int count= 0;
      char file[80];
      sprintf( file, "image%03d.png",count);
      cvSaveImage( file, _src );
      count++;
    }
  }

  return RTC::RTC_OK;
}

/*
  RTC::ReturnCode_t FireWireCamera::onAborting(RTC::UniqueId ec_id)
  {
  return RTC::RTC_OK;
  }
*/

/*
  RTC::ReturnCode_t FireWireCamera::onError(RTC::UniqueId ec_id)
  {
  return RTC::RTC_OK;
  }
*/

/*
  RTC::ReturnCode_t FireWireCamera::onReset(RTC::UniqueId ec_id)
  {
  return RTC::RTC_OK;
  }
*/

/*
  RTC::ReturnCode_t FireWireCamera::onStateUpdate(RTC::UniqueId ec_id)
  {
  return RTC::RTC_OK;
  }
*/

/*
  RTC::ReturnCode_t FireWireCamera::onRateChanged(RTC::UniqueId ec_id)
  {
  return RTC::RTC_OK;
  }
*/



extern "C"
{
 
  void FireWireCameraInit(RTC::Manager* manager)
  {
    coil::Properties profile(firewirecamera_spec);
    manager->registerFactory(profile,
                             RTC::Create<FireWireCamera>,
                             RTC::Delete<FireWireCamera>);
  }
  
};


