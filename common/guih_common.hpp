
//
// GUIH         : GUI in Header files
// GUIH::Window : The base window (not for direct use)
//
// (c) 2006  Jaeil Choi
// last modified in Sep, 2009
//
// --------------------------------------------------------------------------
//
//   This program is free software; you can redistribute it and/or modify it
//  under the terms of the GNU General Public License as published by the
//  Free Software Foundation; either version 2 of the License or later.
//  This program is distributed in the hope that it will be useful, but 
//  WITHOUT ANY WARRANTY. See the GNU General Public License for more details.
//
// --------------------------------------------------------------------------
// 
// GUIH is a very simple, light-weight and cross-platform Graphic User
//  Interface that can visualize images and OpenGL graphics with
//  keyboard and mouse interactions on MS Windows, Apple OS X, and Linux.
// It is intended as a tool for developers and researchers, who don't want
//  to spend time in designing/creating an user interface, but who want 
//  all the features of a general user interface and the freedom to modify. 
//  If you are not satisfied, like me, with GLUT, GLUI, FLTK or HighGUI,
//  you may find GUIH very useful and convenient.
// Key features include dialog windows, timers, default key bindings, 
//  image capture, captions, mapping of mouse click to world-space coordinates,
//  and the built-in cameras for 2D/3D OpenGL.
// GUIH doesn't support a window with controls, except simple input dialogs.
//  For complex settings/configuration of an application, I recommend using
//  a configuration file, rather than a window full of controls.
//
// GUIH package consists of eight C++ header files:
//   guih_common.hpp     : the base window   GUIH::Window
//   guih_image.hpp      : image window      GUIH::ImageWindow
//   guih_opengl2d.hpp   : 2D OpenGL Window  GUIH::OpenGL2DWindow
//   guih_opengl3d.hpp   : 3D OpenGL Window  GUIH::OpenGL3DWindow
//   guih_camera2d.hpp   : 2D camera for     GUIH::OpenGL2DWindow
//   guih_camera3d.hpp   : 3D camera for     GUIH::OpenGL3DWindow
//   guih_args.hpp       : command-line argument processing tool (GUIH::Args)
//   guih_config.hpp     : configuration file processing tool    (GUIH::Config)
// GUIH package defines four window classes:
//   Window              : the base window (Not for direct use!)
//    +- ImageWindow     : for image visualization
//    +- OpenGL2DWindow  : for 2D OpenGL rendering (with built-in 2D camera)
//    +- OpenGL3DWindow  : for 3D OpenGL rendering (with built-in 3D camera)
// GUIH package defines two extra classes:
//   Args                : command-line argument processor (in guih_args.hpp)
//   Config              : configuration file processor (in guih_config.hpp)
//
// To compile on Windows, 
//   - use Visual Studio 8.0 or higher (for GDI+)
//   - create an empty 'Win32 Console Application' project
//   - add 'src' directory to Project->Properties->C/C++->General->AdditionalIncludeDirectories
//   - add '_CRT_SECURE_NO_DEPRECATE' to Project->Properties->C/C++->Preprocessor->PreprocessorDefinitions
//   - add 'opengl32.lib glu32.lib gdiplus.lib' to Project->Properties->Linker->Input->AddtionalDependencies
// To compile on Apple OS X,
//   - install 'GLUT for Mac OS X' (for rendering caption on OpenGL windows), and
//   - use g++ compiler flags : '-framework Carbon -framework AGL -framework OpenGL -framework GLUT'
// To compile on Linux,
//   - install GTKGLExt, and add `pkg-config --cflag --libs gtkglext' to Makefile flags
//
// --------------------------------------------------------------------------
//
// GUIH::Window
//   This class serves as the parent class of all the other window classes,
//   and it's not meant to be used directly by user.
//   It provides common functionalities such as basic window operations, 
//   timer, callback setting, screen capture, input/message dialog, etc.
//
//   Member functions:
//	void openWindow(int w, int h, char *title = NULL, 
// 			void (*init) (Window *w) = NULL,
// 			void (*draw) (Window *w) = NULL,
// 			void (*keys) (Window *w, int k) = NULL, 
// 			void (*mouse)(Window *w, int button, int state, int xy[]) = NULL,
// 			void (*close)(Window *w) = NULL);
//	void closeWindow (void);
//	bool isClosed(void);
//	void setCallbacks(draw, keys, mouse, close );
//	void setTitle(const char *title);
//	void move   (int  x, int  y);
//	void move   (Window *w, bool below);
//	void resize (int  w, int  h, bool interior);
//	void getPos (int *x, int *y);
//	void getSize(int *w, int *h, bool interior);
//	void redraw  (void);
//	bool isActive(void);
//	void present (void);
//	void show    (void);
//	void hide    (void);
//	void iconify (void);
//	void maximize(void);
//	void unmaximize(void);
//	void captureScreenInRGBBuffer(unsigned char buffer[]);
//	void captureScreen(char *filename);
//	void setCaption (char *str=NULL, int r=210, int g=210, int b=255);
//	void createTimer(double interval_sec, void (*ftn)(Window *w), bool activate_now=true);
//	void removeTimer(void);
//	bool toggleTimer(void);
//      void openMessageDialog(char *title, char *message);
//      bool openInputDialog  (char *title, UserInput *u0, UserInput *u1, UserInput *u2);
//      bool openInputDialogScale(char *title, UserInput *u);
//      void setColorForDrawing(float r, float g, float b);
//      void drawLine(int x0, int y0, int x1, int y1);
//      void drawCross(int x, int y, int size);
//      void drawRectangle(int x, int y, int size, bool filled=false);  // centered at (x,y)
//      void drawRectangle(int x, int y, int w, int h, bool filled=false); // starting at (x,y)
// 
//      void runMainLoop(void);		// start GUI main loop
//
// For API example, read the example code at the bottom.
//


#ifndef GUIH_COMMON_HPP
#define GUIH_COMMON_HPP


namespace GUIH {
  
typedef enum { ESC=27, F1=201, F2, F3, F4, F5, F6, F7, F8, F9, F10, F11, F12, 
	       ARROW_UP, ARROW_DOWN, ARROW_LEFT, ARROW_RIGHT, PAGE_UP, PAGE_DOWN,
	       DELETE_KEY, INSERT_KEY, HOME_KEY, END_KEY, ESCAPE_KEY } key_t;
typedef enum { MOUSE_NONE=0, MOUSE_LEFT, MOUSE_MIDDLE, MOUSE_RIGHT } mouse_button_t;  
typedef enum { MOUSE_CLICKED, MOUSE_DRAGGED } mouse_state_t;  
typedef enum { PIXEL_UNKNOWN=0,	// default value (not meant for use)
	       PIXEL_GRAY,	// grayscale		(1 channel  per pixel)
	       PIXEL_GRAYA, 	// grayscale + alpha	(2 channels per pixel)
	       PIXEL_RGB, 	// RGB			(3 channels per pixel)
	       PIXEL_RGBA,  	// RGB + alpha		(4 channels per pixel)
	       PIXEL_BGR, 	// BGR			(3 channels per pixel)
	       PIXEL_YUV   	// YUV			(3 channels per pixel)
}  pixel_t;
  
  
// ===================================================================
// data structure for user-input dialog functions -- openInputDialog()
// ===================================================================
  
class UserInput 
{
public:  // NOT FOR DIRECT ACCESS FROM USER PROGRAMS
  char     type;                // data type ('t','f','d','i');
  char     label[80];		// label for this input entry
  union {
    char   text[80];
    float  fvalues[3];		// [0] input_value, [1] min, [2] max
    double dvalues[3];		// [0] input_value, [1] min, [2] max
    int    ivalues[3];		// [0] input_value, [1] min, [2] max
  };
  int      len;			// size of the input field(entry)
public:  // FOR USER PROGRAMS
  UserInput()  { memset(this, 0, sizeof(UserInput)); }
  ~UserInput() {}
  char* getText(void) { return text; }
  void  setText(const char *label, const char *default_str) {
    strncpy(this->label, (label ? label : ""), 79);
    type = 't';
    if (default_str) strcpy( text, default_str ); else text[0]='\0';
  }
  int*  getInteger(void) { return ivalues; }
  void  setInteger(const char *label, int value, int min, int max) {
    strncpy(this->label, (label ? label : ""), 79);
    type = 'i';  ivalues[0] = value;  ivalues[1] = min;  ivalues[2] = max;
  }
  float* getFloat(void) { return fvalues; }
  void   setFloat(const char *label, float value, float min, float max) {
    strncpy(this->label, (label ? label : ""), 79);
    type = 'f';  fvalues[0] = value;  fvalues[1] = min;  fvalues[2] = max;
  }
  double* getDouble(void) { return dvalues; }
  void    setDouble(const char *label, double value, double min, double max) {
    strncpy(this->label, (label ? label : ""), 79);
    type = 'd';  dvalues[0] = value;  dvalues[1] = min;  dvalues[2] = max;
  }
};
  

}	// end of GUIH namespace


// ===================================================================
#ifdef  WIN32   // beginning of Window version
// ===================================================================

#include <iostream>
#include <windows.h>
#include <windowsx.h>
#include <gdiplus.h>

#ifndef GDIPLUS_DATA
#define GDIPLUS_DATA
static ULONG_PTR			gdiplusToken;
static Gdiplus::GdiplusStartupInput	gdiplusStartupInput;
static bool				gdiplusStarted = false;
static bool getEncoderClsid(char* fname, CLSID *pClsid)
{
	if (fname == NULL) return false;
	int len = (int)strlen(fname);
	if (len < 4 || fname[len-4] != '.') return false;
	wchar_t *format;
	char *ext = fname + strlen(fname)-3;
	if      (strcmp(ext, "png")==0 || strcmp(ext, "PNG")==0) format = L"image/png";
	else if (strcmp(ext, "jpg")==0 || strcmp(ext, "JPG")==0) format = L"image/jpg";
	else if (strcmp(ext, "bmp")==0 || strcmp(ext, "BMP")==0) format = L"image/bmp";
	else return false;
	UINT  num = 0;          // number of image encoders
	UINT  size = 0;         // size of the image encoder array in bytes
	Gdiplus::ImageCodecInfo* pImageCodecInfo = NULL;
	Gdiplus::GetImageEncodersSize(&num, &size);
	if(size == 0) return false;
	pImageCodecInfo = (Gdiplus::ImageCodecInfo*)(malloc(size));
	if(pImageCodecInfo == NULL) return false;
	GetImageEncoders(num, size, pImageCodecInfo);
	for(UINT j = 0; j < num; ++j) 
	  if( wcscmp(pImageCodecInfo[j].MimeType, format) == 0 ) {
		 *pClsid = pImageCodecInfo[j].Clsid;
		 free(pImageCodecInfo);  return true;
	  }    
	free(pImageCodecInfo);
	return false;
}
#endif //GDIPLUS_DATA


namespace GUIH {
  
static LRESULT CALLBACK GUIHWinProc(HWND hwnd, UINT uMsg, WPARAM wParam, LPARAM lParam);
static int  guih_window_count = 0;	// total number of opened windows
static int  guih_serial_number = -1;
static void (*guih_idle_ftn)(void *data) = NULL;
static void *guih_idle_data = NULL;
  
#if defined WIN64 || defined EM64T // --------------------------------
#define get_window_long_ptr	GetWindowLongPtr
#define set_window_long_ptr(hwnd, id, ptr) setWindowLongPtr(hwnd, id, (LONG_PTR)(ptr))
#define get_class_long_ptr	GetClassLongPtr
#define GUIH_USERDATA		GWLP_USERDATA
#define GUIH_WNDPROC		GWLP_WNDPROC
#define GUIH_HCURSOR		GCLP_HCURSOR
#define GUIH_HBRBACKGROUND	GCLP_HBRBACKGROUND
#else // -------------------------------------------------------------
#define get_window_long_ptr	GetWindowLong
#define set_window_long_ptr(hwnd, id, ptr) SetWindowLong(hwnd, id, (LONG)(size_t)ptr)
#define get_class_long_ptr	GetClassLong
#define GUIH_USERDATA		GWL_USERDATA
#define GUIH_WNDPROC		GWL_WNDPROC
#define GUIH_HCURSOR		GCL_HCURSOR
#define GUIH_HBRBACKGROUND	GCL_HBRBACKGROUND
#endif // ------------------------------------------------------------

// -------------------------------------------------------------------

class Window
{
public:
  HWND     hwnd;	// window frame
  wchar_t  caption_text[512];
  unsigned char   caption_rgb[3];
  bool   timer_paused;
  double timer_interval; // in seconds
  int    serial_number;	// unique window index, which begins from 1
  int    penR, penG, penB;
  int    *window_count;

public:
  bool  shift, ctrl, alt, skey_by_user;
  void  (*cb_init)   (Window *w);
  void  (*cb_draw)   (Window *w);
  void  (*cb_keys)   (Window *w, int k);
  void  (*cb_mouse)  (Window *w, int button, int state, int xy[]);
  void  (*cb_timer)  (Window *w);
  void  (*cb_close)  (Window *w);
  
public:
  // Window specific callback functions. These virtual functions should be replaced.
  virtual bool runOnInit(void) { if (cb_init) cb_init(this); return true; }
  virtual void runOnCloseReq(void) { }  // this can be cancelled
  virtual void runOnDestroy(void) { if (cb_close) cb_close(this); }  // This function is called automatically by a event, and cannot be cancelled
  virtual bool runOnDraw(void) { if (cb_draw) cb_draw(this); return true; }
  virtual void runOnKeyboard(int key) { if (key > 0 && cb_keys) cb_keys( this, key ); }
  virtual void runOnMouse(int button, int state, int xy[]) { if (cb_mouse) cb_mouse( this, button, state, xy ); }
  bool runOnTimer(void) { if (cb_timer == NULL || timer_paused) return false;  cb_timer(this);  return true; }
  
public:
  Window () { 
    hwnd = 0;   window_count = &guih_window_count;
    cb_init=NULL;     cb_draw=NULL;   cb_keys=NULL;   
    cb_mouse = NULL;  cb_timer=NULL;  cb_close=NULL;
    timer_paused = true;  shift = ctrl = alt = false;  skey_by_user = true;
    serial_number = ++guih_serial_number;
	penR = penG = penB = 0;
	if (!gdiplusStarted) {
		Gdiplus::GdiplusStartup( &gdiplusToken, &gdiplusStartupInput, NULL );
		gdiplusStarted = true;
	}
	caption_text[0] = '\0';
  }
  Window (int w, int h, const char *title = NULL, 
	  void (*init) (Window *w) = NULL,
	  void (*draw) (Window *w) = NULL,
	  void (*keys) (Window *w, int k) = NULL, 
	  void (*mouse)(Window *w, int button, int state, int xy[]) = NULL,
	  void (*close)(Window *w) = NULL) { 
    hwnd = 0;   window_count = &guih_window_count;
    cb_init=NULL;  cb_draw=NULL;     cb_keys=NULL;   cb_mouse = NULL;
    cb_close=NULL;   cb_timer=NULL;  timer_paused = true; 
    shift = ctrl = alt = 0;  skey_by_user = true;
    serial_number = ++guih_serial_number;
	if (!gdiplusStarted) {
		Gdiplus::GdiplusStartup( &gdiplusToken, &gdiplusStartupInput, NULL );
		gdiplusStarted = true;
	}
    openWindow(w, h, title, init, draw, keys, mouse, close); 
	caption_text[0] = '\0';
  }
  ~Window() { closeWindow(); }
public:
  virtual void openWindow(int w, int h, const char *title = NULL, 
			  void (*init) (Window *w) = NULL,
			  void (*draw) (Window *w) = NULL,
			  void (*keys) (Window *w, int k) = NULL, 
			  void (*mouse)(Window *w, int button, int state, int xy[]) = NULL,
			  void (*close)(Window *w) = NULL) {
    if (w <= 0) w = 200;  if (h <= 0) h = 200;
    closeWindow();
    setCallbacks(draw, keys, mouse, close);
    createBasicObjects( w, h, title );
  }
  virtual void closeWindow(void) { 
    // Note that this is just one of several ways of closing window.
    // Therefore, actual clean-up should be done in 'runOnDestroy()'.
    if (hwnd) { DestroyWindow( hwnd ); hwnd = NULL; }
  }
  
  void createBasicObjects(int w, int h, const char *title) {
    if (hwnd) closeWindow();
    if (!title) title = "Untitled";
    static bool bWCRegistered = false;
    if (!bWCRegistered) {
      WNDCLASS wndc;
      wndc.style = CS_OWNDC | CS_VREDRAW | CS_HREDRAW;
      wndc.lpfnWndProc   = GUIHWinProc;
      wndc.cbClsExtra    = wndc.cbWndExtra   = 0;
      wndc.hInstance     = 0;
      wndc.hIcon         = LoadIcon  (NULL, IDI_APPLICATION);
      wndc.hCursor       = LoadCursor(NULL, IDC_CROSS);
      wndc.hbrBackground = (HBRUSH)GetStockObject(GRAY_BRUSH);
      wndc.lpszClassName = L"GUIH_MainWindow";
      wndc.lpszMenuName  = NULL;
      RegisterClass(&wndc);
      bWCRegistered = true;
    }
    // create the window frame
    DWORD defStyle = WS_VISIBLE | WS_MINIMIZEBOX | WS_MAXIMIZEBOX | WS_SYSMENU;
	wchar_t wbuf[512];
	mbstowcs( wbuf, title, (int)strlen(title)+1 );
    hwnd = CreateWindow( L"GUIH_MainWindow", (LPCWSTR)wbuf, defStyle | WS_OVERLAPPEDWINDOW,
			 CW_USEDEFAULT, CW_USEDEFAULT, w, h, NULL, NULL, NULL, NULL );
    ShowWindow( hwnd, SW_SHOW );
    // set (Window*) pointer in the window handle
    set_window_long_ptr( hwnd, GUIH_USERDATA, this );
  }
  // register user callback functions --------------------------------
  void setCallbacks( void (*draw) (Window *w),
		     void (*keys) (Window *w, int k), 
		     void (*mouse)(Window *w, int button, int state, int xy[]) = NULL,
		     void (*close)(Window *w) = NULL ) {
    cb_draw = draw;  cb_keys = keys;  cb_mouse = mouse;  cb_close = close;
  }
  
public:
  int  getSerialNumber(void) { return serial_number; }
  bool isClosed(void)      { return (hwnd == NULL); }
  
  void setTitle(const char *title) { 
    if (!hwnd) return;
    wchar_t wbuf[512];
    mbstowcs( wbuf, title, (int)strlen(title)+1 );
    SetWindowText( hwnd, wbuf); // (LPCWSTR)title );
  }
  void move   (int x, int y)   { 
    if (!hwnd) return;  RECT rect;  GetWindowRect( hwnd, &rect );
    MoveWindow( hwnd, x, y, rect.right - rect.left, rect.bottom - rect.top, TRUE);
  }
  void move   (Window *win, bool below) {
    if (!win || win->isClosed()) return;
    int x,y,w,h; win->getPos(&x,&y); win->getSize(&w,&h,false); 
    if (below) move(x,y+h); else move(x+w,y); 
  }
  void resize (int w, int h, bool interior=false) {
    if (!hwnd) return; 
    RECT rect, client; 
    GetWindowRect( hwnd, &rect );
    if (!interior) {	// set the size of the window
      MoveWindow( hwnd, rect.left, rect.top, w, h, TRUE);
    } else {		// set the size of the client area
      GetClientRect( hwnd, &client );
      int side_margin   = (rect.right - rect.left) - (client.right - client.left);
      int title_margin  = client.top;
      int bottom_margin = (rect.bottom - rect.top) - client.bottom;
      MoveWindow( hwnd, rect.left, rect.top, w + side_margin, h + title_margin + bottom_margin, TRUE);
    }
  }
  void getPos (int *x, int *y) { 
    if (!hwnd) return;  RECT rect;  GetWindowRect( hwnd, &rect );
    *x = rect.left;  *y = rect.top;
  }
  void getSize(int *w, int *h, bool interior=true) {
    if (!hwnd) return;  RECT rect;  
    if (!interior) GetWindowRect( hwnd, &rect );  // get the size of the window
    else           GetClientRect( hwnd, &rect );  // get the size of the client area
    *w = rect.right - rect.left;  *h = rect.bottom - rect.top;
  }
  void redraw  (void) { if (hwnd) InvalidateRect( hwnd, 0, 0 ); }
  void present (void) { if (hwnd) SetFocus(hwnd); }
  //   bool isActive(void) { return (window && (bool)gtk_window_is_active( GTK_WINDOW(window) )); }
  //   void show    (void) { if (window) gtk_widget_show( window ); }
  //   void hide    (void) { if (window) gtk_widget_hide( window ); }
  //   void iconify (void) { if (window) gtk_window_iconify( GTK_WINDOW(window) ); }
  //   void maximize(void)   { if (window) gtk_window_maximize  ( GTK_WINDOW(window) ); }
  //   void unmaximize(void) { if (window) gtk_window_unmaximize( GTK_WINDOW(window) ); }
  bool captureScreenInRGBBuffer(unsigned char *buffer) { 
    if (!hwnd) return false;
	RECT rect; GetClientRect(hwnd, &rect);
	Gdiplus::Bitmap bitmap(abs(rect.right-rect.left), abs(rect.top-rect.bottom));
	Gdiplus::Graphics *g = Gdiplus::Graphics::FromImage( &bitmap );
	HDC hdc = g->GetHDC();
	HDC sdc = GetDC(hwnd);
	BitBlt( hdc, 0, 0, bitmap.GetWidth(), bitmap.GetHeight(), sdc, 0, 0, SRCCOPY );
	ReleaseDC( hwnd, sdc );
	g->ReleaseHDC( hdc );
	delete( g );
	Gdiplus::Rect rct(0, 0, bitmap.GetWidth(), bitmap.GetHeight());
	Gdiplus::BitmapData bitmapData;
	bitmap.LockBits( &rct, Gdiplus::ImageLockModeRead, PixelFormat24bppRGB, &bitmapData );
	int   row, rowstride = bitmap.GetWidth() * sizeof(unsigned char) * 3;
	char *src = (char*)bitmapData.Scan0, *dst = (char*)buffer;
	for (row = 0; row < (int)bitmap.GetHeight(); row++) {
		memcpy( dst, src, rowstride );
		src += bitmapData.Stride;  dst += rowstride;
	}
	bitmap.UnlockBits( &bitmapData );
	return true;
  }
  bool captureScreen(char *fname=NULL) { 
    if (!hwnd) return false;
	if (fname == NULL) fname = "guih_capture.png";
	RECT rect; GetClientRect(hwnd, &rect);
	Gdiplus::Bitmap bitmap(abs(rect.right-rect.left), abs(rect.top-rect.bottom));
	Gdiplus::Graphics *g = Gdiplus::Graphics::FromImage( &bitmap );
	HDC hdc = g->GetHDC();
	HDC src = GetDC(hwnd);
	BitBlt( hdc, 0, 0, bitmap.GetWidth(), bitmap.GetHeight(), src, 0, 0, SRCCOPY );
	ReleaseDC( hwnd, src );
	g->ReleaseHDC( hdc );
	delete( g );
    wchar_t wfname[512]; mbstowcs( wfname, fname, (int)strlen(fname)+1 );
    CLSID encoderClsid;  if (!getEncoderClsid( fname, &encoderClsid )) return false;
    Gdiplus::Status stat = bitmap.Save( wfname, &encoderClsid, NULL );
	return (stat==0 ? true : false);
  }
  // caption ---------------------------------------------------------
  void setCaption (char *str=NULL, int r=210, int g=210, int b=255) { 
    mbstowcs( caption_text, (str ? str : (char*)""), (int)(str ? strlen(str)+1 : 1) );
    caption_rgb[0]=(unsigned char)r; 
    caption_rgb[1]=(unsigned char)g; 
    caption_rgb[2]=(unsigned char)b;
  }
  // timer -----------------------------------------------------------
  void createTimer(double interval_sec, void (*ftn)(Window *w), bool activate_now=true) {
    cb_timer = ftn;  timer_interval = interval_sec;
    if (activate_now) {
      timer_paused = false;
      SetTimer( hwnd, 100+serial_number, (int)(interval_sec*1000), NULL );
    }
  }
  void removeTimer(void) { timer_paused = true;  KillTimer(hwnd, 100+serial_number); }
  bool toggleTimer(double interval_sec=0, void (*ftn)(Window *w)=NULL) { 
    if (interval_sec > 0) timer_interval = interval_sec;
    if (ftn) cb_timer = ftn;
    if (timer_paused && cb_timer) {
      timer_paused = false; 
      SetTimer( hwnd, 100+serial_number, (int)(timer_interval*1000), NULL );
      return true;
    } else {
      timer_paused = true;  KillTimer(hwnd, 100+serial_number); 
      return false;
    }
  }
  bool isTimerPaused(void) { return timer_paused; }
  
  // -----------------------------------------------------------------
  int remapKeyValue(WPARAM oldkey, bool *skey=NULL) {
    int  key;
    if      (oldkey >= (WPARAM)'A'  && oldkey <= (WPARAM)'Z') { if (shift) key = 'A' + (int)(oldkey - (WPARAM)'A'); else key = 'a' + (int)(oldkey - (WPARAM)'A'); }
    else if (oldkey >= VK_F1 && oldkey <= VK_F12) key = F1 + (int)(oldkey - (WPARAM)112);
    else if (oldkey == VK_ESCAPE) { if (skey) *skey=true; key = ESCAPE_KEY; }
    else if (oldkey == VK_UP   )  { if (skey) *skey=true; key = ARROW_UP; }
    else if (oldkey == VK_DOWN )  { if (skey) *skey=true; key = ARROW_DOWN; }
    else if (oldkey == VK_LEFT )  { if (skey) *skey=true; key = ARROW_LEFT; }
    else if (oldkey == VK_RIGHT)  { if (skey) *skey=true; key = ARROW_RIGHT; }
    else if (oldkey == VK_HOME)   { if (skey) *skey=true; key = HOME_KEY; }
    else if (oldkey == VK_END)    { if (skey) *skey=true; key = END_KEY; }
    else if (oldkey == VK_PRIOR)  { if (skey) *skey=true; key = PAGE_UP; }
    else if (oldkey == VK_NEXT)   { if (skey) *skey=true; key = PAGE_DOWN; }
    else key = 0;  // key input with zero value will be ignored
    return key;
  }
  int remapNumericKey2Integer(char key, bool add10=false) {
    // Remap numeric keys "0123456789" and symbols ")!@#$%^&*(", so that we get [0,1,...19].
    int val = -1;
    if (key >= '0' && key <= '9') val = key - '0';
    else if (key==')') val = 10;  else if (key=='!') val = 11;
    else if (key=='@') val = 12;  else if (key=='#') val = 13;
    else if (key=='$') val = 14;  else if (key=='%') val = 15;
    else if (key=='^') val = 16;  else if (key=='&') val = 17;
    else if (key=='*') val = 18;  else if (key=='(') val = 19;
    return (val<0 ? -1 : (val + (add10 ? 10 : 0)));
  }
  
  // user input dialog -----------------------------------------------
  bool openMessageDialog(char *title, char *message, bool yesno=false) { return true; }
  bool openInputDialog(char *title, UserInput *ui0, UserInput *ui1=NULL, UserInput *ui2=NULL) { return true; }
  bool openInputDialogScale(char *title, UserInput *ui) { return true; }
  
  // drawing ---------------------------------------------------------
  void setColorForDrawing(float r, float g, float b) { 
    penR = (int)(255*r);  penG = (int)(255*g);  penB = (int)(255*b);
  }
  void drawLine(int x0, int y0, int x1, int y1) { 
    Gdiplus::Graphics graphics(hwnd);
    Gdiplus::Pen pen(Gdiplus::Color(penR, penG, penB));
    graphics.DrawLine( &pen, x0, y0, x1, y1 );
  }
  void drawRectCenteredAt(int x, int y, int size, bool filled=false) { 
    // Draw a rectangle centered at (x,y)
    Gdiplus::Graphics graphics(hwnd);
    Gdiplus::Pen pen(Gdiplus::Color(penR, penG, penB));
    graphics.DrawRectangle( &pen, x-size/2, y-size/2, size, size );
  }
  void drawRect(int x, int y, int w, int h, bool filled=false) {
    // Draw a rectangle starting at (x,y) 
    Gdiplus::Graphics graphics(hwnd);
    Gdiplus::Pen pen(Gdiplus::Color(penR, penG, penB));
    graphics.DrawRectangle( &pen, x, y, w, h );
  }
  void drawCross(int x, int y, int size, char type='+') { 
    if        (type=='+') {
      drawLine(x-size, y, x+size, y); drawLine(x, y-size, x, y+size); 
    } else if (type=='x') {
      drawLine(x-size, y-size, x+size, y+size); drawLine(x+size, y-size, x-size, y+size); 
    } else if (type=='|') {
      int hs = (size/2);
      drawLine(x-hs, y, x+hs, y); drawLine(x, y-size, x, y+size); 
    } else if (type=='-') {
      int hs = (size/2);
      drawLine(x-size, y, x+size, y); drawLine(x, y-hs, x, y+hs); 
    } else {
      int dx=0, dy=0;
      switch (type) {
      case '0':  dx = +size;  dy = -size;  break;
      case '1':  dx = -size;  dy = -size;  break;
      case '2':  dx = -size;  dy = +size;  break;
      case '3':  dx = +size;  dy = +size;  break;
      }
      drawLine( x, y,  x+dx, y );   drawLine( x, y,  x, y+dy );
    }
  }
  
  // -----------------------------------------------------------------
  // Non window-specific functions
  // -----------------------------------------------------------------
  void runMainLoop(void) {
    //for (int i=0; i<win_moves; i++) // arrange windows in relative position
    //win_moves[i].tgt->move( win_moves[i].ref, win_moves[i].below );
    int *window_count = this->window_count;  // points to 'guih_window_count'
    this->present();
    MSG msg;
    while (true) {
      if (GetMessage( &msg, NULL, 0, 0 )) {
	TranslateMessage(&msg);	//  Translate keyboard-specific messages
	DispatchMessage(&msg);	//  Send it off to the window procedure
      }
      if (*window_count <= 0) break;
    }
    Gdiplus::GdiplusShutdown( gdiplusToken );
  }

  // idle function
  void addIdleFunction(void (*idle)(void* data), void* data=NULL) { guih_idle_ftn = idle; guih_idle_data = data; }
  void removeIdleFunction(void* data=NULL) { guih_idle_ftn = NULL; guih_idle_data = NULL; }

};

// -------------------------------------------------------------------
// The main message handler
// -------------------------------------------------------------------

static LRESULT CALLBACK GUIHWinProc( HWND hwnd, UINT uMsg, WPARAM wParam, LPARAM lParam )
{
  Window* w = (Window*)(size_t)get_window_long_ptr( hwnd, GUIH_USERDATA );
  //if( !window ) return DefWindowProc(hwnd, uMsg, wParam, lParam);
  switch(uMsg) {
  case WM_CREATE:
    guih_window_count++; 
    PostMessage( hwnd, WM_USER+1, wParam, lParam );
    break;
  case WM_USER+1:  w->runOnInit();  break;	// WM_CREATE
  case WM_DESTROY:
    w->runOnDestroy();
    guih_window_count--;
    w->hwnd = NULL;  
    break;
    //case WM_CLOSE:   break;  // close request ?
  case WM_PAINT: // ==================================================
    w->runOnDraw();
    {
      Gdiplus::FontFamily fontFamily(L"Times New Roman");
      Gdiplus::Font       font(&fontFamily, 14, Gdiplus::FontStyleRegular, Gdiplus::UnitPixel);
      RECT rect; GetClientRect(hwnd, &rect);
      Gdiplus::PointF     fpoint(5.0f, rect.bottom - 20.0f), bpoint(6.0f, rect.bottom - 19.0f);
      Gdiplus::SolidBrush whiteBrush(Gdiplus::Color(255, 255, 255, 255));
      Gdiplus::SolidBrush blackBrush(Gdiplus::Color(255, 0, 0, 255));
      Gdiplus::Graphics graphics(hwnd);
      graphics.DrawString( w->caption_text, -1, &font, bpoint, &blackBrush );
      graphics.DrawString( w->caption_text, -1, &font, fpoint, &whiteBrush );
    }
    return DefWindowProc(hwnd, uMsg, wParam, lParam);
    break;
  case WM_LBUTTONDOWN:  case WM_RBUTTONDOWN:  case WM_MOUSEMOVE: 
  case WM_LBUTTONUP:    case WM_RBUTTONUP:   
    {
      static int mhxy[802], mhpos=0;
      int button = ((uMsg==WM_LBUTTONDOWN || uMsg==WM_LBUTTONUP) ? MOUSE_LEFT :
		    (uMsg==WM_RBUTTONDOWN || uMsg==WM_RBUTTONUP) ? MOUSE_RIGHT : MOUSE_NONE);
      int x = (int)GET_X_LPARAM(lParam);
      int y = (int)GET_Y_LPARAM(lParam);
      if        (uMsg==WM_LBUTTONDOWN || uMsg==WM_RBUTTONDOWN) {
	mhxy[0] = x;  mhxy[1] = y;  mhpos = 2;
      } else if (uMsg==WM_MOUSEMOVE && mhpos>0) {
	if (mhpos<800-1) { mhxy[mhpos++] = x;  mhxy[mhpos++] = y; }
      } else if (uMsg==WM_LBUTTONUP   || uMsg==WM_RBUTTONUP) {
	mhxy[mhpos+0] = -1;  mhxy[mhpos+1] = -1;
	if (mhpos==2) w->runOnMouse( button, MOUSE_CLICKED, mhxy );
	else          w->runOnMouse( button, MOUSE_DRAGGED, mhxy );
	mhpos = 0;  w->redraw();
      }
    }
    break;
  case WM_CHAR: 
    w->runOnKeyboard( (int)wParam ); w->redraw(); 
    break;
  case WM_KEYDOWN:  
    if      (wParam == VK_CONTROL) w->ctrl  = true;
    else if (wParam == VK_SHIFT)   w->shift = true;
    else if (wParam == VK_MENU)    w->alt   = true;
    else if (wParam >= 'A' && wParam <= 'Z') break;
    else if (wParam >= '0' && wParam <= '9') break;
    else {
      int key = remapKeyValue(wParam);
      w->runOnKeyboard( key );
      switch (key) {
      case ESCAPE_KEY:  w->closeWindow();  return TRUE;
      case INSERT_KEY:  w->captureScreen(); break;
      case DELETE_KEY: case ' ':  w->setCaption(NULL); break;
      default: break;
      }
      w->redraw();
    }
    break; 
  case WM_KEYUP:  // =================================================
    if      (wParam == VK_CONTROL) w->ctrl  = false; 
    else if (wParam == VK_SHIFT)   w->shift = false; 
    else if (wParam == VK_MENU)    w->alt   = false; 
    break;
  case WM_SIZE: PostMessage( hwnd, WM_USER+2, wParam, lParam ); break;
  case WM_TIMER:  w->runOnTimer();   break;
  //case WM_IDLE: if (guih_idle_ftn) guih_idle_ftn( guih_idle_data ); break;
  }
  return DefWindowProc(hwnd, uMsg, wParam, lParam);
}
  
}  // end of namespace GUIH

// ===================================================================
#elif defined(__APPLE__) || defined(MACOSX)	// beginning of OS X version
// ===================================================================

#include <iostream>
#include <Carbon/Carbon.h>
#include <ApplicationServices/ApplicationServices.h>

namespace GUIH {

static int      guih_window_count = 0;	// total number of opened windows
static int      guih_serial_number = -1;
static bool     guih_app_ready = false;  
static bool     guih_mbar_ready = false;
static EventLoopTimerRef  idle_ref=NULL;  
static EventLoopIdleTimerUPP  idle_upp=NULL;  
static void (*guih_idle_ftn)(void *data) = NULL;
  
static pascal void TimerAction(EventLoopTimerRef  theTimer, void* userData);
static pascal void IdleTimerAction (EventLoopTimerRef  theTimer, EventLoopIdleTimerMessage idleMsg, void* userData);
static OSStatus AppEventHandler(EventHandlerCallRef nextHandler, EventRef event, void *userData);
static OSStatus WindowEventHandler(EventHandlerCallRef nextHandler, EventRef event, void *userData);
  
class Window
{
public:
  WindowRef	window;
  CGContextRef	windowContext;	
  int		serial_number;	// unique window index, which begins from 1
  bool		timer_paused;
  double	timer_interval; // in seconds
  EventLoopTimerRef  timer_ref;
  EventLoopTimerUPP  timer_upp;
  
  char		caption[256];	// caption
  // user callbacks
  void  (*cb_init)   (Window *w);
  void  (*cb_draw)   (Window *w);
  void  (*cb_keys)   (Window *w, int k);
  void  (*cb_mouse)  (Window *w, int button, int state, int xy[]);
  void  (*cb_close)  (Window *w);
  void  (*cb_timer)  (Window *w);
  
public:
  bool  shift, ctrl, alt, cmmd, skey_by_user;  // user can access these values
  
  Window *move_ref;
  bool   move_below;
  
public:
  Window () : window(NULL), windowContext(NULL), timer_paused(true), timer_ref(NULL), timer_upp(NULL), 
	      cb_init(NULL), cb_draw(NULL), cb_keys(NULL), cb_mouse(NULL), cb_close(NULL), cb_timer(NULL), 
	      shift(0), ctrl(0), alt(0), cmmd(0), skey_by_user(false) {
    caption[0] = '\0'; 
    serial_number = ++guih_serial_number;
    move_ref=NULL; 
    if (!guih_app_ready) {
      // There is a strange bug in Carbon and CoreGraphics in OS X, when used by multithread.
      // You must execute 'setupOSXAppEventHandler()' before any 'CGCreatImage..()' function.
      // Otherwise, 'RunApplicationEventLoop()' may fail and falls to a halt.
      setupOSXAppEventHandler();
    }
  }
  Window (int w, int h, const char *title = NULL, 
	  void (*init) (Window *w) = NULL,
	  void (*draw) (Window *w) = NULL,
	  void (*keys) (Window *w, int key) = NULL, 
	  void (*mouse)(Window *w, int button, int state, int xy[]) = NULL,
	  void (*close)(Window *w) = NULL) 
    : window(NULL), windowContext(NULL), timer_paused(true), timer_ref(NULL), timer_upp(NULL), 
      cb_init(NULL), cb_draw(NULL), cb_keys(NULL), cb_mouse(NULL), cb_close(NULL), cb_timer(NULL), 
      shift(0), ctrl(0), alt(0), cmmd(0), skey_by_user(false) { 
    caption[0] = '\0'; 
    serial_number = ++guih_serial_number;
    openWindow(w, h, title, init, draw, keys, mouse, close); 
  }
  virtual ~Window() { closeWindow(); }
  
public:
  // Window specific callback functions. These virtual functions should be replaced.
  virtual bool runOnInit(WindowRef w) { if (cb_init) cb_init(this); return true; }
  virtual void runOnCloseReq(void) { }  // this can be cancelled
  virtual void runOnDestroy (void) { if (cb_close) cb_close(this); }  // This function is called automatically by a event, and cannot be cancelled
  virtual bool runOnDraw    (void) { if (cb_draw) cb_draw(this); return true; }
  virtual void runOnKeyboard(int key) { if (key > 0 && cb_keys) cb_keys( this, key ); };
  virtual void runOnMouse(int button, int state, int xy[]) { if (cb_mouse) cb_mouse( this, button, state, xy ); }
  bool runOnTimer(void) { if (!cb_timer || timer_paused) return false;  cb_timer(this);  return true; }
  
public:
  virtual void openWindow(int w, int h, const char *title = NULL, 
			  void (*init) (Window *w) = NULL,
			  void (*draw) (Window *w) = NULL,
			  void (*keys) (Window *w, int k) = NULL, 
			  void (*mouse)(Window *w, int button, int state, int xy[]) = NULL,
			  void (*close)(Window *w) = NULL) {
    if (w <= 0) w = 200;  if (h <= 0) h = 200;
    closeWindow();
    setCallbacks(draw, keys, mouse, close);
    createBasicWindow( w, h, title );	// create the window
    ShowWindow( window );
    SelectWindow( window );
  }
  virtual void closeWindow(void) { 
    // Note that this is just one of several ways of closing window.
    // Therefore, actual clean-up should be done in 'runOnDestroy()'.
    if (window) { 
      DisposeWindow( window );  window = NULL; 
    }
  }
  
  void setupOSXMenubar(void) {
    ProcessSerialNumber psn = { 0, kCurrentProcess };
    TransformProcessType (&psn, kProcessTransformToForegroundApplication);
    SetFrontProcess(&psn);
    ClearMenuBar();				//Clear Menu Bar
    MenuRef windMenu;
    CreateStandardWindowMenu(0, &windMenu);	//Create Window Menu
    InsertMenu(windMenu, 0);
    ShowMenuBar();
    guih_mbar_ready = true;
  }
  
  void setupOSXAppEventHandler(void) {
    EventTypeSpec app_events[] = {{ kEventClassCommand, kEventCommandProcess }};
    OSStatus err;
    err = InstallApplicationEventHandler(NewEventHandlerUPP(AppEventHandler), 
					 GetEventTypeCount(app_events), app_events, NULL, NULL);
    guih_app_ready = true;
  }
  
  void createBasicWindow(int w, int h, const char *title, bool resizable=true) {
    if (!guih_app_ready)  setupOSXAppEventHandler();	// set OS X foreground process
    if (!guih_mbar_ready) setupOSXMenubar();
    if (w < 100) w = 100;  if (h < 100) h = 100;
    closeWindow();
    Rect rect;  int  mbar=44, offset = guih_window_count * 22;
    SetRect( &rect, offset+0, mbar+offset+0, offset+w, mbar+offset+h );
    OSStatus err;
    // http://developer.apple.com/documentation/Carbon/Reference/Window_Manager/Reference/reference.html#//apple_ref/c/tdef/WindowAttributes
    err = CreateNewWindow( kDocumentWindowClass, 
			   (kWindowStandardHandlerAttribute | kWindowCloseBoxAttribute | 
			    kWindowFullZoomAttribute | kWindowCollapseBoxAttribute | 
			    (resizable ? kWindowResizableAttribute : kWindowNoAttributes)),
			   &rect, &window );
    resize( w, h, true );  // set the content region to the given size
    setTitle( title );
    EventTypeSpec win_events[] = {{ kEventClassCommand, kEventCommandProcess },
				  { kEventClassWindow, kEventWindowActivated },
				  { kEventClassWindow, kEventWindowDeactivated },
				  { kEventClassWindow, kEventWindowShown },
				  { kEventClassWindow, kEventWindowClose },
				  { kEventClassWindow, kEventWindowClosed },
				  { kEventClassWindow, kEventWindowDrawContent },
				  { kEventClassWindow, kEventWindowZoomed },
				  { kEventClassWindow, kEventWindowBoundsChanged },
				  { kEventClassKeyboard, kEventRawKeyDown },
				  { kEventClassKeyboard, kEventRawKeyRepeat },
				  { kEventClassKeyboard, kEventRawKeyModifiersChanged },
				  { kEventClassMouse, kEventMouseMoved },
				  { kEventClassMouse, kEventMouseWheelMoved },
				  { kEventClassMouse, kEventMouseDown },
				  { kEventClassMouse, kEventMouseUp },
				  { kEventClassMouse, kEventMouseDragged }};
    err = InstallWindowEventHandler( window, NewEventHandlerUPP(WindowEventHandler),
				     GetEventTypeCount(win_events), win_events, this, NULL );
    guih_window_count++;
    //std::cout << guih_window_count << std::endl;
  }
  // register user callback functions --------------------------------
  void setCallbacks( void (*draw) (Window *w),
		     void (*keys) (Window *w, int k), 
		     void (*mouse)(Window *w, int button, int state, int xy[]) = NULL,
		     void (*close)(Window *w) = NULL ) {
    cb_draw = draw;  cb_keys = keys;  cb_mouse = mouse;  cb_close = close;
  }
  
public:
  int  getSerialNumber(void) { return serial_number; }
  bool isClosed(void)      { return (window == NULL); }
  
  void setTitle(const char *title) {
    if (!window) return;
    CFStringRef titleKey = CFStringCreateWithCString(NULL, title, kCFStringEncodingASCII);
    SetWindowTitleWithCFString( window, titleKey );
    CFRelease( titleKey );
  }
  void move (int x, int y)   { 
    Rect rect;
    GetWindowBounds(window, kWindowTitleBarRgn, &rect);
    MoveWindow(window, x, y + (rect.bottom-rect.top), false); 
  }
  void move (Window *win, bool below) { 
    if (!win || win->isClosed()) return;
    int x=0,y=0,w,h;  win->getPos(&x,&y); win->getSize(&w,&h,false); 
    //printf("x=%d y=%d w=%d h=%d\n", x, y, w, h );
    if (below) move(x,y+h); else move(x+w,y); 
    move_ref = win;  move_below = below;
  }
  void resize (int w, int h, bool interior=false)  { 
    Rect in_rect, ex_rect; 
    if (interior) {
      GetWindowBounds(window, kWindowContentRgn, &in_rect); 
      GetWindowBounds(window, kWindowStructureRgn, &ex_rect); 
      SetRect( &ex_rect, ex_rect.left, ex_rect.top, 
	       ex_rect.right  + (w - (in_rect.right - in_rect.left)),
	       ex_rect.bottom + (h - (in_rect.bottom - in_rect.top)) );
      SetWindowBounds(window, kWindowStructureRgn, &ex_rect); 
    } else {
      GetWindowBounds(window, kWindowContentRgn, &ex_rect); 
      SetRect( &ex_rect, ex_rect.left, ex_rect.top, ex_rect.left + w, ex_rect.top + h );
      SetWindowBounds(window, kWindowStructureRgn, &ex_rect); 
    }
    redraw();
  }
  void getPos (int *x, int *y) { Rect rect; GetWindowBounds(window, kWindowStructureRgn, &rect); *x=rect.left; *y=rect.top; }
  void getSize(int *w, int *h, bool interior) { 
    Rect rect; 
    GetWindowBounds(window, (interior ? kWindowContentRgn:kWindowStructureRgn), &rect); 
    *w=(rect.right-rect.left); *h=(rect.bottom-rect.top); 
  }
  virtual void redraw  (void) {
    Rect rect; 
    GetWindowBounds(window, kWindowContentRgn, &rect); 
    rect.left = 0;   rect.right = rect.right - rect.left;
    rect.top = 0;    rect.bottom = rect.bottom - rect.top;
    InvalWindowRect(window, &rect); 
  }
  void present (void) { SelectWindow( window ); }
  bool isActive(void) { return IsWindowActive(window); }
  void show    (void) { ShowWindow(window); }
  void hide    (void) { HideWindow(window); }
  void iconify (void) { }
  void maximize(void)   { }
  void unmaximize(void) { }
  virtual bool captureScreenInRGBBuffer(unsigned char *buffer) { return false; }
  bool captureScreen(char *fname=NULL) { 
    int  w, h;  getSize( &w, &h, true );  // interior only
    unsigned char *imgbuf=(unsigned char*)malloc( w * h * 3 * sizeof(unsigned char) );
    if (!imgbuf) return false;
    bool ret = (captureScreenInRGBBuffer(imgbuf) &&
		saveRGBBufferInFile( (fname ? fname : (char*)"guih_capture.png"), imgbuf, w, h ));
    if (imgbuf) { free(imgbuf); imgbuf = NULL; }
    return ret;
  }
  bool saveRGBBufferInFile(char *fname, void *buffer, int w, int h) {
    unsigned char *imgbuf = (unsigned char *)buffer;
    if (!imgbuf) return false;
    CGColorSpaceRef   cs = CGColorSpaceCreateWithName( kCGColorSpaceGenericRGB );
    CGDataProviderRef dp = CGDataProviderCreateWithData( NULL, imgbuf, 3*w*h, NULL );
    CGImageRef imgref = CGImageCreate( w, h, 8, 8*3, 3*w, cs, kCGBitmapByteOrderDefault, 
				       dp, NULL, false, kCGRenderingIntentDefault );
    CGColorSpaceRelease( cs );  CGDataProviderRelease( dp );
    // save the image using CGImageDestination
    // (http://developer.apple.com/documentation/Carbon/Conceptual/understanding_utis/understand_utis_intro/chapter_1_section_1.html)
    CFStringRef type;  char *ext = getFileExtension( fname );
    if      (strcmp(ext, "png")==0 || strcmp(ext, "PNG")==0)  type = kUTTypePNG;
    else if (strcmp(ext, "jpg")==0 || strcmp(ext, "JPG")==0)  type = kUTTypeJPEG;
    else if (strcmp(ext, "bmp")==0 || strcmp(ext, "BMP")==0)  type = kUTTypeBMP;
    else if (strcmp(ext, "gif")==0 || strcmp(ext, "GIF")==0)  type = kUTTypeGIF;
    else if (strcmp(ext, "tif")==0 || strcmp(ext, "tiff")==0) type = kUTTypeTIFF;
    else { CGImageRelease( imgref );  return false; }
    CFURLRef url = CFURLCreateFromFileSystemRepresentation( NULL, (const UInt8*)fname, strlen(fname), false );
    CGImageDestinationRef idst = CGImageDestinationCreateWithURL( url, type, 1, NULL );
    CFRelease( url );    if (idst == NULL) return false; 
    CGImageDestinationAddImage( idst, imgref, NULL );
    bool status = CGImageDestinationFinalize( idst );
    CGImageRelease( imgref );
    return status;
  }
  char* getFileExtension(char *filename) {
    int  i;  // find the beginning of the file extension of 'filename'
    for (i=strlen(filename); i>0 && filename[i]!='.'; i--);
    return (i>0 ? filename+i+1 : (char*)"");
  }
  
  // caption ---------------------------------------------------------
  void setCaption (char *str=NULL) { 
    if (str==NULL) caption[0] = '\0'; 
    else { 
      unsigned int i;
      for (i=0; str[i] && i< sizeof(caption)-2; i++) caption[i] = str[i];
      caption[i] = '\0';
    }
  }
  virtual void draw_caption(void) {   // This function shouldn't be called by user 
    if (!window || !windowContext || !caption[0]) return;
    // fonts : "Times", "Times-Bold", "Courier", "Courier New", "Arial", etc.
    CGContextSelectFont( windowContext, "Arial", 12, kCGEncodingMacRoman ); // kCGEncodingFontSpecific
    CGContextSetCharacterSpacing( windowContext, 1 );
    CGContextSetTextDrawingMode( windowContext, kCGTextFill );  // kCGTextFillStroke
    //CGAffineTransform xform =  CGAffineTransformMakeRotation  (MyRadians (45));
    //CGContextSetTextMatrix( windowContext, xform );
    //CGContextSetRGBStrokeColor( windowContext, 1, 1, 1, 1 );
    CGContextSetRGBFillColor( windowContext, 0, 0, 0, 1 );
    CGContextShowTextAtPoint( windowContext, 11,  9, caption, strlen(caption) );
    CGContextSetRGBFillColor( windowContext, 1, 0, 0, 1 );
    CGContextShowTextAtPoint( windowContext, 10, 10, caption, strlen(caption) );
  }
  
  // timer -----------------------------------------------------------
  void createTimer(double interval_sec, void (*ftn)(Window *w), bool activate_now=true) {
    cb_timer = ftn;  timer_interval = interval_sec;
    if (activate_now) {
      EventLoopRef       mainLoop = GetMainEventLoop();
      timer_upp = NewEventLoopTimerUPP(TimerAction);
      InstallEventLoopTimer( mainLoop, (EventTimerInterval)interval_sec, 1, timer_upp, this, &timer_ref );
      timer_paused = false;
    }
  }
  void removeTimer(void) { 
    if (timer_ref) { RemoveEventLoopTimer( timer_ref );  timer_ref = NULL; }
    if (timer_upp) { DisposeEventLoopTimerUPP( timer_upp );  timer_upp = NULL; }
    timer_paused = true; 
  }
  bool toggleTimer(double interval_sec=-1, void (*ftn)(Window *w)=NULL) { 
    if (interval_sec >= 0) timer_interval = interval_sec;
    if (ftn) cb_timer = ftn;
    if (timer_paused && cb_timer) {
      createTimer( timer_interval, cb_timer );
      return true;
    } else { 
      removeTimer();
      return false;
    }
  }
  bool isTimerPaused(void) { return timer_paused; }
  
  // keyboard interpretation -----------------------------------------
  int remapKeyValue(EventRef event, bool *skey=NULL) {
    UInt32 macKeyCode;  char macCharCode;  UniChar uc;
    GetEventParameter(event, kEventParamKeyMacCharCodes, typeChar, NULL, sizeof(macCharCode), NULL, &macCharCode);
    GetEventParameter(event, kEventParamKeyCode, typeUInt32, NULL, sizeof(macKeyCode), NULL, &macKeyCode);
    GetEventParameter(event, kEventParamKeyUnicodes, typeUnicodeText, NULL, sizeof(uc), NULL, &uc);
    //printf("macKeyCode=%ld  macCharCode=%d(x:%x)  skey_by_user=%d  uniChar=%04x=%c\n", (long)macKeyCode, macCharCode, macCharCode, skey_by_user, uc, (char)uc);
    //For more information, read  http://developer.apple.com/qa/qa2005/qa1446.html
    if      (ctrl && macCharCode >  0 && macCharCode < 27) macCharCode += 96;  // ^a ~ ^z
    else if (ctrl && macCharCode > 26 && macCharCode < 32) macCharCode += 64;  // ^[ ~ ^_
    else if (alt) 
      switch (macKeyCode) {  // alt + NumericKeys
      case 18: macCharCode = 49; break;   case 19: macCharCode = 50; break;
      case 20: macCharCode = 51; break;   case 21: macCharCode = 52; break;
      case 23: macCharCode = 53; break;   case 22: macCharCode = 54; break;
      case 26: macCharCode = 55; break;   case 28: macCharCode = 56; break;
      case 25: macCharCode = 57; break;   case 29: macCharCode = 48; break;
      }
    if      (macKeyCode == 126) { if (skey) *skey = true; return ARROW_UP; }
    else if (macKeyCode == 125) { if (skey) *skey = true; return ARROW_DOWN; }
    else if (macKeyCode == 123) { if (skey) *skey = true; return ARROW_LEFT; }
    else if (macKeyCode == 124) { if (skey) *skey = true; return ARROW_RIGHT; }
    else if (macKeyCode ==  53) { if (skey) *skey = true; return ESCAPE_KEY; }
    else if (macCharCode ==  8) { if (skey) *skey = false; return '\b'; } //DELETE
    else if (macKeyCode == 122) { if (skey) *skey = false; return F1; }
    else if (macKeyCode == 120) { if (skey) *skey = false; return F2; }
    else if (macKeyCode ==  99) { if (skey) *skey = false; return F3; }
    else if (macKeyCode == 118) { if (skey) *skey = false; return F4; }
    else if (macKeyCode ==  96) { if (skey) *skey = false; return F5; }
    else if (macKeyCode ==  97) { if (skey) *skey = false; return F6; }
    else if (macKeyCode ==  98) { if (skey) *skey = false; return F7; }
    else if (macKeyCode == 100) { if (skey) *skey = false; return F8; }
    else if (macKeyCode == 101) { if (skey) *skey = false; return F9; }
    else if (macKeyCode == 109) { if (skey) *skey = false; return F10; }
    else if (macKeyCode == 103) { if (skey) *skey = false; return F11; }
    else if (macKeyCode == 111) { if (skey) *skey = false; return F12; }
    else if (cmmd && (macCharCode=='i' || macCharCode=='I')) { if (skey) *skey = true; return INSERT_KEY; }
    else if (cmmd && (macCharCode=='-' || macCharCode=='_')) { if (skey) *skey = true; return HOME_KEY; }
    else if (cmmd && (macCharCode=='=' || macCharCode=='+')) { if (skey) *skey = true; return END_KEY; }
    else if (cmmd && (macCharCode=='[' || macCharCode=='{')) { if (skey) *skey = true; return PAGE_UP; }
    else if (cmmd && (macCharCode==']' || macCharCode=='}')) { if (skey) *skey = true; return PAGE_DOWN; }
    else { if (skey) *skey = false; return macCharCode; }
  }
  int remapNumericKey2Integer(char key, bool add10=false) {
    // Remap numeric keys "0123456789" and symbols ")!@#$%^&*(", so that we get [0,1,...19].
    int val = -1;
    if (key >= '0' && key <= '9') val = key - '0';
    else if (key==')') val = 10;  else if (key=='!') val = 11;
    else if (key=='@') val = 12;  else if (key=='#') val = 13;
    else if (key=='$') val = 14;  else if (key=='%') val = 15;
    else if (key=='^') val = 16;  else if (key=='&') val = 17;
    else if (key=='*') val = 18;  else if (key=='(') val = 19;
    return (val<0 ? -1 : (val + (add10 ? 10 : 0)));
  }
  
  // user input dialog -----------------------------------------------
  bool openMessageDialog(char *title, char *message, bool yesno=false) { 
    DialogRef theAlert;
    CFStringRef msg = ::CFStringCreateWithBytes(NULL, (UInt8 *)message, strlen(message), CFStringGetSystemEncoding(), false);
    CreateStandardAlert(kAlertStopAlert, msg, NULL, NULL, &theAlert);
    RunStandardAlert(theAlert, NULL, NULL);
//     ExitToShell();
    // 
//     alertStringExplanation = SomeUTF String;
//     alertStringTitle = SomeUTF String ;
//     l1 = alertStringExplanation.length();
//     l2 = alertStringTitle.length();
//     UInt32 textLength = l2;
//     CFStringRef errMsg = ::CFStringCreateWithBytes(NULL, (UInt8 *) alertStringTitle.get(), textLength, CFStringGetSystemEncoding(), false);
//     textLength = l1;
//     CFStringRef errorExplanation = ::CFStringCreateWithBytes(NULL, (UInt8 *) alertStringExplanation.get(), textLength, CFStringGetSystemEncoding(), false);
//     UInt32 numUnicodeChars = ::CFStringGetMaximumSizeForEncoding(textLength, kCFStringEncodingUnicode);

//     DialogRef alert;
//     DialogItemIndex outHit;
//     CreateStandardAlert(kAlertCautionAlert, errMsg, errorExplanation, NULL, &alert);
//     RunStandardAlert(alert, NULL, &outHit);
    return true; 
  }
  bool openInputDialog(char *title, UserInput *ui0, UserInput *ui1=NULL, UserInput *ui2=NULL) { return true; }
  bool openInputDialogScale(char *title, UserInput *ui) { return true; }
  
  // drawing ---------------------------------------------------------
  void setColorForDrawing(float r, float g, float b) { 
    CGColorRef color = CGColorCreateGenericRGB( r, g, b, 1.0 );
    CGContextSetStrokeColorWithColor( windowContext, color );
    CGColorRelease (color);
  }
  void drawLine(int x0, int y0, int x1, int y1) { 
    int ww, hh;  getSize( &ww, &hh, true );
    CGContextMoveToPoint( windowContext, x0, hh-y0 );
    CGContextAddLineToPoint( windowContext, x1, hh-y1 );
  }
  void drawRectCenteredAt(int x, int y, int size, bool filled=false) { 
    // Draw a rectangle centered at (x,y)
    int ww, hh;  getSize( &ww, &hh, true );
    CGRect rect = CGRectMake( x, hh-(y-size/2), size, -size );
    CGContextAddRect( windowContext, rect );
  }
  void drawRect(int x, int y, int w, int h, bool filled=false) { 
    // Draw a rectangle starting at (x,y) 
    int ww, hh;  getSize( &ww, &hh, true );
    CGRect rect = CGRectMake( x, hh-y, w, -h );
    CGContextAddRect( windowContext, rect );
  }
  void drawCross(int x, int y, int size, char type='+') { 
    if        (type=='+') {
      drawLine(x-size, y, x+size, y); drawLine(x, y-size, x, y+size); 
    } else if (type=='x') {
      drawLine(x-size, y-size, x+size, y+size); drawLine(x+size, y-size, x-size, y+size); 
    } else if (type=='|') {
      int hs = (size/2);
      drawLine(x-hs, y, x+hs, y); drawLine(x, y-size, x, y+size); 
    } else if (type=='-') {
      int hs = (size/2);
      drawLine(x-size, y, x+size, y); drawLine(x, y-hs, x, y+hs); 
    } else {
      int dx=0, dy=0;
      switch (type) {
      case '0':  dx = +size;  dy = -size;  break;
      case '1':  dx = -size;  dy = -size;  break;
      case '2':  dx = -size;  dy = +size;  break;
      case '3':  dx = +size;  dy = +size;  break;
      }
      drawLine( x, y,  x+dx, y );   drawLine( x, y,  x, y+dy );
    }
  }
  
  // -----------------------------------------------------------------
  // Non window-specific functions
  // -----------------------------------------------------------------
  void runMainLoop(void) { 
    // start main loop
    SelectWindow( window );
    RunApplicationEventLoop();
    // There is a strange bug in Carbon and CoreGraphics in OS X, when used by multithread.
    // You must execute 'setupOSXAppEventHandler()' before any 'CGCreatImage..()' function.
    // Otherwise, 'RunApplicationEventLoop()' may fail and falls to a halt.
  }
  void addIdleFunction(void (*idle_ftn)(void* data), void* data=NULL) { 
    // idle  function
    guih_idle_ftn = idle_ftn;
    EventLoopRef mainLoop = GetMainEventLoop();
    idle_upp = NewEventLoopIdleTimerUPP(IdleTimerAction);
    InstallEventLoopIdleTimer( mainLoop, (EventTimerInterval)0.2, (EventTimerInterval)0.01, 
			       idle_upp, this, &idle_ref );
  }
  void removeIdleFunction(void* data=NULL) { 
    if (idle_ref) { RemoveEventLoopTimer(idle_ref); idle_ref = NULL; }
    if (idle_upp) { DisposeEventLoopIdleTimerUPP( idle_upp );  idle_upp = NULL; }
  }

};
  
static pascal void TimerAction(EventLoopTimerRef theTimer, void* userData) 
{
  Window *w = (Window*)userData;
  w->runOnTimer();
}
  
static pascal void IdleTimerAction (EventLoopTimerRef theTimer, EventLoopIdleTimerMessage idleMsg, void* userData)
{
  Window *w = (Window*)userData;
  // idleMsg: kEventLoopIdleTimerStarted kEventLoopIdleTimerIdling kEventLoopIdleTimerStopped
  if (guih_idle_ftn) guih_idle_ftn(w);
  else w->removeIdleFunction();
}
  
static OSStatus AppEventHandler(EventHandlerCallRef nextHandler, EventRef event, void *userData)
{
  UInt32   clss   = GetEventClass (event);
  //UInt32   kind   = GetEventKind (event); 
  OSStatus result = CallNextEventHandler(nextHandler, event);
	
  if (clss == kEventClassCommand) {
    HICommand theHICommand;
    GetEventParameter( event, kEventParamDirectObject, typeHICommand, NULL, sizeof(HICommand), NULL, &theHICommand );
    switch ( theHICommand.commandID )	{
    case kHICommandQuit:			// Quit (Command-Q)
      break;
    default:  result = eventNotHandledErr;  break;
    }
  } else if (clss == kEventClassKeyboard) {
  } else if (clss == kEventClassMouse) {
  }
	
  return result;
}

static OSStatus WindowEventHandler(EventHandlerCallRef nextHandler, EventRef event, void *userData)
{
  Window*   w   = (Window*)userData;
  UInt32   clss   = GetEventClass (event);
  UInt32   kind   = GetEventKind (event); 
  OSStatus result = CallNextEventHandler(nextHandler, event);
	
  if (clss == kEventClassCommand) {
    HICommand theHICommand;
    GetEventParameter( event, kEventParamDirectObject, typeHICommand, NULL, sizeof(HICommand), NULL, &theHICommand );
    switch ( theHICommand.commandID ) {
    case kHICommandQuit:  w->closeWindow(); break;	// Quit (Command-Q)
    default:  result = eventNotHandledErr;  break;
    }
  } else if(clss == kEventClassWindow) {
    WindowRef     window;
    GetEventParameter(event, kEventParamDirectObject, typeWindowRef, NULL, sizeof(WindowRef), NULL, &window);
    switch (kind) {
    case kEventWindowActivated:     w->runOnInit(w->window);            break;  ////// init after draw ?
    case kEventWindowDrawContent:
      if (w->window) {
	QDBeginCGContext( GetWindowPort(w->window), &w->windowContext );	// begin CGContext
	CGContextBeginPath( w->windowContext );
	w->runOnDraw();
	if (!CGContextIsPathEmpty( w->windowContext )) {
	  CGContextClosePath( w->windowContext );
	  //CGContextFillPath( w->windowContext );
	  CGContextStrokePath( w->windowContext );
	}
	w->draw_caption();  // draw caption
	QDEndCGContext( GetWindowPort(w->window), &w->windowContext );	// end   CGContext
      }
      break;
    case kEventWindowClose:  break;   // this event will be skipped by 'Quit' (Command-Q)
    case kEventWindowClosed: 
      w->runOnDestroy(); 
      if (--guih_window_count <= 0) {
	QuitApplicationEventLoop();
	ProcessSerialNumber psn = { 0, kCurrentProcess };
	GetCurrentProcess( &psn );
	GetNextProcess( &psn );
	SetFrontProcess( &psn );
      }
      //std::cout << guih_window_count << std::endl;
      break;
    case kEventWindowZoomed:        break;
    case kEventWindowBoundsChanged: break;
    default:  result = eventNotHandledErr;  break;
    }
  } else if (clss == kEventClassKeyboard) {
    UInt32 macKeyCode, macKeyModifiers;
    GetEventParameter(event, kEventParamKeyModifiers, typeUInt32, NULL, sizeof(macKeyModifiers), NULL, &macKeyModifiers);
    GetEventParameter(event, kEventParamKeyCode, typeUInt32, NULL, sizeof(macKeyCode), NULL, &macKeyCode);
    if (kind == kEventRawKeyModifiersChanged) {
      w->alt   = (macKeyModifiers & optionKey);
      w->ctrl  = (macKeyModifiers & controlKey);
      w->shift = (macKeyModifiers & shiftKey);
      w->cmmd  = (macKeyModifiers & cmdKey);
    } else if (kind == kEventRawKeyDown || kind == kEventRawKeyRepeat) {
      // window-type specific callback tasks
      int key = w->remapKeyValue(event);
      w->runOnKeyboard( key );
      switch (key) {
      case ESCAPE_KEY:  w->closeWindow();  break;
      case INSERT_KEY:  w->captureScreen(); break;
      case DELETE_KEY: case ' ':  w->setCaption(NULL); break;
      default: break;
      }
      w->redraw();
    }
  } else if (clss == kEventClassMouse) {
    Point pos;  int button = MOUSE_LEFT;
    static int mhpos, mhxy[802];
    GetEventParameter(event, kEventParamMouseLocation, typeQDPoint, NULL, sizeof(Point), NULL, &pos);
    int x, y;  w->getPos( &x, &y );
    Rect in_rect, ex_rect; // for the calculation of window frame size
    GetWindowBounds(w->window, kWindowContentRgn, &in_rect); 
    GetWindowBounds(w->window, kWindowStructureRgn, &ex_rect); 
    x = pos.h - ( x + in_rect.left - ex_rect.left);
    y = pos.v - ( y + in_rect.top  - ex_rect.top);
    //printf("pos=(%d %d)  xy=(%d %d)  in=(%d %d) ex=(%d %d)\n", pos.h, pos.v, x, y, in_rect.left, in_rect.top, ex_rect.left, ex_rect.top);
    if      (kind==kEventMouseDown   ) { mhxy[0]=x; mhxy[1]=y; mhpos=2; }
    else if (kind==kEventMouseDragged) { if (mhpos<800-1) { mhxy[mhpos++]=x; mhxy[mhpos++]=y; } }
    else if (kind==kEventMouseUp     ) { 
      mhxy[mhpos+0] = -1;  mhxy[mhpos+1] = -1;
      if (mhxy[mhpos-1]<0) ; // do nothing
      else if (mhpos==2) w->runOnMouse( button, MOUSE_CLICKED, mhxy );
      else               w->runOnMouse( button, MOUSE_DRAGGED, mhxy );
      mhpos = 0;  w->redraw();
    }
  }
	
  return result;
}

}	// end of GUIH namespace

// ===================================================================
#else		// beginning of Linux/Unix version
// ===================================================================

#include <iostream>
#include <gtk/gtk.h>
#include <gdk/gdkkeysyms.h>
#include <gtk/gtkgl.h>
#include <string.h>
#include <stdlib.h>

namespace GUIH {

static void     on_destroy (GtkWidget *widget, gpointer udata);
static gboolean on_close_rq(GtkWidget *widget, GdkEvent *event, gpointer udata);
static gboolean on_init    (GtkWidget *widget, gpointer udata);
static gboolean on_draw    (GtkWidget *widget, GdkEventExpose *event, gpointer udata);
static gboolean on_keypress(GtkWidget *widget, GdkEventKey *event, gpointer udata);
static gboolean on_button_press  (GtkWidget *widget, GdkEventButton *event, gpointer data);
static gboolean on_button_release(GtkWidget *widget, GdkEventButton *event, gpointer data);
static gboolean on_motion_notify (GtkWidget *widget, GdkEventMotion *event, gpointer data);
static gboolean on_timer   (void *data);
static gboolean on_idle    (gpointer data);
static int      guih_window_count = 0;	// total number of opened windows
static int      guih_serial_number = -1;
static void (*guih_idle_ftn)(void *data) = NULL;
  
// -------------------------------------------------------------------

class Window
{
protected:
  GtkWidget	*window;	// GtkWindow *
  GtkWidget	*vbox;	// GtkVBox *
  GtkWidget	*area;	// GtkDrawingArea *
  int		serial_number;	// unique window index, which begins from 1
  bool		timer_paused;
  double	timer_interval; // in seconds
  char		caption[256];	// caption
  // user callbacks
  void  (*cb_init)   (Window *w);
  void  (*cb_draw)   (Window *w);
  void  (*cb_keys)   (Window *w, int k);
  void  (*cb_mouse)  (Window *w, int button, int state, int xy[]);
  void  (*cb_close)  (Window *w);
  void  (*cb_timer)  (Window *w);
  
public:
  // Window specific callback functions. These virtual functions should be replaced.
  virtual bool runOnInit(GtkWidget *widget) { if (cb_init) cb_init(this); return true; }
  virtual void runOnCloseReq(void) { }  // this can be cancelled
  virtual void runOnDestroy (void) { if (cb_close) cb_close(this); }  // This function is called automatically by a event, and cannot be cancelled
  virtual bool runOnDraw(GtkWidget *widget, GdkEventExpose *event) { if (cb_draw) cb_draw(this); return true; }
  virtual void runOnKeyboard(int key) { if (key > 0 && cb_keys) cb_keys( this, key ); }
  virtual void runOnMouse(int button, int state, int xy[]) { if (cb_mouse) cb_mouse( this, button, state, xy ); }
  bool runOnTimer(void) { if (cb_timer == NULL || timer_paused) return false;  cb_timer(this);  return true; }
  int  mhxy[800], mhpos;
  
public:
  bool  shift, ctrl, alt, skey_by_user;  // user can access these values
  
  Window *move_ref;
  bool   move_below;
  
public:
  Window () : window(NULL), vbox(NULL), area(NULL), timer_paused(true), 
	      cb_init(NULL), cb_draw(NULL), cb_keys(NULL),
	      cb_mouse(NULL), cb_close(NULL), cb_timer(NULL), 
	      mhpos(0), shift(0), ctrl(0), alt(0), skey_by_user(false) {
    caption[0] = '\0'; 
    serial_number = ++guih_serial_number;
    move_ref=NULL; 
  }
  Window (int w, int h, const char *title = NULL, 
	  void (*init) (Window *w) = NULL,
	  void (*draw) (Window *w) = NULL,
	  void (*keys) (Window *w, int key) = NULL, 
	  void (*mouse)(Window *w, int button, int state, int xy[]) = NULL,
	  void (*close)(Window *w) = NULL) 
    : window(NULL), vbox(NULL), area(NULL), timer_paused(true), 
      cb_init(NULL), cb_draw(NULL), cb_keys(NULL),
      cb_mouse(NULL), cb_close(NULL), cb_timer(NULL), 
      mhpos(0), shift(0), ctrl(0), alt(0), skey_by_user(false) { 
    caption[0] = '\0'; 
    serial_number = ++guih_serial_number;
    openWindow(w, h, title, init, draw, keys, mouse, close); 
  }
  virtual ~Window() { closeWindow(); }
public:
  virtual void openWindow(int w, int h, const char *title = NULL, 
			  void (*init) (Window *w) = NULL,
			  void (*draw) (Window *w) = NULL,
			  void (*keys) (Window *w, int k) = NULL, 
			  void (*mouse)(Window *w, int button, int state, int xy[]) = NULL,
			  void (*close)(Window *w) = NULL) {
    if (w <= 0) w = 200;  if (h <= 0) h = 200;
    closeWindow();
    setCallbacks(draw, keys, mouse, close);
    createBasicWidgets( w, h, title );	// create the window and other widgets
    RegisterBasicEvents( init );	// setup default event handlers
    gtk_widget_show( window );		// show the window
  }
  virtual void closeWindow(void) { 
    // Note that this is just one of several ways of closing window.
    // Therefore, actual clean-up should be done in 'runOnDestroy()'.
    if (window) { gtk_widget_destroy( window ); window = NULL; }
  }
  //void closeWindow(void) { bool ret;  if (window) g_signal_emit_by_name( window, "delete-event", this, &ret ); }
  
  void createBasicWidgets(int w, int h, const char *title) {
    if (w < 100) w = 100;  if (h < 100) h = 100;
    closeWindow();
    if (++guih_window_count == 1) {
      int argc=0;  char **argv=NULL;
      gtk_init( &argc, &argv );
    }
    window = gtk_window_new( GTK_WINDOW_TOPLEVEL );
    vbox   = gtk_vbox_new( FALSE, 0 );
    area   = gtk_drawing_area_new();
    gtk_window_set_title( GTK_WINDOW(window), (title ? title : "Untitled") );
    gtk_box_pack_end( GTK_BOX(vbox), area, TRUE, TRUE, 0 );
    gtk_widget_show( area );
    gtk_container_add( GTK_CONTAINER(window), vbox );
    gtk_widget_show( vbox );
    gtk_window_resize(GTK_WINDOW(window), w, h );
    //gtk_widget_set_size_request( window, w, h );
  }
  // register user callback functions --------------------------------
  void setCallbacks( void (*draw) (Window *w),
		     void (*keys) (Window *w, int k), 
		     void (*mouse)(Window *w, int button, int state, int xy[]) = NULL,
		     void (*close)(Window *w) = NULL ) {
    if (!window) return;
    cb_draw = draw;  cb_keys = keys;  cb_mouse = mouse;  cb_close = close;
    // setup addtional event handlers
    gtk_widget_add_events (area,
			   GDK_BUTTON1_MOTION_MASK | // dragging with left button
			   GDK_BUTTON_PRESS_MASK   |
			   GDK_BUTTON_RELEASE_MASK );
    g_signal_connect(G_OBJECT(area), "motion_notify_event",  G_CALLBACK(on_motion_notify), this);
    g_signal_connect(G_OBJECT(area), "button_press_event",   G_CALLBACK(on_button_press), this);
    g_signal_connect(G_OBJECT(area), "button_release_event", G_CALLBACK(on_button_release), this);
  }
  void RegisterBasicEvents( void(*init)(Window *w)=NULL ) {
    cb_init = init;
    g_signal_connect_after (G_OBJECT(area), "realize",   G_CALLBACK(on_init),    this);
    g_signal_connect(G_OBJECT(window), "destroy",        G_CALLBACK(on_destroy), this);
    g_signal_connect(G_OBJECT(window), "delete-event",   G_CALLBACK(on_close_rq),this);
    g_signal_connect(G_OBJECT(window), "key-press-event",G_CALLBACK(on_keypress),this);
    g_signal_connect(G_OBJECT(area),   "expose-event",   G_CALLBACK(on_draw),    this);
  }
  
public:
  int  getSerialNumber(void) { return serial_number; }
  bool isClosed(void)      { return (window == NULL); }
  
  void setTitle(const char *title) { if (window) gtk_window_set_title(GTK_WINDOW(window), title); }
  void move (int x, int y)   { if (window) gtk_window_move     (GTK_WINDOW(window), x, y ); move_ref=NULL; }
  void move (Window *win, bool below) { 
    if (!win || win->isClosed()) return;
    struct timespec req, rem;
    req.tv_sec = 0;  req.tv_nsec = 200000; nanosleep( &req, &rem );
    int x,y,w,h;  win->getPos(&x,&y); win->getSize(&w,&h,false); 
    if (below) move(x,y+h); else move(x+w,y); 
    move_ref = win;  move_below = below;
  }
  void resize (int w, int h, bool interior=false)   { 
    if (!window) return;
    if (!interior) gtk_window_resize (GTK_WINDOW(window), w, h ); 
    else           gtk_window_resize (GTK_WINDOW(window), w, h );
//     { gtk_widget_set_size_request( area, w, h ); sleep(1); }
  }
  void getPos (int *x, int *y) { if (window) gtk_window_get_position(GTK_WINDOW(window), x, y ); }
  void getSize(int *w, int *h, bool interior) { 
    if (!window) return;
    if (!interior) {		// get the size of the window (including decorations)
      gtk_window_get_size(GTK_WINDOW(window), w, h ); 
      *w += 10;  *h += 27;
    } else { 			// get the size of the client area
      *w = area->allocation.width; *h = area->allocation.height; 
    }
  }
  void redraw  (void) { if (window) gdk_window_invalidate_rect(area->window, &area->allocation, FALSE); }
  void present (void) { if (window) gtk_window_present( GTK_WINDOW(window) ); }
  bool isActive(void) { return (window && (bool)gtk_window_is_active( GTK_WINDOW(window) )); }
  void show    (void) { if (window) gtk_widget_show( window ); }
  void hide    (void) { if (window) gtk_widget_hide( window ); }
  void iconify (void) { if (window) gtk_window_iconify( GTK_WINDOW(window) ); }
  void maximize(void)   { if (window) gtk_window_maximize  ( GTK_WINDOW(window) ); }
  void unmaximize(void) { if (window) gtk_window_unmaximize( GTK_WINDOW(window) ); }
  bool captureScreenInRGBBuffer(unsigned char *buffer) { 
    if (!window || !area) return false;
    GdkPixbuf *pb = gdk_pixbuf_get_from_drawable(NULL, area->window, 
						 gdk_colormap_get_system(), 0, 0, 0, 0, 
						 area->allocation.width, 
						 area->allocation.height);
    int w = gdk_pixbuf_get_width(pb);
    int h = gdk_pixbuf_get_height(pb);
    void *pixels = gdk_pixbuf_get_pixels(pb);
    memcpy( buffer, pixels, 3 * w * h );
    g_object_unref( pb );    return true;
  }
  bool captureScreen(char *fname=NULL) { 
    if (!window || !area) return false;
    GdkPixbuf *pb = gdk_pixbuf_get_from_drawable(NULL, area->window, 
						 gdk_colormap_get_system(), 0, 0, 0, 0, 
						 area->allocation.width, 
						 area->allocation.height);
    char *fn  = (fname ? fname : (char*)"guih_capture.png");
    char *ext = fn + strlen(fn)-3;
    if (strcmp(ext, "jpg")==0 || strcmp(ext, "JPG")==0)  gdk_pixbuf_save( pb, fn, "jpeg", NULL, NULL );
    else  gdk_pixbuf_save( pb, fn, "png", NULL, NULL );
    if (fname==NULL) printf("catured image saved in 'guih_capture.png'\n");
    g_object_unref( pb );    return true;
  }
  // caption ---------------------------------------------------------
  void setCaption (char *str=NULL) { 
    if (str) { strncpy(caption, str, 255);  str[255] = '\0'; } 
    else caption[0] = '\0'; 
  }
  void draw_caption(void) {   // This function shouldn't be called by user
    if (!caption[0]) return;
    int x = 10, y = area->allocation.height-20;
    GdkColor  color;
    PangoLayout *layout = gtk_widget_create_pango_layout( (GtkWidget*)area, caption );
    color.red = color.green = color.blue = guint16(65535 * 0.15);
    gdk_draw_layout_with_colors( area->window, area->style->fg_gc[GTK_STATE_NORMAL], x, y, layout, &color, NULL );
    color.red = 65535; color.green = color.blue = 0;
    gdk_draw_layout_with_colors( area->window, area->style->fg_gc[GTK_STATE_NORMAL], x-1, y-1, layout, &color, NULL );
  }
  
  // timer -----------------------------------------------------------
  void createTimer(double interval_sec, void (*ftn)(Window *w), bool activate_now=true) {
    cb_timer = ftn;  timer_interval = interval_sec;
    if (activate_now) {
      timer_paused = false;
      g_timeout_add( (guint)(interval_sec*1000), on_timer, this );
    }
  }
  void removeTimer(void) { timer_paused = true; }
  bool toggleTimer(double interval_sec=0, void (*ftn)(Window *w)=NULL) { 
    if (interval_sec > 0) timer_interval = interval_sec;
    if (ftn) cb_timer = ftn;
    if (timer_paused && cb_timer) {
      timer_paused = false; 
      g_timeout_add( (guint)(timer_interval*1000), on_timer, this );
      return true;
    } else { 
      timer_paused = true;  return false; 
    }
  }
  bool isTimerPaused(void) { return timer_paused; }
  
  // -----------------------------------------------------------------
  int remapKeyValue(int oldkey, bool *skey=NULL) {
    int  key;
    if (key==' ') setCaption();
    if      (oldkey >= GDK_a  && oldkey <= GDK_z)   key = 'a' + (oldkey - GDK_a);
    else if (oldkey >= GDK_A  && oldkey <= GDK_Z)   key = 'A' + (oldkey - GDK_A);
    else if (oldkey >= GDK_0  && oldkey <= GDK_9)   key = '0' + (oldkey - GDK_0);
    else if (oldkey >= GDK_F1 && oldkey <= GDK_F12) key = F1 + (oldkey - GDK_F1);
    else if (oldkey >= GDK_space && oldkey <= GDK_slash) key = 32 + (oldkey - GDK_space);
    else if (oldkey == GDK_Return) key = '\n';
    else if (oldkey == GDK_Tab) key = '\t';
    else if (oldkey == GDK_BackSpace) key = '\b';
    else if (oldkey == GDK_backslash) key = '\\';
    else if (oldkey < 255) key = oldkey;
    else if (oldkey == GDK_Up   ) { if (skey) *skey=true; key = ARROW_UP; }
    else if (oldkey == GDK_Down ) { if (skey) *skey=true; key = ARROW_DOWN; }
    else if (oldkey == GDK_Left ) { if (skey) *skey=true; key = ARROW_LEFT; }
    else if (oldkey == GDK_Right) { if (skey) *skey=true; key = ARROW_RIGHT; }
    else if (oldkey == GDK_Home)  { if (skey) *skey=true; key = HOME_KEY; }
    else if (oldkey == GDK_End)   { if (skey) *skey=true; key = END_KEY; }
    else if (oldkey == GDK_Page_Up)   { if (skey) *skey=true; key = PAGE_UP; }
    else if (oldkey == GDK_Page_Down) { if (skey) *skey=true; key = PAGE_DOWN; }
    else key = 0;
    return key;
  }
  int remapNumericKey2Integer(char key, bool add10=false) {
    // Remap numeric keys "0123456789" and symbols ")!@#$%^&*(", so that we get [0,1,...19].
    int val = -1;
    if (key >= '0' && key <= '9') val = key - '0';
    else if (key==')') val = 10;  else if (key=='!') val = 11;
    else if (key=='@') val = 12;  else if (key=='#') val = 13;
    else if (key=='$') val = 14;  else if (key=='%') val = 15;
    else if (key=='^') val = 16;  else if (key=='&') val = 17;
    else if (key=='*') val = 18;  else if (key=='(') val = 19;
    return (val<0 ? -1 : (val + (add10 ? 10 : 0)));
  }
  
  // user input dialog -----------------------------------------------
  bool openMessageDialog(char *title, char *message, bool yesno=false) {
    GtkWidget *dlg, *msg;		// just show the given message
    if (yesno)
      dlg = gtk_dialog_new_with_buttons(title, GTK_WINDOW(window), GTK_DIALOG_MODAL,
					GTK_STOCK_YES, GTK_RESPONSE_ACCEPT,
					GTK_STOCK_NO,  GTK_RESPONSE_CANCEL, NULL);
    else
      dlg = gtk_dialog_new_with_buttons(title, GTK_WINDOW(window), GTK_DIALOG_MODAL,
					GTK_STOCK_OK, GTK_RESPONSE_CANCEL, NULL);
    msg = gtk_label_new(message);
    gtk_container_add(GTK_CONTAINER (GTK_DIALOG(dlg)->vbox), msg);
    gtk_widget_show_all(dlg);
    return (gtk_dialog_run(GTK_DIALOG(dlg)) == GTK_RESPONSE_ACCEPT);
  }  
  bool openInputDialog(char *title, UserInput *ui0, UserInput *ui1=NULL, UserInput *ui2=NULL) {
    UserInput *ui[3] = { ui0, ui1, ui2 };  int i;
    GtkWidget *dlg, *label[3], *entry[3];
    dlg = gtk_dialog_new_with_buttons(title, GTK_WINDOW(window), GTK_DIALOG_MODAL,
				      GTK_STOCK_OK, GTK_RESPONSE_ACCEPT, 
				      GTK_STOCK_CANCEL, GTK_RESPONSE_CANCEL, NULL);
    GtkWidget *hbox, *lbox, *rbox;
    hbox = gtk_hbox_new(false,10);  gtk_container_add(GTK_CONTAINER(GTK_DIALOG(dlg)->vbox), hbox);
    lbox = gtk_vbox_new(true, 10);  gtk_container_add(GTK_CONTAINER(hbox), lbox);
    rbox = gtk_vbox_new(true, 10);  gtk_container_add(GTK_CONTAINER(hbox), rbox);
    for (i = 0; i < 3 && ui[i]; i++) {		// show current value
      label[i] = gtk_label_new( ui[i]->label );
      char buffer[80];
      switch (ui[i]->type) {
      case 't':  
	entry[i] = gtk_entry_new_with_max_length(ui[i]->len ? ui[i]->len : 40);
	sprintf(buffer, "%s", ui[i]->text);  gtk_entry_set_text(GTK_ENTRY(entry[i]), buffer);
	break;
      case 'f':  
	entry[i] = gtk_entry_new_with_max_length(ui[i]->len ? ui[i]->len : 12);
	sprintf(buffer, "%g", ui[i]->fvalues[0]);  gtk_entry_set_text(GTK_ENTRY(entry[i]), buffer);
	break;
      case 'd':  
	entry[i] = gtk_entry_new_with_max_length(ui[i]->len ? ui[i]->len : 12);
	sprintf(buffer, "%g", ui[i]->dvalues[0]);  gtk_entry_set_text(GTK_ENTRY(entry[i]), buffer);
	break;
      case 'i':  
	entry[i] = gtk_entry_new_with_max_length(ui[i]->len ? ui[i]->len : 12);
	sprintf(buffer, "%d", ui[i]->ivalues[0]);  gtk_entry_set_text(GTK_ENTRY(entry[i]), buffer);
	break;
      default: 
	openMessageDialog("Error", "Invalid type for openInputDialog() -- ['t','f','d','i']");
	gtk_widget_destroy(dlg);  return false;
      }
      gtk_container_add(GTK_CONTAINER(lbox), label[i]);
      gtk_container_add(GTK_CONTAINER(rbox), entry[i]);
    }
    gtk_widget_show_all(dlg);			// show the dialog window
    if (gtk_dialog_run(GTK_DIALOG(dlg)) == GTK_RESPONSE_ACCEPT) {
      for (i = 0; i < 3 && ui[i]; i++) {	// update with new values
	const char *text = gtk_entry_get_text(GTK_ENTRY(entry[i]));
	switch (ui[i]->type) {
	case 't':  strcpy( ui[i]->text, text );  break;
	case 'f': 
	  ui[i]->fvalues[0] = atof(text); 
	  if (ui[i]->fvalues[1] != 0 || ui[i]->fvalues[2] != 0) {
	    if (ui[i]->fvalues[0] < ui[i]->fvalues[1]) ui[i]->fvalues[0] = ui[i]->fvalues[1];
	    if (ui[i]->fvalues[0] > ui[i]->fvalues[2]) ui[i]->fvalues[0] = ui[i]->fvalues[2];
	  }
	  break;
	case 'd':  
	  ui[i]->dvalues[0] = atof(text); 
	  if (ui[i]->dvalues[1] != 0 || ui[i]->dvalues[2] != 0) {
	    if (ui[i]->dvalues[0] < ui[i]->dvalues[1]) ui[i]->dvalues[0] = ui[i]->dvalues[1];
	    if (ui[i]->dvalues[0] > ui[i]->dvalues[2]) ui[i]->dvalues[0] = ui[i]->dvalues[2];
	  }
	  break;
	case 'i': 
	  ui[i]->ivalues[0] = atoi(text); 
	  if (ui[i]->ivalues[1] != 0 || ui[i]->ivalues[2] != 0) {
	    if (ui[i]->ivalues[0] < ui[i]->ivalues[1]) ui[i]->ivalues[0] = ui[i]->ivalues[1];
	    if (ui[i]->ivalues[0] > ui[i]->ivalues[2]) ui[i]->ivalues[0] = ui[i]->ivalues[2];
	  }
	  break;
	}
      }
      gtk_widget_destroy (dlg);  redraw();  return true;
    } else {					// dialog window cancelled
      gtk_widget_destroy (dlg);  redraw();  return false;
    }
  }
  bool openInputDialogScale(char *title, UserInput *ui) {
    // get user input in hscale
    GtkWidget *dlg=NULL, *label, *scale, *range;
    char minmax[80];  float min, max;
    switch (ui->type) {
    case 'f': min = (float) ui->fvalues[1]; max = (float) ui->fvalues[2]; break;
    case 'd': min = (double)ui->dvalues[1]; max = (double)ui->dvalues[2]; break;
    case 'i': min = (int)   ui->ivalues[1]; max = (int)   ui->ivalues[2]; break;
    default: 
      openMessageDialog("Error", "Invalid type for openInputDialogScale() -- ['f'/'d'/'i']");
      gtk_widget_destroy(dlg);  return false;
    }
    sprintf(minmax, "Min = %g        Max = %g", min, max);
    dlg = gtk_dialog_new_with_buttons(title, GTK_WINDOW(window), GTK_DIALOG_MODAL,
				      GTK_STOCK_OK, GTK_RESPONSE_ACCEPT, 
				      GTK_STOCK_CANCEL, GTK_RESPONSE_CANCEL, NULL);
    label = gtk_label_new(ui->label);
    range = gtk_label_new(minmax);
    scale = gtk_hscale_new_with_range( min, max, 1000 );
    gtk_scale_set_draw_value(GTK_SCALE(scale), true);
    gtk_scale_set_value_pos (GTK_SCALE(scale), GTK_POS_TOP);
    switch (ui->type) {				// show current value
    case 'f': gtk_range_set_value(GTK_RANGE(scale), ui->fvalues[0]);  break;
    case 'd': gtk_range_set_value(GTK_RANGE(scale), ui->dvalues[0]);  break;
    case 'i': gtk_range_set_value(GTK_RANGE(scale), ui->ivalues[0]);  break;
    }
    gtk_container_add(GTK_CONTAINER(GTK_DIALOG(dlg)->vbox), label);
    gtk_container_add(GTK_CONTAINER(GTK_DIALOG(dlg)->vbox), scale);
    gtk_container_add(GTK_CONTAINER(GTK_DIALOG(dlg)->vbox), range);
    label = gtk_label_new("");  // empty space
    gtk_container_add(GTK_CONTAINER(GTK_DIALOG(dlg)->vbox), label);
    gtk_widget_show_all(dlg);			// show the dialog window
    if (gtk_dialog_run(GTK_DIALOG(dlg)) == GTK_RESPONSE_ACCEPT) {
      switch (ui->type) {			// update with new value
      case 'f': ui->fvalues[0] = (float) gtk_range_get_value(GTK_RANGE(scale));  break;
      case 'd': ui->dvalues[0] = (double)gtk_range_get_value(GTK_RANGE(scale));  break;
      case 'i': ui->ivalues[0] = (int)   gtk_range_get_value(GTK_RANGE(scale));  break;
      }
      gtk_widget_destroy (dlg);  return true;
    } else {					// dialog window cancelled
      gtk_widget_destroy (dlg);  return false;
    }
  }  
  
  // drawing ---------------------------------------------------------
  void setColorForDrawing(float r, float g, float b) { 
    GdkColor color; 
    color.red   = (guint)(65535 * r);
    color.green = (guint)(65535 * g);
    color.blue  = (guint)(65535 * b);
    gdk_gc_set_rgb_fg_color( area->style->fg_gc[GTK_STATE_NORMAL], &color); 
  }
  void drawLine(int x0, int y0, int x1, int y1) { gdk_draw_line(area->window, area->style->fg_gc[GTK_STATE_NORMAL], x0, y0, x1, y1); }
  void drawRectCenteredAt(int x, int y, int size, bool filled=false) { 
    // Draw a rectangle centered at (x,y)
    gdk_draw_rectangle(area->window, area->style->fg_gc[GTK_STATE_NORMAL], filled, x-size, y-size, 2*size+1, 2*size+1); 
  }
  void drawRect(int x, int y, int w, int h, bool filled=false) { 
    // Draw a rectangle starting at (x,y) 
    gdk_draw_rectangle(area->window, area->style->fg_gc[GTK_STATE_NORMAL], filled, x, y, w, h); 
  }
  void drawCross(int x, int y, int size, char type='+') { 
    if        (type=='+') {
      drawLine(x-size, y, x+size, y); drawLine(x, y-size, x, y+size); 
    } else if (type=='x') {
      drawLine(x-size, y-size, x+size, y+size); drawLine(x+size, y-size, x-size, y+size); 
    } else if (type=='|') {
      int hs = (size/2);
      drawLine(x-hs, y, x+hs, y); drawLine(x, y-size, x, y+size); 
    } else if (type=='-') {
      int hs = (size/2);
      drawLine(x-size, y, x+size, y); drawLine(x, y-hs, x, y+hs); 
    } else {
      int dx=0, dy=0;
      switch (type) {
      case '0':  dx = +size;  dy = -size;  break;
      case '1':  dx = -size;  dy = -size;  break;
      case '2':  dx = -size;  dy = +size;  break;
      case '3':  dx = +size;  dy = +size;  break;
      }
      drawLine( x, y,  x+dx, y );   drawLine( x, y,  x, y+dy );
    }
  }
  
  // -----------------------------------------------------------------
  // Non window-specific functions
  // -----------------------------------------------------------------
  // start main loop
  void runMainLoop(void) { 
    this->present();
    gtk_main(); 
  }
  // idle  function
  void addIdleFunction(void (*idle)(void* data), void* data=NULL) { guih_idle_ftn = idle; g_idle_add(on_idle, data); }
  void removeIdleFunction(void* data) { g_idle_remove_by_data(data); }

};


// -------------------------------------------------------------------
// default callback functions
// -------------------------------------------------------------------

static gboolean on_init( GtkWidget *widget, gpointer udata )
{
  Window* w = (Window*)udata;
  if (w->runOnInit( widget ) == false) return false;	// window-type specific callback tasks
  return TRUE;
}

static gboolean on_close_rq( GtkWidget *widget, GdkEvent* event, gpointer udata )
{
  // return TRUE to cancel, FALSE to destroy
  Window* w = (Window*)udata;
  w->runOnCloseReq();
  return FALSE;
}

static void on_destroy( GtkWidget *widget, gpointer udata )
{
  Window* w = (Window*)udata;
  if (--guih_window_count <= 0) gtk_main_quit(); 
  w->runOnDestroy();		// window-type specific callback tasks
}

static gboolean on_draw( GtkWidget *widget, GdkEventExpose *event, gpointer udata )
{
  Window* w = (Window*)udata;
  if (w->move_ref) w->move(w->move_ref, w->move_below);
  // window-type specific callback tasks (defined for each child class)
  w->runOnDraw( widget, event );  
  // draw caption
  w->draw_caption();
  return TRUE;
}

static gboolean on_keypress( GtkWidget *widget, GdkEventKey* event, gpointer udata )
{
  Window* w = (Window*)udata;
  w->shift = (event->state & GDK_SHIFT_MASK);
  w->ctrl  = (event->state & GDK_CONTROL_MASK);
  w->alt   = (event->state & GDK_MOD1_MASK);
  if      (event->keyval == GDK_Escape) { w->closeWindow(); return TRUE; }
  else if (event->keyval == GDK_Insert) { w->captureScreen();  return TRUE; }
  else if (event->keyval == GDK_Delete) w->setCaption(NULL);
  else {
    int key = remapKeyValue(event->keyval); 
    w->runOnKeyboard( key );  // window-type specific callback tasks
  }
  // redraw the scene
  w->redraw();
  return TRUE;
}

static gboolean on_button_press(GtkWidget *widget, GdkEventButton *event, gpointer udata)
{
  Window *w = (Window*)udata;
  w->mhxy[0] = event->x;  w->mhxy[1] = event->y;  w->mhpos = 2;
  return TRUE;
}

static gboolean on_button_release(GtkWidget *widget, GdkEventButton *event, gpointer udata)
{
  Window *w = (Window*)udata;
  mhxy[mhpos+0] = -1;  mhxy[mhpos+1] = -1;
  if (w->mhpos==2) w->runOnMouse( event->button, MOUSE_CLICKED, w->mhxy );
  else             w->runOnMouse( event->button, MOUSE_DRAGGED, w->mhxy );
  w->mhpos = 0;  w->redraw();
  return TRUE;
}

static gboolean on_motion_notify(GtkWidget *widget, GdkEventMotion *event, gpointer udata)
{
  Window *w = (Window*)udata;
  if (w->mhpos > 0 && w->mhpos < 800-1) {
    w->mhx[w->mhpos++] = event->x;  w->mhy[w->mhpos++] = event->y;
  }
  return TRUE;
}

static gboolean on_idle (gpointer data)
{
  if (guih_idle_ftn) guih_idle_ftn(data);
  return true;
}

static gboolean on_timer(void *data)
{
  Window *w = (Window*)data;
  return w->runOnTimer();
}


}	// end of GUIH namespace


// ===================================================================
#endif	// end of linux version
// ===================================================================


#endif	// GUIH_COMMON_HPP


// ===================================================================
#if 0	// Example code starts here
// ===================================================================
// compile on Linux:  g++ -o guih guih_example.cpp `pkg-config --cflags gtkglext-1.0 --libs gtkglext-1.0` -Wl,-rpath,/usr/local/lib
#include <iostream>
#include "guih_common.hpp"
void init(GUIH::Window *w) { /* add your code */ }
void draw(GUIH::Window *w) 
{
  w->setColorForDrawing( 1.0, 0.0, 0.0 );
  w->drawLine( 50, 50, 200, 100 );	// x0, y0, x1, y1
  w->drawRect( 50, 50, 200, 100 );	// xmin, ymin xmax ymax
}
void keys(GUIH::Window *w, int key)
{
  GUIH::UserInput uia, uib;
  printf("Key pressed with ctrl=%d, alt=%d and shift=%d\n", w->ctrl, w->alt, w->shift);
  switch (key) {
  case GUIH::F1:
    cerr << "=====================================================================" << endl;
    cerr << "    F12 : on/off                                              " << endl;
    cerr << "=====================================================================" << endl;
    break;
  case GUIH::F12: break;
  case 'a': case 'b': case 'A': case '0': case '9':  break;
  case GUIH::ARROW_UP: case GUIH::ARROW_LEFT: break;
  case ' ': case '\t': case '\b': case '?': break;
  case 'i':
    uia.setInteger("integer input", 27, 1, 99);  // default value 24, min 1, max 99
    uib.setText   ("string  input", "default text string to be edited");
    bool ret = w->openInputDialog("Test Dialog", &uia, &uib);
    if (ret) printf("User input : %d %s \n", uia.ivalues[0], uib.text);
    else     printf("User input cancelled \n");
    break;
  default: break;
  }
}
void mouse(GUIH::Window *w, int button, int state, int xy[])
{
  if        (button==GUIH::MOUSE_LEFT && state==GUIH::MOUSE_CLICKED) {
    printf(" clicked  at (%d,%d)\n", xy[0], xy[1] );
  } else if (button==GUIH::MOUSE_LEFT && state==GUIH::MOUSE_DRAGGED) {
    int mpos, mcnt;  for(mpos=mcnt=0; xy[mpos]>=0; mpos+=2) mcnt++;
    printf(" dragged from (%d,%d) to (%d,%d)\n", xy[0], xy[1], xy[mcnt*2-2], xy[mcnt*2-1]);
  }
}
void close(GUIH::Window *w) { /* add your code */ }
void timer(GUIH::Window *w)
{
  static int countdown = 5;
  std::cout << "timer countdown " << countdown << std::endl;
  if (--countdown <= 0) w->removeTimer();
}
void idle(void *data) 
{
  char *str = (char*)data;
  std::cout << str << std::endl;
  //GUIH::removeIdleFunction( data );
}
int  main(void)
{
  // You can create windows as many as you want, and
  // you don't have to specify all the callback functions.
  GUIH::Window w0(500, 400, "Untitled", init, draw, keys, mouse, close);
  GUIH::Window w1(300, 300, "GUIH Window 1", NULL, NULL, NULL, NULL, NULL);
  GUIH::Window w2(300, 300, "GUIH Window 2"); 
  GUIH::Window w3;  w3.openWindow(300, 300, "GUIH Window 3"); 
  w3.setCallbacks(draw, keys, mouse, close);
  // timer function (you can set a timer for each window)
  w0.createTimer( 1.0, timer );  // interval in second
  // idle function (only one idle function per program)
  // (It doesn't matter which window you call addIdleFunction() from.)
  w0.addIdleFunction( idle, (void*)"Nothing to do" );
  // Start the main loop -- press ESC to close each window. 
  // (It doesn't matter which window you call runMainLoop() from.)
  w0.runMainLoop();
  // Close each window by pressing ESC on keyboard.
  // The program will quit when no window is left.
  return EXIT_SUCCESS;
}
// ===================================================================
#endif	// Example code ends here
// ===================================================================
