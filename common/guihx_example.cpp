
// 
// Example program for using GUIH
// 
// Jaeil Choi, 2008
// 
// How to build this program:
//   On Mac OSX, link it with some frameworks: Carbon, AGL, OpenGL, GLUT
//     g++ -o guih guihx_example.cpp -framework Carbon -framework AGL -framework OpenGL -framework GLUT
//   On Linux, link it with a libraries: GTKGLExt
//     g++ -o guih guihx_example.cpp `pkg-config --cflags gtkglext-1.0 --libs gtkglext-1.0` -Wl,-rpath,/usr/local/lib
// 

#include <iostream>
#include "guih_image.hpp"
#include "guih_opengl2d.hpp"
#include "guih_opengl3d.hpp"
using namespace std;

void mouse(GUIH::Window *w, int button, int state, int xy[])
{
  // Callback function for handling mouse events (CLICKED or DRAGGED)
  if (button == GUIH::MOUSE_LEFT && state == GUIH::MOUSE_CLICKED)
    printf(" mouse %d event %d at (%d %d)\n", button, state, xy[0], xy[1] );
}
void draw(GUIH::Window *w)
{
  // Callback function for drawing events
  GUIH::OpenGL3DWindow *w3d = (GUIH::OpenGL3DWindow *) w;
  glColor3ub(255, 100, 100);
  w3d->drawTorus( 5, 2, 20, 12 );  // out_radius, in_radius, slices, steps
  glColor3ub(0, 255, 100);
  glBegin(GL_LINES);
  glVertex3f(0, 0, 0);
  glVertex3f(5, 5, 5);
  glEnd();
  glColor3ub(0, 0, 255);
  glPushMatrix();
  glTranslatef( 1, 1, 5 );
  w3d->drawCube( 0.5 );
  glPopMatrix();
}
void keys(GUIH::Window *w, int key)
{
  // Callback function for handling keyboard events
  printf("key pressed: keyAsValue=%d  keyAsChar='%c', (ctrl=%c alt=%c shift=%c)\n",
	 key, (char)key, (w->ctrl ? 'Y':'N'), (w->alt ? 'Y':'N'), (w->shift ? 'Y':'N'));
}
void timer(GUIH::Window *w)
{
  // Callback function for timer events
  GUIH::OpenGL3DWindow *w3 = (GUIH::OpenGL3DWindow*) w;
  w3->camera.circle( +2.0, 0 );
  w3->redraw();
  static int count= 0;
  if (++count >= 90) w->removeTimer();
}

void draw_mandelbrot(int w, int h, unsigned char img[], double cx, double cy, double scale, int max_iteration=0)
{
  if (max_iteration <= 0) max_iteration = (int)(scale * 10);
  for (int row=0; row < h; row++) {
    for (int col=0,iter=0; col < w; col++, img++) {
      double x=0, x0 = cx + (col - w/2) / scale;
      double y=0, y0 = cy + (row - h/2) / scale;
      for (iter=0; x*x+y*y <= (2*2) && iter < max_iteration; iter++) {
	double xtemp = x*x - y*y + x0;
	y = 2*x*y + y0;
	x = xtemp;
      }
      *img = (iter == max_iteration ? 0 : 255);
    }
  }
}

int main(int argc, char **argv)
{
  // A window for displaying an image --------------------------------
  GUIH::ImageWindow wa(600, 300, "Image Window", NULL, NULL, keys, mouse);
  if (argc >= 2) {
    wa.showImage( argv[1] );
  } else {
    printf("To open an existing image, use \"%s <filename>\"\n", argv[0]);
    unsigned char *image = (unsigned char*)malloc(1*600*300*sizeof(unsigned char));
    draw_mandelbrot( 600, 300, image, -0.5, 0, 100 );
    wa.showImage( 600, 300, GUIH::PIXEL_GRAY, image );
  }
  // 2D OpenGL window ------------------------------------------------
  GUIH::OpenGL2DWindow w2(300, 300, "OpenGL2D Window");
  w2.initCamera(-100, -100, 200, 200);
  w2.show_axes = true;
  // 3D OpenGL window ------------------------------------------------
  GUIH::OpenGL3DWindow w3(300, 300, "OpenGL3D Window", NULL, draw, keys, mouse);
  w3.initCamera( 0.0,5.0,30.0, 0,0,0, 0,1,0 );
  w3.show_axes = true;
  w3.createTimer( 0.1, timer );
  // start the main loop ---------------------------------------------
  w2.move( &wa, true );
  w3.move( &w2, false );
  wa.runMainLoop();
  return EXIT_SUCCESS;
}
