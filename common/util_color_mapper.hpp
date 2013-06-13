
//
// UTIL::ColorMapper<> class template
//
// Jaeil Choi
// last modified in June, 2007
//
//
// This class is for the mapping a float value [min~max] to various colors.
// Type of template (type of r/g/b) is either 'unsigned char' or 'int'.
// RGB values have range of [0~255], and they are predefined by macro.
//
// Example1: UTIL::ColorMapper<int> cmap("rgb", 0, 1);
//           cmap.getColor( 0.5, &r, &g, &b );
// Example2: UTIL::ColorMapper<int> cmap;
//           cmap.setPalette("wyrgb", 0, 100, true);  // colors,min,max,repeat
//           cmap.getColor( 50, rgb );
//


#ifndef UTIL_COLOR_MAPPER_HPP
#define UTIL_COLOR_MAPPER_HPP

#include <iostream>

namespace UTIL {
  

template <class T>
class ColorMapper {	// <T> is the type of color
public:
  T      colors[21][3];
  int    ncolors;
  double min, max;
  bool   repeat;
public:
  ColorMapper() { setPalette("bCgYr", 0.0, 1.0, false); }
  ColorMapper(char *cstr, double min=0, double max=1, bool repeat=false) { setPalette(cstr, min, max, repeat); }
  ~ColorMapper() {}
  
  // -----------------------------------------------------------------
  // Setting a color palette
  // -----------------------------------------------------------------
public:
  void setPalette(char *cstr, double min, double max, bool repeat=false) {
    // Set a color palette for a value between 'min' and 'max',
    //   using a sequence of characters that represent a color.
    if (!cstr || min >= max) { printf("Error(ColorMapper::setPalette): invalid arguments\n"); return; }
    this->min = min;		this->max = max;
    this->repeat   = repeat;
    if      (strcmp(cstr, "rainbow"    )==0) cstr = "r{}g[]bp";   // rainbow
    else if (strcmp(cstr, "temp"       )==0) cstr = "bMr";        // blue-magenta-red
    else if (strcmp(cstr, "temperature")==0) cstr = "bcfyr";      // blue-Dcyan-Fgreen-Dyellow-red
    else if (strcmp(cstr, "heatsource" )==0) cstr = "bgryw";      // blue-green-red-yellow-white
    else if (strcmp(cstr, "grayscale"  )==0) cstr = "0123456789"; // black-white grayscale
    int i, len=strlen(cstr);
    for (i = 0; i < len && i < 20; i++) {
      switch (cstr[i]) {
      case '0':  setPaletteColor(i,   0,   0,   0);  break; // black
      case '1':  setPaletteColor(i,  28,  28,  28);  break; // grayscale
      case '2':  setPaletteColor(i,  57,  57,  57);  break; 
      case '3':  setPaletteColor(i,  85,  85,  85);  break;
      case '4':  setPaletteColor(i, 113, 113, 113);  break;
      case '5':  setPaletteColor(i, 142, 142, 142);  break;
      case '6':  setPaletteColor(i, 170, 170, 170);  break;
      case '7':  setPaletteColor(i, 198, 198, 198);  break;
      case '8':  setPaletteColor(i, 227, 227, 227);  break;
      case '9':  setPaletteColor(i, 255, 255, 255);  break; // white
      case 'b':  setPaletteColor(i,   0,   0, 255);  break; // blue
      case 'B':  setPaletteColor(i,  50,  50, 255);  break; // blue bright
      case 'c':  setPaletteColor(i,   0, 255, 255);  break; // cyan (aqua)
      case 'C':  setPaletteColor(i,   0, 150, 150);  break; // cyan dark
      case 'e':  setPaletteColor(i, 245, 245, 220);  break; // beige
      case 'f':  setPaletteColor(i,  70, 160,  70);  break; // forest green
      case 'g':  setPaletteColor(i,   0, 255,   0);  break; // green
      case 'G':  setPaletteColor(i,  50, 255,  50);  break; // green bright
      case 'k':  setPaletteColor(i, 255, 192, 203);  break; // pink
      case 'K':  setPaletteColor(i, 135, 206, 235);  break; // sky blue
      case 'l':  setPaletteColor(i, 230, 230, 250);  break; // lavender
      case 'm':  setPaletteColor(i, 255,   0, 255);  break; // magenta
      case 'M':  setPaletteColor(i, 150,   0, 150);  break; // magenta dark
      case 'n':  setPaletteColor(i,   0,   0, 128);  break; // navy
      case 'N':  setPaletteColor(i, 165,  42,  42);  break; // brown
      case 'o':  setPaletteColor(i, 255, 165,   0);  break; // orange
      case 'p':  setPaletteColor(i, 128,   0, 128);  break; // purple
      case 'r':  setPaletteColor(i, 255,   0,   0);  break; // red
      case 'R':  setPaletteColor(i, 255,  50,  50);  break; // red bright
      case 's':  setPaletteColor(i, 192, 192, 192);  break; // silver
      case 't':  setPaletteColor(i,  64, 224, 208);  break; // turquoise
      case 'T':  setPaletteColor(i,   0, 128, 128);  break; // teal
      case 'v':  setPaletteColor(i, 238, 130, 238);  break; // violet
      case 'V':  setPaletteColor(i, 128, 128,   0);  break; // olive
      case 'w':  setPaletteColor(i, 255, 255, 255);  break; // white
      case 'y':  setPaletteColor(i, 255, 255,   0);  break; // yellow
      case 'Y':  setPaletteColor(i, 150, 150,   0);  break; // yellow dark
      case '{':  setPaletteColor(i, 200, 128,   0);  break; // 
      case '}':  setPaletteColor(i, 128, 200,   0);  break; // 
      case '[':  setPaletteColor(i,   0, 200, 128);  break; // 
      case ']':  setPaletteColor(i,   0, 128, 200);  break; // 
      default :  setPaletteColor(i,   0,   0,   0);  break;
      }
    }
    ncolors = (len < 20 ? len : 20);
    memcpy( colors+ncolors, colors, 3*sizeof(T) );  // for smooth repetition
  }
  inline void setPaletteColor(int i, T r, T g, T b) {
    colors[i][0] = r;  colors[i][1] = g;  colors[i][2] = b;
  }
  
  // -----------------------------------------------------------------
  // getting color values
  // -----------------------------------------------------------------
public:  
  void getColor(double value, T rgb[]) { getColor(value, rgb+0, rgb+1, rgb+2); }
  void getColor(double value, T *r, T *g, T *b) {
    // Get the color that correspond to the value in the range of [ min ~ max ].
    int nc = (repeat ? ncolors + 1 : ncolors);
    if (min  >= max) { min = 0; max = 1; }
    if (repeat) {
      while (value < min) value += (max - min);
      while (value > max) value -= (max - min);
    } else {
      if (value < min) value = min;
      if (value > max) value = max;
    }
    double t = (value - min) / (max - min);
    if (nc > 1) {
      int    idx, i, j;
      double off, step = 1.0 / (nc-1);
      idx = (int)(t / step);
      i = idx; j = i + 1;
      off = (t - (idx * step)) / step;
      *r = (T)((1-off) * colors[i][0] + (off) * colors[j][0]);
      *g = (T)((1-off) * colors[i][1] + (off) * colors[j][1]);
      *b = (T)((1-off) * colors[i][2] + (off) * colors[j][2]);
    } else {
      *r = *g = *b = (T)(255 * t);
    }
  }
  
  void getRandomColorWithSeed(int seed, T rgb[3], int luminance=200) { 
    // Get a random color with specifed luminance, using the provided seed.
    float c[3], avg, ratio;
    srand( seed );
    c[0] = rand() % 256;
    c[1] = rand() % 256;
    c[2] = rand() % 256;
    avg = (c[0] + c[1] + c[2])/3;
    if (avg == 0) rgb[0] = rgb[1] = rgb[2] = (T)luminance;
    else {
      ratio = luminance / avg;
      c[0] *= ratio;  c[1] *= ratio;  c[2] *= ratio;
      rgb[0] = (c[0] > 255 ? 255 : (T)c[0]);
      rgb[1] = (c[1] > 255 ? 255 : (T)c[1]);
      rgb[2] = (c[2] > 255 ? 255 : (T)c[2]);
    }
  }
};

}	// end of namespace UTIL

#endif // UTIL_COLOR_MAPPER_HPP
