
//
// GUIH       : GUI in Header files
// GUIH::Args : command-line argument processor
//
// (c) 2006  Jaeil Choi
// last modified in May, 2007
//
// --------------------------------------------------------------------------
//  This program is free software; you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation; either version 2 of the License, or
//  (at your option) any later version.
//
//  This program is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
//
//  You should have received a copy of the GNU General Public License
//  along with this program; if not, write to the Free Software
//  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
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
//   Args                : command-line argument processor
//   Config              : configuration file processor
//
// For usage, read the example code at the bottom.
//
  

#ifndef GUIH_ARGUMENTS_HPP
#define GUIH_ARGUMENTS_HPP

#include <iostream>
#include <iomanip>
#include <cfloat>

namespace GUIH {
  
  
#define ARG_CHAR	1
#define ARG_INT		2
#define ARG_FLOAT	3
#define ARG_DOUBLE	4
#define ARG_STRING	5
#define ARG_BOOL	6

#define ARG_UNSET (-FLT_MIN)
#define ARG_DSTR  (-FLT_MIN*2)
#define ARG_NAME_MAX_LEN 20

typedef struct {
  char  name[ARG_NAME_MAX_LEN+1];
  bool  processed;
  bool  range_set;
  int   type;
  int   count;
  void  *p;
  double range[2];
} arg_t;


class Args {
public:
  int  argc;		// number of inputs
  char **argv;		// list of input (after args processed)
  int  idx;		// index of current input file (initially 0)
  char *prg;		// program (argv[0])
private:
  arg_t *args;		// list of user-specified arguments
  int   nargs;		// number of user-specified arguments
  char  *default_str;
public:
  bool cfg;		
public:
  Args()  { 
    args = NULL;  nargs = 0; 
    argv = NULL;  argc = idx = 0;
    prg  = NULL;  cfg = false;
  }
  ~Args() { if (args) free(args); }
  
public:
  int setOpt(char *name, int type, void *p, double defaultv=ARG_UNSET) {
    return setOpt(name, type, p, 1, defaultv, ARG_UNSET, ARG_UNSET);
  }
  int setOpt(char *name, int type, void *p, int defaultv) {
    return setOpt(name, type, p, 1, (double)defaultv, ARG_UNSET, ARG_UNSET);
  }
  int setOpt(char *name, int type, void *p, char *defaultv) {
    default_str = defaultv;
    return setOpt(name, type, p, 1, ARG_DSTR, ARG_UNSET, ARG_UNSET);
  }
  int setOpt(char *name, int type, void *p, double defaultv, double min, double max) {
    return setOpt(name, type, p, 1, defaultv, min, max);
  }
  int setOpt(char *name, int type, void *p, int count, double defaultv) {
    return setOpt(name, type, p, count, defaultv, ARG_UNSET, ARG_UNSET);
  }
  int setOpt(char *name, int type, void *p, int count, double defaultv, double min, double max) {
    int idx = findFlag(name), i;
    if (idx < 0) {
      args = (arg_t*)realloc(args, (nargs+1) * sizeof(arg_t));
      idx = nargs++;
      args[idx].processed = false;
    }
    strncpy(args[idx].name, name, ARG_NAME_MAX_LEN);
    args[idx].name[ARG_NAME_MAX_LEN] = '\0';
    args[idx].type = type;
    args[idx].count = count;
    args[idx].p = p;
    args[idx].range_set = false;
    args[idx].range[0] = min;  args[idx].range[1] = max;
    // set default values
    if        (type == ARG_BOOL)
      for (i = 0; i < args[idx].count; i++) ((bool*)(args[idx].p))[i] = (defaultv == 0 || defaultv == ARG_UNSET ? false : true);
    else if (type == ARG_STRING)
      for (i = 0; i < args[idx].count; i++) ((char**)(args[idx].p))[i] = (defaultv == ARG_DSTR ? default_str : NULL);
    else {
      for (i = 0; i < args[idx].count; i++) {
	switch (args[idx].type) {
	case ARG_INT:    *((int*)   (args[idx].p)+i) = (int)   defaultv; break;
	case ARG_FLOAT:  *((float*) (args[idx].p)+i) = (float) defaultv; break;
	case ARG_DOUBLE: *((double*)(args[idx].p)+i) = (double)defaultv; break;
	case ARG_CHAR:   *((char*)  (args[idx].p)+i) = (char)  defaultv; break;
	}
      }
      if (min != ARG_UNSET || max != ARG_UNSET) args[idx].range_set = true;
    }
    return idx;
  }
  
  bool processed(char *name) {
    int idx = findFlag(name);
    if (idx >= 0 && args[idx].processed) return true;
    return false;
  }
  
  bool process(int argc, char **argv, int min_argc=0, bool leave_unknowns_untouched=false) {
    int  pos = 1, idx, i;
    bool show = false;
    prg = argv[0];
    while (pos < argc) {
      if (argv[pos][0] != '-') break;
      idx = findFlag( &(argv[pos][1]) );
      if (idx < 0) {
	if      (strcmp(argv[pos], "-ARG")==0) { show = true; pos++; continue; }
	else if (strcmp(argv[pos], "-INI")==0) { cfg  = true; pos++; continue; }
	else if (leave_unknowns_untouched) {
	  this->argv = (pos >= argc ? NULL : argv + pos);
	  this->argc = argc - pos;  idx = 0;
	  break;
	} else std::cout << "Invalid argument '" << argv[pos] << "'" << std::endl;
	return false;
      }
      //std::cout << "pos = " << pos << "  argv[pos] = '" << argv[pos] << "'" << std::endl;
      for (i = 0; i < args[idx].count; i++) {
	// check 
	if (pos >= argc-1 && 
	    (args[idx].type == ARG_CHAR || args[idx].type == ARG_INT ||
	     args[idx].type == ARG_FLOAT || args[idx].type == ARG_DOUBLE ||
	     args[idx].type == ARG_STRING)) return false;
	// read values
	int    *ip = ((int*)   (args[idx].p)) + i;  
	float  *fp = ((float*) (args[idx].p)) + i;
	double *dp = ((double*)(args[idx].p)) + i;
	char   *cp = ((char*)  (args[idx].p)) + i;
	bool   *bp = ((bool*)  (args[idx].p)) + i;
	char  **sp = ((char**) (args[idx].p)) + i;
	//std::cout << "reading args : '" << argv[pos+1] << "'" << std::endl;
	switch (args[idx].type) {
	case ARG_CHAR  : *cp = argv[++pos][0];     break;
	case ARG_INT   : *ip = atoi(argv[++pos]);  break;
	case ARG_FLOAT : *fp = (float)atof(argv[++pos]);  break;
	case ARG_DOUBLE: *dp = atof(argv[++pos]);  break;
	case ARG_STRING: *sp = argv[++pos]; break;
	case ARG_BOOL  : 
	  if (pos < argc && argv[pos][1] == '\0') {
	    if      (argv[pos][0] == 't' || argv[pos][0] == 'T') { *bp = true; pos++; }
	    else if (argv[pos][0] == 'f' || argv[pos][0] == 'F') { *bp = true; pos++; }
	    else *bp = true; 
	  } else *bp = true; 
	  break;
	}
	// check value range
	if (args[idx].range_set) {
	  switch (args[idx].type) {
	  case ARG_INT:
	    if (*ip < (int)(args[idx].range[0])) *ip = (int)(args[idx].range[0]);
	    if (*ip > (int)(args[idx].range[1])) *ip = (int)(args[idx].range[1]);
	    break;
	  case ARG_FLOAT:
	    if (*fp < args[idx].range[0]) *fp = (float)args[idx].range[0];
	    if (*fp > args[idx].range[1]) *fp = (float)args[idx].range[1];
	    break;
	  case ARG_DOUBLE:
	    if (*dp < args[idx].range[0]) *dp = args[idx].range[0];
	    if (*dp > args[idx].range[1]) *dp = args[idx].range[1];
	    break;
	  case ARG_CHAR:
	    if (*cp < (char)(args[idx].range[0])) *cp = (char)(args[idx].range[0]);
	    if (*cp > (char)(args[idx].range[1])) *cp = (char)(args[idx].range[1]);
	    break;
	  }
	}
      }
      args[idx].processed = true;
      pos++;
    }
    if (leave_unknowns_untouched) {
      this->argv = argv + pos - 1;
      this->argc = argc - pos + 1;
    } else {
      this->argv = (pos >= argc ? NULL : argv + pos);
      this->argc = argc - pos;  idx = 0;
    }
    if (show) printInfo();
    if (min_argc>0 && this->argc<min_argc) return false;
    return true;
  }
  
  char* getFileExtension(char *filename) {
    int  i;  // find the beginning of the file extension of 'filename'
    for (i=strlen(filename); i>0 && filename[i]!='.'; i--);
    return (i>0 ? filename+i+1 : (char*)"");
  }
  
  void printInfo(void) {
    int  idx, i;
    char buf[256];
    void *v;
    std::cout << "GUIH::Args Information for " << prg << std::endl;
    // args
    std::cout << "  User specified args" << std::endl;
    std::cout.setf(std::ios::left);
    for (idx = 0; idx < nargs; idx++) {
      v = args[idx].p;
      if (args[idx].count == 1)  sprintf(buf, "%s", args[idx].name);
      else  sprintf(buf, "%s[%d]", args[idx].name, args[idx].count);
      std::cout << "    -" << std::setw(ARG_NAME_MAX_LEN) << buf << " ";	// names
      buf[0] = '\0';
      for (i = 0; i < args[idx].count; i++) {
	switch (args[idx].type) {
	case ARG_CHAR:   
	  if (((char*)(v))[i])  sprintf(buf+strlen(buf), "%c ", ((char*)(v))[i]); 
	  else                  sprintf(buf+strlen(buf), "%s ", "(0)"); 
	  break;
	case ARG_INT:    sprintf(buf+strlen(buf), "%d ", ((int*)(v))[i]); break;
	case ARG_FLOAT:  sprintf(buf+strlen(buf), "%g ", (((float*)(v))[i] == ARG_UNSET ? 0 : ((float*)(v))[i])); break;
	case ARG_DOUBLE: sprintf(buf+strlen(buf), "%g ", (((double*)(v))[i] == ARG_UNSET ? 0 : ((double*)(v))[i])); break;
	case ARG_STRING: sprintf(buf+strlen(buf), "%s ", ((char**)(v))[i]); break;
	case ARG_BOOL:   sprintf(buf+strlen(buf), "%s ", (((bool*)(v))[i] ? "true" : "false")); break;
	default: buf[strlen(buf)] = '\0'; break;
	}
      }
      std::cout << std::setw(28) << buf << " ";				// value(s)
      std::cout << (args[idx].processed ? "(Y)  " : "( )  ");	// processed ?
      switch (args[idx].type) {				// type
      case ARG_CHAR:   std::cout << std::setw(6) << "char"   << " "; break;
      case ARG_INT:    std::cout << std::setw(6) << "int"    << " "; break;
      case ARG_FLOAT:  std::cout << std::setw(6) << "float"  << " "; break;
      case ARG_DOUBLE: std::cout << std::setw(6) << "double" << " "; break;
      case ARG_STRING: std::cout << std::setw(6) << "string" << " "; break;
      case ARG_BOOL:   std::cout << std::setw(6) << "bool"   << " "; break;
      default:         std::cout << std::setw(6) << ""       << " "; break;
      }
      if (args[idx].range_set) {				// Range[min~max]
	double *r = args[idx].range;
	std::cout << "Range[";
	switch (args[idx].type) {
	case ARG_CHAR:   std::cout << (char)  r[0] << " ~ " << (char)  r[1]; break;
	case ARG_INT:    std::cout << (int)   r[0] << " ~ " << (int)   r[1]; break;
	case ARG_FLOAT:  std::cout << (float) r[0] << " ~ " << (float) r[1]; break;
	case ARG_DOUBLE: std::cout << (double)r[0] << " ~ " << (double)r[1]; break;
	}
	std::cout << "]";
      }
      std::cout << std::endl;
    }
    // print inputs
    std::cout << "  Input strings : total " << this->argc << "  (current=" << this->idx << ")" << std::endl;
    for (i = 0; i < this->argc; i++) {
      std::cout << "    " << this->argv[i];
      if (i % 5 == 4) std::cout << std::endl;
    }
    if (i%5 != 0) std::cout << std::endl;
  }
  
private:
  int  findFlag(char *name) {
    int i;
    for (i = 0; i < nargs; i++) 
      if (strncmp(args[i].name, name, ARG_NAME_MAX_LEN) == 0) break;
    return (i == nargs ? -1 : i);
  }
  
};


}	// end of GUIH namespace


#endif	// GUIH_ARGUMENTS_HPP


// ===================================================================
#if 0	// Example code starts
// ===================================================================
#include <iostream>
#include "guih_args.hpp"
using namespace std;
void usage(char *program, bool bExit)
{
  cerr << program << ": Example of using GUIH::Args" << endl;
  cerr << "Usage:  " << program << " [OPTIONS] input_files" << endl;
  cerr << "        -i0 n     : integer value [0]" << endl;
  cerr << "        -i1 n     : integer value [2]" << endl;
  cerr << "        -i2 n     : integer value [50]" << endl;
  cerr << "        -i3 a b c : integer values [50]" << endl;
  cerr << "        -f value  : float value [27]" << endl;
  cerr << "        -d value  : double value [6]" << endl;
  cerr << "        -s string : string value [NULL]" << endl;
  cerr << "        -v        : boolean value [false]" << endl;
  if (bExit) exit(EXIT_FAILURE);
}
int main(int argc, char **argv) 
{
  GUIH::Args arg;
  int    vi0, vi1, vi2, vi3[3];
  float  vf;
  double vd;
  char*  str;
  bool   verbose;
  // Types can be either ARG_BOOL, ARG_CHAR, ARG_INT, ARG_FLOAT, ARG_DOUBLE, or ARG_STRING
  // If default value is not given, the variable is set to 0/NULL/false.
  arg.setOpt("i0", ARG_INT, &vi0);              // single value ('vi0' is set to 0)
  arg.setOpt("i1", ARG_INT, &vi1, 2);           // single value with default
  arg.setOpt("i2", ARG_INT, &vi2, 50, 0, 50);   // single value with default and min/max range
  arg.setOpt("i3", ARG_INT, vi3, 3, 50, 0, 50); // multiple values with default and range
  arg.setOpt("f", ARG_FLOAT, &vf, 27);          // single value with default
  arg.setOpt("d", ARG_DOUBLE, &vd, 6, 0, 10);   // single value with default and min/max range
  arg.setOpt("s", ARG_STRING, &str);     // string ('str' is set to NULL)
  arg.setOpt("v", ARG_BOOL, &verbose);   // boolean variable ('verbose' is set to false)
  if (!arg.process(argc, argv)) usage(argv[0], true);  // without restriction on the number of arguments
  //if (!arg.process(argc, argv, 1)) usage(argv[0], true);  // minimum number of arguments == 1
  arg.printInfo();  // show how the arguments were processed
  // You can use the rest of the arguments. (Those which do not begin with '-')
  for (int i = 0; i < arg.argc; i++)
    printf("Input file %d : %s \n", i, arg.argv[i]);
  // To tell whether it was actually given the default value or 
  // it wasn't found at all in the arguments, use processed().
  if (arg.processed("i3")) printf("'-i3' flag was present\n");
  return EXIT_SUCCESS;
}
// ===================================================================
#endif	// Example code ends
// ===================================================================
