
//
// UTIL::Clock class object
//
// Jaeil Choi
// last modified in Jan, 2004
//
//
// Example:
//   UTIL::Clock clock;
//   cout << "DateTime : " << clock.getLocalTime() << endl;    
//   // for comparison, 'ctime()' returns "Mon Mar 10 19:22:37 2003\n"
//


#ifndef UTIL_CLOCK_HPP
#define UTIL_CLOCK_HPP


#include <iostream>
#include <iomanip>
#include <ctime>
#include <cstdarg>

// ===================================================================

#ifdef WIN32
#include <Windows.h>
#define sleepd(sec)  Sleep((int)(sec*1000))
#else
#include <sys/time.h>
#define sleepd(sec)  \
  do { \
    struct timespec req, rem; \
    req.tv_sec = (time_t)((long)sec); \
    req.tv_nsec = (long)((sec - (double)req.tv_sec) * 1000000000); \
    nanosleep( &req, &rem ); \
  } while (0)
#endif

namespace UTIL {
  
class Clock 
{
private:
  char    bufdummy[40];
public:
  Clock() {}
  ~Clock() {}
  
  // -----------------------------------------------------------------
  // Local Date & Time
  // -----------------------------------------------------------------
public:  
  inline int Y (struct tm *n) { return n->tm_year + 1900; }
  inline int YY(struct tm *n) { return n->tm_year %  100; }
  inline int M (struct tm *n) { return n->tm_mon+1; }
  inline int D (struct tm *n) { return n->tm_mday;  }
  inline int WD(struct tm *n) { return n->tm_wday;  }  // day in the week
  inline int YD(struct tm *n) { return n->tm_yday;  }  // day in the year
  inline int h (struct tm *n) { return n->tm_hour;  }
  inline int m (struct tm *n) { return n->tm_min;   }
  inline int s (struct tm *n) { return n->tm_sec;   }
  inline int ms(void) { 
#ifdef WIN32
    return 0;
#else
    struct timeval tv; struct timezone tz; gettimeofday(&tv, &tz); return (int)((tv.tv_usec%1000000)/1000); 
#endif
  }

  char *getLocalTime(int format=0, char *buf=NULL) {
    time_t     t = time(NULL);
    struct tm *n = localtime(&t);
    if (buf == NULL) buf = this->bufdummy;
    switch (format) {
    case 0: sprintf(buf,"%02d/%02d/%04d %02d:%02d:%02d", M(n),D(n),Y(n),h(n),m(n),s(n)); break; // MM/DD/YYYY hh:mm:ss
    case 1: sprintf(buf,"%02d/%02d/%04d", M(n), D(n), Y(n)); break; // MM/DD/YYYY
    case 2: sprintf(buf,"%02d:%02d:%02d", h(n), m(n), s(n)); break; // hh:mm:ss
    case 3: sprintf(buf,"%02d:%02d:%02d %02d/%02d/%04d", h(n),m(n),s(n),M(n),D(n),Y(n)); break; // hh:mm:ss MM/DD/YYYY
    case 4: sprintf(buf,"%02d:%02d, %02d/%02d/%04d", h(n),m(n),M(n),D(n),Y(n)); break;   // hh:mm, MM/DD/YYYY
    default:sprintf(buf,"%02d/%02d %02d:%02d:%02d", M(n),D(n),h(n),m(n),s(n)); break; // MM/DD hh:mm:mm
    }
    return buf;
  }
  
  char* getDateTimeString(char *buf, char *head=NULL, char *tail=NULL) {
    // Get the data-and-time string such as "20091231_2359_5900".
    time_t t = time(NULL);
    struct tm *n = localtime(&t);
    if (buf == NULL) buf = this->bufdummy;
    sprintf(buf, "%s%04d%02d%02d_%02d%02d_%02d%02d%s", (head ? head : ""),
	    Y(n), M(n), D(n), h(n), m(n), s(n), ms()/10, 
	    (tail ? tail : ""));
    return buf;
  }
  
  // -----------------------------------------------------------------
  // Etc.
  // -----------------------------------------------------------------

  double getTimeInSec(void) {
    // Get current time in seconds (from the Epoch)
#ifdef WIN32
    FILETIME ft;
    unsigned __int64 tmpres = 0;
    static int tzflag;
    GetSystemTimeAsFileTime(&ft);
    tmpres |= ft.dwHighDateTime;
    tmpres <<= 32;
    tmpres |= ft.dwLowDateTime;
    //converting file time to unix epoch
    tmpres /= 10;  //convert into microseconds
    tmpres -= DELTA_EPOCH_IN_MICROSECS; 
    //   tv->tv_sec = (long)(tmpres / 1000000UL);
    //   tv->tv_usec = (long)(tmpres % 1000000UL);
    return (double)((long)(tmpres / 1000000UL) + (long)(tmpres % 1000000UL)*1.0e-6);
#else
    struct timeval tv;  struct timezone tz; 
    gettimeofday(&tv, &tz);
    return (tv.tv_sec + tv.tv_usec * 1.0e-6);
#endif
  }
  
};


}	// end of namespace UTIL

#endif  // UTIL_CLOCK_HPP

