
//
// UTIL::Timer class object
//
// Jaeil Choi
// last modified in Jan, 2004
//
//
// Example:
//   UTIL::Timer timer;
//   timer.start();	
//   double t  = timer.pause();	
//   timer.resume();	
//   double tt = timer.stop();	
//


#ifndef UTIL_TIMER_HPP
#define UTIL_TIMER_HPP


#include <iostream>
#include <iomanip>
#include <ctime>
#include <cstdarg>

// ===================================================================

#ifdef WIN32
#include <Windows.h>
#define sleepd(sec)  Sleep((int)(sec*1000))
#if defined(_MSC_VER) || defined(_MSC_EXTENSIONS) 
#define DELTA_EPOCH_IN_MICROSECS  11644473600000000Ui64 
#else 
#define DELTA_EPOCH_IN_MICROSECS  11644473600000000ULL 
#endif
#else
#include <sys/time.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#define sleepd(sec)  \
  do { \
    struct timespec req, rem; \
    req.tv_sec = (time_t)((long)sec); \
    req.tv_nsec = (long)((sec - (double)req.tv_sec) * 1000000000); \
    nanosleep( &req, &rem ); \
  } while (0)
#endif

// ===================================================================

namespace UTIL {
  
  
// ===================================================================
  
typedef struct { char name[10]; double mx; double tt; int cnt; } sw_t;
  
class Timer 
{
private:
  bool    running;
  double  accumulated_time;		// in sec
  clock_t last_start_time;
  
  int     sw_size;
  sw_t   *sw_data;
  
public:
  Timer() : running(false), accumulated_time(0), last_start_time(0),
	    sw_size(0), sw_data(NULL) {}
  ~Timer() { if (sw_data) free(sw_data); }
  
  // -----------------------------------------------------------------
  // stopwatch (processor time in seconds)
  // -----------------------------------------------------------------
public:  
  void start(const char* tname = NULL)  {
    if (tname) fprintf(stdout, "Timer Start   :       [%s] \n", tname);  // in sec
    running = true;
    accumulated_time = 0;
    last_start_time  = clock();
    //printf("started last_start_time=%g\n", (double)last_start_time);
  }
  double stop(const char* tname = NULL)  {
    // return the total accumulated time (sec) the timer has run
    if (running) {
      double elapsed = getElapsedTime();
      accumulated_time += elapsed;
      //printf("stopped accumulated=%g\n", accumulated_time);
      if (tname) fprintf(stdout, "Timer Stop    : %7.4f [%s]  total=%7.4f\n", elapsed, tname, accumulated_time);
    }
    running = false;       // change timer status to stopped
    return accumulated_time;
  }
  double pause(const char* tname = NULL) {
    // return elaspsed time (sec) since the the last resume/check/start
    if (!running) return 0;   
    double elapsed = getElapsedTime();
    accumulated_time += elapsed;
    running = false;
    if (tname) fprintf(stdout, "Timer Paused  : %7.4f [%s] \n", elapsed, tname);
    return elapsed;
  }
  void resume(const char* tname = NULL)  {
    if (tname) fprintf(stdout, "Timer Resumed :       [%s] \n", tname);
    running = true;        // change timer status to running
    last_start_time = clock();  // Determine the start time;
  }
  double checkTime(const char* tname = NULL) {
    // return elaspsed time (sec) since the the last resume/check/start
    if (!running) return 0;   
    double elapsed = getElapsedTime();
    accumulated_time += elapsed;
    //printf("  checked elapsed = %g  accumulated=%g\n", elapsed, accumulated_time);
    if (tname) fprintf(stdout, "Timer Checked : %7.4f [%s] \n", elapsed, tname);
    last_start_time = clock();  // Determine the start time;
    return elapsed;
  }
  
  inline double getElapsedTime(void)  {
    // Determine how much time (sec) has passed since the last resume/check/start
    clock_t curr_time = clock();
    double  elapsed = double(curr_time - last_start_time) / CLOCKS_PER_SEC;
    //printf("  elapsed = %g  curr=%g  last=%g\n", elapsed, (double)curr_time, (double)last_start_time);
    return elapsed;
  }
  double getTotalTime(void)  {
    if (running) return accumulated_time + getElapsedTime();
    else         return accumulated_time;
  }
  
  // -----------------------------------------------------------------
  // stopwatch with internal storage,
  //   in order to measure the speed of multiple tasks REPEATEDLY
  // -----------------------------------------------------------------
public:
  void setupSWatch(int n, ... ) {
    // Set up a stopwatch to measure performance of parts of a program,
    //   with the number ('n') and the names (char*) of all the fields.
    if (sw_data) free(sw_data);
    sw_size = n;				// 
    sw_data = (sw_t*)calloc( n, sizeof(sw_t) );	// 
    va_list valist;				// copy names
    va_start( valist, n );
    for (int i=0; i<n; i++) {
      char *name = va_arg( valist, char* );  if (!name) break;
      strncpy( sw_data[i].name, name, 9 );
    }
    va_end( valist );
  }
  void startSWatch(void)  { if (!sw_data) return; start(); }
  void checkSWatch(int i) { if (!sw_data||i<0||i>=sw_size) return;   double et = checkTime(); sw_data[i].tt += et; /*if (et>sw_data[i].mx) sw_data[i].mx = et;*/ sw_data[i].cnt++; }
  void stopSWatch (int i=-1) { if (!sw_data||i<0||i>=sw_size) stop(); else { double et = stop(); sw_data[i].tt += et; /*if (et>sw_data[i].mx) sw_data[i].mx = et;*/ sw_data[i].cnt++; } }
  void clearSWatch(void)  { for (int i=0; i<sw_size; i++) sw_data[i].mx = sw_data[i].tt = sw_data[i].cnt = 0; }
  void printSWatch(int prec=3) { 
    int  i, size=prec+4;
    char format[80];
    printf("SWatch ");		// field titles
    sprintf(format, "%%%ds ", size);
    for (i=0; i<sw_size; i++)  printf( format, sw_data[i].name );
    sprintf(format, "%%%d.%df ", size, prec);
//     printf("\n   Max ");		// maximum execution time
//     for (i=0; i<sw_size; i++)  printf( format, sw_data[i].mx );
    printf("\n   Avg ");		// average execution time
    for (i=0; i<sw_size; i++)  printf( format, (sw_data[i].cnt>0 ? sw_data[i].tt/sw_data[i].cnt : 0) );
    printf("\n   Sum ");		// total execution time
    for (i=0; i<sw_size; i++)  printf( format, sw_data[i].tt );
    printf("\n   Cnt "); 		// number of executions
    sprintf(format, "%%%dd ", size);
    for (i=0; i<sw_size; i++)  printf( format, sw_data[i].cnt );
    printf("\n");
    if (accumulated_time==0)
      printf("Warning: The measured time may not be accurate. (last accumulated_time==0)\n");
  }
  
  // -----------------------------------------------------------------
  // Etc.
  // -----------------------------------------------------------------
public:  
  double getFPS(int n = 10) {
    // Measuring FPS (Frame Per Second) in iterative function calls
    static double history[50];
    static int    h_pos = 50;
    double sum = 0;
    if (--h_pos < 0) h_pos = 49;
    history[h_pos] = stop();
    // average over 'n' history records
    for (int i = 0; i < n; i++) sum += history[ (h_pos+i)%50 ];
    start();
    // if (true) {
    //   printf("n=%d sum=%.3f : ", n, sum);
    //   for (int i=0; i < n; i++) printf("%.3f ", history[ (h_pos+i)%50 ]);
    //   printf("\n");
    // }
    return (sum > 0 ? n/sum : 0.0);
  }
  
  double getTimeInSec(void) {
    // Get current time in seconds
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

#endif  // UTIL_TIMER_HPP


// ===================================================================
// ===================================================================

// // #define iX86_CLOCK
// // #define iX86_CLOCK

// #if defined(iX86_CLOCK)

// // intel X86 CPU clock time (64 bits)
// // This part of code is originally from the following site:
// // http://www.cs.utexas.edu/users/sawada/cs352/hw1/timing.html

// typedef unsigned long long int  timer_value_t;

// typedef struct {
//   union {
//     volatile unsigned long long int timerValue;
//     struct {
//       volatile unsigned long int timerValueLow;
//       volatile unsigned long int timerValueHigh;
//     } separated;
//   } u;
// } x86timer_t;

// static timer_value_t getTimerValue(void) {
//   x86timer_t Timer;
//   asm volatile ("rdtsc; movl %%edx, %0; movl %%eax, %1"
// 		: "=r" (Timer.u.separated.timerValueHigh),
// 		"=r" (Timer.u.separated.timerValueLow)
// 		: 
// 		: "%edx", "%eax");
//   return Timer.u.timerValue;
// }
  
// // -------------------------------------------------------------------
// #else

// // Coordinated Universal Time (UTC)   
// // (seconds since the Epoch, 00:00:00 on Jan 1, 1970)

// typedef time_t timer_value_t;
// #define getTimerValue(void)  time(NULL)

// #endif


/* ================================================================ */
