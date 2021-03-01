#include <string.h>
#include <errno.h>
#include <iostream>
#include <cstdlib>
#include <math.h>  // must not be <cmath> because this breaks isnan.

#ifndef __ERROR_H__
#define __ERROR_H__

#ifdef DEBUGGING
#ifndef __gnu_linux__  
#define FATAL_ERROR(MSG) do{					\
  std::cout << __FILE__ << ':' << __LINE__ << ':' << MSG << std::endl;	\
  std::cout << "aborting." << std::endl ; \
  abort(); \
}while(0)
#elif defined(SX)
#define FATAL_ERROR(MSG) do{					\
  std::cout << __FILE__ << ':' << __LINE__ << ':' << MSG << std::endl;	\
  std::cout << "exiting." << std::endl; \
  exit(1); \
}while(0)
#else
#define FATAL_ERROR(MSG) do{					\
  std::cout << __FILE__ << ':' << __LINE__ << ':' << MSG << std::endl;	\
  std::cout << "provoking a ... " ; \
  std::cout << *(((int *) 0) +3) + 4;\
}while(0)
#endif
#else
#ifdef __gnu_linux__
#define FATAL_ERROR(MSG) do{					\
  std::cout << __FILE__ << ':' << __LINE__ << ':' << MSG << std::endl;	\
  std::cout << "provoking a ... " ; \
  std::cout << *(((int *) 0) +3) + 4;\
}while(0)
#else
#define FATAL_ERROR(MSG) do{					\
  std::cout << __FILE__ << ':' << __LINE__ << ':' << MSG << std::endl;	\
  exit(1);                                                      \
}while(0)
#endif
#endif

#define SYS_ERROR() FATAL_ERROR(strerror(errno));

/* #ifdef DEBUGGING */
#define WARNING(MSG) do{					\
  std::cout << __FILE__ << ':' << __LINE__ << ":WARNING:" << MSG << std::endl;	\
}while(0)
/* #else */
/* #define WARNING(MSG) */
/* #endif */

#ifdef DEBUGGING
#define DEBUG(MSG) do{					\
  std::cout << __FILE__ << ':' << __LINE__ << ':' << MSG << std::endl;	\
}while(0)
#else
#define DEBUG(MSG) 
#endif

#if defined(DEBUGGING) && ! (defined(SX) && defined(PARALLEL))
#define ASSERT(X) do{if(!(X)) FATAL_ERROR("Assertation " #X " failed");}while(0)
#else
#define ASSERT(X)
#endif

#define WARN_IF(X,Y)						\
do{if((X)) std::cout << __FILE__ << ':' << __LINE__		\
		     << ": WARNING: " << (#X) << std::endl	\
		     << "************* " << Y << std::endl;	\
 }while(0)
  
#define ALWAYS_ASSERT(X) do{if(!(X)) FATAL_ERROR("Assertation " #X " failed");}while(0)

#define SYSCALL(X) do{int _MYERROR_;if((_MYERROR_=(X))) FATAL_ERROR("Syscall" << std::endl << "   "<< #X << std::endl << "failed with return value " << _MYERROR_ );}while(0)

#define REPORT(X) std::cout << #X << " = " << (X) << std::endl

#define REPORT_ONCE(X) do{static bool REPORTED=false; if(!REPORTED){REPORTED=true; std::cout << #X << " = " << (X) << std::endl;}}while(0);

#define ERROR_TEST(X,Y) do{if((X)==(Y))WARNING(#Y<<" detected");}while(0);

extern int TRACEFLAG;

#define TRACE_RANDOM   0x01
#define TRACE_ODE      0x02
#define TRACE_FLOWS    0x04
#define TRACE_DYNAMICS 0x08
#define TRACE_LOOPS    0x10
#define TRACE_MAIN     0x20
#define TRACE_SPECIES  0x40
#define TRACE_LINKS  0x80

#ifdef DEBUGGING
#define TRACE(X,FLG) do{if(TRACE_##FLG&TRACEFLAG){std::cout << __FILE__ << ':' << __LINE__ << ":"  << "TRACE(" << TRACE_##FLG << "): ";  REPORT(X);}}while(0)
#else
#define TRACE(X,FLG)
#endif

void outOfMemory();

void signal_handling();
extern int exit_now;
extern int save_now;

class terminal_condition {
  const char * message;
public:
  terminal_condition(const char * m):message(m){};
  terminal_condition():message("no message"){};
  operator const char *(){return message;};
};

class AnalysisBug{};

#define my_isnan(X) std::isnan(X)
#define my_isinf(X) std::isinf(X)

bool test_my_isnan(double f);///< evaluates my_isnan at run time
bool test_my_isinf(double f);///< evaluates my_isinf at run time

int cache_mark(char * begin,char * end);/// try to prevent cache losses

#endif // __ERROR_H__
