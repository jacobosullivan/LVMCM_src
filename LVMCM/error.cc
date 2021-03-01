//$Id: error.cc 2166 2011-05-24 14:49:20Z axel $
// no code in this file
// it is just required for the .h -> .cc systematic

#include "error.h"
#include <signal.h>
#include <stdlib.h>

int TRACEFLAG=0; //report nothing

// Manage adjustable parameters:
#include "cfgList.h"
static cfgStruct cfg[] = 
{
  CFGINT(TRACEFLAG),
  {0, CFG_END, 0}
};
static cfg_add dummy(cfg);

void outOfMemory()
{
  FATAL_ERROR("Probably out of memory");
}

int exit_now=0;
int save_now=0;

void exiter(int i){
  // To make signaling of running jobs work in torque, you need to
  // make sure that the shells running the job to not catch the
  // signal.  For this, put this lines into the file ~/.bash_profile
  // AND into the job execution script:
  //
  // trap "" SIGUSR1 SIGUSR2 SIGFPE
  //
  switch(i){
  case SIGUSR1:
    if(get_cfg_parameter("record_relaxation_dynamics")){
      set_cfg_parameter("record_relaxation_dynamics","0");
    }else{
      set_cfg_parameter("record_relaxation_dynamics","1");
    }
    break;
  case SIGUSR2:
    // This is a mechanism to get diagnostic outputs from running jobs
    // via qsig -s USR2.  To make this work, please activate email
    // forwarding by putting your email address into the file
    // ~/.forward .  
    system("(sleep 1;(echo $PBS_JOBID on $HOSTNAME;echo;top -bn 2 -d 5|grep -B 6 -A15 COMMAND)|mail -s top $USER)&");
    break;
  case SIGFPE:
    abort();
    break;
  default:
    exit_now=1;
    break;
  }
  signal_handling();
}

void signal_handling(){
  signal(SIGUSR1,&exiter);
  signal(SIGUSR2,&exiter);
  signal(SIGXCPU,&exiter);
  signal(SIGHUP,&exiter);
  signal(SIGTERM,&exiter);
  signal(SIGFPE,&exiter);
}

#ifdef __DARWIN_10_6_AND_LATER
#include <cmath>
#endif

bool test_my_isnan(double f){
  return my_isnan(f);
}
bool test_my_isinf(double f){
  return my_isinf(f);
}

int cache_mark(char * begin,char * end){
  int sum=0;
  const int cache_line_size=64; //should be adjustable

  for(char * p=end;p>=begin;p-=cache_line_size){
    sum+=*p;
  }

  return sum;
}
