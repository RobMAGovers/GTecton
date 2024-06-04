#include <time.h>

int smallSleep_(long *nano_seconds, int *ierr){
    struct timespec sleeptime ;
    sleeptime.tv_sec = 0 ;
    sleeptime.tv_nsec = *nano_seconds ;
    *ierr = nanosleep(&sleeptime, NULL);
}
