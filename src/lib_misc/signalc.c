#include <signal.h>
  
#ifdef FORTRAN_UNDERSCORE
  void signalc_(SigNum, ProcPtr, Flag)
#else
  void signalc (SigNum, ProcPtr, Flag)
#endif
int *SigNum, *Flag;
void (*ProcPtr) ();
{
  if ( *Flag == 0 ) 
    signal(*SigNum, SIG_DFL);
  else if ( *Flag > 0 ) 
    signal(*SigNum, SIG_IGN);
  else
    signal(*SigNum, ProcPtr);
  return;
}
