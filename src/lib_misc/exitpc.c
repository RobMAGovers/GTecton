#include <stdlib.h>
  
#ifdef FORTRAN_UNDERSCORE
  void exitpc_(Status)
#else
  void exitpc (Status)
#endif
int *Status;
{
  exit (*Status);
}
