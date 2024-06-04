#include <unistd.h>
#include <stdio.h>
#include <errno.h>
#include <string.h>
  
#define FORTRAN_CHAR_LEN size_t

#ifdef FORTRAN_UNDERSCORE
  void rm_(Path,PathLen)
#else
  void rm(Path,PathLen)
#endif

char *Path;
FORTRAN_CHAR_LEN PathLen;
{
  char lpath[96];
  int i,f,l,n;
  i=0; while (Path[i] == 32 && i<PathLen) i++; f = i;
  i=PathLen-1; while (Path[i] == 32 && i>=0) i--; l = i + 1;
  n=l-f; if (n <= 0) return;
  
  strncpy(lpath,Path+f,n); lpath[n] = '\0';
  if (unlink(lpath)) perror(" rm");
  return;
}
