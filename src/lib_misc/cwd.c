#include <stdio.h>
#include <string.h>

#define FORTRAN_CHAR_LEN size_t

#ifdef FORTRAN_UNDERSCORE
  void cwd_(char *Path, FORTRAN_CHAR_LEN PathLen)
#else
  void cwd(char *Path, FORTRAN_CHAR_LEN PathLen)
#endif
{
  int l;
  char *getcwd(char *, size_t);
  size_t Path_Len;
  Path_Len = (size_t) PathLen;
  if ( getcwd(Path,Path_Len) == NULL ) fprintf(stderr,"getcwd failed\n");
  for (l = strlen(Path); l < PathLen; l++) Path[l] = 32;
  return;
}
