#include <stdio.h>

#ifdef FORTRAN_UNDERSCORE
  int flushc_(mode)
#else
  int flushc (mode)
#endif
int *mode;
{
	int ierr;
	ierr = 0;
	if (*mode == 0) ierr = fflush(stderr);
	if (*mode == 1) ierr = fflush(stdin);
	if (*mode == 2) ierr = fflush(stdout);
	return (ierr);
}
