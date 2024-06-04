/* FORTRAN callable STAT, LSTAT and FSTAT routines */
#include <sys/stat.h>
#include <stdio.h>
#include <string.h>

#define FORTRAN_CHAR_LEN size_t

#ifdef FORTRAN_UNDERSCORE
  int cstat_(fname, statb, fname_len)
#else
  int cstat (fname, statb, fname_len)
#endif

char *fname;
long int statb[13];
FORTRAN_CHAR_LEN fname_len;
{
  struct stat Buffer;
  char file[256];
  int ierr,i,f,l,n;

  /* copy string up 'til lnblnk */
  i=0; while (fname[i] == 32 && i<fname_len) i++; f = i;
  i=fname_len-1; while (fname[i] == 32 && i>=0) i--; l = i + 1;
  n=l-f; if (n <= 0) return (1);
  strncpy(file,fname+f,n); file[n] = '\0';

  /* Get stat-structure */
  ierr = stat(file,&Buffer);
  if (ierr != 0) perror(" stat");

  /* Copy buffer into integer array */
  statb[ 0] = Buffer.st_dev;
  statb[ 1] = Buffer.st_ino;
  statb[ 2] = Buffer.st_mode;
  statb[ 3] = Buffer.st_nlink;
  statb[ 4] = Buffer.st_uid;
  statb[ 5] = Buffer.st_gid;
  statb[ 6] = Buffer.st_rdev;
  statb[ 7] = Buffer.st_size;
  statb[ 8] = Buffer.st_atime;
  statb[ 9] = Buffer.st_mtime;
  statb[10] = Buffer.st_ctime;
  statb[11] = Buffer.st_blksize;
  statb[12] = Buffer.st_blocks;
  return (ierr);
}
