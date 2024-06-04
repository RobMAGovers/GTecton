/*
 * Routine for writing real or integer data to a binary file.
 *
 * char *fname		file name
 * char *mode		"w" for (over)write, or "a" for append.
 * float or int buf[]	data array.
 * int *n		number of data.
 * (size_t fname_len	needs only be supplied if routine is called from another
 * 			C-routine, not if called from FORTRAN)
 * RETURN VALUES	0 OK, -1 error occurred
 */

#include <fcntl.h>	  /* for creat() */
#include <unistd.h>	  /* for close() and lseek() */
#include <errno.h>	  /* for perror() */
#include <stdio.h>	  /* for BUFSIZ and fprintf() */
#include <string.h>   /* for strncpy() */
#define PERMS 0644    /* file access permissions: -rw-r--r-- */
#define NBYTESINT   4 /* number of bytes per integer */
#define NBYTESFLOAT 8 /* number of bytes per Fortran double precision */


#ifdef FORTRAN_UNDERSCORE
int wrdxint_(char *fname,char *mode,int *buf,int *n,unsigned int fname_len)
#else
int wrdxint (char *fname,char *mode,int *buf,int *n,unsigned int fname_len)
#endif
{
    char name[256];
    int fd, i, m, l, f, numbytes, length;
    size_t nbytes;

    /* copy fname(fnblk(fname):lnblnk(fname)) to name */
    length = (int)(fname_len);
    if (length <= 0) {
        fprintf(stderr,"wrdx: no file argument error\n"); return -1;
    }

    i=0; while (fname[i] == 32 && i<length) i++; f = i;
    i=length-1; while (fname[i] == 32 && i>=0) i--; l = i + 1;
    m=l-f; 
    if (m <= 0) {
        fprintf(stderr,"wrdx: no file argument error\n"); return -1;
    }

    strncpy(name,fname+f,(size_t)m); name[m] = '\0';

    /* open the file */
    if (*mode == 'w') {
        if ( (fd = creat(name,PERMS)) == -1) {
            perror("wrdx"); return -1;
        }
    } else if (*mode == 'a') {
        if ((fd = open(name, O_WRONLY, 0)) == -1) {
             if ( (fd = creat(name,PERMS)) == -1) {
                 perror("wrdx"); return -1;
             }
        }
        if (lseek(fd, 0L, 2) == -1) {
            perror("wrdx"); return -1;
        }
    } else {
        fprintf(stderr,"wrdx: mode %c not recognized\n",*mode);
        return -1;
    }

    /* Write buf[n] in portions of BUFSIZ bytes to file */

    numbytes = NBYTESINT*(*n); nbytes = write(fd,buf,numbytes);
    if (nbytes != numbytes) {
        fprintf(stderr,"wrdx: write incomplete\n"); return -1;
    }

    /* close file */
    if (close(fd) == -1) {
        perror("wrdx"); return -1;
    }
    return 0;
}

#ifdef FORTRAN_UNDERSCORE
int wrdxfloat_(char *fname,char *mode,double *buf,int *n,unsigned int fname_len)
#else
int wrdxfloat (char *fname,char *mode,double *buf,int *n,unsigned int fname_len)
#endif
{
    char name[256];
    int fd, i, m, l, f, numbytes, length;
    size_t nbytes;

    /* copy fname(fnblk(fname):lnblnk(fname)) to name */
    length = (int)(fname_len);
    if (length <= 0) {
        fprintf(stderr,"wrdx: no file argument error\n"); return -1;
    }

    i=0; while (fname[i] == 32 && i<length) i++; f = i;
    i=length-1; while (fname[i] == 32 && i>=0) i--; l = i + 1;
    m=l-f; 
    if (m <= 0) {
        fprintf(stderr,"wrdx: no file argument error\n"); return -1;
    }

    strncpy(name,fname+f,(size_t)m); name[m] = '\0';

    /* open the file */
    if (*mode == 'w') {
        if ( (fd = creat(name,PERMS)) == -1) {
            perror("wrdx"); return -1;
        }
    } else if (*mode == 'a') {
        if ((fd = open(name, O_WRONLY, 0)) == -1) {
             if ( (fd = creat(name,PERMS)) == -1) {
                 perror("wrdx"); return -1;
             }
        }
        if (lseek(fd, 0L, 2) == -1) {
            perror("wrdx"); return -1;
        }
    } else {
        fprintf(stderr,"wrdx: mode %c not recognized\n",*mode);
        return -1;
    }

    /* Write buf[n] in portions of BUFSIZ bytes to file */

    numbytes = NBYTESFLOAT*(*n); nbytes = write(fd,buf,numbytes);
    if (nbytes != numbytes) {
        fprintf(stderr,"wrdx: write incomplete\n"); return -1;
    }

    /* close file */
    if (close(fd) == -1) {
        perror("wrdx"); return -1;
    }
    return 0;
}

