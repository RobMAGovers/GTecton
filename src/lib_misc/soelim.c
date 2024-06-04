#include <stdio.h>
#include <string.h>

int resolve_file(source, dest)
FILE	*source, *dest;
{
	FILE	*insert;
	char 	line[BUFSIZ], file[BUFSIZ];
	int	ierr = 0;

	while (fgets(line, BUFSIZ, source)) {
		if (strncmp(line, ".so", 3))
			fputs(line, dest);
		else {
			sscanf(line, "%*s %s", file);
			if ((insert = fopen(file, "r")) == NULL) {
				perror(file);
				ierr = 1;
			} else {
				ierr = resolve_file(insert, dest);
				fclose(insert);
			}
		}
	}
	return ierr;
}

void copy_it (fo,fi,len)
char *fi, *fo;
long len;
{
	int i,f,l,n;

	i=0; while (fi[i] == 32 && i<len) i++; f = i;
	i=len-1; while (fi[i] == 32 && i>=0) i--; l = i + 1;
	n=l-f; if (n <= 0) return;

	strncpy(fo,fi+f,n); fo[n] = '\0';

	return;
}

#ifdef FORTRAN_UNDERSCORE
  int soelim_(fin, fout, fin_len, fout_len)
#else
  int soelim (fin, fout, fin_len, fout_len)
#endif
char	*fin, *fout;
long	fin_len, fout_len;
{
	FILE 	*fp_in, *fp_out;
	char    FIN[BUFSIZ], FOUT[BUFSIZ];
	int	ierr;

	if (fin_len <= 0) {
		fprintf(stderr,"soelim: input file length error\n");
		return 1;
	}
	if (fout_len <= 0) {
		fprintf(stderr,"soelim: output file length error\n");
		return 1;
	}
	copy_it (FIN, fin, fin_len); copy_it (FOUT, fout, fout_len);
	if ((fp_in = fopen(FIN,"r")) == NULL) {
		perror(FIN);
		return 1;
	}
	if ((fp_out = fopen(FOUT,"w")) == NULL) {
		perror(FOUT);
		return 1;
	}
	ierr = resolve_file(fp_in, fp_out);
	if (ierr == 0) {
	    fclose(fp_in);
	    fclose(fp_out);
	}
	return (ierr);
}
