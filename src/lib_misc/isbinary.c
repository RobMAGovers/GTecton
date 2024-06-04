/* 
 * determine binary or ascii data format
 * return -1 on error, 1 if binary and 0 if ascii
*/

#include <stdio.h>
#include <fcntl.h>
#include <ctype.h>
#include <unistd.h>
#define MAX_RD_BUF 1048576

#ifdef FORTRAN_UNDERSCORE
int isbinary_(s, ls, len_s)
#else
int isbinary (s, ls, len_s)
#endif
	char s[];
	int *ls, len_s;
{
	int c,cnt,fd,only_ascii = 0;
	int i = 0;
	unsigned char buf[MAX_RD_BUF];
	s[*ls] = '\0';
	if ((fd = open(s, O_RDONLY)) < 0 ||
	(cnt = read(fd, &buf[0], MAX_RD_BUF)) <= 0) {
		perror(s);
		return -1;
	}
//	do {only_ascii = isascii(buf[i]); i++; }
	do {/*printf("%d\n",(int)(buf[i]));*/only_ascii = ((int)(buf[i]) <= 127 && (int)(buf[i])>0); i++; }
	while (only_ascii && i < cnt);
	close(fd); s[*ls] = 32;
	return ! only_ascii;
}
