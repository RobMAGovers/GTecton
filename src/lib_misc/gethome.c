#include <pwd.h>
#include <sys/types.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define FORTRAN_CHAR_LEN int
#define EOL '\0'
#define SPACE 32

#ifdef FORTRAN_UNDERSCORE
    void gethome_(char *, FORTRAN_CHAR_LEN);
#else
    void gethome (char *, FORTRAN_CHAR_LEN);
#endif
void strtF (char *, FORTRAN_CHAR_LEN);


#ifdef FORTRAN_UNDERSCORE
    void gethome_(ccharfunc_ptr, ccharfunc_len)
#else
    void gethome (ccharfunc_ptr, ccharfunc_len)
#endif
	char *ccharfunc_ptr;
        FORTRAN_CHAR_LEN ccharfunc_len;
{
	int i;
    struct passwd *pd;//*getpwnam();
	char *homedir;
	char *username;
	const char *uservar = "USER";
	for (i=0; i<ccharfunc_len; i++) ccharfunc_ptr[i] = SPACE;

	username = getenv(uservar);
	if (username == NULL)
	{
		printf("gethome: getenv(USER) failed\n");
	}
	else
	{
		pd = getpwnam(username);
		if ( (pd == NULL) || (homedir = pd->pw_dir) == NULL )
		{
			printf("gethome: getpwnam failed\n");
		}
		else
		{
			strcpy(ccharfunc_ptr, homedir);
			strtF(ccharfunc_ptr,ccharfunc_len);
		}
	}
	return;
}

void strtF(string,len)
char *string;
FORTRAN_CHAR_LEN len;
{
    int i,j;
    /* find lowest adress containing an EOL */
    for (j=0; j<len; j++)
        if ( *(string+j) == EOL || *(string+j) < 6 || *(string+j) > 126 )
            break;
    if (j < len) {
        for (i=j; i<len; i++) *(string+i) = 32;
    }
    /* *(string+len) = EOL; */
    return;
}
