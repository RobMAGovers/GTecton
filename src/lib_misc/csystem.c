#include <stdio.h>
#include <stdlib.h>
#include <sys/wait.h>

#ifdef FORTRAN_UNDERSCORE
  int csystem_(cmd, idum)
#else
  int csystem (cmd, idum)
#endif

char *cmd;
int  idum;
{
	int status;
	*(cmd+idum) = '\0';
	status = system(cmd);

    /* status = (WIFEXITED(status) == 0); */
       status = WEXITSTATUS(status);

	/* if (status != 0) perror("csystem"); */
	*(cmd+idum) = ' ';
	return status;
}
