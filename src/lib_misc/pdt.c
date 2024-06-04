#include <utime.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <dirent.h>
#include <stdio.h>
#include <string.h>
#include <time.h>

static char            *usage = "[-rvw] [-d #] file [file ...]";
static struct stat      info;
static time_t           NOW;
static time_t           MAX = 3;
static int              VERBOSE = 0;
static int              RECURSE = 0;
static int              W[7] = { 0, 1, 1, 1, 1, 1, 0 };
extern char            *optarg;
extern int              optind;
#define RANGE 		0x7FFF
#define NEXT 		continue
#define ISDIR   	((info.st_mode & S_IFDIR) == S_IFDIR)
#define MIN_H   	9
#define MAX_H   	17

void 
traverse(path)
	char            path[];
{
	DIR            *loc;
	struct dirent  *entry;
	int             i;
	char            buf[BUFSIZ];

	if ((loc = opendir(path)) == NULL) {
		perror(path);
		return;
	}
	while ((entry = readdir(loc)) != NULL) {
		if (strcmp(entry->d_name, "..") == 0
		||  strcmp(entry->d_name,  ".") == 0 ) NEXT;
		sprintf(buf, "%s/%s", path, entry->d_name);
		if (stat(buf, &info)) {
			perror(buf);
			NEXT;
		}
		ISDIR ? traverse(buf) : (void) update(buf);
	}
	closedir(loc);
}

int 
update(path)
	char            path[];
{
	time_t          DELTA, new[2];
	struct tm      *date;
	do {
		DELTA = (time_t) ((((double) rand()) / RANGE) * MAX);
		new[0] = new[1] = NOW - DELTA;
		date = localtime(new);
	} while (date->tm_hour < MIN_H
	    ||   date->tm_hour > MAX_H
	    ||   W[date->tm_wday] == 0);
	if (utime(path, new)) {
		perror(path);
		return (1);
	}
	if (VERBOSE)
		printf("%-50s%24s", path, asctime(date));
	return (0);
}

main(argc, argv)
	int             argc;
	char           *argv[];
{
	int             flag;

	while ((flag = getopt(argc, argv, "d:nrvw")) != EOF) {
		switch (flag) {
		case 'd':
			MAX = atol(optarg);
			break;
		case 'r':
			RECURSE++;
			break;
		case 'v':
			VERBOSE++;
			break;
		case 'w':
			W[0] = W[6] = 1;
			break;
		case '?':
			printf("usage: %s %s\n", argv[0], usage);
			exit(1);
		}
	}
	if ((MAX *= 24 * 3600) == 0) {
		printf("usage: %s %s\n", argv[0], usage);
		exit(1);
	}
	srand((unsigned) (NOW = time((time_t *) 0)));
	for (; optind < argc; optind++) {
		if (stat(argv[optind], &info)) {
			perror(argv[optind]);
			NEXT;
		}
		if (ISDIR)
			if (RECURSE)
				traverse(argv[optind]);
			else
				NEXT;
		else
			(void) update(argv[optind]);
	}
	exit(0);
}
