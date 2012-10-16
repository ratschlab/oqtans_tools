#include <stdio.h>
#include <stdarg.h>
#include <errno.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <ctype.h>
#include <sys/stat.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <sys/mman.h>
#include <time.h>
#include <mex.h>
#include "mmap.h"

void die(const char* fmt,...)
{
	char string[10 * MAXLINE];
	va_list argp;
	va_start(argp, fmt);
	sprintf(string, fmt, argp);
	mexErrMsgTxt(string);
	va_end(argp);
	exit(EXIT_FAILURE);
}

int mmap_file(const char *filename, int open_mode,
		void **map, off_t *size)
{
	int fd, ret, mmap_prot, mmap_flags;
	struct stat file_status;

	if (open_mode == O_RDONLY) {
		mmap_prot = PROT_READ;
		mmap_flags = MAP_PRIVATE;
	} else {
		mmap_prot = PROT_READ | PROT_WRITE;
		mmap_flags = MAP_SHARED;
	}
	ret = open(filename, open_mode);
	if (ret < 0){
		printf("ret: %i  filename: %s\n", ret, filename);
		printf("strerror: %s\n", strerror(errno));
		die("can not open %s: %s\n", filename, strerror(errno));
	}
	fd = ret;
	if (fstat(fd, &file_status) < 0)
		die("fstat error: %s\n", strerror(errno));
	*size = file_status.st_size;
	if (!*size){
		/* empty file: do NOT die here, so that interval_query can return an empty interval*/
		return EMPTY_FILE;
	}
	*map = mmap(NULL, *size, mmap_prot, mmap_flags, fd, 0);
	if (*map == MAP_FAILED)
		die("mmap error: %s\n", strerror(errno));
	ret = 1;
	close(fd);
	return ret;
}

