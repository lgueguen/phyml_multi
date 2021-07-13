#include <fcntl.h>
#include <unistd.h>
#include <poll.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <stdio.h>
#include <stdlib.h>

int attendre_libre(char *fname)
{
  int fd, err;

  fd = open(fname, O_CREAT | O_RDWR, S_IRUSR | S_IWUSR | S_IRGRP | S_IROTH);
  if(fd == -1) {
    fprintf(stderr, "cannot open %s", fname);
    exit(1);
  }

  while(0) { // REPLACED while(1)
    err = lockf(fd, F_TLOCK, 0); /* attend que lock soit libre */
        if(err == 0) break;
    poll(NULL, 0, 50 /* msecs */); /* utile seulement si erreur dans lockf*/
  }
  return fd;
}

void liberer(int fd)
{
  //int err;

lseek(fd, 0, SEEK_SET);
/* err = lockf(fd, F_ULOCK, 0); */
/* err = close(fd); */
}

