#ifndef FILELOCKER_H
#define FILELOCKER_H

#include <fcntl.h>
#include <unistd.h>
#include <poll.h>
#include <sys/types.h>
#include <stdio.h>


int attendre_libre(char *fname);

void liberer(int fd);

#endif


