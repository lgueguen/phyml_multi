#include <signal.h>
 
#ifndef SIGANALYSE_H
#define SIGANALYSE_H

int check_sigusr1(void);
void do_on_sigusr1(int signum); 
int check_sigusr2(void);
void do_on_sigusr2(int signum); 


#endif
