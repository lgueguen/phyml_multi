#include "utilities.h"
#include "sigAnalyse.h"


/* globals */
int stop_next_point = FALSE; 
int final_stop = FALSE; 


/*********************************************************/

int check_sigusr2(void)
{

  //signal(SIGUSR1, do_on_sigusr1);
 signal(SIGUSR2, do_on_sigusr2);
 return stop_next_point;
} 


/*********************************************************/
void do_on_sigusr2(int signum)
{
  if(signum == SIGUSR2) {
    stop_next_point = TRUE;
    printf("Receiving SIGUSR2\n");      
  }
}


/*********************************************************/
int check_sigusr1(void)
{
  signal(SIGUSR1, do_on_sigusr1);
  return final_stop;
} 


/*********************************************************/
void do_on_sigusr1(int signum)
{
  if (signum == SIGUSR1){
    final_stop = TRUE;
    printf("Receiving SIGUSR1\n");      
  }

}


/*********************************************************/
