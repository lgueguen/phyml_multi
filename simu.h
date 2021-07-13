/*

PHYML :  a program that  computes maximum likelihood  phylogenies from
DNA or AA homologous sequences 

Copyright (C) Stephane Guindon. Oct 2003 onward

All parts of  the source except where indicated  are distributed under
the GNU public licence.  See http://www.opensource.org for details.

*/
#include <unistd.h>
#include <fcntl.h> /* for making non-blocking pipes with  F_SETFL */

#ifndef CURR_H
#define CURR_H

void Simu(arbre *tree,int n_step_max, int write_int, int read_int, int maximum);//, int read_int_end);
void Select_Edges_To_Swap(arbre *tree,edge **sorted_b,int *n_neg);
void Fix_All(arbre *tree);
void Update_Bl(arbre *tree,double fact);
void Make_N_Swap(arbre *tree,edge **b,int beg,int end);
int Make_Best_Swap(arbre *tree);
int Mov_Backward_Topo_Bl(arbre *tree,double lk_old,edge **tested_b,int n_tested, double * lk_other_vec, int maximum);
void Unswap_N_Branch(arbre *tree,edge **b,int beg,int end);
void Swap_N_Branch(arbre *tree,edge **b,int beg,int end);

#endif
