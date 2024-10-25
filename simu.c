/*

  PHYML :  a program that  computes maximum likelihood  phylogenies from
  DNA or AA homologous sequences 

  Copyright (C) Stephane Guindon. Oct 2003 onward

  All parts of  the source except where indicated  are distributed under
  the GNU public licence.  See http://www.opensource.org for details.

*/

#include "utilities.h"
#include "ml.h"
#include "optimiz.h"
#include "models.h"
#include "free.h"
#include "simu.h"
#include "sigAnalyse.h"

double PROP_STEP;


/*********************************************************/
extern int stop_next_point; 


void Simu(arbre *tree, int n_step_max, int write_int, int read_int, int maximum)//, int read_int_end)
{
  double old_loglk,n_iter,lambda;//,diff_lk;
  int i,n_neg,n_tested,n_without_swap,step;//,it_lim_without_swap;
  edge **sorted_b,**tested_b;
  int each,each_invar;
  //    int opt_free_param;
  
  //listening to the user signal 2 in case the server tells to stop the computations
  check_sigusr2();
  
  printf ("Client %d is in charge of tree %d\n", getpid(), tree->input->current_tree);

  sorted_b = (edge **)mCalloc(tree->n_otu-3,sizeof(edge *));
  tested_b = (edge **)mCalloc(tree->n_otu-3,sizeof(edge *));
      
  old_loglk = tree->tot_loglk = UNLIKELY;
  n_iter              = 1.0;
  //    it_lim_without_swap = (tree->mod->invar)?(8):(5);
  n_tested            = 0;
  n_without_swap      = 0;
  step                = 0;
  each                = 4;
  lambda              = 0.75;
  each_invar          = 2;
  old_loglk           =  tree->tot_loglk;
  //    opt_free_param      = 0;

  double *lk_vec; //The likelihood vector of the tree
  double *lk_other_vec; //The likelihood vector which is the sum of the other tree likelihoods
  
  //Allocating the space for the matrix containing the tree likelihoods   
  tree->mod->all_lk_vec = (double **)mCalloc(tree->input->n_trees,sizeof(double *));
  for (i=0;i<tree->input->n_trees;i++) {
    tree->mod->all_lk_vec[i] = (double *)mCalloc(tree->data->crunch_len,sizeof(double));
  }

  lk_vec=(double *)mCalloc(tree->data->crunch_len,sizeof(double));
  lk_other_vec=(double *)mCalloc(tree->data->crunch_len,sizeof(double));
    
  //Zeroing the auxiliary likelihood vector
  For(i,tree->data->crunch_len) {
    lk_other_vec[i]=0.0;
  }
    
  Lk(tree, tree->data, lk_other_vec, maximum);
  printf ("Client %d : Initial tree->tot_loglk : %f\n", getpid(), tree->tot_loglk);
    
  For (i, tree->data->crunch_len) {
    lk_vec[i]=tree->true_site_lk[i];
  }

  //Sending the initial likelihood vector
  int wok = write(write_int,lk_vec, tree->data->crunch_len * sizeof(double));
  if (wok == -1)
    {}

  // printf("Client %d has sent its vector 0\n", getpid());

  //We wait for the "go" from the server : 
  if (read(read_int, &i, sizeof(int))<sizeof(int)) {
    printf("PROBLEM READ in Simu\n");
  }//It is better if the read is blocking here !

  //USING THE FILE TO COMMUNICATE
  if(tree->input->HMM) {
    Read_From_Memory_With_Lambda (tree->input->fname, tree->mod->all_lk_vec, &tree->mod->hmm_lambda, tree->input->n_trees, tree->data->crunch_len);
    //printf("For client %d, lambda is %f\n", getpid(), tree->mod->hmm_lambda);
  }
  else {
    //printf ("Reading without lambda\n");
    Read_From_Memory (tree->input->fname, tree->mod->all_lk_vec, tree->input->n_trees, tree->data->crunch_len);
  }

  Update_Lk_Vector_In_Lk_Table (tree);
    
  lk_other_vec = Compute_Lk_Summary (tree, tree->mod->all_lk_vec, tree->input->n_trees, tree->input->maxi, tree->input->current_tree); 

  Lk(tree,tree->data, lk_other_vec, maximum);

  old_loglk           =  tree->tot_loglk;
  //printf("Client has read its vector :  new loglk : %f\n", tree->tot_loglk);
  //fcntl(read_int, F_SETFL, fcntl(read_int, F_GETFL) | O_NONBLOCK);//Read non blocking
  //  fcntl(write_int, F_SETFL, fcntl(write_int, F_GETFL) | O_NONBLOCK);//write non blocking
  do
    {
      ++step;
      each--;
      each_invar--;
	  
      tree->mod->s_opt->opt_bl = 0;
      tree->both_sides    = 1;
	
      if(tree->mod->s_opt->print)
        {
          /* if(old_loglk < UNLIKELY+1)
             printf("\n. Process %d : Log(lk) :               * -> %15.6f ",getpid(), tree->tot_loglk);*/
          if (step==1) {}
          else {
            /* printf("\n. Process %d : Log(lk) : %15.6f -> %15.6f ",getpid(),old_loglk,tree->tot_loglk);*/
            if(old_loglk > UNLIKELY+1)
              {
                if(n_tested > 1) printf("Process %d : Log(lk) : %15.6f -> %15.6f,  %3d swaps done\n",getpid(),old_loglk,tree->tot_loglk, n_tested);
                else             printf("Process %d : Log(lk) : %15.6f -> %15.6f,  %3d swap  done\n",getpid(), old_loglk,tree->tot_loglk, n_tested);
              }
          }
        }
            
      //fflush(NULL);
            
      //  if((fabs(old_loglk-tree->tot_loglk) < 1.E-03) || (n_without_swap > it_lim_without_swap)) break;
      //Not acceptable in a multiprocess framework
      //  if(n_without_swap > it_lim_without_swap) break;
         
            
      if(tree->tot_loglk < old_loglk) //If the new lk is worse than the old one
        {
          if(tree->mod->s_opt->print)
            printf("\n\n. Moving backward (topology + branch lengths) \n");
          fflush(NULL);
          if(!Mov_Backward_Topo_Bl(tree,old_loglk,tested_b,n_tested, lk_other_vec, maximum))
            Exit("\n. Err: mov_back failed\n");
          if(!tree->n_swap) n_neg = 0;
                    
          For(i,2*tree->n_otu-3) tree->t_edges[i]->l_old = tree->t_edges[i]->l;
          Optimiz_All_Free_Param(tree,tree->mod->s_opt->print, lk_other_vec, maximum);
          Lk(tree,tree->data, lk_other_vec, maximum);
        }
      else 
        { //The new lk is better than the old one
			
          if(!each)
            {
              //                            opt_free_param = 1;
              each = 4;
              if(tree->mod->s_opt->print) printf("\n");
              Optimiz_All_Free_Param(tree,tree->mod->s_opt->print, lk_other_vec, maximum);
              tree->mod->s_opt->opt_bl = 0;
              tree->both_sides    = 1;
              Lk(tree,tree->data, lk_other_vec, maximum);
            }
				
          //We could send our new likelihood vector here
          //As well as receive the new one
          For (i, tree->data->crunch_len) {
            lk_vec[i]=tree->true_site_lk[i];
          }
          if (!stop_next_point)
            wok = write(write_int, lk_vec, tree->data->crunch_len*sizeof(double));
          //printf("Client %d has sent its vector 1\n", getpid());
	
          //USING THE FILE TO COMMUNICATE
          if (!stop_next_point)
            {
              if(tree->input->HMM) {
                Read_From_Memory_With_Lambda (tree->input->fname, tree->mod->all_lk_vec, &tree->mod->hmm_lambda, tree->input->n_trees, tree->data->crunch_len); 
                //printf("For client %d, lambda is %f\n", getpid(), tree->mod->hmm_lambda);
              }
              else {
                Read_From_Memory (tree->input->fname, tree->mod->all_lk_vec, tree->input->n_trees, tree->data->crunch_len);
              }
		 
              Update_Lk_Vector_In_Lk_Table (tree);
		
              lk_other_vec = Compute_Lk_Summary (tree, tree->mod->all_lk_vec, tree->input->n_trees, tree->input->maxi, tree->input->current_tree); 

              Lk(tree,tree->data, lk_other_vec, maximum);
              //printf("Client has read its vector 2 : new loglk : %f\n", tree->tot_loglk);
	
              old_loglk = tree->tot_loglk;
 
              Fix_All(tree);
                
              n_neg = 0;
              For(i,2*tree->n_otu-3)
                if((!tree->t_edges[i]->left->tax) && 
                   (!tree->t_edges[i]->rght->tax)) 
                  NNI(tree,tree->t_edges[i],0, lk_other_vec, maximum);
		
              Select_Edges_To_Swap(tree,sorted_b,&n_neg);
                
              Sort_Edges_Diff_Lk(tree,sorted_b,n_neg);
                
              Optimiz_Ext_Br(tree, lk_other_vec, maximum);	  
                
              Update_Bl(tree,lambda);
		
              n_tested = 0;
              For(i,(int)ceil((double)n_neg*(lambda)))
                tested_b[n_tested++] = sorted_b[i];
		
              Make_N_Swap(tree,tested_b,0,n_tested);
                
              if(n_tested > 0) n_without_swap = 0;
              else             n_without_swap++;
                
              fflush(NULL);
              Lk(tree,tree->data, lk_other_vec, maximum);
            }
        }
      n_iter+=1.0;
    }
  while(!stop_next_point);
    
  Free(sorted_b);
  Free(tested_b);
 
  /*     For (i, tree->data->crunch_len) { */
  /*       lk_vec[i]=tree->true_site_lk[i]; */
  /*     } */
  //write(write_int, lk_vec, tree->data->crunch_len*sizeof(double));
  //printf("Client %d has sent its vector 3\n", getpid());

  free(lk_vec);
  free(lk_other_vec);
  
  return;

}

/*********************************************************/

void Select_Edges_To_Swap(arbre *tree, edge **sorted_b, int *n_neg)
{
  int i;
  edge *b;
  //  int min;


  *n_neg = 0;
  tree->min_diff_lk = .0;
  //  min = 0;
  For(i,2*tree->n_otu-3)
    {
      b = tree->t_edges[i];

      if((!b->left->tax) 
         && (!b->rght->tax) 
         && (b->diff_lk < 0.0-MDBL_MIN)) 
	{

	  if((b->left->b[b->l_v1]->diff_lk < b->diff_lk) ||
	     (b->left->b[b->l_v2]->diff_lk < b->diff_lk) ||
	     (b->rght->b[b->r_v1]->diff_lk < b->diff_lk) ||
	     (b->rght->b[b->r_v2]->diff_lk < b->diff_lk)) continue;

	  if(b->diff_lk < tree->min_diff_lk) 
	    {
	      tree->min_diff_lk = b->diff_lk;
              //	      min = i;
	    }

	  sorted_b[*n_neg] = b;
	  (*n_neg)++;
	}
    }
}

/*********************************************************/

void Fix_All(arbre *tree)
{
  int i;
  edge *b;

  tree->mod->pinvar_old = tree->mod->pinvar;
  tree->mod->alpha_old = tree->mod->alpha;
  tree->mod->kappa_old = tree->mod->kappa;
  tree->mod->lambda_old = tree->mod->lambda;
  
  For(i,2*tree->n_otu-3)
    {
      b = tree->t_edges[i];
      b->l_old = b->l;
    }
}

/*********************************************************/

void Update_Bl(arbre *tree, double fact)
{
  int i;
  edge *b;

  For(i,2*tree->n_otu-3)
    {
      b = tree->t_edges[i];
      b->l = b->l_old + (b->ql[0]-b->l_old)*fact;
    }
}

/*********************************************************/

void Make_N_Swap(arbre *tree,edge **b, int beg, int end)
{
  int i;

  tree->n_swap = 0;
  for(i=beg;i<end;i++)
    {
      /*       printf("make swap on %3d d->%10f\n",b[i]->num,b[i]->diff_lk); */
      /*       if(drand48()>0.75) */
      /* 	{ */
      (b[i]->best_conf == 2)?
        (Swap(b[i]->left->v[b[i]->l_v2],b[i]->left,b[i]->rght,b[i]->rght->v[b[i]->r_v1],tree)):
        (Swap(b[i]->left->v[b[i]->l_v2],b[i]->left,b[i]->rght,b[i]->rght->v[b[i]->r_v2],tree));
	
      b[i]->l = b[i]->ql[b[i]->best_conf-1];
      tree->n_swap++;
      /* 	} */
    }
}

/*********************************************************/

int Make_Best_Swap(arbre *tree)
{
  int i,j,return_value;
  edge *b,**sorted_b;
  

  sorted_b = (edge **)mCalloc(tree->n_otu-3,sizeof(edge *));
  
  j=0;
  For(i,2*tree->n_otu-3) if((!tree->t_edges[i]->left->tax) &&
			    (!tree->t_edges[i]->rght->tax))
    sorted_b[j++] = tree->t_edges[i];

  Sort_Edges_Diff_Lk(tree,sorted_b,tree->n_otu-3);

  if(sorted_b[0]->diff_lk < -0.0)
    {
      b = sorted_b[0];
      return_value = 1;
      (b->best_conf == 2)?
	(Swap(b->left->v[b->l_v2],b->left,b->rght,b->rght->v[b->r_v1],tree)):
	(Swap(b->left->v[b->l_v2],b->left,b->rght,b->rght->v[b->r_v2],tree));
      
      b->l = b->ql[b->best_conf-1];
    }
  else return_value = 0;

  Free(sorted_b);

  return return_value;
}

/*********************************************************/

int Mov_Backward_Topo_Bl(arbre *tree, double lk_old, edge **tested_b, int n_tested, double * lk_other_vec, int maximum)
{
  double *l_init;
  int i,step,n_swp,beg,end;
  edge *b,**swp;


  l_init = (double *)mCalloc(2*tree->n_otu-3,sizeof(double));
  swp = (edge **)mCalloc(tree->n_otu-3,sizeof(edge *));

  For(i,2*tree->n_otu-3) l_init[i] = tree->t_edges[i]->l;
  
  step = 2;
  tree->both_sides = 0;
  do
    {
      n_swp = 0;
      For(i,2*tree->n_otu-3) 
	{
	  b = tree->t_edges[i];
	  b->l = b->l_old + (1./step) * (l_init[i] - b->l_old);
	}

      beg = (int)floor((double)n_tested/(step-1));
      end = 0;
      Unswap_N_Branch(tree,tested_b,beg,end);
      beg = 0;
      end = (int)floor((double)n_tested/step);
      Swap_N_Branch(tree,tested_b,beg,end);
      
      if(end == n_swp) tree->n_swap = 0;
      
      tree->mod->s_opt->opt_bl = 0;
      tree->both_sides    = 0;
      Lk(tree,tree->data, lk_other_vec, maximum);
      
      step++;

    }while((tree->tot_loglk < lk_old) && (step < 100));


  if(step == 100)
    {
      For(i,2*tree->n_otu-3) 
	{
	  b = tree->t_edges[i];
	  b->l = b->l_old;
	}

      tree->mod->s_opt->opt_bl = 0;
      tree->both_sides    = 0;
      Lk(tree,tree->data, lk_other_vec, maximum );
    }

  Free(l_init);
  Free(swp);

  tree->n_swap = 0;
  For(i,2*tree->n_otu-3) 
    {
      if(tree->t_edges[i]->diff_lk < 0.0) tree->n_swap++;
      tree->t_edges[i]->diff_lk = +1.0;
    }

  if(tree->tot_loglk > lk_old)                 return  1;
  else if((tree->tot_loglk > lk_old-MIN_DIFF_LK) && 
	  (tree->tot_loglk < lk_old+MIN_DIFF_LK)) return -1;
  else                                         return  0;
}

/*********************************************************/

void Unswap_N_Branch(arbre *tree, edge **b, int beg, int end)
{
  int i;
 
  if(end>beg)
    {
      for(i=beg;i<end;i++)
	{
	  (b[i]->best_conf == 2)?
	    (Swap(b[i]->left->v[b[i]->l_v2],b[i]->left,b[i]->rght,b[i]->rght->v[b[i]->r_v1],tree)):
	    (Swap(b[i]->left->v[b[i]->l_v2],b[i]->left,b[i]->rght,b[i]->rght->v[b[i]->r_v2],tree));
	  b[i]->l = b[i]->l_old;
	}
    }
  else
    {
      for(i=beg-1;i>=end;i--)
	{
	  (b[i]->best_conf == 2)?
	    (Swap(b[i]->left->v[b[i]->l_v2],b[i]->left,b[i]->rght,b[i]->rght->v[b[i]->r_v1],tree)):
	    (Swap(b[i]->left->v[b[i]->l_v2],b[i]->left,b[i]->rght,b[i]->rght->v[b[i]->r_v2],tree));
	  b[i]->l = b[i]->l_old;
	}
    }
}

/*********************************************************/

void Swap_N_Branch(arbre *tree,edge **b, int beg, int end)
{
  int i;
  
  if(end>beg)
    {
      for(i=beg;i<end;i++)
	{
	  (b[i]->best_conf == 2)?
	    (Swap(b[i]->left->v[b[i]->l_v2],b[i]->left,b[i]->rght,b[i]->rght->v[b[i]->r_v1],tree)):
	    (Swap(b[i]->left->v[b[i]->l_v2],b[i]->left,b[i]->rght,b[i]->rght->v[b[i]->r_v2],tree));
	  b[i]->l = b[i]->ql[b[i]->best_conf-1];
	}
    }
  else
    {
      for(i=beg-1;i>=end;i--)
	{
	  (b[i]->best_conf == 2)?
	    (Swap(b[i]->left->v[b[i]->l_v2],b[i]->left,b[i]->rght,b[i]->rght->v[b[i]->r_v1],tree)):
	    (Swap(b[i]->left->v[b[i]->l_v2],b[i]->left,b[i]->rght,b[i]->rght->v[b[i]->r_v2],tree));
	  b[i]->l = b[i]->ql[b[i]->best_conf-1];
	}

    }
}

/*********************************************************/
