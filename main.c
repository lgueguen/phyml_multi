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
#include "bionj.h"
#include "models.h"
#include "free.h"
#include "options.h"
#include "simu.h"
#include "unistd.h"
#include <sys/types.h>	/* for pid_t, fd_set, FD_ZERO, FD_SET, FD_ISSET */
#include <unistd.h>	/* for pipe, getpid, fork, close, select, read, write */
#include <fcntl.h> /* for making non-blocking pipes with  F_SETFL*/
//#include <signal.h> /* For signal handling*/
#include <sys/wait.h>
#include <errno.h>
#include "sigAnalyse.h"
#include "filelocker.h"

#ifdef PHYML

extern int final_stop; 

int main(int argc, char **argv)
{
  seq **data;
  allseq *alldata;
  seq **tmpData;
  allseq *tmpAlldata;
  option *input;
  char *s_tree, *s_any;
  FILE *fp_phyml_tree,*fp_phyml_stats,*fp_phyml_lk, *fp_phyml_siteLks;//,*fp_best_tree,*fp_best_tree_stats;
  arbre *tree;
  int n_data_sets;//, n_otu;
  matrix *mat;
  matrix *tmpMat;
  model *mod;
  time_t t_beg,t_end;
  div_t hour,min;
  int num_tree;
  //  double best_loglk;
  int *pid, i;  
  char numchar[12];
  int** serv_cli; /* file descriptors for pipes: server speaks to clients */
  int** cli_serv; /* file descriptors for pipes: clients speak to server*/
  // int** serv_cli_end; /* file descriptors for pipes: server speaks to clients to tell them it is soon over*/
  
  fd_set lire; /* descriptor set for select */
  double ** lk_vectors, ** lk_vectors_copy;	/* likelihood vectors sent by clients */
  double * lk_summary; /* A likelihood vector containing the sum of the likelihoods of all the trees but the one the client this vector is sent to works upon */ 
  double mean_loglk; /* the likelihood of the multitrees model */
  double *lk_vec;
  double *temp_other_vec;

  int *remainingClients; /* vector that keeps which clients are still in activity (have not converged yet)*/
  int numberOfRemainingClients;
  //  int improved=0;
  int maxi, site;//, maximum;
  int retour;
  //int t;
  double old_loglk, current_loglk;
  // int oneToSend=1;
  int status;
  int step=0;
  int temp;
  char *siteLksFile;
  double ** siteLks;
  double * failureCount; //To count the number of times a client does not succeed in improving its tree so that the likelihood improves
  int failureLimit = 20;
  int *listeTrees, *listeModels;
  int num_tree2;
  int ind=0;
  arbre **treeListe;
  int *ClientsToRemove;
  int newn_tree;
  //To deal with the file used for communications between the clients and the server.
  char* fname; 
  char* backup_starting_tree;
  char* starting_tree;
  
  /* int colNum, lineNum, j;
     double **lkmat;
  */
  fname = (char *)mCalloc(14,sizeof(char));
  backup_starting_tree= (char *)mCalloc(14,sizeof(char));
  starting_tree= (char *)mCalloc(14,sizeof(char));
  backup_starting_tree = "startingTree";
  //Modified random seed
  srand(getpid()*time(NULL));
  printf("PhyML_Multi : this program has been modified compared to the usual PhyML to try to detect recombination events.\n\t\n\n");
  // printf("Random seed : %d\n",rand()); 
  
  
  tree = NULL;
  mod  = NULL;

  Init_Constant();

  s_any = (char *)mCalloc(T_MAX_FILE,sizeof(char));

  fflush(stdout);


  
  input = (option *)Get_Input(argc,argv);

  Make_Model_Complete(input->mod);

  mod = input->mod;

 

  n_data_sets = 0;


  if(input->inputtree) Test_Multiple_Data_Set_Format(input);
  //else input->n_trees = 1;
  
  //if(input->n_data_sets > 1) input->n_trees = 1;

  //  best_loglk = UNLIKELY;

  printf("Number of trees expected : %d\n", input->n_trees );


  //TEMPORARY : we suppose at least 3 trees in the alignment
  /*  if (!input->inputtree) {
      input->n_trees = 3;
      }
  */


  treeListe=(arbre**)mCalloc(input->n_trees,sizeof(arbre*));


  //Allocates the space for the pipes
  serv_cli = (int **)mCalloc(input->n_trees,sizeof(int*)); //The server speaks to the clients
  cli_serv = (int **)mCalloc(input->n_trees,sizeof(int*)); //The clients speak to the server
  //  serv_cli_end = (int **)mCalloc(input->n_trees,sizeof(int*)); //The server speaks to the clients to tell them to go to the next step
  ///TEMPORARY
  //	  maximum=FALSE;
	 

  lk_summary=NULL;
	  
 
  do
    {

      n_data_sets++;

      time(&t_beg);

      //      n_otu = 0;

      if(n_data_sets > input->n_data_sets) 
        {
          data = NULL;	  
        }
      else
        {
          data = Get_Seq(input,0);
          /* 	  Print_Seq(data,n_otu); */
          /* 	  Exit(""); */
        }

      if(data)
        {
          if(n_data_sets > 1) printf("\n. Data set [#%d]\n",n_data_sets);

          printf("\n. Compressing sequences...\n");

          alldata = Compact_Seq(data,input);

          Free_Seq(data,alldata->n_otu);

          Init_Model(alldata,mod);
	  
          Check_Ambiguities(alldata,input->mod->datatype,input->mod->stepsize);

	
          fname=input->fname;

          lk_vectors = (double **)mCalloc(input->n_trees,sizeof(double*)); 
          For(num_tree,input->n_trees){
            lk_vectors[num_tree]=(double *)mCalloc(alldata->crunch_len,sizeof(double));
          }


          //	  mat = ML_Dist(alldata,mod);
	 
          //For each tree, a client is created whose job is to optimize it
          For(num_tree,input->n_trees)
            {
              input->current_tree = num_tree;
              if(!input->inputtree)
                {
                  //We read again the matrix, but we only use segments of length (total length/input->n_trees)
                  //This way we obtain input->n_trees different trees.
	
                  tmpData = Get_Seq_MultiTrees(input,  0);
                  if(tmpData)
                    {
                      tmpAlldata = Compact_Seq(tmpData,input);
		      
                      Free_Seq(tmpData,tmpAlldata->n_otu);
		       
                      Init_Model(tmpAlldata,mod);
		       
                      Check_Ambiguities(tmpAlldata,input->mod->datatype,input->mod->stepsize);
		       
                      printf("\n. Computing pairwise distances...\n");
		       
                      tmpMat = ML_Dist(tmpAlldata,mod);
		       
                      printf("\n. Building BIONJ tree...\n");
		       
                      tmpMat->tree = Make_Tree(tmpAlldata);
                      //printf ("After Make_Tree\n");
                      Bionj(tmpMat); //printf ("After bionj\n");
                      //Bionj_Br_Length(tmpMat); printf ("After Bionj_Br_Length\n");

		       
                      //The branch lengths are modified so that they are computed using the whole alignment, not a few sites. 
                      //This should reduce the risks of predicting a recombination where there are only rate differences
		      
                      tmpMat->method = 0;
		       
                      tree = tmpMat->tree;
                      mat = ML_Dist(alldata,mod);
                      mat->tree = tree;
                      mat->method = 0;
                      //Bionj_Br_Length(mat);
		      		       
                      Free_Mat(tmpMat);
                      Free_Mat(mat);
		       
                      /* 		       Order_Tree_CSeq(tree,alldata); */
                      /* 		       mat = ML_Dist(alldata,mod); */
		       
                      /* 		       mat->tree = tree; */
		       
                      /* 		       mat->method = 0; */
		       
                      /* 		       Bionj_Br_Length(mat); */
		       
                      /* 		       Free_Mat(mat); */


		       
                    }
                }
              else // We have many trees in the input file
                {
                  if(input->n_trees > 1) printf("\n. Reading user tree [#%d]\n",num_tree+1);
                  else printf("\n. Reading user tree...\n");
		  
                  if(input->n_trees == 1) rewind(input->fp_input_tree);

                  tree = Read_Tree_File(input->fp_input_tree);

                  if(!tree) 
                    {
                      printf("\n. Missing tree for data set #%d\n",n_data_sets);
                      printf("  This data set is not analyzed.\n");
                      data = NULL;
                    }
 
                  if(!tree->has_branch_lengths)
                    {
                      printf("\n. Computing branch length estimates...\n");
		      
                      Order_Tree_CSeq(tree,alldata);
		      
                      mat = ML_Dist(alldata,mod);
		      
                      mat->tree = tree;
		      
                      mat->method = 0;

                      Bionj_Br_Length(mat);

                      Free_Mat(mat);
                    }
                }
	      
              if(!tree) continue;
	      
	  
              tree->mod        = mod;
              tree->input      = input;
              tree->data       = alldata;
              tree->both_sides = 1;
              tree->n_pattern  = tree->data->crunch_len/tree->mod->stepsize;

              Order_Tree_CSeq(tree,alldata);

              Make_Tree_4_Lk(tree,alldata,alldata->init_len);

              if (tree->input->HMM) {
                tree->mod->message_size = tree->data->crunch_len * sizeof(double) + sizeof(double);  //we send the lambda parameter
              }
              else {
                tree->mod->message_size = tree->data->crunch_len * sizeof(double);
              }

              if (tree->input->HMM) {
                lk_summary = (double *)mCalloc(alldata->crunch_len+1,sizeof(double));  
                //We initiate lk_summary with 0s.
                For(i,alldata->crunch_len+1) {
                  lk_summary[i]=0.0;
                }
              }
              else {
                lk_summary = (double *)mCalloc(alldata->crunch_len,sizeof(double));  
                //We initiate lk_summary with 0s.
                For(i,alldata->crunch_len) {
                  lk_summary[i]=0.0;
                }
              }

	 
              //Allocating the space for the matrix containing the tree likelihoods   
              tree->mod->all_lk_vec = (double **)mCalloc(tree->input->n_trees,sizeof(double *));
              for (i=0;i<tree->input->n_trees;i++) {
                tree->mod->all_lk_vec[i] = (double *)mCalloc(tree->data->crunch_len,sizeof(double));
              }

              Lk(tree,tree->data, lk_summary, input->maxi);
              printf("Tree log-likelihood : %f\n", tree->tot_loglk);
              For (i, tree->data->crunch_len) {
                lk_vectors[num_tree][i]=tree->true_site_lk[i];
              }

              treeListe[num_tree] = tree;
	      

            }//We have computed the likelihoods for all trees; we can think of discarding some of these
	  
          //freeing the space for the matrix containing the tree likelihoods   
          for (i=0;i<tree->input->n_trees;i++) {
            free(tree->mod->all_lk_vec[i]);
          }
          free(tree->mod->all_lk_vec);
	  
          ClientsToRemove = (int *)mCalloc(tree->input->n_trees,sizeof(int));
          For(num_tree,tree->input->n_trees){
            ClientsToRemove[num_tree]=1; //0 when we want to remove
          }

          //Here we can examine the likelihoods and give them to the segmentation algorithm; this should limit the number of trees that are going to be optimized afterwards.
          //Checking whether we have the same tree (=same likelihoods) more than once
          For(num_tree,input->n_trees){
            for(num_tree2=num_tree+1; num_tree2<input->n_trees; num_tree2++){
              if (ClientsToRemove[num_tree2]==1) {
                ind=0;
                For(i,tree->data->crunch_len){
                  if (lk_vectors[num_tree][i]!=lk_vectors[num_tree2][i]) {
                    break;
                  }
                  else {
                    ind=ind+1;
                  }
                }
              }
              if (ind==tree->data->crunch_len) {
                printf ("Two trees among the first %d trees are identical, one is discarded\n", input->n_trees);
                ClientsToRemove[num_tree2]=0;
              }
            }
          }
	  
          //Now we should have removed all the trees that were identical. We set their values in the matrix to 0.0
	  
          For(num_tree,input->n_trees){
            if (ClientsToRemove[num_tree]==0) {
              printf("ClientsToRemove[num_tree]==0 : %d\n", num_tree);
              For(i,tree->data->crunch_len){
                lk_vectors[num_tree][i]=0.0;
              }
            }
          }
	  
          //We output the matrix and ask Sarment to segment it or we use our in-built phylo-HMM ?
          //The phylo-HMM is more secure, because with MPP we might have the problem of not finding a good optimum...
          //So we try the phylo-HMM
	  
          Optimize_Hmm_Lambda(tree, lk_vectors, &(tree->mod->hmm_lambda), 
                              tree->mod->hmm_lambda, 
                              0.01, 0.9999999999,
                              100, tree->input->n_trees);


          printf("New value of the auto-correlation parameter : %f\n", tree->mod->hmm_lambda);
	  
          listeModels = (int *)mCalloc(tree->data->crunch_len,sizeof(int));
          listeModels = Compute_Viterbi (tree, lk_vectors, tree->input->n_trees);
	  
          listeTrees = (int *)mCalloc(tree->input->n_trees,sizeof(int));

          //Now we analyse listeModels so that we can remove useless models
	  
          For(num_tree,tree->input->n_trees){
            listeTrees[num_tree]=0;
          }
          For(i,tree->data->crunch_len){
            listeTrees[listeModels[i]] = 1;
          }
	  
          For(num_tree,tree->input->n_trees){
            if (listeTrees[num_tree]==0) {
              ClientsToRemove[num_tree]=0;
            }
          }

          newn_tree=0;
          For(num_tree,tree->input->n_trees){
            if (ClientsToRemove[num_tree]==1) {//we count the trees that we keep
              newn_tree=newn_tree+1;
            }
          }

          printf("Tree log-likelihood : %f\n", tree->tot_loglk);
          printf("New Number of input trees after the removal of apparently useless ones : %d\n",newn_tree); 



          //Allocates the space for the clients pid
          pid = (int *)mCalloc(newn_tree,sizeof(int));
	  
          /* Creates the pipes */
          For(num_tree,newn_tree){
            serv_cli[num_tree] = (int*)mCalloc(2,sizeof(int));
            cli_serv[num_tree] = (int*)mCalloc(2,sizeof(int));
            if ( pipe(serv_cli[num_tree]) == -1)
              {
                printf("Error while creating pipe serv_cli\n");
                exit(EXIT_FAILURE);
              }
            if ( pipe(cli_serv[num_tree]) == -1)
              {
                printf("Error while creating pipe cli_serv\n");
                exit(EXIT_FAILURE);
              }
          }
          //Now we should have enough pipes
	  
          //Because clients won't need it, we can free it
          For(num_tree,input->n_trees){
            free(lk_vectors[num_tree]);
          }
          free(lk_vectors);

          num_tree=-1;
	  
          For(num_tree2,input->n_trees)
            {
              if (ClientsToRemove[num_tree2]==0) {
                printf("We have broken : num_tree : %d, num_tree2 : %d\n", num_tree,num_tree2);
                //break;
              }
              else {
                num_tree=num_tree+1;
                printf ("We have not broken: num_tree : %d, num_tree2 : %d\n", num_tree,num_tree2);
                //Outputting the starting trees
		
                sprintf(numchar, "%d", num_tree);
		
                strcpy(starting_tree,backup_starting_tree);
                strcat(starting_tree, numchar);
                fp_phyml_tree = Openfile(starting_tree,input->phyml_tree_file_open_mode);
                s_tree = Write_Tree(treeListe[num_tree]);
                fprintf(fp_phyml_tree,"%s\n",s_tree);
                Free(s_tree);
                fclose(fp_phyml_tree);



                //I guess it would be right to fork right now

                pid[num_tree]=fork();

                if (pid[num_tree]!=0) {//Server
                  //We do nothing here, so that we open the next tree and give it to the next client.
                  /*	if (close(cli_serv[num_tree][1]) == -1)
                      {
                      printf("Error while closing cli_serv[1]\n");
                      }
                      if (close(serv_cli[num_tree][0]) == -1)
                      {
                      printf("Error while closing serv_cli[0]\n");
                      }*/
                  /*	if (close(serv_cli_end[num_tree][0]) == -1)
                      {
                      printf("Error while closing serv_cli[0]\n");
                      } 
                  */
                  //fcntl(serv_cli[num_tree][1], F_SETFL, fcntl(serv_cli[num_tree][1], F_GETFL) | O_NONBLOCK);//write non blocking
                  //fcntl(cli_serv[num_tree][0], F_SETFL, fcntl(cli_serv[num_tree][0], F_GETFL) | O_NONBLOCK);//read non blocking
                }



                else { //Client
                  free(ClientsToRemove);
		
                  input->n_trees = newn_tree;
                  input->current_tree = num_tree;

                  tree = treeListe[num_tree];

                  tree->mod        = mod;
                  tree->input      = input;
                  tree->data       = alldata;
                  tree->both_sides = 1;
                  tree->n_pattern  = tree->data->crunch_len/tree->mod->stepsize;

                  //We close the right ends of the pipes
                  For (i, input->n_trees){
                    if (i==num_tree){// if we deal with the pipes that concern the right tree
                      if (close(cli_serv[num_tree][0]) == -1)
                        {
                          printf("Error while closing cli_serv[0]\n");
                        }
                      if (close(serv_cli[num_tree][1]) == -1)
                        {
                          printf("Error while closing serv_cli[1]\n");
                        }
                      //We make the read and the write non-blocking
                      // fcntl(serv_cli[num_tree][0], F_SETFL, fcntl(serv_cli[num_tree][0], F_GETFL) | O_NONBLOCK);//Read non blocking
                      //fcntl(cli_serv[num_tree][1], F_SETFL, fcntl(cli_serv[num_tree][1], F_GETFL) | O_NONBLOCK);//write non blocking
                      // fcntl(serv_cli_end[num_tree][0], F_SETFL, fcntl(serv_cli_end[num_tree][0], F_GETFL) | O_NONBLOCK);
		    
                    }
                    else { // these pipes do not concern the current client, so we close them all
                      if (close(serv_cli[i][0]) == -1)
                        {
                          printf("Error while closing serv_cli[0], %d ERRNO=%s\n", i, strerror(errno));
                        }
                      if (close(serv_cli[i][1]) == -1)
                        {
                          printf("Error while closing serv_cli[1], %d\n", i);
                        }
                      if (close(cli_serv[i][0]) == -1)
                        {
                          printf("Error while closing cli_serv[0], %d ERRNO=%s\n", i, strerror(errno));
                        }
                      if (close(cli_serv[i][1]) == -1)
                        {
                          printf("Error while closing cli_serv[1], %d ERRNO=%s\n", i, strerror(errno));
                        }
                    }
                  }
		
                  //This version is not able to do bootstrap
                  if (! tree->mod->bootstrap) {
                    if (tree->mod->s_opt->opt_topo){
                      Simu(tree,1000, cli_serv[num_tree][1], serv_cli[num_tree][0], input->maxi);//, serv_cli_end[num_tree][0]);
                      /* check_sigusr1();*/
                      //  pause(); 
                    }
                    else
                      {
                        if(tree->mod->s_opt->opt_free_param){
                          Round_Optimize(tree,tree->data, cli_serv[num_tree][1], serv_cli[num_tree][0], input->maxi);//, serv_cli_end[num_tree][0]);
                          /*	check_sigusr1();
                          //	do {}//printf("Waiting for the end\n");}
                          //	while(stop_next_point);*/
                          //	pause(); 
                        }
                        else
                          {   
                            //Allocating the space for the matrix containing the tree likelihoods   
                            tree->mod->all_lk_vec = (double **)mCalloc(tree->input->n_trees,sizeof(double *));
                            for (i=0;i<tree->input->n_trees;i++) {
                              tree->mod->all_lk_vec[i] = (double *)mCalloc(tree->data->crunch_len,sizeof(double));
                            }
                            lk_vec = (double *)mCalloc(tree->data->crunch_len,sizeof(double)); 
                            //listening to the user signal 2 in case the server tells to stop the computations
                            //	  check_sigusr2();
                            //Computing the likelihood
                            Lk(tree,tree->data, lk_summary, input->maxi);
                            //printf("\n. Log(lk) :               * -> %15.6f ",tree->tot_loglk);
                            For (i, tree->data->crunch_len) {
                              lk_vec[i]=tree->true_site_lk[i];
                            }
                            //Sending the likelihood vector
                            if(write(cli_serv[num_tree][1], lk_vec, tree->data->crunch_len*sizeof(double))!=tree->data->crunch_len*sizeof(double)) {
                              printf("Client %d write error\n", num_tree);
                              Exit("");
                            }
                            // do {printf("Waiting for the end\n");}
                            //while(!stop_next_point);
                            //pause();
                            //  while(!stop_next_point){};
                          }
                      }
                  }
                  else {
                    printf("This program does not do bootstrap. Sorry.\n");
                    Exit("");
                    //Bootstrap(tree);// I guess bootstrap is totally useless here
                  }

                  printf("Optimization finished for client %d", getpid());
                  //printf("HERERE\n");
                  s_tree = Write_Tree(tree);

                  //Creating the output file for this specific tree
                  sprintf(numchar, "%d", num_tree);
                  strcat(input->phyml_stat_file, numchar);
                  fp_phyml_stats = Openfile(input->phyml_stat_file,input->phyml_stat_file_open_mode);
		
                  fprintf(fp_phyml_stats,"\n- PHYML %s -\n\n", VERSION);
		
                  strcat(input->phyml_tree_file, numchar);
                  fp_phyml_tree = Openfile(input->phyml_tree_file,input->phyml_tree_file_open_mode);
		
                  strcat(input->phyml_lk_file, numchar);
                  fp_phyml_lk = fopen(input->phyml_lk_file,"w");
		
	
                  fprintf(fp_phyml_tree,"%s\n",s_tree);
                  Free(s_tree);
		
                  Unconstraint_Lk(tree);
	
                  if (input->n_data_sets==1)
                    Print_Fp_Out(fp_phyml_stats, t_beg, t_end, tree, input, n_data_sets);
                  else
                    Print_Fp_Out_Lines(fp_phyml_stats, t_beg, t_end, tree, input, n_data_sets);
	
		
                  fprintf(fp_phyml_lk,"%f\n",tree->tot_loglk);

                  fclose(fp_phyml_lk);
		
                  fclose(fp_phyml_tree);
		
                  fclose(fp_phyml_stats);
	 
		
                  Free_Tree_Lk(tree);
		
                  Free_Tree(tree);
		
		
                  Free_Cseq(alldata);
	  
                  Free_Model(mod);
		
                  if(input->fp_seq ) fclose(input->fp_seq );
                  if(input->fp_input_tree) fclose(input->fp_input_tree);
		
                  Free_Input(input);    
		
                  Free(s_any);
	
                  //	while (final_stop==FALSE) {
                  printf("Client pauses\n");
                  /*	check_sigusr1();
                      if (!final_stop)
                      pause();*/
                  printf("Job done for client %d\n", getpid());
                  //The client has finished its job
                  //	}

                  return 0; //The client has returned ! 
                }


              }
            } //End of the loop on the trees
	  
          /*//////////////////////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\*/
            /***********************THE SERVER INSTRUCTIONS********************************/
          /*\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\///////////////////////////////////////////*/
          printf("\n\tThe server has finished distributing the trees\n\n");

          //TEST : HMM only in the server
          //  input->HMM = TRUE;

          input->n_trees = newn_tree;
          //Closing all the unused extremities of the pipes
	     
          For(num_tree,input->n_trees){
	   
            if (close(cli_serv[num_tree][1]) == -1)
              {
                printf("Error while closing cli_serv[%d][1]\n", num_tree);
              }
            if (close(serv_cli[num_tree][0]) == -1)
              {
                printf("Error while closing serv_cli[%d][0]\n", num_tree);
              }
          }
	 

          //Allocating the likelihood vectors
          lk_vectors = (double **)mCalloc(input->n_trees,sizeof(double*)); 
          lk_vectors_copy = (double **)mCalloc(input->n_trees,sizeof(double*)); 
          temp_other_vec = (double *)mCalloc(tree->data->crunch_len,sizeof(double)); 
	  
          For(num_tree,input->n_trees){
            lk_vectors[num_tree]=(double *)mCalloc(tree->data->crunch_len,sizeof(double));
          }
          For(num_tree,input->n_trees){
            lk_vectors_copy[num_tree]=(double *)mCalloc(tree->data->crunch_len,sizeof(double));
          }
	  
          failureCount = (double *)mCalloc(tree->input->n_trees,sizeof(double));


          //If we do not want to optimize the trees
          if ((!(tree->mod->s_opt->opt_topo))&&(!tree->mod->s_opt->opt_free_param)) {
            // Now the server should receive the initial likelihoods from the clients	
            For(num_tree,input->n_trees){
              /*  temp=0;
                  while (temp<tree->data->crunch_len * sizeof(double)) {
                  temp+=read(cli_serv[num_tree][0], lk_vectors[num_tree], tree->data->crunch_len * sizeof(double));
                  }*/
              if(read(cli_serv[num_tree][0], lk_vectors[num_tree], tree->data->crunch_len * sizeof(double))< tree->data->crunch_len * sizeof(double)) {
                printf("Problem read in main l711\n"); 
              } 
              else{
                printf("Initial likelihood vector read from tree %d\n", num_tree);
              }
            }

            if (input->HMM){
              mean_loglk=Compute_HMM_Multi_Likelihood (tree, lk_vectors, input->n_trees);
            }
            else{
              mean_loglk=Compute_Simplified_Multi_Likelihood(tree, lk_vectors, input->n_trees, input->maxi);
            }
            //printf("mean_loglk : %f\n", mean_loglk);
            For(num_tree, input->n_trees) {
              kill(pid[num_tree], SIGUSR2);
            }
          }//END WE DO NOT WANT TO OPTIMIZE THE TREEs
	  
          else //If we do  want to optimize the trees
            {
              //Allocating the remaining clients vector
              remainingClients = (int *)mCalloc(input->n_trees,sizeof(int)); 
              For(num_tree,input->n_trees){
                remainingClients[num_tree]=1;
              }
              numberOfRemainingClients = input->n_trees;


              // Now the server should receive the initial likelihoods from the clients	
              For(num_tree,input->n_trees){
                if((temp=read(cli_serv[num_tree][0], lk_vectors[num_tree], tree->data->crunch_len * sizeof(double)))< tree->data->crunch_len * sizeof(double)) {
                  printf("Problem read in main l744 :  read %d\n", temp);
                }
                else {
                  printf("Initial likelihood vector read from tree %d\n", num_tree);
                }
              }
	      
	     

              //Writing into the file used for communications
              if (input->HMM) { 
                Write_To_Memory_With_Lambda (fname, lk_vectors, &(tree->mod->hmm_lambda), input->n_trees, tree->data->crunch_len);
              }
              else {
                //We've got all the likelihood vectors for all the trees
                //Building and sending the lk_summary to each client
                /*	For(num_tree,input->n_trees){
                    lk_summary = Compute_Lk_Summary(tree, lk_vectors, input->n_trees, input->maxi, num_tree);
                    write(serv_cli[num_tree][1],lk_summary,  tree->mod->message_size);
                    printf("Supplementary vector sent to client %d\n", num_tree);
                    }*/
                Write_To_Memory (fname, lk_vectors, input->n_trees, tree->data->crunch_len);
              }

              //We tell the clients that they can read the file now
              For(num_tree,input->n_trees){ 
                i=1;
                temp=write(serv_cli[num_tree][1],&i, sizeof(int));
                if (temp!=sizeof(int)) {
                  printf("ERROR WRITE MAIN l602 : %d\n", temp);
                } 
              }
              /////////////////////Computing the initial MultiTree Likelihood/////////////////////////
              if (input->HMM) {
                mean_loglk=Compute_HMM_Multi_Likelihood (tree, lk_vectors, input->n_trees);
              }
              else {
                mean_loglk=Compute_Simplified_Multi_Likelihood(tree, lk_vectors, input->n_trees, input->maxi);
              }
              printf("\t\t\tInitial Log-likelihood of the multiTrees model : %f\n",mean_loglk);
		
              For(num_tree,input->n_trees){
                failureCount[num_tree]=0;
              }


              /***************************************************************************/
              ////////////////////////////OPTIMIZATION LOOP////////////////////////////////
              /***************************************************************************/
              do
                {
                  step++;
                  if (tree->input->HMM) {
                    if (step%10==0){
                      printf ("\nstep : %d : Autocorrelation parameter before optimization : %f",step, tree->mod->hmm_lambda);
                      Optimize_Hmm_Lambda(tree, lk_vectors, &(tree->mod->hmm_lambda), 
                                          tree->mod->hmm_lambda, 
                                          0.01, 0.9999999999,
                                          100, tree->input->n_trees);
                      printf (" and then after : %f\n\n",tree->mod->hmm_lambda);
                    }
                  }
                  //improved=0; 
                  if(numberOfRemainingClients==0){
                    numberOfRemainingClients=-1;
                    break;
                  }
		   
                  old_loglk=current_loglk=mean_loglk;
                  printf("current_loglk : %f\n", current_loglk);
		    
                  ///////////////////Listening to the clients////////////////////
                  /* Reset the descriptors set for select */
                  FD_ZERO(&lire);
                  For(num_tree,input->n_trees){
                    if (remainingClients[num_tree]!=0) {
                      FD_SET(cli_serv[num_tree][0], &lire);
                    }
                  }
                  maxi=-1;
                  For(num_tree,input->n_trees){
                    if (remainingClients[num_tree]!=0) {
                      if (cli_serv[num_tree][0]>maxi) {
                        maxi=cli_serv[num_tree][0];
                      }
                    }
                  }
                  
                  //Building a copy to safely learn whether it is wise to accept the new vector from a client
                  For(num_tree,input->n_trees){
                    For(site, tree->data->crunch_len) {
                      lk_vectors_copy[num_tree][site]=lk_vectors[num_tree][site];
                    }
                  }
                  retour=-1;
                  do {
                    //The server is awaken if a client has something to say.
                    retour = select(maxi+1, &lire, NULL, NULL, NULL);
                  } while (retour==-1);

                  printf("stuck here %d\n",step);

                  if (retour <0) {
                    perror("select\n");
                    return -1;
                  }
                  if (retour==0) {
                    fprintf(stderr, "too long a wait for clients\n");
                    return -1;
                  }
                  ///////////////////Loop on the clients : who is talking ?////////////////////
                  For(num_tree,input->n_trees){
                    if (remainingClients[num_tree]!=0) {
                      if (FD_ISSET(cli_serv[num_tree][0],&lire))
                        {
                          if ((temp=read(cli_serv[num_tree][0], temp_other_vec, tree->data->crunch_len*sizeof(double)))!=tree->data->crunch_len*sizeof(double)) {
                            printf ("problem read in main : temp : %d, tree->data->crunch_len*sizeof(double) : %ld\n",temp, tree->data->crunch_len*sizeof(double));
                          }
                          else {
                            For (i, tree->data->crunch_len) {
                              lk_vectors[num_tree][i]=temp_other_vec[i];
                            }
                          }
                          //Need to test if the new vector provided by our client does provide 
                          //an improved likelihood or not. If it does, we take it ; if not, we don't
                          if (input->HMM) {
                            mean_loglk=Compute_HMM_Multi_Likelihood (tree, lk_vectors, input->n_trees);
                          }
                          else {
                            mean_loglk=Compute_Simplified_Multi_Likelihood(tree, lk_vectors, input->n_trees, input->maxi);
                          }
                          if(mean_loglk<=current_loglk) {//The likelihood is not good enough: the vector is discarded
                            printf("\n\t\t\tThe likelihood has not increased after reading a new likelihood vector\n");
                            For(site, tree->data->crunch_len) {
                              lk_vectors[num_tree][site]=lk_vectors_copy[num_tree][site];
                            }
                            mean_loglk=current_loglk;
                            failureCount[num_tree]=failureCount[num_tree]+1;
                            if (failureCount[num_tree]==failureLimit) {
                              if (remainingClients[num_tree]!=0) {
                                numberOfRemainingClients=numberOfRemainingClients-1;
                                printf("We remove an unefficient client !\n");    
                                kill(pid[num_tree], SIGUSR2);				    
                              }
                              remainingClients[num_tree]=0;
                            }
                          }
                          else {//The likelihood is better: we update the saving copy and the current lk
                            failureCount[num_tree]=0; //We reset the counter
                            //improved=1;
                            For(site, tree->data->crunch_len) {
                              lk_vectors_copy[num_tree][site]=lk_vectors[num_tree][site];
                            }
                            printf("\n\t\t\tNew log-likelihood of the multiTrees model :%f -> %f\n",current_loglk, mean_loglk);
                            current_loglk=mean_loglk;
                            /*  For(t, input->n_trees){//Building and sending the lk_summary to each client
                                if (remainingClients[t]!=0) {
                                lk_summary = Compute_Lk_Summary(tree, lk_vectors, input->n_trees, input->maxi, t);
                                temp=write(serv_cli[t][1],lk_summary, tree->mod->message_size);
                                if (temp!=tree->mod->message_size) {
                                printf("ERROR WRITE MAIN l594 : %d\n", temp);
                                } 
                                } 
                                }*/
                            //Writing into the file used for communications
                            if (input->HMM) {
                              Write_To_Memory_With_Lambda (fname, lk_vectors, &(tree->mod->hmm_lambda), input->n_trees, tree->data->crunch_len);
                            }
                            else {
                              Write_To_Memory (fname, lk_vectors, input->n_trees, tree->data->crunch_len);
                            }
                            if ((mean_loglk<old_loglk+MIN_DIFF_LK) && (step >10)){
                              if (remainingClients[num_tree]!=0) {
                                numberOfRemainingClients=numberOfRemainingClients-1;
                                printf("A client has converged !\n");    
                                kill(pid[num_tree], SIGUSR2);				    
                              }
                              remainingClients[num_tree]=0;
                              // For (i, input->n_trees) {printf("numtree : %d, remainingClients[i] : %d, number : %d\n",num_tree,remainingClients[i], numberOfRemainingClients);}
                            }
                          }
                        }//END if (FD_ISSET
                    }
                  }//END For(num_tree,input->n_trees){
                  if(mean_loglk>=old_loglk+MIN_DIFF_LK) {
                    /*  For(t, input->n_trees){//Building and sending the lk_summary to each client
                        lk_summary = Compute_Lk_Summary(tree, lk_vectors, input->n_trees, input->maxi, t);
                        temp=write(serv_cli[t][1],lk_summary, tree->mod->message_size);
                        if (temp!=tree->mod->message_size) {
                        // printf("ERROR WRITE MAIN l594 : %d\n", temp);
                        }
                        }*/
                    //Writing into the file used for communications
                    if (input->HMM) {
                      Write_To_Memory_With_Lambda (fname, lk_vectors, &(tree->mod->hmm_lambda), input->n_trees, tree->data->crunch_len);
                    }
                    else {
                      Write_To_Memory (fname, lk_vectors, input->n_trees, tree->data->crunch_len);
                    }
                  }
                  //This loop ends when the likelihood function does not increase significantly	  
                } while (numberOfRemainingClients>=0);//&&((mean_loglk>old_loglk+MIN_DIFF_LK/100)&&(improved==1)&&(step>10)));
              //while (((improved==1)&&(mean_loglk>old_loglk+MIN_DIFF_LK))||(improved==0)||(step<=10));
	  	      
              //We tell the clients that we enter the last round of optimizations
              printf("\n\n\t\t\tEntering the last round of optimizations, if any\n\n");
              /*     For(num_tree, input->n_trees) {
                     kill(pid[num_tree], SIGUSR1);
                     }
              */
              old_loglk=current_loglk=mean_loglk;
              /* 
              //We make the read non-blocking
              For(num_tree, input->n_trees) {
              fcntl(cli_serv[num_tree][0], F_SETFL, fcntl(cli_serv[num_tree][0], F_GETFL) | O_NONBLOCK);//read non blocking
              }
              //We read their new vectors, if any.
              For(num_tree, input->n_trees) {
              if ((temp=read(cli_serv[num_tree][0], lk_vectors[num_tree], tree->data->crunch_len * sizeof(double)))&&(temp<tree->data->crunch_len * sizeof(double))) {
              //printf("Error final read in main : read %d bytes instead of %d\n", temp, tree->data->crunch_len * sizeof(double));
              }
              else {
              //Need to test if the new vector provided by our client does provide an improved likelihood or not
              //If it does, we take it ; if not, we don't
              if (input->HMM) {
              mean_loglk=Compute_HMM_Multi_Likelihood (tree, lk_vectors, input->n_trees);
              }
              else {
              mean_loglk=Compute_Simplified_Multi_Likelihood(tree, lk_vectors, input->n_trees, input->maxi);
              }
              if((mean_loglk<old_loglk)||(isnan(mean_loglk))) {//The likelihood is not good enough : the vector is discarded
              mean_loglk=old_loglk;
              For(site, tree->data->crunch_len) {
              lk_vectors[num_tree][site]=lk_vectors_copy[num_tree][site];
              }
              }
              else {//We update the saving copy
              For(site, tree->data->crunch_len) {
              lk_vectors_copy[num_tree][site]=lk_vectors[num_tree][site];
              }
              printf("\t\t\tNew log-likelihood of the multiTrees model :%f -> %f\n\n",old_loglk, mean_loglk);
              old_loglk=mean_loglk;
              }
              }
              }//END For(num_tree, input->n_trees) {
              */
            }//END if we optimize the trees
          ///////////////////////////Waiting for the clients to stop to avoid zombis//////////
          /*	    For(num_tree, input->n_trees) {
                  kill (pid[num_tree], SIGUSR1);
                  } */
          /* For(num_tree, input->n_trees) {
             kill(pid[num_tree], SIGUSR1);
             }*/
          printf("WAITING\n");
          /* For(num_tree, input->n_trees) {
             while(wait(& status)>=0) {
             printf ("\nSTATUS : %d\n",status); 
             For(num_tree, input->n_trees) {
             kill(pid[num_tree], SIGUSR1);
             } 
             }
             }*/ 
		  
          for (i = 0; i<input->n_trees; i++) {
            //kill(pid[i], SIGUSR1);
            waitpid(pid[i], & status, 0);
            if (WIFEXITED(status)) {
              printf("Over, code=%d (normal behaviour)\n", WEXITSTATUS(status));
            } else if (WIFSIGNALED(status)) {
              printf("Killed by signal %d (normal behaviour)\n", WTERMSIG(status));
            } else if (WIFSTOPPED(status)) {
              printf("Killed by signal %d (normal behaviour)\n", WSTOPSIG(status));
            } else if (WIFCONTINUED(status)) {
              printf("Continued (normal behaviour)\n");
            }



            // printf ("\nPID : %d, STATUS : %d\n",pid[i], status);
          } 
		 

          /*
/////////////////////////////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\The END///////////////////////////////////////
          */

          printf("\t\t\tFinal Log-likelihood of the multiTrees model : %f\n",mean_loglk);

          //Creating and filling the vectors containing the real position site likelihoods.

          siteLks=(double **)mCalloc(input->n_trees,sizeof(double*));
          For(num_tree,input->n_trees){
            siteLks[num_tree]=(double *)mCalloc(tree->data->clean_len+1,sizeof(double));
          }
          For(site, tree->data->crunch_len) {
            For(i, tree->data->wght[site]) {
              For(num_tree,input->n_trees){
                siteLks[num_tree][tree->data->pospatt[site][i]]=lk_vectors[num_tree][site];
              }
            }
          }
	   
          //Outputting the site likelihoods
          //Open the file that is going to be used to output the results 
          siteLksFile = (char *)mCalloc(T_MAX_FILE,sizeof(char));
          strcpy(siteLksFile,input->seqfile);
          strcat(siteLksFile,"_phyml_siteLks.txt");
          fp_phyml_siteLks = Openfile(siteLksFile,1);
          //Outputting the header
          fprintf(fp_phyml_siteLks,"Site\t");
          For(i,input->n_trees){
            fprintf(fp_phyml_siteLks,"TreeLk%d\t", i);
          }
          For(i,input->n_trees){
            fprintf(fp_phyml_siteLks,"TreelogLk%d\t", i);
          }
          fprintf(fp_phyml_siteLks,"\n");
          //Outputting the site Lks
          For(site, tree->data->clean_len) {
            fprintf(fp_phyml_siteLks,"%d\t", site);
            For(num_tree,input->n_trees){
              fprintf(fp_phyml_siteLks,"%E\t", siteLks[num_tree][site]);
            }
            For(num_tree,input->n_trees){
              fprintf(fp_phyml_siteLks,"%f\t", log(siteLks[num_tree][site]));
            }
            fprintf(fp_phyml_siteLks,"\n");
          }

          fp_phyml_lk = fopen(input->phyml_lk_file,"w");
          fprintf(fp_phyml_lk,"logLikelihood : %f\n",mean_loglk);	 
          if (input->HMM) {
            fprintf(fp_phyml_lk,"Autocorrelation parameter : %f\n",tree->mod->hmm_lambda);
          }
	   
          time(&t_end);
	   
          hour = div(t_end-t_beg,3600);
          min  = div(t_end-t_beg,60  );
	   
          min.quot -= hour.quot*60;
          printf("\n\n. Time used %dh%dm%ds\n", hour.quot,min.quot,(int)(t_end-t_beg)%60);
          printf("\n%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n\n");
          fprintf(fp_phyml_lk,"\n\n. Time used %dh%dm%ds\n", hour.quot,min.quot,(int)(t_end-t_beg)%60);
          fclose(fp_phyml_lk);

          fclose(fp_phyml_siteLks);
          Free_Tree(tree);
	   
          Free_Cseq(alldata);
        }
    }while(data);
  Free_Model(mod);
  
  if(input->fp_seq ) fclose(input->fp_seq );
  if(input->fp_input_tree) fclose(input->fp_input_tree);
  
  Free_Input(input);    
  
  Free(s_any);
  remove(fname);

  return 0;
}

#endif
/*********************************************************/

/*	For(t, input->n_trees){//Building and sending the lk_summary to each client
    lk_summary = Compute_Lk_Summary(tree, lk_vectors, input->n_trees, input->maxi, t);
    if (write(serv_cli[t][1],lk_summary, tree->data->crunch_len * sizeof(double))!=tree->data->crunch_len * sizeof(double)) {
    printf("ERROR WRITE MAIN l540");
    }
    }*/

/* if ((temp=read(cli_serv[num_tree][0], lk_vectors[num_tree], tree->data->crunch_len * sizeof(double)))&&(temp<tree->data->crunch_len * sizeof(double))) {
   printf("Error read in main : read %d bytes instead of %d\n", temp, tree->data->crunch_len * sizeof(double));
   }
   else {*/
