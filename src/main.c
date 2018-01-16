/*
IgPhyML: a program that computes maximum likelihood phylogenies under
non-reversible codon models designed for antibody lineages.

Copyright (C) Kenneth B Hoehn. Sept 2016 onward.

built upon

codonPHYML: a program that  computes maximum likelihood phylogenies from
CODON homologous sequences.

Copyright (C) Marcelo Serrano Zanetti. Oct 2010 onward.

built upon

PhyML:  a program that  computes maximum likelihood phylogenies from
DNA or AA homologous sequences.

Copyright (C) Stephane Guindon. Oct 2003 onward.

All parts of the source except where indicated are distributed under
the GNU public licence. See http://www.opensource.org for details.

*/
/*! \file main.c
    \brief The main file
*/ 
#include "spr.h"
#include "utilities.h"
#include "lk.h"
#include "optimiz.h"
#include "bionj.h" 
#include "models.h"
#include "free.h" 
#include "help.h"
#include "simu.h"
#include "eigen.h"
#include "pars.h" 
#include "alrt.h"
#include "io.h"


#ifdef PHYML

int main(int argc, char **argv){
  
  calign *cdata; //!< Pointer that will hold the input sequences.
  option *io; //!< Pointer for the simulation options.
  t_tree *tree; //!< Pointer for a tree
  int n_otu,/*!< number of taxa.*/ num_data_set; /*!< for multiple data sets.*/ 
  int num_tree,tree_line_number,num_rand_tree; 
  matrix *mat;
  time_t t_beg,	/*!< Start Time.*/ t_end; /*!< Stop Time.*/ 
  phydbl best_lnL,most_likely_size,tree_size;
  int r_seed;
  char *most_likely_tree=NULL;

  tree             = NULL;
  //mod              = NULL;
  best_lnL         = UNLIKELY;
  most_likely_size = -1.0;
  tree_size        = -1.0;
  
  r_seed = abs(4*(int)time(NULL)*(int)time(NULL)+4*(int)time(NULL)+1); //!< Modified by Marcelo
  r_seed=1234;
  srand(r_seed);
  SetSeed(r_seed);
  
  io = (option *)Get_Input(argc,argv); //!< Read the simulation options from interface or command line.
  //io->r_seed = (io->r_seed<0)?r_seed:io->r_seed;
  io->r_seed=1234;
  if(io->mod->whichrealmodel != HLP17 && io->mod->partfilespec != 0){
	  printf("\n. Site-partitioned omega only available with HLP17 right now. Sorry.");
	  printf("\n. This is mostly due to laziness, so feel free to complain to Ken about this.");
	  printf("\n. In the meantime, try running -m HLP17 --hotness 0 instead of -m GY\n\n");
	  exit(EXIT_FAILURE);
  }
  
  mat = NULL;
  tree_line_number = 0;
  
  //ADDED BY KEN
  //TURNS OFF ALRT
  io->ratio_test = 0;
  io->n_trees=1;
  if((io->n_data_sets > 1) && (io->n_trees > 1)){
    io->n_data_sets = MIN(io->n_trees,io->n_data_sets);
    io->n_trees     = MIN(io->n_trees,io->n_data_sets);
  }

  Make_Model_Complete(io->mod);

  For(num_data_set,io->ntrees){
    n_otu = 0;
    best_lnL = UNLIKELY;
    printf("On data %d\n",num_data_set);
    printf("copying model\n");
    model *mod = Copy_Partial_Model(io->mod); //!< Pointer that will hold the model applied.
    printf("copied model\n");

    strcpy(mod->in_tree_file,io->treefs[num_data_set]); //copy input tree to model
    strcpy(mod->in_align_file,io->datafs[num_data_set]); //copy input data to model
    strcpy(mod->rootname,io->rootids[num_data_set]); //copy root name

    printf("%s\n",mod->in_tree_file);
    printf("%s\n",mod->in_align_file);
    printf("%s\n",mod->rootname);

    if(io->in_tree == 2)mod->fp_in_tree = Openfile(mod->in_tree_file,0);
    mod->fp_in_align = Openfile(mod->in_align_file,0);

    if(mod->print_trace) {
    	strcpy(mod->out_trace_stats_file,mod->in_align_file);
    	strcat(mod->out_trace_stats_file,"_igphyml_stats_trace");
    	/*strcat(mod->out_trace_stats_file,num_data_set);
    	strcat(mod->out_trace_stats_file,"_");*/
    	strcpy(mod->out_trace_tree_file,mod->in_align_file);
    	strcat(mod->out_trace_tree_file,"_igphyml_tree_trace");
    	/*strcat(mod->out_trace_tree_file,"_");
    	strcat(mod->out_trace_tree_file,num_data_set);
    	strcat(mod->out_trace_tree_file,"_");*/
    	printf("%s %s\n",mod->out_trace_stats_file,mod->out_trace_tree_file);
        mod->fp_out_tree_trace = Openfile(mod->out_trace_tree_file, 1 );//openOutputFile(mod->out_trace_tree_file, "_igphyml_tree_trace", ".txt", io);
        mod->fp_out_stats_trace = Openfile(mod->out_trace_stats_file, 1 );//openOutputFile(mod->out_trace_stats_file, "_igphyml_stats_trace", ".txt", io);
        printf("made output files\n");
    }

    printf("about to test data format\n");
    if(io->in_tree == 2) Test_Multiple_Data_Set_Format(io,mod);
    else io->n_trees = 1;
    printf("tested data format\n");


    Get_Seq(io,mod);
    //if(io->convert_NT_to_AA) Conv_NT_seq_to_AA_seq(io); unused ken 5/1

    if(io->threshold_exp) {if(num_data_set<io->dataset-1) {Free_Seq(mod->data,mod->n_otu);continue;}}//!< Added by Marcelo. Run arbitrary data set.

    Make_Model_Complete(mod);

    Set_Model_Name(mod);

    Print_Settings(io,mod);
	  printf("printed settings\n");

    //mod = io->mod;
    if(mod->data){
      if(io->n_data_sets > 1) PhyML_Printf("\n. Data set [#%d]\n",num_data_set+1);
      //Added by Ken 18/8
      //Don't compress data if doing a model with multiple partitions
      if(io->mod->partfilespec==1){
      	  io->colalias = 0;
      }else{
    	  io->colalias = 0;//for now, don't compress sequences at all
      }
      printf("compacting data\n");
      cdata = Compact_Data(mod->data,io,mod);
      printf("compacted data\n");
      Free_Seq(mod->data,cdata->n_otu);
	
      if(cdata) Check_Ambiguities(cdata,mod->datatype,mod->state_len);
      else{
    	  PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
	  	  Warn_And_Exit("");
      }
      io->n_trees=1;

      printf("printed ambigs\n");
      printf("%d\n",io->n_trees);
      printf("%d\n",io->mod->s_opt->random_input_tree);
      printf("%d\n",io->mod->s_opt->n_rand_starts);
      for(num_tree=(io->n_trees == 1)?(0):(num_data_set);num_tree < io->n_trees;num_tree++){

    	  if(!io->mod->s_opt->random_input_tree) io->mod->s_opt->n_rand_starts = 1;

    	 printf("random starts: %d\n",io->mod->s_opt->n_rand_starts);
	For(num_rand_tree,io->mod->s_opt->n_rand_starts){
	  if((io->mod->s_opt->random_input_tree) && (io->mod->s_opt->topo_search != NNI_MOVE))
	  if(!io->quiet) PhyML_Printf("\n. [Random start %3d/%3d]\n",num_rand_tree+1,io->mod->s_opt->n_rand_starts);

	  io->init_run=1; //! Added by Marcelo.
	  printf("initing model\n");
	  Init_Model(cdata,mod,io);
	  printf("initing model\n");

	  io->init_run=0;

	  switch(io->in_tree){
	    case 0 : case 1 : { tree = Dist_And_BioNJ(cdata,mod,io); break; }
	    case 2 :          { tree = Read_User_Tree(cdata,mod,io); break; }
	  }
	  printf("read tree\n");
	 /* int bri=0;
	  For(bri,2*tree->n_otu-3) printf("in main: %lf\n",tree->t_edges[bri]->l);*/
	  if(!tree) continue;
	    
	  #if defined OMP || defined BLAS_OMP

	  t_beg=omp_get_wtime();
	  tree->t_beg=t_beg;
	    
	  #else
	    
	  time(&t_beg);
	  time(&(tree->t_beg));
	   
	  #endif

	  io->tree          = tree;  //will need to change this
	  tree->mod         = mod;
	  tree->io          = io;
	  tree->data        = cdata;
	  tree->both_sides  = 1;
	  tree->n_pattern   = tree->data->crunch_len;

	  //previously global variables in Optimiz, now tree variables to ease in analyzing multiple data sets
	  tree->SIZEp=0;
	  tree->noisy=0;
	  tree->Iround=0;
	  tree->NFunCall=0;
	  tree->AlwaysCenter=0;
	  tree->gemin=1e-6;
	  tree->Small_Diff=.5e-6;

	  //added by Ken
	  //Find location of root node if HLP17
	  if(mod->whichrealmodel == HLP17){
		  printf("rootname %s\n",mod->rootname);
	  	  mod->startnode = -1;
	      int nodepos;
	      for(nodepos=0;nodepos<((tree->n_otu-1)*2);nodepos++){
	      	if(strcmp(io->tree->noeud[nodepos]->name,mod->rootname)==0){
	      		mod->startnode=nodepos;
	      		//update ancestors of tree, now that root node is found
	      		Update_Ancestors_Edge(io->tree->noeud[nodepos],io->tree->noeud[nodepos]->v[0],io->tree->noeud[nodepos]->b[0],tree);
	      		PhyML_Printf("\n. Start node found: %d %s\n",nodepos,mod->rootname);
	      	}
	      }

	     if(mod->startnode==-1){
	    	 PhyML_Printf("\n\nRoot sequence ID not found in data file!\n");
	    	 exit(EXIT_FAILURE);
	     }
	  }
	  printf("setting up partition model\n");
	  //Set up default partition model if necessary
	  if(tree->mod->partfilespec==0){
	     tree->mod->partIndex = (int *)mCalloc(tree->n_pattern,sizeof(int));
	     int indexi;
	     for(indexi=0;indexi<tree->n_pattern;indexi++){
	      	 mod->partIndex[indexi]=0;
	     }
	     mod->nparts=1;
	     mod->nomega_part=mod->nparts;
    	 mod->partNames = (char**)mCalloc(mod->nparts,sizeof(char*));
    	 mod->partNames[0]=(char*)mCalloc(T_MAX_OPTION,sizeof(char));
    	 strcpy(mod->partNames[0],"SINGLE");
	  }

	  if(mod->s_opt->random_input_tree) Random_Tree(tree);
	    
	  if((!num_data_set) && (!num_tree) && (!num_rand_tree)) Check_Memory_Amount(tree);
	  //For(bri,2*tree->n_otu-3) printf("in main2: %lf\n",tree->t_edges[bri]->l);
	  printf("prepping tree for lhood\n");
	  Prepare_Tree_For_Lk(tree);
	  printf("prepping tree for lhood\n");


	  if(tree->mod->ambigprint && tree->mod->whichrealmodel == HLP17){
		  FILE *ambigfile = fopen(tree->mod->ambigfile, "w");
		  if (ambigfile == NULL){
		      printf("Error opening ambig file!\n");
		      exit(0);
		  }
		  Print_Ambig_States(tree->noeud[tree->mod->startnode],tree,ambigfile);
		  fclose(ambigfile);
		  printf("\n. Printed out ambiguous states to %s\n",tree->mod->ambigfile);
	  }else if(tree->mod->ambigprint){
		  printf("\n. Can only print ambiguous characters with HLP17 model\n");
	  }

	  io->mod_s[num_data_set]=mod;
	  io->tree_s[num_data_set]=tree;

      /*tree->br_len_invar_set = NO;
	    
	  if(io->in_tree == 1) Spr_Pars(tree);
	 
	  if(io->testInitTree){ //!< Added by Marcelo.
	    //!< Do nothing!	    
	  }else if(io->lkExperiment){ //!< Added by Marcelo.
	    lkExperiment(tree,num_tree);
	  }
	  else if(tree->mod->s_opt->opt_topo){
	    if(tree->mod->s_opt->topo_search      == NNI_MOVE) Simu_Loop(tree);
	    else if(tree->mod->s_opt->topo_search == SPR_MOVE){
	  	  Speed_Spr_Loop(tree);
	    }else{
	    	Best_Of_NNI_And_SPR(tree);
	    }
	  }else{
	    if(tree->mod->s_opt->opt_subst_param || tree->mod->s_opt->opt_bl){
	    	Round_Optimize(tree,tree->data,ROUND_MAX*2); //! *2 Added by Marcelo to match codeml values.
	    }else{
	    	printf("doing lhood\n");
	    	Lk(tree);
	    }
	  }
	  printf("Finished lhood\n");
	  

	  tree->both_sides = 1;
	  Lk(tree);
	  Pars(tree);
	  Get_Tree_Size(tree);
	  PhyML_Printf("\n. Log-likelihood of the current tree:\t\t\t\t%.2f\n",tree->c_lnL);

	  if(tree->mod->whichrealmodel==HLP17){
		  Update_Ancestors_Edge(io->tree->noeud[io->tree->mod->startnode],io->tree->noeud[io->tree->mod->startnode]->v[0],io->tree->noeud[io->tree->mod->startnode]->b[0],tree);
	  }
	  
	  // Print the tree estimated using the current random (or BioNJ) starting tree
	  if(io->mod->s_opt->n_rand_starts > 1){
	    Br_Len_Involving_Invar(tree);
	    Print_Tree(io->fp_out_trees,tree);
	    fflush(NULL);
	  }
	    
	  // Record the most likely tree in a string of characters
	  if(tree->c_lnL > best_lnL){
	    best_lnL = tree->c_lnL;
	    Br_Len_Involving_Invar(tree);
	    if(most_likely_tree) free(most_likely_tree);
	    most_likely_tree = Write_Tree(tree);
	    most_likely_size = Get_Tree_Size(tree);
	  }
	    
	  #if defined OMP || defined BLAS_OMP
	    
	  t_end=omp_get_wtime();
	    
	  #else
	    
	  time(&t_end);
	    
	  #endif
	    printf("printing output\n");
	  Print_Fp_Out(io->fp_out_stats,t_beg,t_end,tree,
	   io,num_data_set+1,
	   (tree->mod->s_opt->n_rand_starts > 1)?
	   (num_rand_tree):(num_tree),mod);//was io->mod
	  printf("printed to output\n");
			 
	  if(tree->io->print_site_lnl) Print_Site_Lk(tree,io->fp_out_lk);
			 
	  // Start from BioNJ tree
	  if((num_rand_tree == io->mod->s_opt->n_rand_starts-1) && (tree->mod->s_opt->random_input_tree)){
	    // Do one more iteration in the loop, but don't randomize the tree
	    num_rand_tree--;
	    tree->mod->s_opt->random_input_tree = 0;
	  }
			 
	  Free_Spr_List(tree);
	  Free_One_Spr(tree->best_spr);
	  if(tree->mat) Free_Mat(tree->mat);
	  Free_Triplet(tree->triplet_struct);
	  Free_Tree_Pars(tree);
	  Free_Tree_Lk(tree);
	  Free_Tree(tree);*/
	}
	/*printf("Ratio test? %d\n",io->ratio_test);

	if(io->testInitTree){ //!< Added by Marcelo
	  //!< Do nothing!
	}
	else
	// Launch bootstrap analysis
	if(mod->bootstrap){
	  if(!io->quiet) PhyML_Printf("\n. Launch bootstrap analysis on the most likely tree...\n");
	  most_likely_tree = Bootstrap_From_String(most_likely_tree,cdata,mod,io);
	}
	else if(io->ratio_test){
	  // Launch aLRT
	  if(!io->quiet) PhyML_Printf("\n. Compute aLRT branch supports on the most likely tree...\n");
	  most_likely_tree = aLRT_From_String(most_likely_tree,cdata,mod,io);
	}
	  
	// Print the most likely tree in the output file
	if(!io->quiet) PhyML_Printf("\n. Printing the most likely tree in file:\n" ".... '%s'\n", Basename(io->out_tree_file));
	if(io->n_data_sets == 1) rewind(io->fp_out_tree);
 
	PhyML_Fprintf(io->fp_out_tree,"%s\n",most_likely_tree);
	  
	if(io->fp_out_compare){ //!< Added by Marcelo.
	  fprintf(io->fp_out_compare,"%s\n",most_likely_tree);
	}
	  
	if(io->n_trees > 1 && io->n_data_sets > 1) break;*/


      }//for(num_tree=(io->n_trees == 1)?(0):(num_data_set);num_tree < io->n_trees;num_tree++)


   /*   Free_Cseq(cdata);*/
    }else{
      PhyML_Printf("\n. No data was found.\n");
      PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
      Warn_And_Exit("");
    }
    //Free_Model_Complete(mod);
  } //For(num_data_sets

  printf("about to do stuff\n");
 /*Originally had these declared as private, however, they will be uninitialized in the openmp loop :-/
  *  int testInitTree = io->testInitTree;
  int lkExperimentv = io->lkExperiment;
  int in_tree = io->in_tree;
  int ntrees = io->ntrees;
  int RM = ROUND_MAX;*/

  //omp_lock_t *writelock;
 // omp_init_lock(writelock);
  printf("about to do stuff\n");

#if defined OMP || defined BLAS_OMP

//#pragma omp parallel for private(in_tree,lkExperimentv,testInitTree) if(io->splitByTree)
#pragma omp parallel for if(io->splitByTree)


#endif
  For(num_data_set,io->ntrees){
	  //omp_set_lock(writelock);
	  printf("DOING LHOODS ON %d\n",num_data_set);
	  model* mod = io->mod_s[num_data_set];
	  t_tree* tree = io->tree_s[num_data_set];
	  mod->num=num_data_set;
	  tree->br_len_invar_set = NO;
	  int testInitTree = io->testInitTree;
	  int lkExperimentv = io->lkExperiment;
	  int in_tree = io->in_tree;
	  int ntrees = io->ntrees;
	  int RM = ROUND_MAX;
	  //omp_unset_lock(writelock);

	  if(in_tree == 1) Spr_Pars(tree);

	  if(testInitTree){ //!< Added by Marcelo.
	    //!< Do nothing!
	  }else if(lkExperimentv){ //!< Added by Marcelo.
		  lkExperiment(tree,num_tree);
	  }
	  else if(tree->mod->s_opt->opt_topo){
	    if(tree->mod->s_opt->topo_search      == NNI_MOVE) Simu_Loop(tree);
	    else if(tree->mod->s_opt->topo_search == SPR_MOVE){
	    	printf("about to do spr moves\n");
	  	  Speed_Spr_Loop(tree);
	    }else{
	    	Best_Of_NNI_And_SPR(tree);
	    }
	  }else{
	    if(tree->mod->s_opt->opt_subst_param || tree->mod->s_opt->opt_bl){
	    	Round_Optimize(tree,tree->data,RM*2); //! *2 Added by Marcelo to match codeml values.
	    }else{
	    	printf("doing lhood\n");
	    	Lk(tree);
	    }
	  }
	  //omp_set_lock(writelock);
	  printf("Finished lhood\n");
	  //omp_unset_lock(writelock);


	  tree->both_sides = 1;
	  Lk(tree);
	  Pars(tree);
	  Get_Tree_Size(tree);
	  //omp_set_lock(writelock);
	  PhyML_Printf("\n. Log-likelihood of the current tree:\t\t\t\t%.2f\t%d\n",tree->c_lnL, mod->num);
	  //omp_unset_lock(writelock);


	  tree->both_sides = 1;
	  Lk(tree);
	  Pars(tree);
	  Get_Tree_Size(tree);
	  PhyML_Printf("\n. Log-likelihood of the current tree:\t\t\t\t%.2f\n",tree->c_lnL);

	  if(tree->mod->whichrealmodel==HLP17){
	  	  Update_Ancestors_Edge(tree->noeud[tree->mod->startnode],tree->noeud[tree->mod->startnode]->v[0],tree->noeud[tree->mod->startnode]->b[0],tree);
	   }

	  if(io->ntrees==1){
	  	  // Print the tree estimated using the current random (or BioNJ) starting tree
	  	  if(mod->s_opt->n_rand_starts > 1){
	  	    Br_Len_Involving_Invar(tree);
	  	    printf("about to print tree\n");
	  	    Print_Tree(io->fp_out_tree,tree);
	  	    printf("about to print tree\n");
	  	    fflush(NULL);
	  	  }

	  	  // Record the most likely tree in a string of characters
	  	  if(tree->c_lnL > best_lnL){
	  	    best_lnL = tree->c_lnL;
	  	    Br_Len_Involving_Invar(tree);
	  	    if(most_likely_tree) free(most_likely_tree);
	  	    most_likely_tree = Write_Tree(tree);
	  	    most_likely_size = Get_Tree_Size(tree);
	  	  }

	  	  #if defined OMP || defined BLAS_OMP

	  	  t_end=omp_get_wtime();

	  	  #else

	  	  time(&t_end);

	  	  #endif
	  	    printf("printing output\n");
	  	  Print_Fp_Out(io->fp_out_stats,t_beg,t_end,tree,
	  	   io,num_data_set+1,
	  	   (tree->mod->s_opt->n_rand_starts > 1)?
	  	   (num_rand_tree):(num_tree),mod);//was io->mod
	  	  printf("printed to output\n");

	  	  //if(tree->io->print_site_lnl) Print_Site_Lk(tree,io->fp_out_lk);
	  }

  }
  


  /*if(most_likely_tree) free(most_likely_tree);
  
  if(io->mod->s_opt->n_rand_starts > 1) PhyML_Printf("\n. Best log likelihood: %.2f\n",best_lnL);
  
  Free_Optimiz(io->mod->s_opt);
  if(io->mod->whichmodel==GTR) Free_Custom_Model(io->mod); //!< Added by Marcelo.
  Free_Model_Basic(io->mod);
  
  fflush(NULL);
  
  //if(io->fp_in_align)    fclose(io->fp_in_align);
  //if(io->mod->fp_in_tree)     fclose(io->mod->fp_in_tree);
  if(io->fp_out_lk)      fclose(io->fp_out_lk);
  if(io->fp_out_tree)    fclose(io->fp_out_tree);
  if(io->fp_out_trees)   fclose(io->fp_out_trees);
  if(io->fp_out_stats)   fclose(io->fp_out_stats);
  if(io->mod->print_trace){
	  fclose(io->fp_out_tree_trace);
	  fclose(io->fp_out_stats_trace);
  }
  if(io->fp_out_ps)      fclose(io->fp_out_ps); //!< Added by Marcelo.
  if(io->fp_out_compare) fclose(io->fp_out_compare); //!< Added by Marcelo.
	
  Free_Input(io);
    */
  #if defined OMP || defined BLAS_OMP
     
  t_end=omp_get_wtime();
    
  #else
    
  time(&t_end);
    
  #endif
    
  Print_Time_Info(t_beg,t_end);
    
  return 0;
}

#endif

