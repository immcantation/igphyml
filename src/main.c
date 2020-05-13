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


int main(int argc, char **argv){
  calign *cdata; //!< Pointer that will hold the input sequences.
  option *io;    //!< Pointer for the simulation options.
  t_tree *tree;   //!< Pointer for a tree
  int n_otu,/*!< number of taxa.*/ num_data_set; /*!< for multiple data sets.*/
  int num_tree,num_rand_tree;
  time_t t_beg,	/*!< Start Time.*/ t_end; /*!< Stop Time.*/ 
  phydbl best_lnL,most_likely_size,tree_size;
  int r_seed;
  char *most_likely_tree=NULL;
  int i,j,k;

  tree             = NULL;
  best_lnL         = UNLIKELY;
  most_likely_size = -1.0;
  tree_size        = -1.0;
  
  /*r_seed = abs(4*(int)time(NULL)*(int)time(NULL)+4*(int)time(NULL)+1); //!< Modified by Marcelo

  //r_seed=1234;
  srand(r_seed);
  SetSeed(r_seed);*/
  
  io = (option *)Get_Input(argc,argv); //!< Read the options from the command line.

  if(io->mod->whichrealmodel > HLP17 && io->mod->partfilespec != 0){
	  Warn_And_Exit("\n. Site-partitioned omega only available with HLP17 right now. Sorry."
	  "\n. This is mostly due to laziness, so feel free to complain to Ken about this."
	  "\n. In the meantime, try running -m HLP17 --hotness 0 -f empirical --rootpi instead of -m GY\n\n");
  }

  //declare data structures for upper model
  Make_Model_Complete(io->mod);
   //   printf("\n state_freq %d",io->mod->s_opt->opt_state_freq);


  //read in each dataset and set up respective model
  int last_otu=0;
  For(num_data_set,io->ntrees){
    n_otu = 0;
    best_lnL = UNLIKELY;
    if(io->mod->optDebug)printf("On data %d\ncopying mode`l\n",num_data_set);
    model *mod = Copy_Partial_Model(io->mod,num_data_set); //!< Pointer that will hold the model applied.
    if(io->mod->optDebug)printf("copied model\n");
    mod->num=num_data_set;
    mod->quiet=YES;
    strcpy(mod->in_tree_file,io->treefs[num_data_set]); //copy input tree to model
    strcpy(mod->in_align_file,io->datafs[num_data_set]); //copy input data to model
    strcpy(mod->rootname,io->rootids[num_data_set]); //copy root name
    if(io->mod->optDebug){ //Check if sequence and tree files exist
    	printf("\n. tree file: %s\t",mod->in_tree_file);
    	printf("align file: %s\t",mod->in_align_file);
    	printf("root name: %s",mod->rootname);
    }
    if(!Filexists(mod->in_align_file)) {
    	char* tmp = (char *) mCalloc (T_MAX_FILE, sizeof(char));
        strcpy(tmp, "\n. The alignment file '");
        strcat(tmp, mod->in_align_file);
        strcat(tmp, "' does not exist.\n");
        Warn_And_Exit(tmp);
    }
    if(strcmp(mod->in_tree_file,"N")!=0 && !Filexists(mod->in_tree_file)){
    	 char* tmp = (char *) mCalloc (T_MAX_FILE, sizeof(char));
    	 strcpy(tmp, "\n. The tree file '");
    	 strcat(tmp, mod->in_tree_file);
    	 strcat(tmp, "' does not exist.\n");
    	 Warn_And_Exit(tmp);
    }
    mod->fp_in_align = Openfile(mod->in_align_file,0);

    Get_Seq(io,mod);
    if(!mod->data){
      PhyML_Printf("\n. No data was found.\n");
      PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
      Warn_And_Exit("");
    }

    if(io->min_otu>0){//check if minimum sequence count has been reached
       	if(num_data_set > 0 && mod->n_otu > last_otu)
       		Warn_And_Exit("\n repfile is not in order! --minseq will not function properly. Exiting\n");
       	last_otu=mod->n_otu;
       	if(mod->n_otu<io->min_otu){
       		printf("\n%s\n..has %d sequences.",mod->in_align_file,mod->n_otu);
       		io->ntrees=num_data_set;
       		continue;
       	}
    }


    //declare data structures for sub-model
    Make_Model_Complete(mod);
    Set_Model_Name(mod);

    //process data and count codons
    cdata = Compact_Data(mod->data,io,mod);
    Free_Seq(mod->data,cdata->n_otu);

    if(cdata) Check_Ambiguities(cdata,mod->datatype,mod->state_len);
    else{
    	PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
	  	Warn_And_Exit("");
    }

   if(!io->mod->s_opt->random_input_tree) io->mod->s_opt->n_rand_starts = 1;
   if((io->mod->s_opt->random_input_tree) && (io->mod->s_opt->topo_search != NNI_MOVE)) io->init_run=1; //! Added by Marcelo.
   if(!io->precon)Init_Model(cdata,mod,io);
   io->init_run=0;
   if(io->in_tree<=1 || strcmp(mod->in_tree_file,"N")==0){
	   tree = Dist_And_BioNJ(cdata,mod,io); //estimate initial tree topology
   }else{
	   mod->fp_in_tree = Openfile(mod->in_tree_file,0);
	   tree = Read_User_Tree(cdata,mod,io); //read in user tree topology
     fclose(mod->fp_in_tree);
   }
   if(io->mod->optDebug)printf("read tree\n");
   if(!tree) continue;

    #if defined OMP || defined BLAS_OMP
   	   t_beg=omp_get_wtime();
   	   tree->t_beg=t_beg;
    #else
   	   time(&t_beg);
   	   time(&(tree->t_beg));
    #endif

   	tree->mod         = mod;
    tree->io          = io;
    tree->data        = cdata;
    tree->both_sides  = 1;
    tree->n_pattern   = tree->data->crunch_len;
    mod->tree         = tree;

    //Find location of root node if HLP
    if(mod->whichrealmodel <= HLP17)Get_Root_Pos(mod,tree,io);
    if(mod->permute_tips == 1)Permute_Tips(tree);
    if(tree->mod->partfilespec==0){ //Set up default partition model if necessary
    	tree->mod->partIndex = (int *)mCalloc(tree->n_pattern,sizeof(int));
    	For(i,tree->n_pattern)mod->partIndex[i]=0;
    }
    if(mod->s_opt->random_input_tree) Random_Tree(tree);
    if(io->mod->s_opt->opt_subst_param || io->mod->s_opt->opt_bl || !io->precon)Prepare_Tree_For_Lk(tree);
    if(tree->mod->ambigprint && tree->mod->whichrealmodel <= HLP17){
    	FILE *ambigfile = fopen(tree->mod->ambigfile, "w");
    	if (ambigfile == NULL)Warn_And_Exit("Error opening ambig file!\n");
    	Print_Ambig_States(tree->noeud[tree->mod->startnode],tree,ambigfile);
    	fclose(ambigfile);
    	printf("\n. Printed out ambiguous states to %s\n",tree->mod->ambigfile);
    }else if(tree->mod->ambigprint){
    	printf("\n. Can only print ambiguous characters with HLP models\n");
    }

    //Store model and data in upper hierarchy
    io->mod_s[num_data_set]=mod;
    io->tree_s[num_data_set]=tree;
    fclose(mod->fp_in_align);
  } //For(num_data_sets

  if(io->mod->permute_tips == 2)Permute_All_MetaData(io, 0);

  //Print user warnings
  if(!io->mod->omega_opt_spec){
      int Npart=0;
      For(i,io->ntrees)if(strcmp(io->partfs[i],"N")==0)Npart++;
      	if(Npart == io->ntrees && io->mod->modeltypeOpt != GY){ //if no partition files specified, use single omega
      		printf("\nDEFAULT: No partition file(s) specified so impossible to partition omega by FWR/CDRs.\n");
      	}else if(io->mod->modeltypeOpt != GY){
      		printf("\nDEFAULT: Partition file(s) specified so partitioning omega by FWR/CDRs."
      				"\n........ Use '--omega e' if you just want one omega.\n");
     }
  }
  if(io->threads > io->ntrees)printf("\nWarning: Number of threads (%d) exceeds number of lineages (%d).\n"
		  "........ This will not speed up computations beyond %d threads.\n",io->threads,io->ntrees,io->ntrees);
  if(io->threads==1 && io->ntrees > 1)printf("\nDEFAULT: Running multiple trees on one thread."
		  "\n........ Use the '--threads' option to specify more (might speed things up).\n");

  //Set up base frequencies and additional data structures
  if(!io->precon)Setup_Repertoire_Models(io);

  //Do topology and parameter estimation if requested
  if(!io->testInitTree && !io->lkExperiment && !io->precon){
  	if(tree->mod->s_opt->opt_topo){ //Estimate topology?
	  if(tree->mod->s_opt->topo_search   == NNI_MOVE)
		  Simu_Loop(io);
  	  else if(tree->mod->s_opt->topo_search == SPR_MOVE){
  		  io->mod->print_trace=0;
  		  io->mod_s[0]->print_trace=0;
  	  	  Speed_Spr_Loop(io);
  	  }else Lazy_Exit("Best of NNI and SPR",__FILE__,__LINE__);
  	}else{
  	  if(io->mod->s_opt->opt_subst_param || io->mod->s_opt->opt_bl)
  	       Round_Optimize(io,ROUND_MAX*2); //*2 Added by Marcelo to match codeml values.
  	}
  }

   //convert base frequencies back to properly output results
  if(io->mod->freq_model<ROOT && io->mod->s_opt->opt_state_freq==YES){
   Freq_to_UnsFreq(io->mod->base_freq,   io->mod->uns_base_freq,   4, 0);
   Freq_to_UnsFreq(io->mod->base_freq+4, io->mod->uns_base_freq+4, 4, 0);
   Freq_to_UnsFreq(io->mod->base_freq+8, io->mod->uns_base_freq+8, 4, 0);
  }

  io->both_sides = 1;
  io->mod->update_eigen=1;
  if(!io->precon){
	  Lk_rep(io);
	  Print_Lk_rep(io,"Final likelihood");
  }

  //Estimate confidence intervals using profile likelihood curves
  if(io->CIest>0)CI_Wrapper(io);


    //ASR, if desired
  if(io->mod->ASR)ASR_Wrapper(io);

  //re-arrange topology based on parsimony model, if desired
  if(io->precon==2 || io->precon==-2 || io->precon==4 ||
		  io->precon==-4 || io->precon==-6|| io->precon==7)
	  Prepars_Wrapper(io);

  t_tree* t = io->tree_s[0];
  /*Init_Class_Tips(tree,io->precon);
  printf("pars: %d\n",Fill_Sankoff(tree->noeud[tree->mod->startnode],tree,1));*/

  #if defined OMP || defined BLAS_OMP
  t_end=omp_get_wtime();
  #else
  time(&t_end);
  #endif
  //output data!
  if(io->out_stats_format != OUTTXT){
	  Print_Tab_Out(io);
  }else{
	  Print_IgPhyML_Out(io);
  }
  //exit(1);
  //parsimony reconstructions, if desired
  if(io->precon)Pars_Reconstructions(io);

  //TODO: Free data structures

  #if defined OMP || defined BLAS_OMP
     
  t_end=omp_get_wtime();
    
  #else
    
  time(&t_end);
    
  #endif
    
  Print_Time_Info(t_beg,t_end);
    
  return 0;
}

