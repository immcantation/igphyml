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
  int i,j,k;

  tree             = NULL;
  //mod              = NULL;
  best_lnL         = UNLIKELY;
  most_likely_size = -1.0;
  tree_size        = -1.0;
  
  r_seed = abs(4*(int)time(NULL)*(int)time(NULL)+4*(int)time(NULL)+1); //!< Modified by Marcelo
  //r_seed=1234;
  srand(r_seed);
  //SetSeed(r_seed);
  
  io = (option *)Get_Input(argc,argv); //!< Read the simulation options from interface or command line.

  io->command=mCalloc(T_MAX_LINE,sizeof(char));
    For(i,argc){
  	  strcat(io->command,argv[i]);
  	  strcat(io->command," ");
    }
    printf("COMMAND: %s\n",io->command);

  //io->r_seed = (io->r_seed<0)?r_seed:io->r_seed;
  io->r_seed=1234;
  if(io->mod->whichrealmodel > HLP17 && io->mod->partfilespec != 0){
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

    int* stopCodons=malloc(64*sizeof(int));
  	int* indexSenseCodons=malloc(64*sizeof(int));
  	int* senseCodons=malloc(64*sizeof(int));
  	 For (i,64){
  	    stopCodons[i]=0;
  	    senseCodons[i]=-1;
  	    indexSenseCodons[i]=-1;
  	  }

  	stopCodons[10]=1; //!< Set stop codons according to the genetic code.
  	stopCodons[11]=1;
  	stopCodons[14]=1;
  	j=0;
  	  For(i,64){ //!< Shift the sense codons down the array.
  	    if(!stopCodons[i]){
  	    	//printf("%d\t%d\n",i,j);
  	      senseCodons[j++]=i;
  	    }
  	  }
  	  For(i,61){ //!< Correct the index for recovery.
  	    indexSenseCodons[senseCodons[i]]=i;
  	  }
  	  io->stopCodons=stopCodons;
  	  io->senseCodons=senseCodons;
  	  io->indexSenseCodons=indexSenseCodons;


  Make_Model_Complete(io->mod);

  //previously global variables in Optimiz, now tree variables to ease in analyzing multiple data sets

  io->SIZEp=0;
  io->noisy=0;
  io->Iround=0;
  io->NFunCall=0;
  io->AlwaysCenter=0;
  io->gemin=1e-6;
  io->Small_Diff=.5e-6;
  io->both_sides=1;

  int last_otu=0;
  For(num_data_set,io->ntrees){
    n_otu = 0;
    best_lnL = UNLIKELY;
    if(io->mod->optDebug)printf("On data %d\n",num_data_set);
    if(io->mod->optDebug)printf("copying mode`l\n");
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

    io->n_trees=1;

    Get_Seq(io,mod);
    if(io->min_otu>0){
        	printf("\n%s\n..has %d sequences.",mod->in_align_file,mod->n_otu);
        	if(num_data_set > 0 && mod->n_otu > last_otu){
        		printf("\n repfile is not in order! --minSeq will not function properly. Exiting\n");
        		exit(EXIT_FAILURE);
        	}
        	last_otu=mod->n_otu;
        	if(mod->n_otu<io->min_otu){
        		io->ntrees=num_data_set;
        		continue;
        	}
    }

    if(io->threshold_exp) {if(num_data_set<io->dataset-1) {Free_Seq(mod->data,mod->n_otu);continue;}}//!< Added by Marcelo. Run arbitrary data set.


    Make_Model_Complete(mod);

    Set_Model_Name(mod);

    mod->nedges=2*mod->n_otu-3;

    //Print_Settings(io,mod);

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
      if(mod->freq_model == ROOT){
    	  	  mod->root_pi = mCalloc(mod->nomega_part,sizeof(phydbl*));
    	  	  For(i,io->mod->nomega_part)mod->root_pi[i]=mCalloc(61,sizeof(phydbl));
      }
      if(io->mod->optDebug)printf("compacting data\n");
      cdata = Compact_Data(mod->data,io,mod);
      if(io->mod->optDebug)printf("compacted data\n");
      Free_Seq(mod->data,cdata->n_otu);
	
      if(cdata) Check_Ambiguities(cdata,mod->datatype,mod->state_len);
      else{
    	  PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
	  	  Warn_And_Exit("");
      }
      io->n_trees=1;

      if(io->mod->optDebug){
    	  printf("%d\n",io->n_trees);
    	  printf("%d\n",io->mod->s_opt->random_input_tree);
    	  printf("%d\n",io->mod->s_opt->n_rand_starts);
      }
      for(num_tree=(io->n_trees == 1)?(0):(num_data_set);num_tree < io->n_trees;num_tree++){

    	  if(!io->mod->s_opt->random_input_tree) io->mod->s_opt->n_rand_starts = 1;

    	if(io->mod->optDebug)printf("random starts: %d\n",io->mod->s_opt->n_rand_starts);
	For(num_rand_tree,io->mod->s_opt->n_rand_starts){
	  if((io->mod->s_opt->random_input_tree) && (io->mod->s_opt->topo_search != NNI_MOVE))
	  if(!io->quiet) PhyML_Printf("\n. [Random start %3d/%3d]\n",num_rand_tree+1,io->mod->s_opt->n_rand_starts);

	  io->init_run=1; //! Added by Marcelo.
	  if(io->mod->optDebug)printf("initing model\n");
	  Init_Model(cdata,mod,io);
	  if(io->mod->optDebug)printf("initing model\n");

	  io->init_run=0;

	  if(io->in_tree<=1 || strcmp(mod->in_tree_file,"N")==0){
	  	  tree = Dist_And_BioNJ(cdata,mod,io);
	   }else{
	  	  mod->fp_in_tree = Openfile(mod->in_tree_file,0);
	  	  tree = Read_User_Tree(cdata,mod,io);
	  }
	  /*switch(io->in_tree){
	    case 0 : case 1 : { tree = Dist_And_BioNJ(cdata,mod,io); break; }
	    case 2 :          { tree = Read_User_Tree(cdata,mod,io); break; }
	  }*/
	  if(io->mod->optDebug)printf("read tree\n");
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

	  //io->tree          = tree;  //will need to change this
	  tree->mod         = mod;
	  tree->io          = io;
	  tree->data        = cdata;
	  tree->both_sides  = 1;
	  tree->n_pattern   = tree->data->crunch_len;

	  //mod->quiet=NO; //turn back on alerts

	  //added by Ken
	  //Find location of root node if HLP17
	  if(mod->whichrealmodel <= HLP17){
		  if(io->mod->optDebug)printf("rootname %s\n",mod->rootname);
	  	  mod->startnode = -1;
	      int nodepos;
	      for(nodepos=0;nodepos<((tree->n_otu-1)*2);nodepos++){
	      	if(strcmp(tree->noeud[nodepos]->name,mod->rootname)==0){
	      		mod->startnode=nodepos;
	      		//update ancestors of tree, now that root node is found
	      		//PhyML_Printf("\n. Start node found: %d %s\n",nodepos,mod->rootname);
	      		Update_Ancestors_Edge(tree->noeud[nodepos],tree->noeud[nodepos]->v[0],tree->noeud[nodepos]->b[0],tree);
	      		//PhyML_Printf("\n. Start node found: %d %s\n",nodepos,mod->rootname);
	      		For(i,mod->nomega_part && mod->freq_model==ROOT){
	      			For(j,mod->ns){
	      				tree->noeud[nodepos]->partfreqs[i][j]=mod->root_pi[i][j];
	      				if(i==0)mod->pi[j]=mod->root_pi[i][j];
	      				//printf("%d\t%lf\n",j,mod->pi[j]);
	      			}
	      		}
	      	}
	      }
	     if(mod->startnode==-1){
	    	 PhyML_Printf("\n\nRoot sequence ID not found in data file! %s %s\n",mod->rootname,mod->in_align_file);
	    	 exit(EXIT_FAILURE);
	     }
	  }
	  if(io->mod->optDebug)printf("setting up partition model\n");
	  //Set up default partition model if necessary
	  if(tree->mod->partfilespec==0){
	     tree->mod->partIndex = (int *)mCalloc(tree->n_pattern,sizeof(int));
	     int indexi;
	     for(indexi=0;indexi<tree->n_pattern;indexi++){
	      	 mod->partIndex[indexi]=0;
	     }
	     mod->nparts=1;
	     mod->nomega_part=mod->nparts;
    	 mod->partNames = (char**)mCalloc(mod->nomega_part,sizeof(char*));
    	 mod->partNames[0]=(char*)mCalloc(T_MAX_OPTION,sizeof(char));
    	 strcpy(mod->partNames[0],"SINGLE");
	  }

	  if(mod->s_opt->random_input_tree) Random_Tree(tree);
	    
	  //if((!num_data_set) && (!num_tree) && (!num_rand_tree)) Check_Memory_Amount(tree);
	  //For(bri,2*tree->n_otu-3) printf("in main2: %lf\n",tree->t_edges[bri]->l);
	  if(io->mod->optDebug)printf("prepping tree for lhood\n");
	  if(io->mod->s_opt->opt_subst_param || io->mod->s_opt->opt_bl || !io->precon){
		  Prepare_Tree_For_Lk(tree);
	  }
	  if(io->mod->optDebug)printf("prepping tree for lhood\n");


	  if(tree->mod->ambigprint && tree->mod->whichrealmodel <= HLP17){
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
	  //For(i,12)io->mod->base_freq[i]=mod->base_freq[i];
	  //For(i,12)io->mod->uns_base_freq[i]=mod->uns_base_freq[i];
     }

      }//for(num_tree=(io->n_trees == 1)?(0):(num_data_set);num_tree < io->n_trees;num_tree++)

   /*   Free_Cseq(cdata);*/

      fclose(mod->fp_in_align);
    }else{
      PhyML_Printf("\n. No data was found.\n");
      PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
      Warn_And_Exit("");
    }
    //Free_Model_Complete(mod);
  } //For(num_data_sets

  if(io->mod->motifstringopt==0 && io->mod->modeltypeOpt <= HLP17){
  	printf("\nDEFAULT: Estimating only symmetric WRC/GYW hotness."
  		   "\n........ Use '--motifs FCH' to model all hot and coldspots (not recommended with small datasets)"
    	   "\n........ Use '--motifs' and '--hotness' for more motif models.\n");
  }
  if(!io->mod->omega_opt_spec){
      int Npart=0;
      	For(i,io->ntrees){
      		if(strcmp(io->partfs[i],"N")==0){
      			Npart++;
      		}
      	}
      	if(Npart == io->ntrees && io->mod->modeltypeOpt != GY){ //if no partition files specified, use single omega
      		printf("\nDEFAULT: No partition file(s) specified so impossible to partition omega by FWR/CDRs.\n");
      	}else if(io->mod->modeltypeOpt != GY){
      		printf("\nDEFAULT: Partition file(s) specified so partitioning omega by FWR/CDRs."
      				"\n........ Use '--omegaOpt e' if you just want one omega.\n");
     }
  }
  if(io->threads > io->ntrees)printf("\nWarning: Number of threads (%d) exceeds number of lineages (%d).\n"
		  "........ This will not speed up computations beyond %d threads.\n",io->threads,io->ntrees,io->ntrees);
  if(io->threads==1 && io->ntrees > 1)printf("\nDEFAULT: Running multiple trees on one thread."
		  "\n........ Use the '--threads' option to specify more (might speed things up).\n");

if(io->mod->freq_model != ROOT){
  //Set up equilibrium base freqs for the repertoire
  if(io->eq_freq_handling != USER){
  For(i,io->ntrees){
	  For(j,12){
		  io->mod->baseCounts[j]+=io->mod_s[i]->baseCounts[j];
		  if(io->mod->optDebug)printf("%d\t%lf\n",j,io->mod->baseCounts[j]);
	  }
  }
  io->mod->base_freq[0]= (io->mod->baseCounts[0]) /(io->mod->baseCounts[0]+io->mod->baseCounts[1]+io->mod->baseCounts[2]+io->mod->baseCounts[3]);
  io->mod->base_freq[1]= (io->mod->baseCounts[1]) /(io->mod->baseCounts[0]+io->mod->baseCounts[1]+io->mod->baseCounts[2]+io->mod->baseCounts[3]);
  io->mod->base_freq[2]= (io->mod->baseCounts[2]) /(io->mod->baseCounts[0]+io->mod->baseCounts[1]+io->mod->baseCounts[2]+io->mod->baseCounts[3]);
  io->mod->base_freq[3]= (io->mod->baseCounts[3]) /(io->mod->baseCounts[0]+io->mod->baseCounts[1]+io->mod->baseCounts[2]+io->mod->baseCounts[3]);
  io->mod->base_freq[4]= (io->mod->baseCounts[4]) /(io->mod->baseCounts[4]+io->mod->baseCounts[5]+io->mod->baseCounts[6]+io->mod->baseCounts[7]);
  io->mod->base_freq[5]= (io->mod->baseCounts[5]) /(io->mod->baseCounts[4]+io->mod->baseCounts[5]+io->mod->baseCounts[6]+io->mod->baseCounts[7]);
  io->mod->base_freq[6]= (io->mod->baseCounts[6]) /(io->mod->baseCounts[4]+io->mod->baseCounts[5]+io->mod->baseCounts[6]+io->mod->baseCounts[7]);
  io->mod->base_freq[7]= (io->mod->baseCounts[7]) /(io->mod->baseCounts[4]+io->mod->baseCounts[5]+io->mod->baseCounts[6]+io->mod->baseCounts[7]);
  io->mod->base_freq[8]= (io->mod->baseCounts[8]) /(io->mod->baseCounts[8]+io->mod->baseCounts[9]+io->mod->baseCounts[10]+io->mod->baseCounts[11]);
  io->mod->base_freq[9]= (io->mod->baseCounts[9]) /(io->mod->baseCounts[8]+io->mod->baseCounts[9]+io->mod->baseCounts[10]+io->mod->baseCounts[11]);
  io->mod->base_freq[10]=(io->mod->baseCounts[10])/(io->mod->baseCounts[8]+io->mod->baseCounts[9]+io->mod->baseCounts[10]+io->mod->baseCounts[11]);
  io->mod->base_freq[11]=(io->mod->baseCounts[11])/(io->mod->baseCounts[8]+io->mod->baseCounts[9]+io->mod->baseCounts[10]+io->mod->baseCounts[11]);
  For(j,12)io->mod->base_freq[j]=roundf(io->mod->base_freq[j]*10000.0f)/10000.0f; //round base freqs to 5 decimal places to help make comparisons to previous versions
  if(io->mod->optDebug)For(j,12){printf("%lf\t%lf\n",io->mod->baseCounts[j],io->mod->base_freq[j]);}
  CF3x4(io->mod->base_freq, io->mod->genetic_code);
  }else{
	  //CF3x4(io->mod->base_freq, io->mod->genetic_code);
	  For(i,12)io->mod->base_freq[i]=io->mod->user_b_freq[i];
  }
  if(io->mod->optDebug)printf("cf3x4\n");


  //JUST TO ELIMINATE SOME VARIABLES
  /*For(j,12)io->mod->base_freq[j]=io->mod_s[0]->base_freq[j];
  For(j,12)io->mod->uns_base_freq[j]=io->mod_s[0]->uns_base_freq[j];*/

  if(io->mod->optDebug)For(j,12){printf("%lf\t%lf\n",io->mod->baseCounts[j],io->mod->base_freq[j]);}
  Freq_to_UnsFreq(io->mod->base_freq,   io->mod->uns_base_freq,   4, 1);
  Freq_to_UnsFreq(io->mod->base_freq+4, io->mod->uns_base_freq+4, 4, 1);
  Freq_to_UnsFreq(io->mod->base_freq+8, io->mod->uns_base_freq+8, 4, 1);
  if(io->mod->optDebug)For(j,12){printf("%lf\t%lf\t%lf\t%lf\n",io->mod->baseCounts[j],io->mod->base_freq[j],io->mod->uns_base_freq[j],io->mod_s[0]->base_freq[j]);}
}
  int nparams=0;
    int nmodparams=2+io->mod->nomega_part+io->mod->nhotness+12;
    nparams += nmodparams*(io->ntrees+1);
    For(j,io->ntrees){
  	  nparams += 2*io->tree_s[j]->n_otu-3;
    }
  if(io->mod->optDebug)printf("\nNeed to store %d parameters",nparams);
  io->paramStore=mCalloc(nparams,sizeof(phydbl));
  io->nparams=nparams;
  io->mod->s_opt->min_diff_lk_global = io->min_diff_lk_global;

#if defined OMP || defined BLAS_OMP
  io->t_beg=omp_get_wtime();
#else
  time(&io->t_beg);
#endif

  if(io->mod->optDebug)printf("about to do stuff\n");
  if(io->testInitTree){ //!< Added by Marcelo.
  	    //!< Do nothing!
  }else if(io->lkExperiment){ //!< Added by Marcelo.
  		  //lkExperiment(tree,num_tree);
  }else if(tree->mod->s_opt->opt_topo){
	  if(tree->mod->s_opt->topo_search      == NNI_MOVE){
		  if(io->mod->optDebug)printf("model type %d\n",io->mod->whichrealmodel);
		  if(io->mod->whichrealmodel<=HLP17){
			  if(io->mod->optDebug)printf("here0\n");
			  For(i,io->ntrees){
				  //Lk(io->tree_s[i]);
				  Get_UPP(io->tree_s[i]->noeud[io->tree_s[i]->mod->startnode], io->tree_s[i]->noeud[io->tree_s[i]->mod->startnode]->v[0], io->tree_s[i]);
			  }
		  }
		  if(io->mod->optDebug)printf("here\n");
		  Simu_Loop(io);
	  }
  	  else if(tree->mod->s_opt->topo_search == SPR_MOVE){
  		  if(io->mod->optDebug)printf("about to do spr moves\n");
  		  io->mod->print_trace=0;
  		  io->mod_s[0]->print_trace=0;
  		  if(io->mod->whichrealmodel<=HLP17){
  			  For(i,io->ntrees){
  				 // Lk(io->tree_s[i]);
  				  Get_UPP(io->tree_s[i]->noeud[io->tree_s[i]->mod->startnode], io->tree_s[i]->noeud[io->tree_s[i]->mod->startnode]->v[0], io->tree_s[i]);
  				 // exit(EXIT_FAILURE);
  			  }
  		  }
  	  	  Speed_Spr_Loop(io);
  		  //Speed_Spr_Loop(io->tree_s[0]);
  	  }else{
  		  Lazy_Exit("Best of NNI and SPR",__FILE__,__LINE__);
  	    //Best_Of_NNI_And_SPR(tree);
  	  }
  }else{
  	  if(io->mod->s_opt->opt_subst_param || io->mod->s_opt->opt_bl){
  	    	Round_Optimize(io,ROUND_MAX*2); //! *2 Added by Marcelo to match codeml values.
  	  }else{
  		if(io->mod->optDebug)printf("doing lhood\n");
  	    io->mod->update_eigen=1;
  	    io->both_sides=1;
  	    if(!io->precon){
  	    	Lk_rep(io);
  	    	Print_Lk_rep(io,"Repertoire likelihood!");
  	    	For(i,io->ntrees){
  	    		Print_Lk(io->tree_s[i],"Subtree");
  	    	}
  	    }
  	  }
  }

  if(io->mod->freq_model!=ROOT){
  //convert base frequencies back to properly output results
   if(io->mod->optDebug)For(j,12){printf("before %lf\t%lf\n",io->mod->baseCounts[j],io->mod->base_freq[j]);}
   Freq_to_UnsFreq(io->mod->base_freq,   io->mod->uns_base_freq,   4, 0);
   Freq_to_UnsFreq(io->mod->base_freq+4, io->mod->uns_base_freq+4, 4, 0);
   Freq_to_UnsFreq(io->mod->base_freq+8, io->mod->uns_base_freq+8, 4, 0);
   if(io->mod->optDebug)For(j,12){printf("after %lf\t%lf\t%lf\t%lf\n",io->mod->baseCounts[j],io->mod->base_freq[j],io->mod->uns_base_freq[j],io->mod_s[0]->base_freq[j]);}
  }
  io->both_sides = 1;
  io->mod->update_eigen=1;
  if(!io->precon){
	  Lk_rep(io);
	  Print_Lk_rep(io,"Final likelihood");
  }

  //Estimate confidence intervals using profile likelihood curves
  if(io->CIest>0){
  	char* CIf = (char*)malloc(T_MAX_FILE*sizeof(char));
  	strcpy(CIf,io->mod->in_align_file);
  	strcat(CIf,"_igphyml_CIlog.txt");
  	if(io->append_run_ID){
  		strcat(CIf,"_");
  		strcat(CIf,io->run_id_string);
  	}
  	FILE* CI = Openfile(CIf, 1 );//openOutputFile(mod->out_trace_tree_file, "_igphyml_tree_trace", ".txt", io);
  io->mod->s_opt->print=1;
  For(i,io->ntrees)io->tree_s[i]->mod->s_opt->print=0;
  findCIs(io->mod,io,CI);
  For(i,io->ntrees){
	  findCIs(io->mod_s[i],io,CI);
  }
  fclose(CI);
  Print_Lk_rep(io,"Final likelihood after CI estimation");
  }
  if(io->mod->ASR){
	  For(j,io->ntrees){
		  t_tree* tree=io->tree_s[j];
		  model* mod = io->mod_s[j];
		  Get_UPP(tree->noeud[tree->mod->startnode],tree->noeud[tree->mod->startnode]->v[0],tree);
		  mod->mlASR=mCalloc(mod->nedges+1,sizeof(int*));
		  mod->probASR=mCalloc(mod->nedges+1,sizeof(phydbl*));
		  mod->mlCodon=mCalloc(mod->nedges+1,sizeof(char*));
  	  	  For(i,mod->nedges+1){
  	  		  mod->mlASR[i]=mCalloc(mod->init_len/3,sizeof(int));
  	  		  mod->probASR[i]=mCalloc(mod->init_len/3,sizeof(phydbl));
  			  mod->mlCodon[i]=mCalloc(mod->init_len,sizeof(char));

  	  		  //printf("LK ON EDGE %d\n",i);
  	  		  //tree->mod->num=-1;
  	  		  if(i<mod->nedges)ASR_At_Given_Edge(tree->t_edges[i],tree,0);
  	  		  else ASR_At_Given_Edge(tree->noeud[tree->mod->startnode]->b[0],tree,1);
  	  		  //printf("\n%s\n",mod->mlCodon[i]);
  	  	  }
	  }
  }
  if(io->precon==2 || io->precon==-2 || io->precon==4 || io->precon==-4 || io->precon==-6){
	  printf("\n");
	  int precon = io->precon;

io->threads=0;
#if defined OMP || defined BLAS_OMP
#pragma omp parallel for if(io->splitByTree)
#endif
	For(i,io->ntrees){ //do likelihood calculations in parallel
		t_tree* tree;
#if defined OMP || defined BLAS_OMP
#pragma omp critical
#endif
		{
			tree= io->tree_s[io->threads++];
		}
	  	  t_node* r = tree->noeud[tree->mod->startnode];
	  	  Init_Class_Tips(tree,precon); //no io
	  	  int pars1 = Fill_Sankoff(r,tree,1); //no io
	  	  //printf("\n\n. %d Initial maximum parsimony score: %d",i,pars1);
	  	  Set_Pars_Counters(r,tree,1); //no io
	  	  //Get_First_Path(r,0,tree,1);
	  	  //printf("\n. Resolving polytomies using isotype information");
	  	  int pars2 = Resolve_Polytomies_Pars(tree,0.001);
	  	  #if defined OMP || defined BLAS_OMP
		  #pragma omp critical
		  #endif
	  	  {
	  		  printf("\n. %d Initial/resolved maximum parsimony score: %d %d %s",tree->mod->num,pars1,pars2,tree->mod->rootname);
	  	  }
	  }
  }

  #if defined OMP || defined BLAS_OMP
  t_end=omp_get_wtime();
  #else
  time(&t_end);
  #endif
  Print_IgPhyML_Out(io);
  //exit(EXIT_FAILURE);

  if(tree->io->print_site_lnl) Print_Site_Lk(tree,io->fp_out_lk);

  if(io->precon){
	  parsReconstructions(io);
  }



  	  /* Start from BioNJ tree */
  	  if((num_rand_tree == io->mod->s_opt->n_rand_starts-1) && (tree->mod->s_opt->random_input_tree)){
  	    /* Do one more iteration in the loop, but don't randomize the tree */
  	    num_rand_tree--;
  	    tree->mod->s_opt->random_input_tree = 0;
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

