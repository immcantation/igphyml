/*
IgPhyML: a program that computes maximum likelihood phylogenies under
non-reversible codon models designed for antibody lineages.

Copyright (C) Kenneth B Hoehn. Sept 2016 onward.

built upon

codonPHYML: a program that  computes maximum likelihood phylogenies from
CODON homologous sequences.

Copyright (C) Marcelo Serrano Zanetti. Oct 2010 onward.

built upon

PHYML :  a program that  computes maximum likelihood  phylogenies from
DNA or AA homologous sequences 

Copyright (C) Stephane Guindon. Oct 2003 onward

All parts of  the source except where indicated  are distributed under
the GNU public licence.  See http://www.opensource.org for details.

*/


#include "pars.h"
#include "io.h"

extern int     stopCodons[64];
extern int   senseCodons[64];
extern char aminoAcidmap[65];
extern int indexSenseCodons[64];




/*********************************************************/

void Prepars_Wrapper(option* io){
	int i;
	int precon = io->precon;

io->threads=0;
#if defined OMP || defined BLAS_OMP
#pragma omp parallel for if(io->splitByTree)
#endif
	For(i,io->ntrees){ //do parsimony-based rearrangements in parallel
		t_tree* tree;
#if defined OMP || defined BLAS_OMP
#pragma omp critical
#endif
		{
			tree= io->tree_s[io->threads++];
		}
	  	  t_node* r = tree->noeud[tree->mod->startnode];
	  	  Init_Class_Tips(tree,precon);
	  	  int pars1 = Fill_Sankoff(r,tree,1); //fill Sankoff matrixes
	  	  //printf("\n\n. %d Initial maximum parsimony score: %d",i,pars1);
	  	  Set_Pars_Counters(r,tree,1); //set initial parsimony counters
	  	  if(tree->mod->optDebug)printf("\n. Resolving polytomies using isotype information");
	  	  int pars2 = Resolve_Polytomies_Pars(tree,0.001);
	  	  #if defined OMP || defined BLAS_OMP
		  #pragma omp critical
		  #endif
	  	  {
	  		  printf("\n. %d Initial/resolved maximum parsimony score: %d %d %s",tree->mod->num,pars1,pars2,tree->mod->rootname);
	  	  }
	  }
}

/*********************************************************/
// TODO: Edit to output and then free each tree at a time
void Pars_Reconstructions(option* io){
	int j,i;
  	int maxtrees=io->maxparstrees;
  	int maxotu=io->maxparsotu;
  	int sample=io->parssample;

  	//Set up stats file
  	char foutp[T_MAX_FILE];
  	strcpy(foutp,io->mod->in_align_file);
  	strcat(foutp,"_igphyml_parstats");
  	if(io->append_run_ID){
  		strcat(foutp, "_");
  		strcat(foutp, io->run_id_string);
  	}
  	strcat(foutp,".txt");
  	FILE* pstatf = Openfile(foutp,1);

  	For(j,io->ntrees){
  		 //set up basic tree stuff
  		 //read in final tree topology sent to outfile
  		 //set up root placement, etc
  		 char foutt[T_MAX_FILE];
		 strcpy(foutt,io->datafs[j]);
		 if(io->mod->ASR)strcat(foutt,"_igphyml_figtree");
		 else strcat(foutt,"_igphyml_tree");
		 if(io->append_run_ID){
			 strcat(foutt, "_");
			 strcat(foutt, io->run_id_string);
		 }
		 strcat(foutt,".txt");
  		  io->mod_s[j]->fp_in_tree = Openfile(foutt,0);
  		  t_tree* tree = Read_User_Tree(io->tree_s[j]->data,io->mod_s[j],io);
  		  model* mod = io->mod_s[j];
  		  tree->mod=mod;
  		  tree->io=io;
  		  tree->data = io->tree_s[j]->data;
  		  io->tree_s[j] = tree;
  		  mod->startnode = -1;
  		  int nodepos;
  		  For(nodepos,((tree->n_otu-1)*2)){
  			 if(strcmp(tree->noeud[nodepos]->name,mod->rootname)==0){
  			      mod->startnode=nodepos;
  			      Update_Ancestors_Edge(tree->noeud[nodepos],tree->noeud[nodepos]->v[0],tree->noeud[nodepos]->b[0],tree);
  			  }
  		  }
  		  if(mod->startnode==-1){
  			 PhyML_Printf("\n\nRoot sequence ID not found in data file! %s %s\n",mod->rootname,mod->in_align_file);
  			 exit(EXIT_FAILURE);
  		  }
  		  //setup joint trees file
  		  char fout[T_MAX_FILE];
  		  strcpy(fout,io->datafs[j]);
  		  strcat(fout,"_igphyml_jointpars");
  		  if(io->append_run_ID){
  		  	 strcat(fout, "_");
  		   	 strcat(fout, io->run_id_string);
  		  }
  		  strcat(fout,".nex");

  		  //array of all possible reconstructions
  		  t_tree** trees = mCalloc(maxtrees,sizeof(t_tree*));
  		  trees[0]=tree;
  		  t_node* r = tree->noeud[tree->mod->startnode]; //root node
  		  Init_Class_Tips(tree,io->precon); //initialize tip states and data structures
  		  int pars = Fill_Sankoff(r,tree,1); //max parsimony score
  		  Set_Pars_Counters(r,tree,1); //set counters to their first minimum position

  		  //Data structures for counting the number of type switches and branch lengths in each state
  	      phydbl* switches = mCalloc(tree->nstate*tree->nstate,sizeof(phydbl));
  	      phydbl* classl = mCalloc(tree->nstate,sizeof(phydbl));

  		  int npars = maxtrees;
  		  int minrootsar[tree->nstate];
  		  int minroots=0;
  		  int minroot=-1;
  		  For(i,tree->nstate){
  			if(r->sroot[i] == pars){
  				minrootsar[minroots] = i;
  				minroots++;
  				minroot=i;
  			}
  		  }
  		  if(tree->n_otu < maxotu){ //if too many taxa, don't bother trying to solve for all reconstructions
  			  npars=0;
  			  For(i,tree->nstate){
  				if(r->sroot[i] == pars){ //recursively solve for all maximum parsimony paths
  					//t_tree* tree2; //read in copy of that tree
					t_tree* tree2 = Read_User_Tree(tree->io->tree_s[j]->data,tree->io->mod_s[j],tree->io);
					tree2->mod=tree->mod;
					tree2->mod->startnode = tree->mod->startnode;
					t_node* r2 = tree2->noeud[tree->mod->startnode];
					tree2->io=tree->io;
					Update_Ancestors_Edge(r2,r2->v[0],r2->b[0],tree);
					Copy_Sankoff_Tree(tree,tree2);
					trees[npars]=tree2;
  					npars = Get_All_Paths(r2,i,tree2,trees,1,maxtrees,npars,j,i)+1;
  				}
  				if(npars >= maxtrees)break; //unless you get too many trees, then just sample
  			  }
  		  	  //printf("\n. %d Found %d maximum parsimony trees\n",j,npars);
  		  	  if(npars < maxtrees){ //if not too many trees found, record stats
  		  	  	  FILE* treeout = Openfile(fout, 1 );
  		  	  	  fprintf(treeout,"#NEXUS\nBegin taxa;\nDimensions ntax=%d;\nTaxlabels\n",tree->n_otu);
  		  	  	  For(i,(tree->n_otu-1)*2){
  		  	  		  if(tree->noeud[i]->tax){
  		  	  			  fprintf(treeout,"%s\n",tree->noeud[i]->name);
  		  	  		  }
  		  	  	  }
  		  	  	  fprintf(treeout,";\nEnd\nBegin trees;\n");
  		  	  	  For(i,npars){
  		  	  		  //printTreeState(trees[i]->noeud[trees[i]->mod->startnode],trees[i],1);printf(" 0\n");
  		  	  		  fprintf(treeout,"Tree TREE%d = [&R] ",i+1);
  		  	  		  Fill_Pars_Stats(trees[i]->noeud[trees[i]->mod->startnode],trees[i],switches,classl,1); //summarize statistics of the reconstruction
 		    		  io->precon *= 10;
  		  	  		  char* ts = Write_Tree(trees[i]);
  		    		  io->precon /= 10;
  		  	  		  ts[strlen(ts)-1] = 0;
  		  	  		  strcat(ts,"[");
  		  	  		  strcat(ts,trees[i]->chars[trees[i]->noeud[tree->mod->startnode]->pstate]);
  		  	  		  strcat(ts,"]:0;");
  		  	  		  fprintf(treeout,"%s\n",ts);
  		  	  		  if(i > 0){
  		  	  			  Clean_Tree(trees[i]);
  		  	  			  Free_Tree(trees[i]);
  		  	  		  }
  		  	  	  }
  		  	  	  fprintf(treeout,"END;\n");
  		  	  	  fclose(treeout);
  		  	  }
  		  }else{
  			printf("\n. Too many sequences for exhaustive parsimony search!");
  		  }
  		  if(npars>=maxtrees){ //if too many trees found, sample!
  			  npars = sample;
  			  /*if(minroots != 1){
  				  printf("NEED TO FIX RAND PATH TO WORK WITH AMBIGUOUS ROOTS %d %d\n",minroots,minroot);
  			  	  exit(EXIT_FAILURE);
  			  }*/
  	  		  FILE* treeout1 = Openfile(fout, 1 );
  			  printf("\n. Sampling %d trees instead for tree %d %s.",sample,j,tree->mod->rootname);
  			  for(i=1;i<maxtrees;i++){ //free trees from attempt to solve for all
  				  if(tree->n_otu < maxotu){
  					  Clean_Tree(trees[i]);
  					  Free_Tree(trees[i]);
  				  }
  			  }
  			  fprintf(treeout1,"#NEXUS\nBegin taxa;\nDimensions ntax=%d;\nTaxlabels\n",tree->n_otu);
  		  	  For(i,(tree->n_otu-1)*2){
  		  	  	  if(tree->noeud[i]->tax){
  		  	  		  fprintf(treeout1,"%s\n",tree->noeud[i]->name);
  		  	  	  }
  		  	  }
  		  	  fprintf(treeout1,";\nEnd\nBegin trees;\n");
  			  For(i,sample){
	  	  		  fprintf(treeout1,"Tree TREE%d = [&R] ",i+1);
	  	  		  int n = rand() % minroots;
	  	  		  int minroot = minrootsar[n];
	  	  		  //printf("root %d %d\n",minroot,minroots);
  			  	  Get_Rand_Path(r,minroot,tree,1);
  			  	  Fill_Pars_Stats(r,tree, switches,classl,1);
  			  	  io->precon *= 10;
  			  	  char* ts = Write_Tree(tree);
  			  	  io->precon /= 10;
  			  	  ts[strlen(ts)-1] = 0;
  			  	  strcat(ts,"[");
  			  	  strcat(ts,tree->chars[tree->noeud[tree->mod->startnode]->pstate]);
  			  	  strcat(ts,"]:0;");
  			  	  fprintf(treeout1,"%s\n",ts);
  			  	  free(ts);
  			  }
  			  fprintf(treeout1,"END;\n");
  			  fclose(treeout1);
  		  }
  		 int tposi, tposj;
  		 For(tposi,tree->nstate){
  		 	classl[tposi] = classl[tposi]/(npars*1.0);
  		 	For(tposj,tree->nstate){
  		 		//printf("%lf\t%d\n",switches[tposi*tree->nstate+tposj],npars);
  		 		switches[tposi*tree->nstate+tposj]=switches[tposi*tree->nstate+tposj]/(npars*1.0);
  		 	}
  		 }
  		 phydbl tswitch, tlen;
  		 tswitch=tlen=0;
  		 //printf("\nMINROOTS: %d",minroots);
  		 For(tposi,tree->nstate){
  		 	fprintf(pstatf,"%d\t%s\tN\t%lf\n",j,tree->chars[tposi],classl[tposi]);
  		 	tlen+=classl[tposi];
  		 }
  		 fprintf(pstatf,"%d\tUCA\t%s",j,tree->chars[minrootsar[0]]);
  		 for(i = 1; i < minroots; i++)fprintf(pstatf,":%s",tree->chars[minrootsar[i]]);
  		 fprintf(pstatf,"\t0.0\n");
  		 	For(tposi,tree->nstate){
  		 		For(tposj,trees[0]->nstate){
  		 			fprintf(pstatf,"%d\t%s\t%s\t%lf\n",j,tree->chars[tposi],tree->chars[tposj],switches[tposi*tree->nstate+tposj]);
  		 				tswitch+=switches[tposi*tree->nstate+tposj];
  		 			}
  		 		}
  		 // printf("\n. Switches: %lf. Length: %lf\n",tswitch,tlen);
  		  fclose(io->mod_s[j]->fp_in_tree);
  		  Clean_Tree(tree);
  		  Free_Tree(tree);
  		  free(trees);
  	}
}

/*********************************************************
 * Initialize Sankoff dynamic programming tables at the tips of the tree
 * -6, -5 : GC/Mem model.
 * 1- 4: Isotype models
 * -1 - -4: Unconstrained isotype models
 */
void Init_Class_Tips(t_tree* tree, int precon){
	 int i,j;
	 char* mtemp;
	 if(precon == 7){
		 Setup_Custom_Pars_Model(tree);
	 }else{
	  if(precon <= -5){ //create model
		 tree->nstate=6;
		 tree->chars = mCalloc(tree->nstate,sizeof(char*));
		 For(i,tree->nstate)tree->chars[i]=mCalloc(10,sizeof(char));
		 strcpy(tree->chars[0],"Naive");
		 strcpy(tree->chars[1],"GC");
		 strcpy(tree->chars[2],"UnMem");
		 strcpy(tree->chars[3],"MemHi");
		 strcpy(tree->chars[4],"MemLo");
		 strcpy(tree->chars[5],"Bmem");
	 }else if(precon < 3 && precon > -3){
		 tree->nstate=7;
		 tree->chars = mCalloc(tree->nstate,sizeof(char*));
		 For(i,tree->nstate)tree->chars[i]=mCalloc(10,sizeof(char));
		 strcpy(tree->chars[0],"M");
		 strcpy(tree->chars[1],"D");
		 strcpy(tree->chars[2],"G31");
		 strcpy(tree->chars[3],"A1");
		 strcpy(tree->chars[4],"G24");
		 strcpy(tree->chars[5],"E");
		 strcpy(tree->chars[6],"A2");
	 }else{
		 tree->nstate=9;
		 tree->chars = mCalloc(tree->nstate,sizeof(char*));
		 For(i,tree->nstate)tree->chars[i]=mCalloc(10,sizeof(char));
		 strcpy(tree->chars[0],"M");
		 strcpy(tree->chars[1],"D");
		 strcpy(tree->chars[2],"G3");
		 strcpy(tree->chars[3],"G1");
		 strcpy(tree->chars[4],"A1");
		 strcpy(tree->chars[5],"G2");
		 strcpy(tree->chars[6],"G4");
		 strcpy(tree->chars[7],"E");
		 strcpy(tree->chars[8],"A2");
	 }
	 For(i,(tree->n_otu-1)*2){ //initialize data structures for each node
	 	tree->noeud[i]->pl = (int *)mCalloc(tree->nstate*tree->nstate,sizeof(int ));//all possible pointer scores, left
	 	tree->noeud[i]->pr = (int *)mCalloc(tree->nstate*tree->nstate,sizeof(int ));//all possible pointer scores, right
	 	tree->noeud[i]->s = (int *)mCalloc(tree->nstate,sizeof(int)); //minimum score of subtrees given each state
	 	tree->noeud[i]->sroot = (int *)mCalloc(tree->nstate,sizeof(int)); //s but for root
	 	tree->noeud[i]->lmin = (int *)mCalloc(tree->nstate,sizeof(int)); //minimum score on left given state at node
	 	tree->noeud[i]->rmin = (int *)mCalloc(tree->nstate,sizeof(int)); //minimum score on right given state at node
	 	tree->noeud[i]->prc = (int *)mCalloc(tree->nstate,sizeof(int));//pointer to the left given state at the current node
	 	tree->noeud[i]->plc = (int *)mCalloc(tree->nstate,sizeof(int));//pointer to the right given state at the current node
	 	tree->noeud[i]->llock = (int *)mCalloc(tree->nstate,sizeof(int));//whether or not pointers on l or r are locked
	 	tree->noeud[i]->rlock = (int *)mCalloc(tree->nstate,sizeof(int));
	 	For(j,tree->nstate){
	 		tree->noeud[i]->s[j]=1000;
	 		tree->noeud[i]->lmin[j]=MAX_PARS;
	 		tree->noeud[i]->rmin[j]=MAX_PARS;
	 	}
	 }
	 tree->step_mat = mCalloc(tree->nstate*tree->nstate,sizeof(int));
	 For(i,tree->nstate){ //set up step mat
	 	 For(j,tree->nstate){
	 			 if(i < j && precon>0)tree->step_mat[j*tree->nstate+i]=10000; //if precon is negative, no costraint
	 			 else if(i == j)tree->step_mat[j*tree->nstate+i]=0; //no same state, no penalty
	 			 else tree->step_mat[j*tree->nstate+i] =1; //otherwise penalty of 1
	 	 }
	 }
	 For(i,tree->n_otu){//read in information from the ends of the sequence names
		 int nelements=0;
		 char* state = mCalloc(T_MAX_OPTION,sizeof(char));
		 char* minfo1 = strdup(tree->noeud[i]->name);
		 char* minfo2 = strdup(tree->noeud[i]->name);
		 while ((mtemp = strsep(&minfo1, "_")) != NULL){nelements++;}
		 For(j,nelements){
			 strcpy(state,strsep(&minfo2, "_"));
		 }
		 if(tree->mod->optDebug)printf("\n%s\t%s",tree->noeud[i]->name,state);
		  if(precon <= -5){
		 	 if(strcmp(state,"Naive")==0||strcmp(state,"GERM")==0)tree->noeud[i]->s[0]=0;
		 	 if(strcmp(state,"GC")==0)   tree->noeud[i]->s[1]=0;
		 	 if(strcmp(state,"UnMem")==0)tree->noeud[i]->s[2]=0;
		 	 if(strcmp(state,"MemHi")==0)tree->noeud[i]->s[3]=0;
		 	 if(strcmp(state,"MemLo")==0)tree->noeud[i]->s[4]=0;
		 	 if(strcmp(state,"Bmem")==0)tree->noeud[i]->s[5]=0;
		 }
		 if(precon < 3 && precon > -3){
		 	 if(strcmp(state,"M")==0)tree->noeud[i]->s[0]=0;
		 	 if(strcmp(state,"D")==0)tree->noeud[i]->s[1]=0;
		 	 if(strcmp(state,"G3")==0||strcmp(state,"G1")==0)tree->noeud[i]->s[2]=0;
		 	 if(strcmp(state,"A1")==0)tree->noeud[i]->s[3]=0;
		 	 if(strcmp(state,"G2")==0||strcmp(state,"G4")==0)tree->noeud[i]->s[4]=0;
		 	 if(strcmp(state,"E")==0)tree->noeud[i]->s[5]=0;
		 	 if(strcmp(state,"A2")==0)tree->noeud[i]->s[6]=0;
		 	 if(strcmp(state,"G")==0){
		 		 tree->noeud[i]->s[2]=0;
		 		 tree->noeud[i]->s[4]=0;
		 	 }
		 }else{
			 if(strcmp(state,"M")==0)tree->noeud[i]->s[0]=0;
			 if(strcmp(state,"D")==0)tree->noeud[i]->s[1]=0;
			 if(strcmp(state,"G3")==0)tree->noeud[i]->s[2]=0;
			 if(strcmp(state,"G1")==0)tree->noeud[i]->s[3]=0;
			 if(strcmp(state,"A1")==0)tree->noeud[i]->s[4]=0;
			 if(strcmp(state,"G2")==0)tree->noeud[i]->s[5]=0;
			 if(strcmp(state,"G4")==0)tree->noeud[i]->s[6]=0;
			 if(strcmp(state,"E")==0)tree->noeud[i]->s[7]=0;
			 if(strcmp(state,"A2")==0)tree->noeud[i]->s[8]=0;
			 if(strcmp(state,"G")==0){
			 	 tree->noeud[i]->s[2]=0;
			 	 tree->noeud[i]->s[3]=0;
			 	 tree->noeud[i]->s[5]=0;
			 	 tree->noeud[i]->s[6]=0;
			  }
		 }
		 For(j,tree->nstate){
			 if(tree->mod->optDebug)printf(" %d",tree->noeud[i]->s[j]);
		 }
		 if(tree->mod->optDebug)printf("\n");
	 }
   }
}

/*********************************************************
* Permute tip names except for germline
*/
void Permute_Tips(t_tree* tree){
	int i, j;
	int indexes[tree->n_otu-1];

	int count = 0;
	//For(i,tree->n_otu)printf("1 %d\t%s\n",i,tree->noeud[i]->name);
	int n = tree->n_otu;
	char temp[T_MAX_NAME];
    For(i,n - 1) {
    		if(strcmp(tree->noeud[i]->name,tree->mod->rootname) != 0){
    			j = i + rand() / (RAND_MAX / (n - i) + 1);
    			while(strcmp(tree->noeud[j]->name,tree->mod->rootname) == 0)
    				j = i + rand() / (RAND_MAX / (n - i) + 1);
    			strcpy(temp,tree->noeud[i]->name);
    			strcpy(tree->noeud[i]->name,tree->noeud[j]->name);
    			strcpy(tree->noeud[j]->name,temp);
    		}
    }
	//For(i,tree->n_otu)printf("2 %d\t%s\n",i,tree->noeud[i]->name);

}

/*********************************************************
* Permute tip names except for germline
*/
void Permute_MetaData(t_tree* tree, int pos){
	int i, j;
	int indexes[tree->n_otu-1];

	//permute last element of tips and re-assign to sequence IDs
	int count = 0;
	//For(i,tree->n_otu)printf("1 %d\t%s\n",i,tree->noeud[i]->name);
	int n = tree->n_otu;
	char temp[T_MAX_NAME];
    For(i,n - 1) {
    		if(strcmp(tree->noeud[i]->name,tree->mod->rootname) != 0){
    			j = i + rand() / (RAND_MAX / (n - i) + 1);
    			while(strcmp(tree->noeud[j]->name,tree->mod->rootname) == 0)
    				j = i + rand() / (RAND_MAX / (n - i) + 1);
    			strcpy(temp,tree->noeud[i]->name);
    			strcpy(tree->noeud[i]->name,tree->noeud[j]->name);
    			strcpy(tree->noeud[j]->name,temp);
    		}
    }
	//For(i,tree->n_otu)printf("2 %d\t%s\n",i,tree->noeud[i]->name);

}


/*********************************************************
* read in parimsony model
*/
void Setup_Custom_Pars_Model(t_tree* tree){
	//printf("Opening %s\n",tree->mod->preconfile);
	int i,j;
	char* mtemp;
	FILE* PARS = Openfile(tree->mod->preconfile,0);
	char* line = mCalloc(T_MAX_LINE,sizeof(char));
	int fscn;
	do{//skip over beginning comments
		fscn = fscanf(PARS, "%s\n",line);
		//printf("LINE %s\n",line);
	}while(strcmp(line,"#BEGIN")!=0 && fscn != EOF);

	//read in states
	fscn = fscanf(PARS, "%d\n",&tree->nstate);
	 tree->chars = mCalloc(tree->nstate,sizeof(char*));

	 //associative arrays holding all possible states mapped back to the states of the model
	 //ambigfrom: nstate x max option
	 char** ambigstatesfrom = mCalloc(tree->nstate*T_MAX_OPTION,sizeof(char*));
	 For(i,tree->nstate*T_MAX_OPTION)ambigstatesfrom[i]=mCalloc(T_MAX_OPTION,sizeof(char*));
	 int* ambigstatesto = mCalloc(tree->nstate*T_MAX_OPTION,sizeof(int));
	 int statecount = 0;//total number of states allowed in data
	 fscn = fscanf(PARS, "%s\n",line);
	 For(i,tree->nstate){
		 tree->chars[i]=mCalloc(T_MAX_OPTION,sizeof(char));
		 fscn=fscanf(PARS, "%s\n",tree->chars[i]);
		 strcpy(ambigstatesfrom[i],tree->chars[i]);
		 ambigstatesto[i]=i;
		 statecount++;
		 //printf("state %d %s\n",i,tree->chars[i]);
	 }
	 For(i,(tree->n_otu-1)*2){ //initialize data structures for each node
	 	tree->noeud[i]->pl = (int *)mCalloc(tree->nstate*tree->nstate,sizeof(int ));//all possible pointer scores, left
	 	tree->noeud[i]->pr = (int *)mCalloc(tree->nstate*tree->nstate,sizeof(int ));//all possible pointer scores, right
	 	tree->noeud[i]->s = (int *)mCalloc(tree->nstate,sizeof(int)); //minimum score of subtrees given each state
	 	tree->noeud[i]->sroot = (int *)mCalloc(tree->nstate,sizeof(int)); //s but for root
	 	tree->noeud[i]->lmin = (int *)mCalloc(tree->nstate,sizeof(int)); //minimum score on left given state at node
	 	tree->noeud[i]->rmin = (int *)mCalloc(tree->nstate,sizeof(int)); //minimum score on right given state at node
	 	tree->noeud[i]->prc = (int *)mCalloc(tree->nstate,sizeof(int));//pointer to the left given state at the current node
	 	tree->noeud[i]->plc = (int *)mCalloc(tree->nstate,sizeof(int));//pointer to the right given state at the current node
	 	tree->noeud[i]->llock = (int *)mCalloc(tree->nstate,sizeof(int));//whether or not pointers on l or r are locked
	 	tree->noeud[i]->rlock = (int *)mCalloc(tree->nstate,sizeof(int));
	 	For(j,tree->nstate){
	 		tree->noeud[i]->s[j]=1000;
	 		tree->noeud[i]->lmin[j]=MAX_PARS;
	 		tree->noeud[i]->rmin[j]=MAX_PARS;
	 	}
	 }
	 tree->step_mat = mCalloc(tree->nstate*tree->nstate,sizeof(int));
	 For(i,tree->nstate){ //set up step mat
	 	 For(j,tree->nstate){
	 		if(i == j)tree->step_mat[i*tree->nstate+j]=0; //no same state, no penalty
	 		else tree->step_mat[i*tree->nstate+j] =1; //otherwise penalty of 1
	 	 }
	 }

	 //read in step constraints
	 fscn = fscanf(PARS, "%s\n",line);
	 char* from = mCalloc(T_MAX_OPTION,sizeof(char));
	 char* to = mCalloc(T_MAX_OPTION,sizeof(char));
	 int val=MAX_PARS;
	 if(fscn != EOF && strcmp(line,"#CONSTRAINTS")==0){
		 do{
			 fscn = fscanf(PARS, "%s %s %d\n",from,to,&val);
			 //printf("LINE %s %s %d\n",from,to,val);
			 if(strcmp(from,"#AMBIGUOUS")!=0){
			 For(i,tree->nstate){ //set up step mat
				 For(j,tree->nstate){
					 if(strcmp(tree->chars[i],from)==0 && strcmp(tree->chars[j],to)==0){
						 //printf("%s\t%s\t%d\n",tree->chars[i],tree->chars[j],tree->step_mat[i*tree->nstate+j]);
						 tree->step_mat[i*tree->nstate+j]=val;
					 }
				 }
			 }
			 }
		 }while(strcmp(from,"#AMBIGUOUS")!=0 && fscn != EOF);
		 rewind(PARS);
	 }
	  For(i,tree->nstate){ //set up step mat
	 	 For(j,tree->nstate){
	 		 //printf("%s\t%s\t%d\n",tree->chars[i],tree->chars[j],tree->step_mat[i*tree->nstate+j]);
	 	 }
	 }
	  //read in ambiguous states
	 if(fscn != EOF){
		 do{
			 fscn = fscanf(PARS,"%s\n",from);
		 }while(strcmp(from,"#AMBIGUOUS")!=0 && fscn != EOF);
		 if(fscn != EOF && strcmp(from,"#AMBIGUOUS")==0){
			 do{
				 fscn = fscanf(PARS, "%s %s\n",from,to);
				 if(fscn==-1)break;
				 //printf("AMB %s %s\n",from,to);
				 strcpy(ambigstatesfrom[statecount],from);
				 For(j,tree->nstate){
					 if(strcmp(tree->chars[j],to)==0){
						 ambigstatesto[statecount]=j;
					 }
				 }
				 statecount++;
			 }while(fscn != EOF);
		 }
	 }
	 For(j,statecount){
		 //printf("state %d\t%s\t%d\n",j,ambigstatesfrom[j],ambigstatesto[j]);
	 }
	 For(i,tree->n_otu){//read in information from the ends of the sequence names
		 int nelements=0;
		 char* state = mCalloc(T_MAX_OPTION,sizeof(char));
		 char* minfo1 = strdup(tree->noeud[i]->name);
		 char* minfo2 = strdup(tree->noeud[i]->name);
		 while ((mtemp = strsep(&minfo1, "_")) != NULL){nelements++;}
		 For(j,nelements-tree->mod->mdpos){
			 strcpy(state,strsep(&minfo2, "_"));
		 }
		 if(tree->mod->optDebug)printf("\n%s\t%s",tree->noeud[i]->name,state);
		 //printf("\n%s\t%s",tree->noeud[i]->name,state);
		 int found = 0;
		 For(j,statecount){
			 if(strcmp(state,ambigstatesfrom[j])==0){
				 tree->noeud[i]->s[ambigstatesto[j]]=0;
				 found++;
			 }
		 }
		 if(!found){
			 printf("\nState %s from sequence %s not found in model.\n",state,tree->noeud[i]->name);
			 Warn_And_Exit("");
		 }
		 For(j,tree->nstate){
			 if(tree->mod->optDebug)printf(" %d",tree->noeud[i]->s[j]);
			 //printf(" %d",tree->noeud[i]->s[j]);
		 }
		 if(tree->mod->optDebug)printf("\n");
		 free(state);
		 //printf("\n");
	 }
	 free(from);
	 free(to);
	 free(line);
	 free(ambigstatesto);
	 For(i,tree->nstate*T_MAX_OPTION)free(ambigstatesfrom[i]);
	 free(ambigstatesfrom);
	 fclose(PARS);
}


/*********************************************************
* Free extraneous data structures
*/
void Clean_Tree(t_tree* tree){
	int i,j;
    For(i,(tree->n_otu-1)*2){
		 free(tree->noeud[i]->pl);
		 free(tree->noeud[i]->pr);
		 free(tree->noeud[i]->s);
		 free(tree->noeud[i]->sroot);
		 free(tree->noeud[i]->lmin);
		 free(tree->noeud[i]->rmin);
		 free(tree->noeud[i]->prc);
		 free(tree->noeud[i]->plc);
		 free(tree->noeud[i]->llock);
		 free(tree->noeud[i]->rlock);
	 }
}


/*********************************************************
* Copy all tree structures and values to a new trees
*/
void Copy_Sankoff_Tree(t_tree* tree1,t_tree* tree2){
	int i,j,k;
	tree2->nstate=tree1->nstate;
    For(i,(tree1->n_otu-1)*2){
    	tree2->noeud[i]->pl = (int *)mCalloc(tree2->nstate*tree2->nstate,sizeof(int ));
    	tree2->noeud[i]->pr = (int *)mCalloc(tree2->nstate*tree2->nstate,sizeof(int ));
    	tree2->noeud[i]->s = (int *)mCalloc(tree2->nstate,sizeof(int));
    	tree2->noeud[i]->sroot = (int *)mCalloc(tree2->nstate,sizeof(int));
    	tree2->noeud[i]->lmin = (int *)mCalloc(tree2->nstate,sizeof(int));
    	tree2->noeud[i]->rmin = (int *)mCalloc(tree2->nstate,sizeof(int));
    	tree2->noeud[i]->prc = (int *)mCalloc(tree2->nstate,sizeof(int));
    	tree2->noeud[i]->plc = (int *)mCalloc(tree2->nstate,sizeof(int));
    	tree2->noeud[i]->llock = (int *)mCalloc(tree2->nstate,sizeof(int));
    	tree2->noeud[i]->rlock = (int *)mCalloc(tree2->nstate,sizeof(int));
    	For(j,tree2->nstate){
    		 tree2->noeud[i]->s[j]=tree1->noeud[i]->s[j];
    		 tree2->noeud[i]->plc[j]=tree1->noeud[i]->plc[j];
    		 tree2->noeud[i]->prc[j]=tree1->noeud[i]->prc[j];
    		 tree2->noeud[i]->lmin[j]=tree1->noeud[i]->lmin[j];
    		 tree2->noeud[i]->rmin[j]=tree1->noeud[i]->rmin[j];
    		 tree2->noeud[i]->llock[j]=tree1->noeud[i]->llock[j];
    		 tree2->noeud[i]->rlock[j]=tree1->noeud[i]->rlock[j];
    		 tree2->noeud[i]->sroot[j]=tree1->noeud[i]->sroot[j];
    		 For(k,tree2->nstate){
    			 tree2->noeud[i]->pl[j*tree2->nstate+k]=tree1->noeud[i]->pl[j*tree2->nstate+k];
    			 tree2->noeud[i]->pr[j*tree2->nstate+k]=tree1->noeud[i]->pr[j*tree2->nstate+k];
    		 }
    	}
	 }
    tree2->step_mat=tree1->step_mat;
    tree2->chars=tree1->chars;
}

/*********************************************************/
// Resolve polytomies based on parsimony score
int Resolve_Polytomies_Pars(t_tree* tree, phydbl thresh){
	int i;
	int nni=1;
	int nnifound=1; //initial parsimony score
	int pars0 = Fill_Sankoff(tree->noeud[tree->mod->startnode],tree,1);
	int iters = 0;
	while(nnifound){ //execute if a rearrangement is found in a tree
		nnifound=0;
		For(i,(tree->n_otu-1)*2){ //if no taxa on either side of edge and length is below threshold, search for NNIs
			if(!tree->noeud[i]->tax && !tree->noeud[i]->anc->tax && tree->noeud[i]->anc_edge->l < thresh){
				pars0 = Fill_Sankoff(tree->noeud[tree->mod->startnode],tree,1);
				t_edge* b = tree->noeud[i]->anc_edge;
				nni = NNI_Pars_Search(tree->noeud[i],tree->noeud[i]->anc,tree->noeud[i]->anc_edge,tree->noeud[i]->anc_edge,pars0,thresh,tree);
				if(!nni)nni = NNI_Pars_Search(tree->noeud[i]->anc,tree->noeud[i],tree->noeud[i]->anc_edge,tree->noeud[i]->anc_edge,pars0,thresh,tree);
				nnifound+=nni;
			}
		}
		iters++;
	}
	int j;
	/*t_node* a;
	t_node* b;
	t_node* c;
	t_node* d;
	t_node* b10;
	t_node* a11;
	t_node* d2;
	For(i,(tree->n_otu-1)*2){
		printf("%d",tree->noeud[i]->num);
		if(tree->noeud[i]->num != tree->mod->startnode)printf("\t%d\n",tree->noeud[i]->anc->num);
		if(tree->noeud[i]->tax)printf("\t%s\n",tree->noeud[i]->name);
		if(tree->noeud[i]->num==9)a=tree->noeud[i];
		if(tree->noeud[i]->num==8)b=tree->noeud[i];
		if(tree->noeud[i]->num==7)c=tree->noeud[i];
		if(tree->noeud[i]->num==0)d=tree->noeud[i];
		if(tree->noeud[i]->num==10)b10=tree->noeud[i];
		if(tree->noeud[i]->num==11)a11=tree->noeud[i];
		if(tree->noeud[i]->num==2)d2=tree->noeud[i];
	}
	printf("resolved tree: %s\n",Write_Tree(tree));
	printf("%lf\n",Score_Polytomy(tree->noeud[11],a,thresh,1,1000,tree));
	Swap(a,b,c,d,tree);
	printf("swapped tree: %s\n",Write_Tree(tree));
	pars0 = Fill_Sankoff(tree->noeud[tree->mod->startnode],tree,1);
	phydbl paths = Score_Polytomy(tree->noeud[11],a,thresh,1,1000,tree);
	printf("PARS: %d PATHS: %lf\n",pars0,paths);
	For(i,(tree->n_otu-1)*2){
		printf("%d",tree->noeud[i]->num);
		if(tree->noeud[i]->num != tree->mod->startnode)printf("\t%d\n",tree->noeud[i]->anc->num);
		if(tree->noeud[i]->tax)printf("\t%s\n",tree->noeud[i]->name);
	}
	Swap(a11,b10,a,d2,tree);
	pars0 = Fill_Sankoff(tree->noeud[tree->mod->startnode],tree,1);
	paths = Score_Polytomy(tree->noeud[11],a,thresh,1,1000,tree);
	printf("PARS: %d PATHS: %lf\n",pars0,paths);
	printf("swapped tree: %s\n",Write_Tree(tree));
	/*exit(1);*/
	nnifound=1; //initial parsimony score
	phydbl unresolved = 1.0;
	iters = 0;
	int maxtrees = 1000;
	t_tree** trees = mCalloc(maxtrees,sizeof(t_tree*));
	while(nnifound){ //execute if a rearrangement is found in a tree
		nnifound=0;
		unresolved=0;
		For(i,(tree->n_otu-1)*2){ //if no taxa on either side of edge and length is below threshold, search for NNIs
			t_node* node = tree->noeud[i];
			if(!tree->noeud[i]->tax && !tree->noeud[i]->anc->tax && tree->noeud[i]->anc_edge->l < thresh){
				pars0 = Fill_Sankoff(tree->noeud[tree->mod->startnode],tree,1);
				phydbl score = Score_Polytomy(node,node->anc,thresh,1,1000,tree);
				/*t_edge* b = node->anc_edge;
				t_node* top = node->anc;
				t_edge* pedge = node->anc_edge;
				while(top->anc->num != tree->mod->startnode && top->anc_edge->l < thresh){
					pedge = top->anc_edge;
					top = top->anc;
				}*/
				        nni = NNI_Paths_Search(node,node->anc,node->anc_edge,node->anc_edge,pars0,thresh,tree,maxtrees,trees,score,iters);
				if(!nni)nni = NNI_Paths_Search(node->anc,node,node->anc_edge,node->anc_edge,pars0,thresh,tree,maxtrees,trees,score,iters);
				//if(nni)printf("%d\t%d\n",iters,node->num);
				nnifound+=nni;
				phydbl resolve = Score_Polytomy(node,node->anc,thresh,0,1000,tree);
				unresolved += resolve;
			}
		}
						//break;
		printf("ITER %d %d %d %lf\n",iters,pars0,nnifound,unresolved);
		//printf("%s\n",Write_Tree(tree));
		if(unresolved < 0.01 && unresolved > -0.01)unresolved = 0;
		iters++;
	}
	return pars0;
}

/*********************************************************
 *
 */
phydbl Score_Polytomy(t_node* d, t_node* b, phydbl thresh, int sumscore, int maxtrees, t_tree* tree){
	int debug = 1;
	maxtrees = 1;
	if(debug)printf("NODE: %d\n",d->num);
	int i,j,k;
	t_node* top = d;
	/*if(d->anc->num == tree->mod->startnode){
		if(debug)printf("switch\n");
		top = b->anc;
	}*/
	int nedges = (tree->n_otu-1)*2;
	while(top->anc->num != tree->mod->startnode && top->anc_edge->l < thresh){
		//printf("trying\n");
		top = top->anc;
	}


	t_tree** btrees = mCalloc(maxtrees,sizeof(t_tree*));
	t_node* r = tree->noeud[tree->mod->startnode];
	//Init_Class_Tips(tree,tree->io->precon); //initialize tip states and data structures
  	int pars = Fill_Sankoff(r,tree,1); //max parsimony score
  	Set_Pars_Counters(r,tree,1);

	t_node* tr2 = tree->noeud[top->num];
	int taxinfo[nedges];
	For(i,nedges)taxinfo[i] = tree->noeud[i]->tax;
	For(i,3){
		if(tr2->anc->num != tr2->v[i]->num)
		Isolate_Polytomy(tr2->v[i],thresh,tree);
	}

	btrees[0]=tree;
	//paths = Get_All_Paths(tree->noeud[top->num], mins, tree, btrees, 0, maxtrees, paths, tree->mod->num, mins)+1;
	int paths = 0;
	int smin = INT_MAX;
	int mins = 0;
	For(i,tree->nstate){
		if(r->sroot[i] < smin){
			smin = r->sroot[i];
			mins = i;
		}
	}
	paths = 1;
	Get_First_Path(tree->noeud[tree->mod->startnode],mins,tree,1);

	tree->io->precon *= 10;
	if(debug)printf("%s\n",Write_Tree(tree));
  	tree->io->precon /= 10;
  	if(debug)printf("%s\n",Write_Tree(tree));

	phydbl mscore = 0;
	int nzero = 0;
	For(i,paths){
		int* scores = mCalloc(tree->nstate,sizeof(int));
		//if(d->num == 1251)debug=1;
		Score_Mono(btrees[i]->noeud[top->num],1,scores,debug,tree);
		For(j,tree->nstate){
			if(debug)printf("%s\t%d\t%d\n",tree->chars[j],scores[j],sumscore);
			if(mscore < scores[j] && !sumscore)mscore = scores[j];
			if(sumscore){
				if(scores[j] > 0)scores[j] *= scores[j];
				mscore += scores[j];
			}
			if(scores[j] > 0)nzero++;
		}
		free(scores);
	}

	phydbl best;
	int ninternal = nzero - 2;
	phydbl worst = 1.0;
	if(nzero > 1)worst = floor(ninternal/2.0) + (ninternal % 2) + 2;
	if(!sumscore){
		//printf("%d\t%lf\t%lf\n",nzero,mscore,worst);
		mscore -= worst;
	}else{
		best = 1 + (nzero-1)*worst;
	}

	For(i,nedges)tree->noeud[i]->tax = taxinfo[i];
	for(i=1;i<paths;i++){
		Clean_Tree(btrees[i]);
		Free_Tree(btrees[i]);
	}
	free(btrees);
	return(mscore);
}

/*********************************************************/

void Score_Mono(t_node* d, int level, int* scores, int debug, t_tree* tree){
	int i;
	if(debug)printf("polytomy\t%d\t%d\t%lf\t%s\n",d->anc->num,d->num,d->anc_edge->l,tree->chars[d->pstate]);
	if(scores[d->pstate] == 0){
		//printf("score %d\t%s\t%d\n",d->num,tree->chars[d->pstate],level);
		scores[d->pstate] = -level;
	}
	if(d->tax){
		if(scores[d->pstate] < 0)scores[d->pstate]=-scores[d->pstate];
	}
	if(!d->tax){
		For(i,3){
			if(d->v[i]->num != d->anc->num){
				Score_Mono(d->v[i],level+1,scores,debug,tree);
			}
		}
	}
}

/*********************************************************
 *
 */
void Isolate_Polytomy(t_node* d, phydbl thresh, t_tree* tree){
	int i;
	//printf("iso %d\n",d->num);
	if(d->tax && d->num != tree->mod->startnode)return;
	if(d->anc_edge->l < thresh){
		For(i,3){
			if(d->v[i]->num != d->anc->num){
				Isolate_Polytomy(d->v[i],thresh,tree);
			}
		}
	}else{
		//printf("%d is now tax\n",d->num);
		d->tax = 1;
	}
}

/*********************************************************
 * Starting from an intitial small length branch, check for NNI moves that would increase parsimony score, then recu333rsively search across
 * all adjacent branches (spreading c and d f) connected with at most thresh length
 * \b           /e
 *  \ c_fcus   /
 *   \c_...__ /d
 *   /  d_fcus\
 *  /          \
 * /a           \f
 *
 * d_fcus does not necessarily connect d and c. c_fcus is not necessarily d_fcus */
int NNI_Paths_Search(t_node *c, t_node *d,t_edge* c_fcus,t_edge* d_fcus, int pars0,phydbl thresh,
		t_tree* tree, int maxtrees,t_tree** trees,phydbl paths0, int iter){

	int dir1,dir2,dir3,dir4,i;
	dir1=dir2=dir3=dir4-1;
	For(i,3) if(d->b[i]->num != d_fcus->num) (dir1<0)?(dir1=i):(dir2=i);
	For(i,3) if(c->b[i]->num != c_fcus->num) (dir3<0)?(dir3=i):(dir4=i);
	t_node* e = d->v[dir1];
	t_node* f = d->v[dir2];
	t_edge* ee = d->b[dir1];
	t_edge* fe = d->b[dir2];

	t_node* a = c->v[dir3];
	t_node* b = c->v[dir4];


	int parsv[3];
	int parsmin[3];
	phydbl pathv[3];

	parsv[0] = pars0;
	parsv[1] = NNI_ParsSwaps(a,c,d,e,tree);
	parsv[2] = NNI_ParsSwaps(b,c,d,e,tree);
	parsmin[0] = 0;
	parsmin[1] = 0;
	parsmin[2] = 0;

	pathv[0] = paths0;
	pathv[1] = NNI_PathsSwaps(a,c,d,e,tree,thresh, maxtrees, trees);
	pathv[2] = NNI_PathsSwaps(b,c,d,e,tree,thresh, maxtrees, trees);

	int minpars = INT_MAX;
	phydbl maxpath = DBL_MAX;
	For(i,3)if(parsv[i] < minpars)minpars = parsv[i];
	For(i,3)if(parsv[i] == minpars)if(pathv[i] < maxpath)maxpath = pathv[i];

	if(parsv[0] == minpars && pathv[0] == maxpath){//Nothing!
		//printf("Not swapping!\n");
	}else if(parsv[1] == minpars && pathv[1] == maxpath){
	 // printf("0 Swapping!\n");
	  Swap(a,c,d,e,tree);
		  return 1;
	}else if(parsv[2] == minpars && pathv[2] == maxpath){
	 // printf("1 Swapping!\n");
	  Swap(b,c,d,e,tree);
		  return 1;
	}else{
		printf("SOMETHING WEIRD\n");
	}
	if(!e->tax){ //search along edge between d and e
		if(ee->l<thresh){
			int pt = NNI_Paths_Search(c,e,c_fcus,ee,pars0,thresh,tree,maxtrees,trees,paths0,iter);
			if(pt)return pt;
		}
	}
	if(!f->tax){ //search along edge between d and f
		if(fe->l<thresh){
			int pt = NNI_Paths_Search(c,f,c_fcus,fe,pars0,thresh,tree,maxtrees,trees,paths0,iter);
			if(pt)return pt;
		}
	}
	//exit(1);
	return 0;
}

/********************************************************
  *Swap specified nodes, swap back, and return parsimony score of the swap
  * \             /d      \             /a
  *  \           /         \           /
  *   \b__...__c/    ->     \b__...__c/
  *   /         \	   		 /		   \
  *  /           \	        /	        \
  * /a            \  	   /d            \
  *
  * nodes b and c are not necessarily on the same branch */
int NNI_PathsSwaps(t_node *a, t_node *b, t_node *c, t_node *d, t_tree *tree, phydbl thresh, int maxtrees, t_tree** trees){
	int l_r, r_l, l_v1, l_v2, r_v3, r_v4;
	int pars0,pars1,pars2;
	phydbl paths0,paths1,paths2;

	  pars0 = Fill_Sankoff(tree->noeud[tree->mod->startnode],tree,1);
	  paths0 = Score_Polytomy(b,c,thresh,1,maxtrees,tree);
	  if(tree->mod->optDebug){
	  	  Get_First_Path(tree->noeud[tree->mod->startnode],1,tree,1);
	  	  printTreeState(tree->noeud[tree->mod->startnode],tree,1);printf(" %lf ",paths0);printf(" 0\n");
	  }
	  //	  exit(1);
	  Swap(a,b,c,d,tree);
	  pars1 = Fill_Sankoff(tree->noeud[tree->mod->startnode],tree,1);
	  paths1 = Score_Polytomy(b,c,thresh,1,maxtrees,tree);
	  if(tree->mod->optDebug){
		  Get_First_Path(tree->noeud[tree->mod->startnode],1,tree,1);
		  printTreeState(tree->noeud[tree->mod->startnode],tree,1);printf(" %lf ",paths1);printf(" 1\n");
	  }
	  Swap(d,b,c,a,tree); //swap nodes
	  pars2 = Fill_Sankoff(tree->noeud[tree->mod->startnode],tree,1);
	  paths2 = Score_Polytomy(b,c,thresh,1,maxtrees,tree);
	  if(tree->mod->optDebug){
	  	  Get_First_Path(tree->noeud[tree->mod->startnode],1,tree,1);
	  	  printTreeState(tree->noeud[tree->mod->startnode],tree,1);printf(" %lf ",paths2);printf(" 2\n");
	  }
	  //	  exit(1);
	  //printf("pars score %d\t%d\t%d\t%lf\t%lf\t%lf\n",pars0,pars1,pars2,paths0,paths1,paths2);
	  if(pars0 != pars2 || paths0 != paths2){
		  printf("pars score %d\t%d\t%d\n",pars0,pars1,pars2);
		  printf("paths score %lf\t%lf\t%lf\n",paths0,paths1,paths2);
		  printf("\n.\tParsimony reconstruction swap inconsistent!\n");
		  exit(EXIT_FAILURE);
	  }
	  //exit(1);
	  return paths1;
}
/*********************************************************
 * Starting from an intitial small length branch, check for NNI moves that would increase parsimony score, then recursively search across
 * all adjacent branches (spreading c and d f) connected with at most thresh length
 * \b           /e
 *  \ c_fcus   /
 *   \c_...__ /d
 *   /  d_fcus\
 *  /          \
 * /a           \f
 *
 * d_fcus does not necessarily connect d and c. c_fcus is not necessarily d_fcus */
int NNI_Pars_Search(t_node *c, t_node *d,t_edge* c_fcus,t_edge* d_fcus, int pars0, phydbl thresh,t_tree* tree){

	int dir1,dir2,dir3,dir4,i;
	dir1=dir2=dir3=dir4-1;
	For(i,3) if(d->b[i]->num != d_fcus->num) (dir1<0)?(dir1=i):(dir2=i);
	For(i,3) if(c->b[i]->num != c_fcus->num) (dir3<0)?(dir3=i):(dir4=i);
	t_node* e = d->v[dir1];
	t_node* f = d->v[dir2];
	t_edge* ee = d->b[dir1];
	t_edge* fe = d->b[dir2];

	t_node* a = c->v[dir3];
	t_node* b = c->v[dir4];

	int pars1 = NNI_ParsSwaps(a,c,d,e,tree);
	int pars2 = NNI_ParsSwaps(b,c,d,e,tree);

	if(tree->mod->optDebug)printf("%d\t%d\t%d\n",pars0,pars1,pars2);
	if(pars0 <= MIN(pars1,pars2)){//Nothing!
		if(tree->mod->optDebug)printf("Not swapping!\n");
	}else if(pars1 < MIN(pars2,pars0)){
	  if(tree->mod->optDebug)printf("0 Swapping!\n");
	  Swap(a,c,d,e,tree);
	  return 1;
	}else if(pars2 < MIN(pars0,pars1)){
	  if(tree->mod->optDebug)printf("1 Swapping!\n");
	  Swap(b,c,d,e,tree);
	  return 1;
	}else if(pars1 == pars2){//if both options are equally better than the original, do the first
		if(tree->mod->optDebug)printf("2 Swapping!\n");
		Swap(a,c,d,e,tree);
		return 1;
	}else{
		printf("SOMETHING WEIRD\n");
	}
	if(!e->tax){ //search along edge between d and e
		if(ee->l<thresh){
			int pt = NNI_Pars_Search(c,e,c_fcus,ee,pars0,thresh,tree);
			if(pt)return pt;
		}
	}
	if(!f->tax){ //search along edge between d and f
		if(fe->l<thresh){
			int pt = NNI_Pars_Search(c,f,c_fcus,fe,pars0,thresh,tree);
			if(pt)return pt;
		}
	}
	return 0;
}

/********************************************************
  *Swap specified nodes, swap back, and return parsimony score of the swap
  * \             /d      \             /a
  *  \           /         \           /
  *   \b__...__c/    ->     \b__...__c/
  *   /         \	   		 /		   \
  *  /           \	        /	        \
  * /a            \  	   /d            \
  *
  * nodes b and c are not necessarily on the same branch */
int NNI_ParsSwaps(t_node *a, t_node *b, t_node *c, t_node *d, t_tree *tree){
	int l_r, r_l, l_v1, l_v2, r_v3, r_v4;
	int pars0,pars1,pars2;

	  pars0 = Fill_Sankoff(tree->noeud[tree->mod->startnode],tree,1);
	  if(tree->mod->optDebug){
	  	  Get_First_Path(tree->noeud[tree->mod->startnode],0,tree,1);
	  	  printTreeState(tree->noeud[tree->mod->startnode],tree,1);printf(" 0\n");
	  }
	  Swap(a,b,c,d,tree);
	  pars1 = Fill_Sankoff(tree->noeud[tree->mod->startnode],tree,1);
	  if(tree->mod->optDebug){
		  Get_First_Path(tree->noeud[tree->mod->startnode],0,tree,1);
		  printTreeState(tree->noeud[tree->mod->startnode],tree,1);printf(" 1\n");
	  }
	  Swap(d,b,c,a,tree); //swap nodes
	  pars2 = Fill_Sankoff(tree->noeud[tree->mod->startnode],tree,1);
	  if(tree->mod->optDebug){
	  	  Get_First_Path(tree->noeud[tree->mod->startnode],0,tree,1);
	  	  printTreeState(tree->noeud[tree->mod->startnode],tree,1);printf(" 2\n");
	  }
	  //printf("pars score %d\t%d\t%d\n",pars0,pars1,pars2);
	  if(pars0 != pars2){
		  printf("pars score %d\t%d\t%d\n",pars0,pars1,pars2);
		  printf("\n.\tParsimony reconstruction swap inconsistent!\n");
		  exit(EXIT_FAILURE);
	  }
	  return pars1;
}

/*********************************************************
 * Fill in dynamic programming tables of the Sankoff algorithm
 * */
int Fill_Sankoff(t_node *d, t_tree *tree, int root){
  int i,j,k,dir1,dir2; //Recurse to a tip!
  if(!d->tax || root){
	  if(!root){ //if at root, only have a single descendant
		  t_node* a = d->anc;
		  dir1=dir2=-1;
		  For(i,3) if(d->v[i]->num != a->num) (dir1<0)?(dir1=i):(dir2=i);
	  }else{
		  dir1=0;
	  }
	  For(j,tree->nstate){
		  if(!root)d->s[j]=1000;
		  else d->sroot[j]=1000;
		  d->lmin[j]=MAX_PARS;
		  d->rmin[j]=MAX_PARS;
	  }
	  if(!root && tree->mod->optDebug)printf("%d\t%d\t%d\n",d->num,d->v[dir1]->num,d->v[dir2]->num);
      Fill_Sankoff(d->v[dir1],tree,0); //left = 1, right = 2
      if(!root)Fill_Sankoff(d->v[dir2],tree,0);
      //fill in pointers and minimums of dynamic programming table
      For(i,tree->nstate){
    		For(j,tree->nstate){
    			if(root){
    				d->pl[i*tree->nstate+j]=d->s[i]+tree->step_mat[i*tree->nstate+j]+d->v[dir1]->s[j];
    				if(d->pl[i*tree->nstate+j] < d->lmin[i])d->lmin[i]=d->pl[i*tree->nstate+j];
    			}else{
    				d->pl[i*tree->nstate+j]=d->v[dir1]->s[j]+tree->step_mat[i*tree->nstate+j];
    				d->pr[i*tree->nstate+j]=d->v[dir2]->s[j]+tree->step_mat[i*tree->nstate+j];
    				if(d->pl[i*tree->nstate+j] < d->lmin[i])d->lmin[i]=d->pl[i*tree->nstate+j];
    				if(d->pr[i*tree->nstate+j] < d->rmin[i])d->rmin[i]=d->pr[i*tree->nstate+j];
    			}
    			if(tree->mod->optDebug)printf("%d\t%d\t%d\t%d\t%d\t%d\n",i,j,d->pl[i*tree->nstate+j],d->pr[i*tree->nstate+j],d->lmin[i],d->rmin[i]);
    		}
    		if(!root)d->s[i]=d->rmin[i]+d->lmin[i];
    		else d->sroot[i]=d->lmin[i];
      }
      if(tree->mod->optDebug)printf("\n");
  }
  int min=MAX_PARS;
  For(i,tree->nstate){
	  if(!root){
		  if(d->s[i]<min)min=d->s[i];
	  }else{
		  if(d->sroot[i]<min)min=d->sroot[i];
	  }
  }
  return min;
}

/********************************************************
 * Get "leftmost" maximally parsimonious labeling of tree
 * */
void Get_First_Path(t_node *d, int index, t_tree *tree,int root){
	int i,j,dir1,dir2;
	int rfound,lfound=0;
	d->pstate=index;
	if(tree->mod->optDebug)printf("%s\n",tree->chars[index]);
	if(!d->tax || root){
		if(!root){
			t_node* a = d->anc;
			dir1=dir2=-1;
			For(i,3) if(d->v[i]->num != a->num) (dir1<0)?(dir1=i):(dir2=i);
		}else{
			dir1=0;
		}
		Get_First_Path(d->v[dir1],d->plc[index],tree,0);
		if(!root)Get_First_Path(d->v[dir2],d->prc[index],tree,0);
	}
}

/********************************************************
 * Recurse down tree, randomly choosing ambiguous pointers
 * */
void Get_Rand_Path(t_node *d, int index, t_tree *tree, int root){
	int i,j,dir1,dir2;
	int ldraw=0;
	int rdraw=0;
	int lfound=0;
	int rfound=0;
	d->pstate=index;
	if(!d->tax&&tree->mod->optDebug)printf("%s\t%d\n",tree->chars[index],index);
	if(!d->tax || root){
		if(!root){
			t_node* a = d->anc;
			dir1=dir2=-1;
			For(i,3) if(d->v[i]->num != a->num) (dir1<0)?(dir1=i):(dir2=i);
		}else{
			dir1=0;
		}
		int lmins=0;int rmins=0; //tally up number of minimal paths
		For(i,tree->nstate){
			if(d->pl[index*tree->nstate+i] == d->lmin[index])lmins++;
			if(!root)if(d->pr[index*tree->nstate+i] == d->rmin[index])rmins++;
		}
		if(lmins>1)ldraw = rand() % lmins; //not perfectly random but probably okay
		if(!root && rmins>1)rdraw = rand() % rmins;
		For(i,tree->nstate){
			if(d->pl[index*tree->nstate+i] == d->lmin[index]){
				if(lfound == ldraw)Get_Rand_Path(d->v[dir1],i,tree,0);
				lfound++;
			}
			if(d->pr[index*tree->nstate+i] == d->rmin[index] && !root ){
				if(rfound == rdraw)Get_Rand_Path(d->v[dir2],i,tree,0);
				rfound++;
			}
		}
	}else{
		if(tree->mod->optDebug)printf("%s\t%s\n",tree->chars[index],d->name);
	}
}

/********************************************************
 * Set counters to "left most" option
 * */
void Set_Pars_Counters(t_node *d, t_tree *tree,int root){
	int i,j,k,dir1,dir2; //Recurse to a tip!
	if(!d->tax || root){
	  if(!root){ //if at root, only have a single descendant
		  t_node* a = d->anc;
		  dir1=dir2=-1;
		  For(i,3) if(d->v[i]->num != a->num) (dir1<0)?(dir1=i):(dir2=i);
	  }else{
		  dir1=0;
	  }
	  int index;
	  int lfound,rfound;
	  For(index,tree->nstate){
		  lfound=rfound=0;
		  For(i,tree->nstate){
			  if(d->pl[index*tree->nstate+i] == d->lmin[index]){
				  if(!lfound)d->plc[index]=i;
				  lfound++;
			  }
			  if(d->pr[index*tree->nstate+i] == d->rmin[index] && !root ){
				  if(!rfound)d->prc[index]=i;
				  rfound++;
			  }
		  }
		  if(lfound==1)d->llock[index]=1;
		  if(!root)if(rfound==1)d->rlock[index]=1;
		  if(!root)if(lfound==0 || rfound==0){
			  printf("No valid pointer found in Sankoff algorithm!\n");
			  exit(EXIT_FAILURE);
		  }
	  }
	  Set_Pars_Counters(d->v[dir1],tree,0);
	  if(!root)Set_Pars_Counters(d->v[dir2],tree,0);
	}
}

/********************************************************
 * Get all possible maximum parsimony labels of internal nodes
 * */
int Get_All_Paths(t_node *d, int index, t_tree *tree, t_tree** btrees, int root,int maxtrees,int treeindex, int repindex, int rootstate){
	//printf("ti %d\n",treeindex);
	int i,j,dir1,dir2;
	int lfound=0;
	int rfound=0;
	d->pstate=index;
	if(!d->tax || root){
		if(!root){
			t_node* a = d->anc;
			dir1=dir2=-1;
			For(i,3) if(d->v[i]->num != a->num) (dir1<0)?(dir1=i):(dir2=i);
		}else{
			dir1=0;
		}
		int lm = d->lmin[index]; //store original minimums and pointers
		int rm = d->rmin[index];
		int lc = d->plc[index];
		int rc = d->prc[index];

		//printf("Left tree! %s\t%d\t%d\t%d\n",tree->chars[lc],treeindex,d->llock[index],d->plc[index]);
		//exit(EXIT_FAILURE);
		if(d->llock[index]){ //if left node pointer is locked, move down with that assignment
			treeindex=Get_All_Paths(d->v[dir1],d->plc[index],tree,btrees,0,maxtrees,treeindex,repindex,rootstate);
		}else{
			For(i,tree->nstate){
				if(d->pl[index*tree->nstate+i] == d->lmin[index]){
					if(lfound>0 && treeindex+1 < maxtrees){ //if alternate left path found
						//t_tree* tree2; //read in copy of that tree
						//tree2 = Read_User_Tree(tree->io->tree_s[repindex]->data,tree->io->mod_s[repindex],tree->io);
						t_tree* tree2 = Make_Tree_From_Scratch(tree->n_otu,tree->data);
						Copy_Tree(tree,tree2);
						tree2->mod=tree->mod;
						tree2->mod->startnode = tree->mod->startnode;
						t_node* r2 = tree2->noeud[tree->mod->startnode];
						tree2->io=tree->io;
						Update_Ancestors_Edge(r2,r2->v[0],r2->b[0],tree);
						d->plc[index]=i;
						d->llock[index]=1; //copy all data structures to new tree, but with this node locked and set to the new pointer
						Copy_Sankoff_Tree(tree,tree2);
						treeindex++;
						btrees[treeindex]=tree2; //start recursion over at the root with this selection locked
						treeindex = Get_All_Paths(r2,rootstate,tree2,btrees,1,maxtrees,treeindex,repindex,rootstate);
					}
					lfound++;
				}
			}
			d->plc[index]=lc;//continue down with original pointer
			d->llock[index]=1; //with path through node fixed
			treeindex=Get_All_Paths(d->v[dir1],d->plc[index],tree,btrees,0,maxtrees,treeindex,repindex,rootstate);
		}
		if(!root){
			//printf("Right tree! %s\t%d\t%d\t%d\n",tree->chars[rc],treeindex,d->rlock[index],d->prc[index]);
			if(d->rlock[index]){ //right node
				treeindex=Get_All_Paths(d->v[dir2],d->prc[index],tree,btrees,0,maxtrees,treeindex,repindex,rootstate);
			}else{
				For(i,tree->nstate){
					if(d->pr[index*tree->nstate+i] == d->rmin[index]){
						if(rfound>0 && treeindex+1 < maxtrees){ //if alternate left path found
							//t_tree* tree2;
							//tree2 = Read_User_Tree(tree->io->tree_s[repindex]->data,tree->io->mod_s[repindex],tree->io);
							t_tree* tree2 = Make_Tree_From_Scratch(tree->n_otu,tree->data);
							Copy_Tree(tree,tree2);
							tree2->mod=tree->mod;
							tree2->mod->startnode = tree->mod->startnode;
							t_node* r2 = tree2->noeud[tree2->mod->startnode];
							tree2->io=tree->io;
							Update_Ancestors_Edge(r2,r2->v[0],r2->b[0],tree);
							d->prc[index]=i;
							d->rlock[index]=1;
							Copy_Sankoff_Tree(tree,tree2);
							treeindex++;
							btrees[treeindex]=tree2;
							treeindex = Get_All_Paths(r2,rootstate,tree2,btrees,1,maxtrees,treeindex,repindex,rootstate);
						}
						rfound++;
					}
				}
				d->prc[index]=rc;//continue down original tree
				d->rlock[index]=1; //with path through node fixed
				treeindex=Get_All_Paths(d->v[dir2],d->prc[index],tree,btrees,0,maxtrees,treeindex,repindex,rootstate);
			}
		}
	}else{
		if(tree->mod->optDebug)printf("%s\t%s\n",tree->chars[index],d->name);
	}
	return treeindex;
}


/********************************************************
 * Recursively add values to switch and length matrices for each type
 * */
void Fill_Pars_Stats(t_node* d,t_tree* tree, phydbl* switches, phydbl* classl, int root){
	int i,j,dir1,dir2;
	dir1=dir2=-1;
	if(!root){
		//printf("%d\t%d\n",d->anc->pstate,d->pstate);
		//if(d->pstate != d->anc->pstate){
			//printf("%d\t%d\n",d->anc->pstate,d->pstate);
			switches[d->anc->pstate*tree->nstate+d->pstate]++;
			classl[d->pstate] += d->anc_edge->l*0.5; //split switched branch lengths
			classl[d->anc->pstate] += d->anc_edge->l*0.5;
		/*}else{
			classl[d->pstate] += d->anc_edge->l;
		}*/
	}
	if(!d->tax){
		t_node* a = d->anc;
		For(i,3) if(d->v[i]->num != a->num) (dir1<0)?(dir1=i):(dir2=i);
	}else{
		dir1=0;
	}

	if(!d->tax||root){
		Fill_Pars_Stats(d->v[dir1],tree,switches,classl,0);
		if(!root)Fill_Pars_Stats(d->v[dir2],tree,switches,classl,0);
	}
}


/********************************************************
 * Print summary of states at each node. Used for debugging
 * */
void printTreeState(t_node *d, t_tree* tree, int root){
	printf("%s,",tree->chars[d->pstate]);
	int dir1,dir2,i;
	if(!d->tax || root){
		if(!root){
			t_node* a = d->anc;
			dir1=dir2=-1;
			For(i,3) if(d->v[i]->num != a->num) (dir1<0)?(dir1=i):(dir2=i);
		}else{
			dir1=0;
		}
		printTreeState(d->v[dir1],tree,0);
		if(!root)printTreeState(d->v[dir2],tree,0);
	}
}


/*********************************************************/

void Make_Tree_4_Pars(t_tree *tree, int n_site)
{
  int i;
  tree->site_pars = (int *)mCalloc(tree->n_pattern, sizeof(int));
  tree->step_mat = (int *)mCalloc(tree->mod->ns * tree->mod->ns, sizeof(int));
  For(i,2*tree->n_otu-3) Make_Edge_Pars(tree->t_edges[i],tree);
  Init_Ui_Tips(tree);
  Init_P_Pars_Tips(tree); /* Must be called after Init_Ui_Tips is called */
  Get_Step_Mat(tree);
}

/*********************************************************/

int Pars(t_tree *tree)
{
  int site,n_patterns;

  n_patterns = tree->n_pattern;

  Post_Order_Pars(tree->noeud[0],tree->noeud[0]->v[0],tree);
  if(tree->both_sides) Pre_Order_Pars(tree->noeud[0],tree->noeud[0]->v[0],tree);
  
  tree->c_pars = 0;
  For(site,n_patterns)
    {
      tree->site_pars[site] = 0;
      tree->curr_site       = site;
      Site_Pars(tree);
      tree->c_pars += tree->site_pars[site] * tree->data->wght[site];
    }

  return tree->c_pars;
}

/*********************************************************/

void Post_Order_Pars(t_node *a, t_node *d, t_tree *tree)
{
  int i,dir;

  dir = -1;

  if(d->tax) return;
  else
    {
      For(i,3)
	{
	  if(d->v[i] != a)
	    Post_Order_Pars(d,d->v[i],tree);
	  else dir = i;
	}
      Get_All_Partial_Pars(tree,d->b[dir],a,d);
    }
}

/*********************************************************/

void Pre_Order_Pars(t_node *a, t_node *d, t_tree *tree)
{
  int i;

  if(d->tax) return;
  else
    {
      For(i,3)
	{
	  if(d->v[i] != a)
	    {
	      Get_All_Partial_Pars(tree,d->b[i],d->v[i],d);
	      Pre_Order_Pars(d,d->v[i],tree);
	    }
	}
    }
}

/*********************************************************/

void Get_All_Partial_Pars(t_tree *tree, t_edge *b_fcus, t_node *a, t_node *d)
{
  if(d->tax) return;
  else Update_P_Pars(tree,b_fcus,d);
}

/*********************************************************/

void Site_Pars(t_tree *tree)
{
  tree->site_pars[tree->curr_site] = Pars_Core(tree->noeud[0]->b[0],tree);
}

/*********************************************************/

void Init_P_Pars_Tips(t_tree *tree)
{
  int curr_site,i,j;
  short int *state_v;
  int dim1;
  short int array[64];
  dim1 = tree->mod->ns;

  state_v = (short int *)mCalloc(tree->mod->ns,sizeof(short int));

  For(curr_site,tree->data->crunch_len)
    {
      For(i,tree->n_otu)
	{
	  if(tree->noeud[i]->b[0]->rght->tax != 1)
	    {
	      PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
	      Warn_And_Exit("\n");	    
	    }

	  if(tree->mod->datatype == NT)
	    {
	      Init_Tips_At_One_Site_Nucleotides_Int(tree->data->c_seq[i]->state[curr_site],
						    0,
						    state_v);	      
	      For(j,tree->mod->ns) tree->noeud[i]->b[0]->p_pars_r[curr_site*dim1+j] = MAX_PARS;
	      For(j,tree->mod->ns) if(state_v[j] > 0.5) tree->noeud[i]->b[0]->p_pars_r[curr_site*dim1+j] =  0;
	    }
	  else if(tree->mod->datatype == AA)
	    {
	      Init_Tips_At_One_Site_AA_Int(tree->data->c_seq[i]->state[curr_site],
					   0,
					   state_v);
	      For(j,tree->mod->ns) tree->noeud[i]->b[0]->p_pars_r[curr_site*dim1+j] = MAX_PARS;
	      For(j,tree->mod->ns) if(state_v[j] > 0.5) tree->noeud[i]->b[0]->p_pars_r[curr_site*dim1+j] =  0;
	    }
	  else if(tree->mod->datatype == GENERIC)
	    {
	      Init_Tips_At_One_Site_Generic_Int(tree->data->c_seq[i]->state+curr_site*tree->mod->state_len,
						tree->mod->ns,
						tree->mod->state_len,
						0,
						state_v);
	      For(j,tree->mod->ns) tree->noeud[i]->b[0]->p_pars_r[curr_site*dim1+j] = MAX_PARS;
	      For(j,tree->mod->ns) if(state_v[j] > 0.5) tree->noeud[i]->b[0]->p_pars_r[curr_site*dim1+j] =  0;
	    }
	    else if(tree->mod->datatype == CODON)                                    //!<Added by Marcelo.
	    {
	      Init_Tips_At_One_Site_Codons_Int(tree->data->c_seq[i]->state[curr_site],
					       0,
					       array, 
					       tree->data->c_seq[i]->alternativeCodons[curr_site]);
	      For(j,tree->mod->ns) tree->noeud[i]->b[0]->p_pars_r[curr_site*dim1+j] = MAX_PARS;
	      For(j,tree->mod->ns) if(array[j] > 0.5) tree->noeud[i]->b[0]->p_pars_r[curr_site*dim1+j] =  0;
	    }
	}
    }
  free(state_v);
}

/*********************************************************/

void Init_Ui_Tips(t_tree *tree)
{  
  int curr_site,i,j,br;
  short int *state_v;
  short int array[64];
  state_v = (short int *)mCalloc(tree->mod->ns,sizeof(short int));

  For(curr_site,tree->data->crunch_len)
    {
      For(i,tree->n_otu)
	{
	  if(tree->mod->datatype == NT)
	    {
	      if(tree->noeud[i]->b[0]->rght->tax != 1)
		{
		  PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
		  Warn_And_Exit("\n");
		}

	      Init_Tips_At_One_Site_Nucleotides_Int(tree->data->c_seq[i]->state[curr_site],
						    0,
						    state_v);	      
	      tree->noeud[i]->b[0]->ui_r[curr_site] = 0;
	      For(j,tree->mod->ns) tree->noeud[i]->b[0]->ui_r[curr_site] += (unsigned int)(state_v[j] * POW(2,j));
	    }
	  else if(tree->mod->datatype == AA)
	    {
	      Init_Tips_At_One_Site_AA_Int(tree->data->c_seq[i]->state[curr_site],
					   0,
					   state_v);
	      tree->noeud[i]->b[0]->ui_r[curr_site] = 0;
	      For(j,tree->mod->ns) tree->noeud[i]->b[0]->ui_r[curr_site] += (unsigned int)(state_v[j] * POW(2,j));
	    }
	  else if(tree->mod->datatype == GENERIC)
	    {
	      Init_Tips_At_One_Site_Generic_Int(tree->data->c_seq[i]->state+curr_site*tree->mod->state_len,
						tree->mod->ns,
						tree->mod->state_len,
						0,
						state_v);
	      tree->noeud[i]->b[0]->ui_r[curr_site] = 0;
	      For(j,tree->mod->ns) tree->noeud[i]->b[0]->ui_r[curr_site] += (unsigned int)(state_v[j] * POW(2,j));
	    }
	    else if(tree->mod->datatype == CODON)                                    //!<Added by Marcelo.
	    {
	      Init_Tips_At_One_Site_Codons_Int(tree->data->c_seq[i]->state[curr_site],
					       0,
					       array, 
					       tree->data->c_seq[i]->alternativeCodons[curr_site]);
	      tree->noeud[i]->b[0]->ui_r[curr_site] = 0;
	      For(j,tree->mod->ns) tree->noeud[i]->b[0]->ui_r[curr_site] += (unsigned int)(array[j] * POW(2,j));
	    }


	}
    }


  For(br,2*tree->n_otu-3)
    {
      For(curr_site,tree->data->crunch_len)
	{
	  tree->t_edges[br]->pars_r[curr_site] = 0;
	  tree->t_edges[br]->pars_l[curr_site] = 0;
	}
    }


  free(state_v);
}

/*********************************************************/

void Update_P_Pars(t_tree *tree, t_edge *b_fcus, t_node *n)
{
/*  
           |
	   |<- b_cus
	   |
	   n
          / \
       	 /   \
       	/     \
*/

  int i,j;
  int site;
  unsigned int *ui, *ui_v1, *ui_v2;
  int *p_pars_v1, *p_pars_v2, *p_pars;
  int *pars, *pars_v1, *pars_v2;
  int n_patterns,matches;
  int min_v1,min_v2;
  int v;
  int dim1;

  dim1 = tree->mod->ns;
  matches = 0;
  ui = ui_v1 = ui_v2 = NULL;
  p_pars = p_pars_v1 = p_pars_v2 = NULL;
  pars = pars_v1 = pars_v2 = NULL;

  n_patterns = tree->n_pattern;

  if(n == b_fcus->left)
    {	     
      ui = b_fcus->ui_l;

      pars = b_fcus->pars_l;
      p_pars = b_fcus->p_pars_l;

      ui_v1 = 
      (n == n->b[b_fcus->l_v1]->left)?
      (n->b[b_fcus->l_v1]->ui_r):
      (n->b[b_fcus->l_v1]->ui_l);

      ui_v2 = 
      (n == n->b[b_fcus->l_v2]->left)?
      (n->b[b_fcus->l_v2]->ui_r):
      (n->b[b_fcus->l_v2]->ui_l);

      p_pars_v1 = 
      (n == n->b[b_fcus->l_v1]->left)?
      (n->b[b_fcus->l_v1]->p_pars_r):
      (n->b[b_fcus->l_v1]->p_pars_l);

      p_pars_v2 = 
      (n == n->b[b_fcus->l_v2]->left)?
      (n->b[b_fcus->l_v2]->p_pars_r):
      (n->b[b_fcus->l_v2]->p_pars_l);

      pars_v1 = 
      (n == n->b[b_fcus->l_v1]->left)?
      (n->b[b_fcus->l_v1]->pars_r):
      (n->b[b_fcus->l_v1]->pars_l);

      pars_v2 = 
      (n == n->b[b_fcus->l_v2]->left)?
      (n->b[b_fcus->l_v2]->pars_r):
      (n->b[b_fcus->l_v2]->pars_l);
    }
  else
    {
      ui = b_fcus->ui_r;
      
      pars = b_fcus->pars_r;
      p_pars = b_fcus->p_pars_r;

      ui_v1 = 
      (n == n->b[b_fcus->r_v1]->left)?
      (n->b[b_fcus->r_v1]->ui_r):
      (n->b[b_fcus->r_v1]->ui_l);

      ui_v2 = 
      (n == n->b[b_fcus->r_v2]->left)?
      (n->b[b_fcus->r_v2]->ui_r):
      (n->b[b_fcus->r_v2]->ui_l);

      p_pars_v1 = 
      (n == n->b[b_fcus->r_v1]->left)?
      (n->b[b_fcus->r_v1]->p_pars_r):
      (n->b[b_fcus->r_v1]->p_pars_l);

      p_pars_v2 = 
      (n == n->b[b_fcus->r_v2]->left)?
      (n->b[b_fcus->r_v2]->p_pars_r):
      (n->b[b_fcus->r_v2]->p_pars_l);

      pars_v1 = 
      (n == n->b[b_fcus->r_v1]->left)?
      (n->b[b_fcus->r_v1]->pars_r):
      (n->b[b_fcus->r_v1]->pars_l);

      pars_v2 = 
      (n == n->b[b_fcus->r_v2]->left)?
      (n->b[b_fcus->r_v2]->pars_r):
      (n->b[b_fcus->r_v2]->pars_l);
    }


  if(tree->mod->s_opt->general_pars)
    {
      For(site,n_patterns)
	{
	  For(i,tree->mod->ns)
	    {
	      min_v1 = MAX_PARS;
	      For(j,tree->mod->ns)
		{
		  v = p_pars_v1[site*dim1+j] + tree->step_mat[i*tree->mod->ns+j]; 
		  if(v < min_v1) min_v1 = v;
		}
	      
	      min_v2 = MAX_PARS;
	      For(j,tree->mod->ns)
		{
		  v = p_pars_v2[site*dim1+j] + tree->step_mat[i*tree->mod->ns+j]; 
		  if(v < min_v2) min_v2 = v;
		}	      
	      p_pars[site*dim1+i] = min_v1 + min_v2;
	    }	  
	}
    }
  else
    {
      For(site,n_patterns)
	{
	  pars[site] = pars_v1[site] + pars_v2[site];

	  ui[site] = ui_v1[site] & ui_v2[site];

	  if(!ui[site])
	    {
	      pars[site]++;
	      ui[site] = ui_v1[site] | ui_v2[site];
	    }
	}
    }
}

/*********************************************************/

int Pars_Core(t_edge *b, t_tree *tree)
{
  int site;
  int i,j;
  int site_pars;
  int min_l,min_r;
  int v;
  int dim1;

  dim1 = tree->mod->ns;
  site = tree->curr_site;
  site_pars = MAX_PARS;

  if(tree->mod->s_opt->general_pars){
      For(i,tree->mod->ns){
    	  	  min_l = MAX_PARS;
    	  	  For(j,tree->mod->ns){
    	  		  v = b->p_pars_l[site*dim1+j] + tree->step_mat[i*tree->mod->ns+j];
    	  		  if(v < min_l) min_l = v;
    	  	  }

    	  	  min_r = MAX_PARS;
    	  	  For(j,tree->mod->ns){
    	  		  v = b->p_pars_r[site*dim1+j] + tree->step_mat[i*tree->mod->ns+j];
    	  		  if(v < min_r) min_r = v;
    	  	  }
    	  	  if((min_l + min_r) < site_pars) site_pars = min_l + min_r;
      }
    }else{
      site_pars = b->pars_l[site] + b->pars_r[site];      
      if(!(b->ui_l[site] & b->ui_r[site])) site_pars++;
    }
  return site_pars;
}

/*********************************************************/
/* Is there one or more parsimoniy step(s) along this t_edge ? 
   0 -> NO; 1 -> YES
*/
int One_Pars_Step(t_edge *b,t_tree *tree)
{
  int site;
  int init_general_pars;

  init_general_pars = tree->mod->s_opt->general_pars;
  
  tree->mod->s_opt->general_pars = 0;
  tree->both_sides   = 1;
  Pars(tree);

  For(site,tree->n_pattern)
    {
      if(!(b->ui_l[site] & b->ui_r[site])) break;
    }
  tree->mod->s_opt->general_pars = init_general_pars;
  if(site == tree->n_pattern) return 0;
  else                        
    {  
      PhyML_Printf("\n. One parsimony step ocurred at site %4d",site);
      return 1;
    }
}

/*********************************************************/
int Pars_At_Given_Edge(t_edge *b, t_tree *tree)
{
  int site,n_patterns;
  
/*   n_patterns = (int)FLOOR(tree->n_pattern*tree->prop_of_sites_to_consider); */
  n_patterns = tree->n_pattern;

  tree->c_pars = .0;
  For(site,n_patterns)
    {
      tree->site_pars[site] = 0;
      tree->curr_site = site;
      tree->site_pars[site] = Pars_Core(b,tree);
      tree->c_pars += tree->site_pars[site] * tree->data->wght[site];
    }
  return tree->c_pars;
}

/*********************************************************/

int Update_Pars_At_Given_Edge(t_edge *b_fcus, t_tree *tree)
{
  Update_P_Pars(tree,b_fcus,b_fcus->left);
  Update_P_Pars(tree,b_fcus,b_fcus->rght);
  tree->c_pars = Pars_At_Given_Edge(b_fcus,tree);
  return tree->c_pars;
}

/*********************************************************/

void Get_Step_Mat(t_tree *tree)
{
  int i,j,k,codoni,codonj,icodon[3],jcodon[3],diff;
  if(tree->mod->datatype== CODON)      //!< Added by Marcelo.
  {
    For(i,tree->mod->ns)
    {
      For(j,tree->mod->ns)
      {
	diff=0;
	codoni=senseCodons[i];
	codonj=senseCodons[j];
	For(k,3)
	{
	  icodon[k]=codoni-((codoni>>2)<<2);
	  codoni=codoni>>2;
	  jcodon[k]=codonj-((codonj>>2)<<2);
	  codonj=codonj>>2;
	  if(icodon[k]!=jcodon[k]) diff++;
	}
	tree->step_mat[ i*tree->mod->ns + j] =    diff ;
      }
    }
  }
  else if(tree->mod->datatype == AA)
    {
      tree->step_mat[ 0*tree->mod->ns+ 0] =    0 ;
      tree->step_mat[ 0*tree->mod->ns+ 1] =    3 ;
      tree->step_mat[ 0*tree->mod->ns+ 2] =    3 ;
      tree->step_mat[ 0*tree->mod->ns+ 3] =    2 ;
      tree->step_mat[ 0*tree->mod->ns+ 4] =    3 ;
      tree->step_mat[ 0*tree->mod->ns+ 5] =    3 ;
      tree->step_mat[ 0*tree->mod->ns+ 6] =    2 ;
      tree->step_mat[ 0*tree->mod->ns+ 7] =    2 ;
      tree->step_mat[ 0*tree->mod->ns+ 8] =    3 ;
      tree->step_mat[ 0*tree->mod->ns+ 9] =    3 ;
      tree->step_mat[ 0*tree->mod->ns+10] =    3 ;
      tree->step_mat[ 0*tree->mod->ns+11] =    3 ;
      tree->step_mat[ 0*tree->mod->ns+12] =    3 ;
      tree->step_mat[ 0*tree->mod->ns+13] =    3 ;
      tree->step_mat[ 0*tree->mod->ns+14] =    2 ;
      tree->step_mat[ 0*tree->mod->ns+15] =    2 ;
      tree->step_mat[ 0*tree->mod->ns+16] =    2 ;
      tree->step_mat[ 0*tree->mod->ns+17] =    3 ;
      tree->step_mat[ 0*tree->mod->ns+18] =    3 ;
      tree->step_mat[ 0*tree->mod->ns+19] =    2 ;
      tree->step_mat[ 1*tree->mod->ns+ 0] =    3 ;
      tree->step_mat[ 1*tree->mod->ns+ 1] =    0 ;
      tree->step_mat[ 1*tree->mod->ns+ 2] =    2 ;
      tree->step_mat[ 1*tree->mod->ns+ 3] =    3 ;
      tree->step_mat[ 1*tree->mod->ns+ 4] =    2 ;
      tree->step_mat[ 1*tree->mod->ns+ 5] =    2 ;
      tree->step_mat[ 1*tree->mod->ns+ 6] =    3 ;
      tree->step_mat[ 1*tree->mod->ns+ 7] =    2 ;
      tree->step_mat[ 1*tree->mod->ns+ 8] =    2 ;
      tree->step_mat[ 1*tree->mod->ns+ 9] =    2 ;
      tree->step_mat[ 1*tree->mod->ns+10] =    2 ;
      tree->step_mat[ 1*tree->mod->ns+11] =    2 ;
      tree->step_mat[ 1*tree->mod->ns+12] =    2 ;
      tree->step_mat[ 1*tree->mod->ns+13] =    3 ;
      tree->step_mat[ 1*tree->mod->ns+14] =    2 ;
      tree->step_mat[ 1*tree->mod->ns+15] =    2 ;
      tree->step_mat[ 1*tree->mod->ns+16] =    2 ;
      tree->step_mat[ 1*tree->mod->ns+17] =    2 ;
      tree->step_mat[ 1*tree->mod->ns+18] =    3 ;
      tree->step_mat[ 1*tree->mod->ns+19] =    3 ;
      tree->step_mat[ 2*tree->mod->ns+ 0] =    3 ;
      tree->step_mat[ 2*tree->mod->ns+ 1] =    2 ;
      tree->step_mat[ 2*tree->mod->ns+ 2] =    0 ;
      tree->step_mat[ 2*tree->mod->ns+ 3] =    2 ;
      tree->step_mat[ 2*tree->mod->ns+ 4] =    2 ;
      tree->step_mat[ 2*tree->mod->ns+ 5] =    2 ;
      tree->step_mat[ 2*tree->mod->ns+ 6] =    2 ;
      tree->step_mat[ 2*tree->mod->ns+ 7] =    3 ;
      tree->step_mat[ 2*tree->mod->ns+ 8] =    2 ;
      tree->step_mat[ 2*tree->mod->ns+ 9] =    2 ;
      tree->step_mat[ 2*tree->mod->ns+10] =    3 ;
      tree->step_mat[ 2*tree->mod->ns+11] =    1 ;
      tree->step_mat[ 2*tree->mod->ns+12] =    2 ;
      tree->step_mat[ 2*tree->mod->ns+13] =    2 ;
      tree->step_mat[ 2*tree->mod->ns+14] =    3 ;
      tree->step_mat[ 2*tree->mod->ns+15] =    2 ;
      tree->step_mat[ 2*tree->mod->ns+16] =    2 ;
      tree->step_mat[ 2*tree->mod->ns+17] =    3 ;
      tree->step_mat[ 2*tree->mod->ns+18] =    2 ;
      tree->step_mat[ 2*tree->mod->ns+19] =    3 ;
      tree->step_mat[ 3*tree->mod->ns+ 0] =    2 ;
      tree->step_mat[ 3*tree->mod->ns+ 1] =    3 ;
      tree->step_mat[ 3*tree->mod->ns+ 2] =    2 ;
      tree->step_mat[ 3*tree->mod->ns+ 3] =    0 ;
      tree->step_mat[ 3*tree->mod->ns+ 4] =    2 ;
      tree->step_mat[ 3*tree->mod->ns+ 5] =    2 ;
      tree->step_mat[ 3*tree->mod->ns+ 6] =    1 ;
      tree->step_mat[ 3*tree->mod->ns+ 7] =    2 ;
      tree->step_mat[ 3*tree->mod->ns+ 8] =    2 ;
      tree->step_mat[ 3*tree->mod->ns+ 9] =    3 ;
      tree->step_mat[ 3*tree->mod->ns+10] =    3 ;
      tree->step_mat[ 3*tree->mod->ns+11] =    2 ;
      tree->step_mat[ 3*tree->mod->ns+12] =    3 ;
      tree->step_mat[ 3*tree->mod->ns+13] =    2 ;
      tree->step_mat[ 3*tree->mod->ns+14] =    3 ;
      tree->step_mat[ 3*tree->mod->ns+15] =    3 ;
      tree->step_mat[ 3*tree->mod->ns+16] =    3 ;
      tree->step_mat[ 3*tree->mod->ns+17] =    3 ;
      tree->step_mat[ 3*tree->mod->ns+18] =    2 ;
      tree->step_mat[ 3*tree->mod->ns+19] =    2 ;
      tree->step_mat[ 4*tree->mod->ns+ 0] =    3 ;
      tree->step_mat[ 4*tree->mod->ns+ 1] =    2 ;
      tree->step_mat[ 4*tree->mod->ns+ 2] =    2 ;
      tree->step_mat[ 4*tree->mod->ns+ 3] =    2 ;
      tree->step_mat[ 4*tree->mod->ns+ 4] =    0 ;
      tree->step_mat[ 4*tree->mod->ns+ 5] =    3 ;
      tree->step_mat[ 4*tree->mod->ns+ 6] =    3 ;
      tree->step_mat[ 4*tree->mod->ns+ 7] =    2 ;
      tree->step_mat[ 4*tree->mod->ns+ 8] =    2 ;
      tree->step_mat[ 4*tree->mod->ns+ 9] =    3 ;
      tree->step_mat[ 4*tree->mod->ns+10] =    2 ;
      tree->step_mat[ 4*tree->mod->ns+11] =    3 ;
      tree->step_mat[ 4*tree->mod->ns+12] =    3 ;
      tree->step_mat[ 4*tree->mod->ns+13] =    2 ;
      tree->step_mat[ 4*tree->mod->ns+14] =    3 ;
      tree->step_mat[ 4*tree->mod->ns+15] =    2 ;
      tree->step_mat[ 4*tree->mod->ns+16] =    3 ;
      tree->step_mat[ 4*tree->mod->ns+17] =    1 ;
      tree->step_mat[ 4*tree->mod->ns+18] =    2 ;
      tree->step_mat[ 4*tree->mod->ns+19] =    3 ;
      tree->step_mat[ 5*tree->mod->ns+ 0] =    3 ;
      tree->step_mat[ 5*tree->mod->ns+ 1] =    2 ;
      tree->step_mat[ 5*tree->mod->ns+ 2] =    2 ;
      tree->step_mat[ 5*tree->mod->ns+ 3] =    2 ;
      tree->step_mat[ 5*tree->mod->ns+ 4] =    3 ;
      tree->step_mat[ 5*tree->mod->ns+ 5] =    0 ;
      tree->step_mat[ 5*tree->mod->ns+ 6] =    2 ;
      tree->step_mat[ 5*tree->mod->ns+ 7] =    3 ;
      tree->step_mat[ 5*tree->mod->ns+ 8] =    1 ;
      tree->step_mat[ 5*tree->mod->ns+ 9] =    3 ;
      tree->step_mat[ 5*tree->mod->ns+10] =    2 ;
      tree->step_mat[ 5*tree->mod->ns+11] =    2 ;
      tree->step_mat[ 5*tree->mod->ns+12] =    2 ;
      tree->step_mat[ 5*tree->mod->ns+13] =    3 ;
      tree->step_mat[ 5*tree->mod->ns+14] =    2 ;
      tree->step_mat[ 5*tree->mod->ns+15] =    3 ;
      tree->step_mat[ 5*tree->mod->ns+16] =    3 ;
      tree->step_mat[ 5*tree->mod->ns+17] =    2 ;
      tree->step_mat[ 5*tree->mod->ns+18] =    2 ;
      tree->step_mat[ 5*tree->mod->ns+19] =    3 ;
      tree->step_mat[ 6*tree->mod->ns+ 0] =    2 ;
      tree->step_mat[ 6*tree->mod->ns+ 1] =    3 ;
      tree->step_mat[ 6*tree->mod->ns+ 2] =    2 ;
      tree->step_mat[ 6*tree->mod->ns+ 3] =    1 ;
      tree->step_mat[ 6*tree->mod->ns+ 4] =    3 ;
      tree->step_mat[ 6*tree->mod->ns+ 5] =    2 ;
      tree->step_mat[ 6*tree->mod->ns+ 6] =    0 ;
      tree->step_mat[ 6*tree->mod->ns+ 7] =    2 ;
      tree->step_mat[ 6*tree->mod->ns+ 8] =    2 ;
      tree->step_mat[ 6*tree->mod->ns+ 9] =    3 ;
      tree->step_mat[ 6*tree->mod->ns+10] =    3 ;
      tree->step_mat[ 6*tree->mod->ns+11] =    2 ;
      tree->step_mat[ 6*tree->mod->ns+12] =    2 ;
      tree->step_mat[ 6*tree->mod->ns+13] =    3 ;
      tree->step_mat[ 6*tree->mod->ns+14] =    3 ;
      tree->step_mat[ 6*tree->mod->ns+15] =    3 ;
      tree->step_mat[ 6*tree->mod->ns+16] =    3 ;
      tree->step_mat[ 6*tree->mod->ns+17] =    2 ;
      tree->step_mat[ 6*tree->mod->ns+18] =    2 ;
      tree->step_mat[ 6*tree->mod->ns+19] =    2 ;
      tree->step_mat[ 7*tree->mod->ns+ 0] =    2 ;
      tree->step_mat[ 7*tree->mod->ns+ 1] =    2 ;
      tree->step_mat[ 7*tree->mod->ns+ 2] =    3 ;
      tree->step_mat[ 7*tree->mod->ns+ 3] =    2 ;
      tree->step_mat[ 7*tree->mod->ns+ 4] =    2 ;
      tree->step_mat[ 7*tree->mod->ns+ 5] =    3 ;
      tree->step_mat[ 7*tree->mod->ns+ 6] =    2 ;
      tree->step_mat[ 7*tree->mod->ns+ 7] =    0 ;
      tree->step_mat[ 7*tree->mod->ns+ 8] =    3 ;
      tree->step_mat[ 7*tree->mod->ns+ 9] =    3 ;
      tree->step_mat[ 7*tree->mod->ns+10] =    3 ;
      tree->step_mat[ 7*tree->mod->ns+11] =    3 ;
      tree->step_mat[ 7*tree->mod->ns+12] =    3 ;
      tree->step_mat[ 7*tree->mod->ns+13] =    3 ;
      tree->step_mat[ 7*tree->mod->ns+14] =    3 ;
      tree->step_mat[ 7*tree->mod->ns+15] =    2 ;
      tree->step_mat[ 7*tree->mod->ns+16] =    3 ;
      tree->step_mat[ 7*tree->mod->ns+17] =    2 ;
      tree->step_mat[ 7*tree->mod->ns+18] =    3 ;
      tree->step_mat[ 7*tree->mod->ns+19] =    2 ;
      tree->step_mat[ 8*tree->mod->ns+ 0] =    3 ;
      tree->step_mat[ 8*tree->mod->ns+ 1] =    2 ;
      tree->step_mat[ 8*tree->mod->ns+ 2] =    2 ;
      tree->step_mat[ 8*tree->mod->ns+ 3] =    2 ;
      tree->step_mat[ 8*tree->mod->ns+ 4] =    2 ;
      tree->step_mat[ 8*tree->mod->ns+ 5] =    1 ;
      tree->step_mat[ 8*tree->mod->ns+ 6] =    2 ;
      tree->step_mat[ 8*tree->mod->ns+ 7] =    3 ;
      tree->step_mat[ 8*tree->mod->ns+ 8] =    0 ;
      tree->step_mat[ 8*tree->mod->ns+ 9] =    3 ;
      tree->step_mat[ 8*tree->mod->ns+10] =    2 ;
      tree->step_mat[ 8*tree->mod->ns+11] =    2 ;
      tree->step_mat[ 8*tree->mod->ns+12] =    3 ;
      tree->step_mat[ 8*tree->mod->ns+13] =    2 ;
      tree->step_mat[ 8*tree->mod->ns+14] =    2 ;
      tree->step_mat[ 8*tree->mod->ns+15] =    3 ;
      tree->step_mat[ 8*tree->mod->ns+16] =    3 ;
      tree->step_mat[ 8*tree->mod->ns+17] =    3 ;
      tree->step_mat[ 8*tree->mod->ns+18] =    2 ;
      tree->step_mat[ 8*tree->mod->ns+19] =    3 ;
      tree->step_mat[ 9*tree->mod->ns+ 0] =    3 ;
      tree->step_mat[ 9*tree->mod->ns+ 1] =    2 ;
      tree->step_mat[ 9*tree->mod->ns+ 2] =    2 ;
      tree->step_mat[ 9*tree->mod->ns+ 3] =    3 ;
      tree->step_mat[ 9*tree->mod->ns+ 4] =    3 ;
      tree->step_mat[ 9*tree->mod->ns+ 5] =    3 ;
      tree->step_mat[ 9*tree->mod->ns+ 6] =    3 ;
      tree->step_mat[ 9*tree->mod->ns+ 7] =    3 ;
      tree->step_mat[ 9*tree->mod->ns+ 8] =    3 ;
      tree->step_mat[ 9*tree->mod->ns+ 9] =    0 ;
      tree->step_mat[ 9*tree->mod->ns+10] =    2 ;
      tree->step_mat[ 9*tree->mod->ns+11] =    2 ;
      tree->step_mat[ 9*tree->mod->ns+12] =    1 ;
      tree->step_mat[ 9*tree->mod->ns+13] =    2 ;
      tree->step_mat[ 9*tree->mod->ns+14] =    3 ;
      tree->step_mat[ 9*tree->mod->ns+15] =    2 ;
      tree->step_mat[ 9*tree->mod->ns+16] =    2 ;
      tree->step_mat[ 9*tree->mod->ns+17] =    3 ;
      tree->step_mat[ 9*tree->mod->ns+18] =    3 ;
      tree->step_mat[ 9*tree->mod->ns+19] =    2 ;
      tree->step_mat[10*tree->mod->ns+ 0] =    3 ;
      tree->step_mat[10*tree->mod->ns+ 1] =    2 ;
      tree->step_mat[10*tree->mod->ns+ 2] =    3 ;
      tree->step_mat[10*tree->mod->ns+ 3] =    3 ;
      tree->step_mat[10*tree->mod->ns+ 4] =    2 ;
      tree->step_mat[10*tree->mod->ns+ 5] =    2 ;
      tree->step_mat[10*tree->mod->ns+ 6] =    3 ;
      tree->step_mat[10*tree->mod->ns+ 7] =    3 ;
      tree->step_mat[10*tree->mod->ns+ 8] =    2 ;
      tree->step_mat[10*tree->mod->ns+ 9] =    2 ;
      tree->step_mat[10*tree->mod->ns+10] =    0 ;
      tree->step_mat[10*tree->mod->ns+11] =    3 ;
      tree->step_mat[10*tree->mod->ns+12] =    2 ;
      tree->step_mat[10*tree->mod->ns+13] =    2 ;
      tree->step_mat[10*tree->mod->ns+14] =    2 ;
      tree->step_mat[10*tree->mod->ns+15] =    3 ;
      tree->step_mat[10*tree->mod->ns+16] =    3 ;
      tree->step_mat[10*tree->mod->ns+17] =    2 ;
      tree->step_mat[10*tree->mod->ns+18] =    2 ;
      tree->step_mat[10*tree->mod->ns+19] =    2 ;
      tree->step_mat[11*tree->mod->ns+ 0] =    3 ;
      tree->step_mat[11*tree->mod->ns+ 1] =    2 ;
      tree->step_mat[11*tree->mod->ns+ 2] =    1 ;
      tree->step_mat[11*tree->mod->ns+ 3] =    2 ;
      tree->step_mat[11*tree->mod->ns+ 4] =    3 ;
      tree->step_mat[11*tree->mod->ns+ 5] =    2 ;
      tree->step_mat[11*tree->mod->ns+ 6] =    2 ;
      tree->step_mat[11*tree->mod->ns+ 7] =    3 ;
      tree->step_mat[11*tree->mod->ns+ 8] =    2 ;
      tree->step_mat[11*tree->mod->ns+ 9] =    2 ;
      tree->step_mat[11*tree->mod->ns+10] =    3 ;
      tree->step_mat[11*tree->mod->ns+11] =    0 ;
      tree->step_mat[11*tree->mod->ns+12] =    2 ;
      tree->step_mat[11*tree->mod->ns+13] =    3 ;
      tree->step_mat[11*tree->mod->ns+14] =    3 ;
      tree->step_mat[11*tree->mod->ns+15] =    2 ;
      tree->step_mat[11*tree->mod->ns+16] =    2 ;
      tree->step_mat[11*tree->mod->ns+17] =    2 ;
      tree->step_mat[11*tree->mod->ns+18] =    2 ;
      tree->step_mat[11*tree->mod->ns+19] =    3 ;
      tree->step_mat[12*tree->mod->ns+ 0] =    3 ;
      tree->step_mat[12*tree->mod->ns+ 1] =    2 ;
      tree->step_mat[12*tree->mod->ns+ 2] =    2 ;
      tree->step_mat[12*tree->mod->ns+ 3] =    3 ;
      tree->step_mat[12*tree->mod->ns+ 4] =    3 ;
      tree->step_mat[12*tree->mod->ns+ 5] =    2 ;
      tree->step_mat[12*tree->mod->ns+ 6] =    2 ;
      tree->step_mat[12*tree->mod->ns+ 7] =    3 ;
      tree->step_mat[12*tree->mod->ns+ 8] =    3 ;
      tree->step_mat[12*tree->mod->ns+ 9] =    1 ;
      tree->step_mat[12*tree->mod->ns+10] =    2 ;
      tree->step_mat[12*tree->mod->ns+11] =    2 ;
      tree->step_mat[12*tree->mod->ns+12] =    0 ;
      tree->step_mat[12*tree->mod->ns+13] =    2 ;
      tree->step_mat[12*tree->mod->ns+14] =    3 ;
      tree->step_mat[12*tree->mod->ns+15] =    2 ;
      tree->step_mat[12*tree->mod->ns+16] =    2 ;
      tree->step_mat[12*tree->mod->ns+17] =    2 ;
      tree->step_mat[12*tree->mod->ns+18] =    3 ;
      tree->step_mat[12*tree->mod->ns+19] =    2 ;
      tree->step_mat[13*tree->mod->ns+ 0] =    3 ;
      tree->step_mat[13*tree->mod->ns+ 1] =    3 ;
      tree->step_mat[13*tree->mod->ns+ 2] =    2 ;
      tree->step_mat[13*tree->mod->ns+ 3] =    2 ;
      tree->step_mat[13*tree->mod->ns+ 4] =    2 ;
      tree->step_mat[13*tree->mod->ns+ 5] =    3 ;
      tree->step_mat[13*tree->mod->ns+ 6] =    3 ;
      tree->step_mat[13*tree->mod->ns+ 7] =    3 ;
      tree->step_mat[13*tree->mod->ns+ 8] =    2 ;
      tree->step_mat[13*tree->mod->ns+ 9] =    2 ;
      tree->step_mat[13*tree->mod->ns+10] =    2 ;
      tree->step_mat[13*tree->mod->ns+11] =    3 ;
      tree->step_mat[13*tree->mod->ns+12] =    2 ;
      tree->step_mat[13*tree->mod->ns+13] =    0 ;
      tree->step_mat[13*tree->mod->ns+14] =    3 ;
      tree->step_mat[13*tree->mod->ns+15] =    2 ;
      tree->step_mat[13*tree->mod->ns+16] =    3 ;
      tree->step_mat[13*tree->mod->ns+17] =    2 ;
      tree->step_mat[13*tree->mod->ns+18] =    2 ;
      tree->step_mat[13*tree->mod->ns+19] =    2 ;
      tree->step_mat[14*tree->mod->ns+ 0] =    2 ;
      tree->step_mat[14*tree->mod->ns+ 1] =    2 ;
      tree->step_mat[14*tree->mod->ns+ 2] =    3 ;
      tree->step_mat[14*tree->mod->ns+ 3] =    3 ;
      tree->step_mat[14*tree->mod->ns+ 4] =    3 ;
      tree->step_mat[14*tree->mod->ns+ 5] =    2 ;
      tree->step_mat[14*tree->mod->ns+ 6] =    3 ;
      tree->step_mat[14*tree->mod->ns+ 7] =    3 ;
      tree->step_mat[14*tree->mod->ns+ 8] =    2 ;
      tree->step_mat[14*tree->mod->ns+ 9] =    3 ;
      tree->step_mat[14*tree->mod->ns+10] =    2 ;
      tree->step_mat[14*tree->mod->ns+11] =    3 ;
      tree->step_mat[14*tree->mod->ns+12] =    3 ;
      tree->step_mat[14*tree->mod->ns+13] =    3 ;
      tree->step_mat[14*tree->mod->ns+14] =    0 ;
      tree->step_mat[14*tree->mod->ns+15] =    2 ;
      tree->step_mat[14*tree->mod->ns+16] =    2 ;
      tree->step_mat[14*tree->mod->ns+17] =    3 ;
      tree->step_mat[14*tree->mod->ns+18] =    3 ;
      tree->step_mat[14*tree->mod->ns+19] =    3 ;
      tree->step_mat[15*tree->mod->ns+ 0] =    2 ;
      tree->step_mat[15*tree->mod->ns+ 1] =    2 ;
      tree->step_mat[15*tree->mod->ns+ 2] =    2 ;
      tree->step_mat[15*tree->mod->ns+ 3] =    3 ;
      tree->step_mat[15*tree->mod->ns+ 4] =    2 ;
      tree->step_mat[15*tree->mod->ns+ 5] =    3 ;
      tree->step_mat[15*tree->mod->ns+ 6] =    3 ;
      tree->step_mat[15*tree->mod->ns+ 7] =    2 ;
      tree->step_mat[15*tree->mod->ns+ 8] =    3 ;
      tree->step_mat[15*tree->mod->ns+ 9] =    2 ;
      tree->step_mat[15*tree->mod->ns+10] =    3 ;
      tree->step_mat[15*tree->mod->ns+11] =    2 ;
      tree->step_mat[15*tree->mod->ns+12] =    2 ;
      tree->step_mat[15*tree->mod->ns+13] =    2 ;
      tree->step_mat[15*tree->mod->ns+14] =    2 ;
      tree->step_mat[15*tree->mod->ns+15] =    0 ;
      tree->step_mat[15*tree->mod->ns+16] =    2 ;
      tree->step_mat[15*tree->mod->ns+17] =    2 ;
      tree->step_mat[15*tree->mod->ns+18] =    2 ;
      tree->step_mat[15*tree->mod->ns+19] =    3 ;
      tree->step_mat[16*tree->mod->ns+ 0] =    2 ;
      tree->step_mat[16*tree->mod->ns+ 1] =    2 ;
      tree->step_mat[16*tree->mod->ns+ 2] =    2 ;
      tree->step_mat[16*tree->mod->ns+ 3] =    3 ;
      tree->step_mat[16*tree->mod->ns+ 4] =    3 ;
      tree->step_mat[16*tree->mod->ns+ 5] =    3 ;
      tree->step_mat[16*tree->mod->ns+ 6] =    3 ;
      tree->step_mat[16*tree->mod->ns+ 7] =    3 ;
      tree->step_mat[16*tree->mod->ns+ 8] =    3 ;
      tree->step_mat[16*tree->mod->ns+ 9] =    2 ;
      tree->step_mat[16*tree->mod->ns+10] =    3 ;
      tree->step_mat[16*tree->mod->ns+11] =    2 ;
      tree->step_mat[16*tree->mod->ns+12] =    2 ;
      tree->step_mat[16*tree->mod->ns+13] =    3 ;
      tree->step_mat[16*tree->mod->ns+14] =    2 ;
      tree->step_mat[16*tree->mod->ns+15] =    2 ;
      tree->step_mat[16*tree->mod->ns+16] =    0 ;
      tree->step_mat[16*tree->mod->ns+17] =    3 ;
      tree->step_mat[16*tree->mod->ns+18] =    3 ;
      tree->step_mat[16*tree->mod->ns+19] =    3 ;
      tree->step_mat[17*tree->mod->ns+ 0] =    3 ;
      tree->step_mat[17*tree->mod->ns+ 1] =    2 ;
      tree->step_mat[17*tree->mod->ns+ 2] =    3 ;
      tree->step_mat[17*tree->mod->ns+ 3] =    3 ;
      tree->step_mat[17*tree->mod->ns+ 4] =    1 ;
      tree->step_mat[17*tree->mod->ns+ 5] =    2 ;
      tree->step_mat[17*tree->mod->ns+ 6] =    2 ;
      tree->step_mat[17*tree->mod->ns+ 7] =    2 ;
      tree->step_mat[17*tree->mod->ns+ 8] =    3 ;
      tree->step_mat[17*tree->mod->ns+ 9] =    3 ;
      tree->step_mat[17*tree->mod->ns+10] =    2 ;
      tree->step_mat[17*tree->mod->ns+11] =    2 ;
      tree->step_mat[17*tree->mod->ns+12] =    2 ;
      tree->step_mat[17*tree->mod->ns+13] =    2 ;
      tree->step_mat[17*tree->mod->ns+14] =    3 ;
      tree->step_mat[17*tree->mod->ns+15] =    2 ;
      tree->step_mat[17*tree->mod->ns+16] =    3 ;
      tree->step_mat[17*tree->mod->ns+17] =    0 ;
      tree->step_mat[17*tree->mod->ns+18] =    2 ;
      tree->step_mat[17*tree->mod->ns+19] =    3 ;
      tree->step_mat[18*tree->mod->ns+ 0] =    3 ;
      tree->step_mat[18*tree->mod->ns+ 1] =    3 ;
      tree->step_mat[18*tree->mod->ns+ 2] =    2 ;
      tree->step_mat[18*tree->mod->ns+ 3] =    2 ;
      tree->step_mat[18*tree->mod->ns+ 4] =    2 ;
      tree->step_mat[18*tree->mod->ns+ 5] =    2 ;
      tree->step_mat[18*tree->mod->ns+ 6] =    2 ;
      tree->step_mat[18*tree->mod->ns+ 7] =    3 ;
      tree->step_mat[18*tree->mod->ns+ 8] =    2 ;
      tree->step_mat[18*tree->mod->ns+ 9] =    3 ;
      tree->step_mat[18*tree->mod->ns+10] =    2 ;
      tree->step_mat[18*tree->mod->ns+11] =    2 ;
      tree->step_mat[18*tree->mod->ns+12] =    3 ;
      tree->step_mat[18*tree->mod->ns+13] =    2 ;
      tree->step_mat[18*tree->mod->ns+14] =    3 ;
      tree->step_mat[18*tree->mod->ns+15] =    2 ;
      tree->step_mat[18*tree->mod->ns+16] =    3 ;
      tree->step_mat[18*tree->mod->ns+17] =    2 ;
      tree->step_mat[18*tree->mod->ns+18] =    0 ;
      tree->step_mat[18*tree->mod->ns+19] =    3 ;
      tree->step_mat[19*tree->mod->ns+ 0] =    2 ;
      tree->step_mat[19*tree->mod->ns+ 1] =    3 ;
      tree->step_mat[19*tree->mod->ns+ 2] =    3 ;
      tree->step_mat[19*tree->mod->ns+ 3] =    2 ;
      tree->step_mat[19*tree->mod->ns+ 4] =    3 ;
      tree->step_mat[19*tree->mod->ns+ 5] =    3 ;
      tree->step_mat[19*tree->mod->ns+ 6] =    2 ;
      tree->step_mat[19*tree->mod->ns+ 7] =    2 ;
      tree->step_mat[19*tree->mod->ns+ 8] =    3 ;
      tree->step_mat[19*tree->mod->ns+ 9] =    2 ;
      tree->step_mat[19*tree->mod->ns+10] =    2 ;
      tree->step_mat[19*tree->mod->ns+11] =    3 ;
      tree->step_mat[19*tree->mod->ns+12] =    2 ;
      tree->step_mat[19*tree->mod->ns+13] =    2 ;
      tree->step_mat[19*tree->mod->ns+14] =    3 ;
      tree->step_mat[19*tree->mod->ns+15] =    3 ;
      tree->step_mat[19*tree->mod->ns+16] =    3 ;
      tree->step_mat[19*tree->mod->ns+17] =    3 ;
      tree->step_mat[19*tree->mod->ns+18] =    3 ;
      tree->step_mat[19*tree->mod->ns+19] =    0 ;
    }
  else
    {
      tree->step_mat[0*tree->mod->ns+0] = 0;
      tree->step_mat[0*tree->mod->ns+1] = 2;
      tree->step_mat[0*tree->mod->ns+2] = 1;
      tree->step_mat[0*tree->mod->ns+3] = 2;

      tree->step_mat[1*tree->mod->ns+0] = 2;
      tree->step_mat[1*tree->mod->ns+1] = 0;
      tree->step_mat[1*tree->mod->ns+2] = 2;
      tree->step_mat[1*tree->mod->ns+3] = 1;

      tree->step_mat[2*tree->mod->ns+0] = 1;
      tree->step_mat[2*tree->mod->ns+1] = 2;
      tree->step_mat[2*tree->mod->ns+2] = 0;
      tree->step_mat[2*tree->mod->ns+3] = 2;

      tree->step_mat[3*tree->mod->ns+0] = 2;
      tree->step_mat[3*tree->mod->ns+1] = 1;
      tree->step_mat[3*tree->mod->ns+2] = 2;
      tree->step_mat[3*tree->mod->ns+3] = 0;

/*       tree->step_mat[0*tree->mod->ns+0] = 0; */
/*       tree->step_mat[0*tree->mod->ns+1] = 1; */
/*       tree->step_mat[0*tree->mod->ns+2] = 1; */
/*       tree->step_mat[0*tree->mod->ns+3] = 1; */

/*       tree->step_mat[1*tree->mod->ns+0] = 1; */
/*       tree->step_mat[1*tree->mod->ns+1] = 0; */
/*       tree->step_mat[1*tree->mod->ns+2] = 1; */
/*       tree->step_mat[1*tree->mod->ns+3] = 1; */

/*       tree->step_mat[2*tree->mod->ns+0] = 1; */
/*       tree->step_mat[2*tree->mod->ns+1] = 1; */
/*       tree->step_mat[2*tree->mod->ns+2] = 0; */
/*       tree->step_mat[2*tree->mod->ns+3] = 1; */

/*       tree->step_mat[3*tree->mod->ns+0] = 1; */
/*       tree->step_mat[3*tree->mod->ns+1] = 1; */
/*       tree->step_mat[3*tree->mod->ns+2] = 1; */
/*       tree->step_mat[3*tree->mod->ns+3] = 0; */

    }
  
  For(i,tree->mod->ns) tree->step_mat[i*tree->mod->ns+i] = 0;
}


/*********************************************************/
