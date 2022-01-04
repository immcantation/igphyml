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
#include "utilities.h"

extern int     stopCodons[64];
extern int   senseCodons[64];
extern char aminoAcidmap[65];
extern int indexSenseCodons[64];



/********************************************************
 * Resolve polytomies within tree if desired
 * */

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
	  	  Init_Class_Tips(tree,precon); //initialize tree structures
	  	  int pars1 = Fill_Sankoff(r,tree,1); //fill Sankoff matrixes
	  	  Set_Pars_Counters(r,tree,1); //set initial parsimony counters
	  	  if(tree->mod->optDebug)printf("\n. Resolving polytomies using isotype information");
	  	  io->thresh = 0.001;
	  	  tree->charindex = mCalloc(tree->nstate,sizeof(int));
	  	  tree->polytomy_states = mCalloc((tree->n_otu-1)*2,sizeof(int));
	  	  tree->polytomy_swaps = mCalloc((tree->n_otu-1)*2,sizeof(int*));
	  	  if(io->mod->polytomyresolve >= 1){
	  		  int pars2 = Resolve_Polytomies_Pars(tree,io->thresh);
	  	  	  #if defined OMP || defined BLAS_OMP
		  	  #pragma omp critical
		  	  #endif
	  		  {
	  			  printf("\n. %d Initial/resolved maximum parsimony score: %d %d %s",tree->mod->num,pars1,pars2,tree->mod->rootname);
	  		  }
        }
	  }
}

/********************************************************
 * Input: full data io object
 * Perform Sankoff parsimony reconstruction on set of trees
 * Calculate the number of switches in each direction for
 * each tree.
 * */
void Pars_Reconstructions(option* io){
	int j,i,k;
  	int sample=io->parssample; // Number of paths to sample down each tree

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
          int nedges = (tree->n_otu-1)*2;
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

  		  //array of all possible reconstructions (cut)
  		  t_node* r = tree->noeud[tree->mod->startnode]; //root node
  		  Init_Class_Tips(tree,io->precon); //initialize tip states and data structures
  		  int pars = Fill_Sankoff(r,tree,1); //max parsimony score
  		  Set_Pars_Counters(r,tree,1); //set counters to their first minimum position, possibly not necessary

  		  //Data structures for counting the number of type switches and branch lengths in each state
          tree->polytomy_states = mCalloc((tree->n_otu-1)*2,sizeof(int)); //document!
          tree->polytomy_swaps = mCalloc((tree->n_otu-1)*2,sizeof(int*));
  	      phydbl* switches = mCalloc(tree->nstate*tree->nstate,sizeof(phydbl));
  	      phydbl* classl = mCalloc(tree->nstate,sizeof(phydbl));
          tree->charindex = mCalloc(tree->nstate,sizeof(int));

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

  		  // set up tree output file
  		  int npars = sample;
  	  	  FILE* treeout1 = Openfile(fout, 1 );
  		  fprintf(treeout1,"#NEXUS\nBegin taxa;\nDimensions ntax=%d;\nTaxlabels\n",tree->n_otu);
  		  For(i,nedges){
  		  	  if(tree->noeud[i]->tax){
  		  		  fprintf(treeout1,"%s\n",tree->noeud[i]->name);
  		  	  }
  		  }
  		  fprintf(treeout1,";\nEnd\nBegin trees;\n");

  		  // make a copy of the original tree, only copies tree structure
          t_tree* tree2 = Make_Tree_From_Scratch(tree->n_otu,tree->data);
          Copy_Tree(tree,tree2);
          tree2->chars = tree->chars;
          tree2->mod = tree->mod;
          tree2->io = tree->io;

          // Set state at each node to its "leftmost" max parsimony state
          Set_Pars_Counters(r,tree,1);
          Get_First_Path(tree->noeud[tree->mod->startnode], minrootsar[0], tree,1);

          // go through all nodes, if an unresolved polytomy is found, resolve it and repeat loop
          // continue until no unresolved polytomies are found
          if(tree->mod->polytomyresolve >= 2){ //if maximum ambiguity resolution desired
              int iters = 0;
              int nnifound = 1;
              while(nnifound){ //execute if a rearrangement is found in a tree
            	  nnifound = 0;
            	  For(k,nedges){ //if no taxa on either side of edge and length is below threshold, search for NNIs
            		  t_node* node = tree->noeud[k];
            		  // if the node is internal, not the MRCA, the ancestral edge is < threshold, and not previously resolved
            		  if(!node->tax && !node->anc->tax && node->anc_edge->l < tree->io->thresh && node->polytomy == 0){
            			  Resolve_Polytomy_Mono(node, tree->io->thresh, 1, 0, tree);
            			  nnifound++;
            			  if(nnifound){
            				  break;
            			  }
            		  }
            	  }
            	  iters++;
              }
          }
          // reset polytomy indicators for all nodes
          For(k,(tree->n_otu-1)*2)tree->noeud[k]->polytomy = 0;
          // Loop through all sampling repetitions
          For(i,sample){
        	if(tree->mod->optDebug)printf("\n Repetition %d",i);
        	// print to tree nexus file
            fprintf(treeout1,"Tree TREE%d = [&R] ",i+1);

        	//randomly choose a root state
	  	  	int n = rand() % minroots;
	  	  	int minroot = minrootsar[n];
	  	  	For(k,(tree->n_otu-1)*2)tree->noeud[k]->polytomy = 0;

	  	  	// if resolved polytomies desired, randomly swap the order of top
	  	  	// nodes in each polytomy
            if(tree->mod->polytomyresolve >= 2){
            	int poly;
            	For(poly,tree->polytomies){
            		int states = tree->polytomy_states[poly]; //number of states in this polytomy
            		if(states > 2){
            			int* swaps = tree->polytomy_swaps[poly];
            			// randomly swap two top non-adjacent nodes representing different states
            			int f = rand() % states;
            			int s;
            			t_node* first = tree->noeud[swaps[f]];
            			t_node* second = first;
            			while(first->anc->num == second->anc->num){
            				s = rand() % states;
            				second = tree->noeud[swaps[s]];
            			}
            			if(tree->mod->optDebug)printf("\nswapping %d %d %d %d",first->num,second->num,f,s);
            			Swap(first,first->anc,second->anc,second,tree);
            		}
            	}
              }
              // Fill in Sankoff structures with newly-swapped trees
              Fill_Sankoff(tree->noeud[tree->mod->startnode],tree,1);

              Get_Rand_Path(r,minroot,tree,1);

  			  Fill_Pars_Stats(r,tree, switches,classl,1);

              For(k,(tree->n_otu-1)*2){
                //non-polytomy nodes keep their numbers after polytomy resolution 
                tree->noeud[k]->polytomy = 0;
                tree2->noeud[k]->pstate = tree->noeud[k]->pstate;
              }
  			  io->precon *= 10;
  			  char* ts = Write_Tree(tree2);
  			  io->precon /= 10;
  			  ts[strlen(ts)-1] = 0;
  			  strcat(ts,"[");
  			  strcat(ts,tree2->chars[tree2->noeud[tree2->mod->startnode]->pstate]);
  			  strcat(ts,"]:0;");
  			  fprintf(treeout1,"%s\n",ts);
  			  free(ts);
          } //For each repetition
  		  fprintf(treeout1,"END;\n");
  		  fclose(treeout1);
          Clean_Tree(tree2);
          Free_Tree(tree2);

          // Divide switch counts by the number of repetitions
  		  int tposi, tposj;
          phydbl total = 0.0;
  		  For(tposi,tree->nstate){
  		 	classl[tposi] = classl[tposi]/(npars*1.0);
  		 	For(tposj,tree->nstate){
  		 		if(tree->mod->optDebug)printf("%lf\t%d\n",switches[tposi*tree->nstate+tposj],npars);
  		 		switches[tposi*tree->nstate+tposj] = switches[tposi*tree->nstate+tposj]/(npars*1.0);
  		 		total += switches[tposi*tree->nstate+tposj];
  		 	}
  		 }
  		 if(tree->mod->optDebug)printf("\nTOTAL %lf\n",total);
  		 phydbl tswitch, tlen;
  		 tswitch=tlen=0;

  		 // Print output file
  		 For(tposi,tree->nstate){ //"time" in each state
  		 	fprintf(pstatf,"%d\t%s\tN\t%lf\n",j,tree->chars[tposi],classl[tposi]);
  		 	tlen+=classl[tposi];
  		 }
  		 fprintf(pstatf,"%d\tUCA\t%s",j,tree->chars[minrootsar[0]]); //what were the possible UCA states
  		 for(i = 1; i < minroots; i++)fprintf(pstatf,":%s",tree->chars[minrootsar[i]]);
  		 fprintf(pstatf,"\t0.0\n");
  		 fprintf(pstatf,"%d\tNTIP\tNTIP\t%d\n",j,tree->n_otu); //number of tips
  		 	For(tposi,tree->nstate){
  		 		For(tposj,tree->nstate){ //swtiches from each state to each other state
  		 			fprintf(pstatf,"%d\t%s\t%s\t%lf\n",j,tree->chars[tposi],tree->chars[tposj],switches[tposi*tree->nstate+tposj]);
  		 				tswitch+=switches[tposi*tree->nstate+tposj];
  		 			}
  		 		}
  		  fclose(io->mod_s[j]->fp_in_tree);
  		  Clean_Tree(tree);
  		  Free_Tree(tree);
  	}
}


/********************************************************
 * Recurse down tree, randomly choosing ambiguous pointers.
 * Quite simple unless polytomy skipping is desired.
 *
 * If polytomy skipping is desired (polytomyresolve=2):
 * If you enter an unmarked polytomy node
 * recurse through polytomy using Get_Rand_Path_Polytomy.
 * Get_Rand_Path_Polytomy will:
 * 1. Assign pstate and pancstate to all nodes in polytomy
 * 2. Mark all polytomy nodes with polytomy=1
 * 3. Record the state of polytomy "tip" nodes
 *
 * Once the above is complete, recurse through polytomy again but don't re-assign states.
 * When a descendant node is not in the polytomy, find all states at the current node
 * that are states represented by the polytomy and have an available path with
 * the same number of switches as the current node state. Randomly choose among these.
 * This effectively allows the algorithm to skip polytomies while maintaining the same
 * number of switches along the tree. This can slightly alter the distribution of switches
 * among states though, as it will disfavor switches from states that are ancestral to, but
 * not part of, a polytomy.
 *
 * d is current node
 * index is index of state to be assigned at that node (pstate)
 * tree is tree object
 * root indicates whether this is the root node of the tree
 * */
void Get_Rand_Path(t_node *d, int index, t_tree *tree, int root){
  int i,j,dir1,dir2;
  int ldraw=0;
  int rdraw=0;
  int lfound=0;
  int rfound=0;
  //if not a polytomy, set state to index
  if(tree->mod->optDebug)printf("\nGet_Rand_Path %d %d",d->num, index);

  if(d->polytomy == 0){ //if node was already visited by Get_Rand_Path_Polytomy, state already assigned
    d->pstate=index;
  }
  if(!d->tax && tree->mod->optDebug)printf("\n%s\t%d",tree->chars[index], index);

  if(!d->tax || root){ //if not at a tip
    if(!root){
    	// find ancestral and descendant directions
    	t_node* a = d->anc;
    	dir1=dir2=-1;
    	For(i,3) if(d->v[i]->num != a->num) (dir1<0)?(dir1=i):(dir2=i);
    	// if d is top of a polytomy, set node states and tally up states in the polytomy
    	if((d->polytomy == 0) && (d->b[dir1]->l < tree->io->thresh ||
          d->b[dir2]->l < tree->io->thresh) && tree->mod->polytomyresolve == 2){
    		int* scores = mCalloc(tree->nstate,sizeof(int));
        	scores[d->pancstate]++; //include ancestral state as a possibility
        	Get_Rand_Path_Polytomy(d, scores, index, tree, 1);
        	d->polytomy = 1;
    	}
    }else{
      dir1=0;
    }

    if(d->polytomy > 0 && tree->mod->polytomyresolve != 2){
    	Warn_And_Exit("Polytomy marked with invalid specification.");
    }
    int* scores;
    if(d->polytomy == 0){ //make scores table with 1 at index..
      scores = mCalloc(tree->nstate,sizeof(int));
      scores[index] = 1;
    }else{ // if this is a polytomy, allow state to be any state in polytomy
      scores = d->polystates;
      if(tree->mod->optDebug){
    	  printf("\npolystates %d %s\n",d->num, tree->chars[d->pstate]);
    	  printf("\n");
    	  For(i,tree->nstate)printf("%d\t",scores[i]);
    	  printf("\n");
    	  For(i,tree->nstate)printf("%s\t",tree->chars[i]);
    	  printf("\n%d mins %d %d",d->num, d->lmin[d->pstate],d->rmin[d->pstate]);
      }
    }
    // Given the state(s) at current node (scores) find all states
    // at the left and right nodes that require the same number of switches
    // as the minimum of the state at the current node
    int lmins=0; int rmins=0;
    For(j,tree->nstate){ //allow for any state in a polytomy to be used as the ancestor of the descendant node
    	if(scores[j] > 0){
    		For(i,tree->nstate){ //pstate used because it can differ from index if polytomy == 0
    			if(d->pl[j*tree->nstate+i] == d->lmin[d->pstate])lmins++;
        		if(!root)if(d->pr[j*tree->nstate+i] == d->rmin[d->pstate])rmins++;
    		}
    	}
    }
    //randomly select left and right states
    if(lmins>1)ldraw = rand() % lmins; //not perfectly random but probably okay
    if(!root && rmins>1)rdraw = rand() % rmins;

    For(j,tree->nstate){
      if(scores[j] > 0){ //if not polytomy, j is just index
        For(i,tree->nstate){
        	//if j is a possible state at this node and j to i is minimal
        	if(d->pl[j*tree->nstate+i] == d->lmin[d->pstate]){
        		if(lfound == ldraw){ // if this corresponds to randomly drawn index
        			//if descendant node not part of polytomy, set its ancestral state to j, so panc can be different from anc->pstate
        			//otherwise if descendant node is in polytomy, this will have already been set
        			if(d->v[dir1]->polytomy==0){
        				d->v[dir1]->pancstate = j;
        			}
        			Get_Rand_Path(d->v[dir1],i,tree,0);
        		}
        		lfound++;
        	}
        	if(d->pr[j*tree->nstate+i] == d->rmin[d->pstate] && !root ){
        		if(rfound == rdraw){
        			if(d->v[dir2]->polytomy==0){
        				d->v[dir2]->pancstate = j;
        			}
        			Get_Rand_Path(d->v[dir2],i,tree,0);
        		}
        	 rfound++;
        	}
        }
      }
    }
  }else{
    if(tree->mod->optDebug)printf("%s\t%s\n",tree->chars[index],d->name);
  }
}

/********************************************************
 * Recurse down polytomy, randomly choosing ambiguous pointers
 * Same as Get_Rand_Path but only operates on nodes within and adjacent
 * to a polytomy. Scores should be an empty vector that is stored in each polytomy node
 * in order to represent the possible states of that polytomy
 * */
void Get_Rand_Path_Polytomy(t_node *d, int* scores, int index, t_tree *tree, int top){
  int i,j,dir1,dir2;
  int ldraw=0;
  int rdraw=0;
  int lfound=0;
  int rfound=0;
  int root = 0;
  if(d->num == tree->mod->startnode)root=1;
  if(tree->mod->optDebug)printf("\nGet_Rand_Path_Polytomy %d %d",d->num, index);
  d->pstate=index; //set state to index
  d->polystates = scores;
  // stop if we reach the end of the polytomy
  if((!d->tax && d->anc_edge->l < tree->io->thresh) || top){
    d->polytomy = 1; //identify d as a polytomy
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
    if(lmins>1)ldraw = rand() % lmins;
    if(!root && rmins>1)rdraw = rand() % rmins;

    For(i,tree->nstate){
        if(d->pl[index*tree->nstate+i] == d->lmin[index]){
          if(lfound == ldraw){
            d->v[dir1]->pancstate = index;
            Get_Rand_Path_Polytomy(d->v[dir1], scores, i, tree, 0);
          }
          lfound++;
        }
        if(d->pr[index*tree->nstate+i] == d->rmin[index] && !root ){
          if(rfound == rdraw){
            d->v[dir2]->pancstate = index;
            Get_Rand_Path_Polytomy(d->v[dir2], scores, i, tree, 0);
          }
          rfound++;
        }
    }
  }else{
      scores[d->pstate]++; //scores contains the number of polytomy bottoms with a particular state
      if(tree->mod->optDebug)printf("%s\t%s\n",tree->chars[index],d->name);
  }
}

/********************************************************
 * Recursively add values to switch and length matrices for each type
 * t_node d is current node
 * tree is tree object
 * switches is an array holding the number of switches from each state to the other
 * classl is the branch lengths of the tree occupied by a given state
 * root indicates whether we're at the root node
 * */
void Fill_Pars_Stats(t_node* d, t_tree* tree, phydbl* switches, phydbl* classl, int root){
  int i,j,dir1,dir2;
  dir1=dir2=-1;
  if(!root){
    if(d->pstate != d->pancstate){
      if(tree->mod->optDebug && d->anc->pstate != d->pancstate){
    	 printf("\nanc switch! %d\t%s\t%s\t%s",d->num,tree->chars[d->anc->pstate],tree->chars[d->pancstate],tree->chars[d->pstate]);
      }
      if(d->anc->pstate != d->pancstate && tree->mod->polytomyresolve != 2){
    	  Warn_And_Exit("Polytomy state inconsistency with invalid specification (polytomyresolve != 2).");
      }
      switches[d->pancstate*tree->nstate+d->pstate]++;
      classl[d->pstate] += d->anc_edge->l*0.5; //split switched branch lengths
      classl[d->pancstate] += d->anc_edge->l*0.5;
    }else{
      classl[d->pstate] += d->anc_edge->l;
    }
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
 * Count all possible switches along polytomies
 *
void Count_Polytomy_Switches(t_node* d, phydbl* switches, phydbl thresh, t_tree* tree){
  int i,j;
  printf("counting polytomy switches\n");
  t_node* top = d;
  while(top->anc_edge->l < thresh && top->anc->num != tree->mod->startnode){
    top = top->anc;
  }
  //printf("top %d\n",top->num);
  int ancstate = top->anc->pstate;
  if(top->anc->polytomy != 0)ancstate = top->pstate;
  //printf("anc state %s\n",tree->chars[ancstate]);
  Count_Polytomy_States(top, NULL, 1, thresh, 0, tree);
  int* scores = d->polystates;
 // For(i,tree->nstate)printf("%d\t",d->polystates[i]);
  //printf("\n");
  phydbl nswitches = 0;
  phydbl nsteps = 0; 
  For(i,tree->nstate){
    if(scores[i] > 0 && i != ancstate)nsteps++;
    For(j,tree->nstate){
      if(scores[i] > 0 && scores[j] > 0 && i != j){
        if(tree->step_mat[i*tree->nstate+j] <= 1 && j != ancstate){
          printf("switch counts %s\t%s\t%lf\n",tree->chars[i],tree->chars[j],nsteps/nswitches);
          nswitches++;
        }
      }
    }
  }
  For(i,tree->nstate){
    For(j,tree->nstate){
      if(scores[i] > 0 && scores[j] > 0 && i != j){
        if(tree->step_mat[i*tree->nstate+j] <= 1 && j != ancstate){
          printf("switch counts %s\t%s\t%lf\n",tree->chars[i],tree->chars[j],nsteps/nswitches);
          switches[i*tree->nstate+j] += nsteps/nswitches;
        }
      }
    }
  }
}
*/

/*********************************************************
 * Initialize Sankoff dynamic programming tables at the tips of the tree
 */
void Init_Class_Tips(t_tree* tree, int precon){
	 int i,j;
	 char* mtemp;
	 if(precon == 7){
		 Setup_Custom_Pars_Model(tree);
	 }else{
		 Warn_And_Exit("Must use precon=7, other options no longer supported after v1.1.3");
	 }
}

/*********************************************************
* Permute tip names except for germline
*
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
    			//printf("%s\n",tree->noeud[i]->name);
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
void Permute_Tips(t_tree* tree){
	int i, j;
	int indexes[tree->n_otu-1];

	int count = 0;
	int n = tree->n_otu;
	char temp[T_MAX_NAME];
    For(i,n - 1) {
    	if(strcmp(tree->noeud[i]->name,tree->mod->rootname) != 0){
    		j = i + rand() / (RAND_MAX / (n - i) + 1);
    		// re-sample if you select the root
    		while(strcmp(tree->noeud[j]->name,tree->mod->rootname) == 0)
    			j = i + rand() / (RAND_MAX / (n - i) + 1);
    		strcpy(temp,tree->noeud[i]->name);
    		strcpy(tree->noeud[i]->name,tree->noeud[j]->name);
    		strcpy(tree->noeud[j]->name,temp);
    	}
    }
}

/*********************************************************
* Permute tip names except for germline
*/
void Permute_All_MetaData(option* io, int pos){
	int i, j;
	int n_otu = 0;
	For(i,io->ntrees){
		n_otu += io->tree_s[i]->n_otu-1;
	}
	int* indexes = mCalloc(n_otu,sizeof(int));
	char** names = mCalloc(n_otu,sizeof(char*));

	//make big list of names with similar list of indexes
	int index = 0;
	For(i,io->ntrees){
		For(j,io->tree_s[i]->n_otu){
			if(strcmp(io->tree_s[i]->noeud[j]->name,io->mod_s[i]->rootname) != 0){
				indexes[index] = index;
				names[index] = mCalloc(T_MAX_NAME,sizeof(char));
				strcpy(names[index],io->tree_s[i]->noeud[j]->name);
				index++;
			}
		}
	}
	//shuffle big list of indexes
	int n = n_otu;
	int temp;
	For(i,n - 1){
    	j = i + rand() / (RAND_MAX / (n - i) + 1);
    	temp = indexes[i];
    	indexes[i] = indexes[j];
    	indexes[j] = temp;
    }

	//re-assign names baed on shuffled indexes
	index = 0;
	For(i,io->ntrees){
		For(j,io->tree_s[i]->n_otu){
			if(strcmp(io->tree_s[i]->noeud[j]->name,io->mod_s[i]->rootname) != 0){
				strcpy(io->tree_s[i]->noeud[j]->name,names[indexes[index++]]);
			}
		}
	}

    free(indexes);
    For(i,n_otu)free(names[i]);
	free(names);
}


/*********************************************************
* read in parimsony model
* Input is a tree object
* File name is stored in preconfile variable
*/
void Setup_Custom_Pars_Model(t_tree* tree){
	int i,j;
	char* mtemp;
	FILE* PARS = Openfile(tree->mod->preconfile,0);
	char* line = mCalloc(T_MAX_LINE,sizeof(char));
	int fscn;
	do{//skip over beginning comments
		fscn = fscanf(PARS, "%s\n",line);
	}while(strcmp(line,"#BEGIN")!=0 && fscn != EOF);

	//read in number of states
	fscn = fscanf(PARS, "%d\n",&tree->nstate);
	tree->chars = mCalloc(tree->nstate,sizeof(char*));

	 //associative arrays holding all possible states mapped back to the states of the model
	 //ambigfrom: nstate x max option
	 char** ambigstatesfrom = mCalloc(tree->nstate*T_MAX_OPTION,sizeof(char*));
	 For(i,tree->nstate*T_MAX_OPTION)ambigstatesfrom[i]=mCalloc(T_MAX_OPTION,sizeof(char*));
	 int* ambigstatesto = mCalloc(tree->nstate*T_MAX_OPTION,sizeof(int));

	 // Read states into tree->chars object
	 // They will be refered to by their index in this array
	 int statecount = 0;
	 fscn = fscanf(PARS, "%s\n",line);
	 For(i,tree->nstate){
		 tree->chars[i]=mCalloc(T_MAX_OPTION,sizeof(char));
		 fscn=fscanf(PARS, "%s\n",tree->chars[i]);
		 strcpy(ambigstatesfrom[i],tree->chars[i]);
		 ambigstatesto[i]=i;
		 statecount++;
		 if(tree->mod->optDebug)printf("state %d %s\n",i,tree->chars[i]);
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
	 	For(j,tree->nstate){ //initialize states and scores to maximum
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
			 if(tree->mod->optDebug)printf("LINE %s %s %d\n",from,to,val);
			 if(strcmp(from,"#AMBIGUOUS")!=0){
			 For(i,tree->nstate){ //set up step mat
				 For(j,tree->nstate){
					 if(strcmp(tree->chars[i],from)==0 && strcmp(tree->chars[j],to)==0){
						 if(tree->mod->optDebug)printf("%s\t%s\t%d\n",tree->chars[i],tree->chars[j],tree->step_mat[i*tree->nstate+j]);
						 tree->step_mat[i*tree->nstate+j]=val;
					 }
				 }
			 }
			 }
		 }while(strcmp(from,"#AMBIGUOUS")!=0 && fscn != EOF);
		 rewind(PARS);
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
				 if(tree->mod->optDebug)printf("ambiguous %s %s\n",from,to);
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
		 int found = 0;
		 For(j,statecount){ //mark tip as also being in the state it's ambiguous with
			 if(strcmp(state,ambigstatesfrom[j])==0){
				 tree->noeud[i]->s[ambigstatesto[j]]=0;
				 found++;
			 }
		 }
		 if(!found){
			 printf("\nState %s at node %d from sequence %s not found in model.\n",state,i,tree->noeud[i]->name);
			 Warn_And_Exit("");
		 }
		 For(j,tree->nstate){
			 if(tree->mod->optDebug)printf(" %d",tree->noeud[i]->s[j]);
		 }
		 if(tree->mod->optDebug)printf("\n");
		 free(state);
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
		 Clean_Sankoff_Node(tree->noeud[i]);
	 }
}

void Clean_Sankoff_Node(t_node* node){
     free(node->pl);
     free(node->pr);
     free(node->s);
     free(node->sroot);
     free(node->lmin);
     free(node->rmin);
     free(node->prc);
     free(node->plc);
     free(node->llock);
     free(node->rlock);
     //free(node->polystates);
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

/********************************************************
 * Rearrange polytomies to a maximum parsimony state using NNIs
 * */
int Resolve_Polytomies_Pars(t_tree* tree, phydbl thresh){
	int i,j;
	int nni=1;
	int nnifound=1; //initial parsimony score
	int pars0 = Fill_Sankoff(tree->noeud[tree->mod->startnode],tree,1); //get starting score
	int iters = 0;
	int nedges = (tree->n_otu-1)*2;
	while(nnifound){ //execute if a rearrangement is found in a tree
		nnifound=0;
		For(i,nedges){ //if no taxa on either side of edge and length is below threshold, search for NNIs
			if(!tree->noeud[i]->tax && !tree->noeud[i]->anc->tax && tree->noeud[i]->anc_edge->l < thresh){
				pars0 = Fill_Sankoff(tree->noeud[tree->mod->startnode],tree,1);
				t_edge* b = tree->noeud[i]->anc_edge;
				// Look for potentially favorable NNIs on either side of the edge
				nni =         NNI_Pars_Search(tree->noeud[i],tree->noeud[i]->anc,
						tree->noeud[i]->anc_edge,tree->noeud[i]->anc_edge,pars0,thresh,tree);
				if(!nni)nni = NNI_Pars_Search(tree->noeud[i]->anc,tree->noeud[i],
						tree->noeud[i]->anc_edge,tree->noeud[i]->anc_edge,pars0,thresh,tree);
				nnifound+=nni;
			}
		}
		iters++;
	}

	// if a maximum ambiguity resolution desired, do that
	if(tree->mod->polytomyresolve >= 2){
		t_node* r = tree->noeud[tree->mod->startnode];
		int smin = INT_MAX;
		int mins = 0;
		For(i,tree->nstate){
			if(r->sroot[i] < smin){
				smin = r->sroot[i];
				mins = i;
			}
		}
		Set_Pars_Counters(r,tree,1);
		Get_First_Path(tree->noeud[tree->mod->startnode],mins,tree,1);
		iters = 0;
		nnifound = 1;
		while(nnifound){ //execute if a rearrangement is found in a tree
			nnifound=0;
			For(i,nedges){ //if no taxa on either side of edge and length is below threshold, search for NNIs
				t_node* node = tree->noeud[i];
				if(!node->tax && !node->anc->tax && node->anc_edge->l < thresh && node->polytomy == 0){
					Resolve_Polytomy_Mono(node,thresh,1,0,tree);
					nnifound++;
					if(nnifound){
						break;
					}
				}
			}
		iters++;
		}
		pars0 = Fill_Sankoff(tree->noeud[tree->mod->startnode],tree,1);
	}
	return pars0;
}


/********************************************************
 * Resolve polytomy to maximum ambiguity configuration
 * input is node of interest, threshold branch length, whether to randomize,
 * the level of the polytomy relative to the polytomy root
 * and the tree object.
 * From Pars_Reconstruction, default is level=0, randomize=1
 * */
t_node* Resolve_Polytomy_Mono(t_node* b, phydbl thresh, int randomize, int level, t_tree* tree){
	int i,j;
	int nedges = (tree->n_otu-1)*2;
	// Recursively move up to the top node of the polytomy
	t_node* top = b;
	while(top->anc->num != tree->mod->startnode && top->anc_edge->l < thresh){
		top = top->anc;
	}
	t_node* anc = top->anc;

	Fill_Sankoff(tree->noeud[tree->mod->startnode], tree, 1); //fill parsimony objects (cut?)
	int pscore = Score_Polytomy(top, thresh, 1000, 0, tree); //get initial score of polytomy (see description)

	if(tree->mod->optDebug){ //may be overkill
		tree->io->precon *= 10;
		printf("%s\n",Write_Tree(tree));
		tree->io->precon /= 10;
	}

	// Get array of nodes at the bottom of the polytomy, collect other info as well
	t_node** nodes = mCalloc(nedges,sizeof(t_node*)); //array of nodes at bottom of polytomy
	t_node** tops = mCalloc(tree->nstate,sizeof(t_node*)); //array of top polytomy nodes
	int* scores = mCalloc(tree->nstate,sizeof(int)); //number of bottom nodes with each state
	int* nums = malloc(nedges*sizeof(int)); //numbers of internal polytomy nodes
    int* edges = malloc(nedges*sizeof(int)); //numbers of internal polytomy edges
	int* numindex = malloc(sizeof(int)); //current index of nums
    int* edgeindex = malloc(sizeof(int)); //current index of edges
	*numindex = 0;
    *edgeindex = 0;
	int nnodes = Prune_Polytomy(top, nodes, scores, nums, edges, 0, numindex, edgeindex, 1, thresh, tree);

	// Count number of unique states in bottom polytomy nodes
	int nzero = 0;
	For(i,tree->nstate)if(scores[i] > 0)nzero++;

	//if polytomy only two nodes, skip rest of function
	if(nnodes < 3 || nzero == 1) return top;

	int n = tree->nstate;
    if(randomize){ //randomize char index order so they are processed randomly below
    	For(i,tree->nstate){ //Need to be out of this if statement?
    		tree->charindex[i]=i;
    	}
    	int temp;
    	j = 0;
    	For(i,n - 1) {
    		j = i + rand() / (RAND_MAX / (n - i) + 1);
    		temp = tree->charindex[i];
    		tree->charindex[i] = tree->charindex[j];
    		tree->charindex[j] = temp;
    	}
    }
    //For(i,tree->nstate)printf("states %d\t%d\n",i,states[i]);
    // Loop through each state (randomly if desired) and join nodes with that
    // state together in an asymmetric manner. Produce array of "top" nodes
    // that lead to lower single-type polytomies
	int node;
	int tindex = -1;
	int nindex = nedges+1; //new nodes get temporary indexes > current index range
	For(n,tree->nstate){ //loop through all states in polytomy
		j = tree->charindex[n];
		if(scores[j] == 0)continue;
		For(node, nnodes){ //loop through all bottom nodes
			if(nodes[node]->pstate == j){ //if node has state j
				nodes[node]->polytomy = 1; //identify node as polytomy node
				if(scores[j] > 0){ //if score positive, first node of this state found
					tops[++tindex] = nodes[node]; //make this node a state top
					scores[j] = -1*scores[j]; //make scores negative
				}else{ //otherwise, add this node to its respective top
					tops[tindex] = Join_Nodes(tops[tindex], nodes[node], nindex+=3);
				}
			}
		}
	}

	// Make a list of node numbers of the top nodes of this polytomy
	int* swaplist = mCalloc(tindex+1,sizeof(int));
	tree->polytomy_states[tree->polytomies] = tindex+1;
	For(n,tindex+1)swaplist[n] = tops[n]->num;
	tree->polytomy_swaps[tree->polytomies] = swaplist;

	// Join top nodes in a balanced manner, get new polytomy top
	Join_Nodes_Balanced(tops,tindex+1,nindex);
	t_node* newtop = tops[0];

	// attach new polytomy top to the anecstor of the old polytomy top
	For(i,3){
		For(j,3){
			if(anc->v[i]->num == top->num && newtop->b[j]->num == newtop->anc_edge->num){
				Attach_Edge(anc, newtop, i, j, anc->b[i]->l);
				i=3;j=3;
				break;
			}
		}
	}
	if(tree->mod->optDebug)printf("\npolytomy top %d\t%d\t%d",top->num,newtop->num,nums[0]);

	// Replace new node and edge numbers with the numbers of the original nodes
	// also delete old nodes to prevent memory leakage
	Fix_Node_Numbers(newtop, nums, 0, 1, thresh, tree);
	Fix_Edge_Numbers(newtop, edges, 0, 1, thresh, tree);
	
	// Get score of this newly resolved polytomy
	int firstscore = Score_Polytomy(newtop, thresh, 1000, 0, tree);

	// Score_Polytomy will reset the states at the nodes of the tree.
	// If all is well, there should be the same number of each type of node at the bottoms
	// of the polytomy as before. Sometimes there can be differences, though, because sometimes
	// resolving the polytomy changes the underlying reconstructions. If this happens, try to fix it
	int redo = 0;
	int* nscores = mCalloc(tree->nstate,sizeof(int));
	Count_Polytomy_States(newtop,nscores,1,thresh,1,tree);
	For(i,tree->nstate)if(nscores[i] != -scores[i])redo=1;

	// If state change detected, try again
	if(redo && level < 10){
		printf("\nstate change detected at tree %s node %d level %d, recursively fixing\n",tree->mod->rootname,newtop->num,level);
		For(i,tree->nstate)printf("%s\t",tree->chars[i]);
		printf("\n");
		For(i,tree->nstate)printf("%d\t",-scores[i]);
		printf("\n");
		For(i,tree->nstate)printf("%d\t",nscores[i]);
		printf("\n");
		For(i,nnodes-1)tree->noeud[i]->polytomy = 0;
		level++;
		newtop = Resolve_Polytomy_Mono(newtop, thresh, 0, level, tree);
	}
	// if tried for more than 10 tries, give up
	if(level >= 10){
		printf("WARNING! Failed to recursively resolve node! Stopping at 10th iteration");
	}

	int score = Score_Polytomy(newtop,thresh,1000,0,tree);
	int rscore = Score_Polytomy(newtop,thresh,1000,1,tree);

	// update edge directions and ancestral edges
	Update_Dirs(tree);
	Update_Ancestors_Edge(tree->noeud[tree->mod->startnode],
			tree->noeud[tree->mod->startnode]->v[0],
			tree->noeud[tree->mod->startnode]->b[0],tree);

	// index the current polytomy
	tree->polytomies++;

	if(rscore > 0){ //if relative score is > 0, poyltomy resolution has failed :-/
		printf("\n!! Polytomy not full resolved! Will continue if using --force_resolve");
		printf("\nPolytomy at %d, PSCORE: %d, SCORE: %d, RSCORE %d",
				newtop->num, pscore, score, rscore);
		if(tree->mod->optDebug){
			tree->io->precon *= 10;
			printf("%s\n",Write_Tree(tree));
			tree->io->precon /= 10;
		}
		if(!tree->mod->force_resolve){
			Warn_And_Exit("Exiting");
		}
	}
	free(nodes);
	free(tops);
	free(nums);
	free(edges);
	free(scores);
	free(nscores);
	free(edgeindex);
	free(numindex);
	return newtop;
}

/*********************************************************
 * Fixes node numbers to those in the original tree
 * nums contains original numbers of nodes within the polytomy
 * index is the index within that array to use, starts at zero
 */
int Fix_Node_Numbers(t_node* d, int* nums, int index, int root, phydbl thresh, t_tree* tree){
	int i;
	if(!(d->tax || d->anc_edge->l > thresh) || root){ //if a node in the polytomy
		if(tree->mod->optDebug)printf("Fix_Node_Numbers replacing %d %d\n",d->num, nums[index]);
		//if current node is in polytomy swap list, replace it with new number
		For(i,tree->polytomy_states[tree->polytomies]){
			if(tree->polytomy_swaps[tree->polytomies][i] == d->num){
				if(tree->mod->optDebug)printf("Fix_Node_Numbers replacing swap %d %d\n",d->num, nums[index]);
				tree->polytomy_swaps[tree->polytomies][i] = nums[index];
			}
		}
		// re-assign node number then post-increment
		d->num = nums[index++];
		Clean_Sankoff_Node(tree->noeud[d->num]); //delete old node
		Free_Node(tree->noeud[d->num]);
		tree->noeud[d->num] = d; //add new node pointer to noeud array

		//declare Sankoff structures for new node
		d->pl = (int *)mCalloc(tree->nstate*tree->nstate,sizeof(int ));
		d->pr = (int *)mCalloc(tree->nstate*tree->nstate,sizeof(int ));
		d->s = (int *)mCalloc(tree->nstate,sizeof(int));
		d->sroot = (int *)mCalloc(tree->nstate,sizeof(int));
		d->lmin = (int *)mCalloc(tree->nstate,sizeof(int));
		d->rmin = (int *)mCalloc(tree->nstate,sizeof(int));
		d->prc = (int *)mCalloc(tree->nstate,sizeof(int));
		d->plc = (int *)mCalloc(tree->nstate,sizeof(int));
		d->llock = (int *)mCalloc(tree->nstate,sizeof(int));
		d->rlock = (int *)mCalloc(tree->nstate,sizeof(int));

		// recurse to each descendant node
		For(i,3){
			if(d->v[i]->num != d->anc->num){
				index = Fix_Node_Numbers(d->v[i], nums, index, 0, thresh, tree);
			}
		}
	}
	return index;
}

/********************************************************
 * Fix edges numbers to those in the original tree
 * Similar to Fix_Node_Numbers but with edges
 */
int Fix_Edge_Numbers(t_node* d, int* edges, int index, int root, phydbl thresh, t_tree* tree){
	int i;
	d->anc_edge->num = edges[index++];
	tree->t_edges[d->anc_edge->num] = d->anc_edge;

	if(!(d->tax || d->anc_edge->l > thresh) || root){
		if(tree->mod->optDebug)printf("replacing edges %d %d\n",d->anc_edge->num, edges[index]);
		For(i,3){
			if(d->v[i]->num != d->anc->num){
				index = Fix_Edge_Numbers(d->v[i], edges, index, 0, thresh, tree);
			}
		}
	}
	return index;
}

/*********************************************************
 * Join list of nodes together in a balanced fashion
 * Assumes all nodes are the same level
 */
void Join_Nodes_Balanced(t_node** nodes, int nnodes, int index){
	int i,j;
	int lowest = INT_MAX;
	int li = -1; //find i and j of the lowest combination
	int lj = -1;
	if(nnodes == 1)return;

	// identify the pair of nodes containing the fewest "top" nodes
	For(i,nnodes){
		for(j=i+1; j<nnodes; j++){
			if(nodes[i]->polytomy + nodes[j]->polytomy < lowest){
				lowest = nodes[i]->polytomy + nodes[j]->polytomy;
				li = i;
				lj = j;
			}
		}
	}
	// join these nodes
	// make the top polytomy score the number of descendant "top" nodes
	t_node* top = Join_Nodes(nodes[li],nodes[lj], index+=3);
	top->polytomy = lowest;

	t_node* temp[nnodes-1]; //make a smaller array of nodes
	j = 0; //add non-joined nodes
	For(i,nnodes){
		if(i != li && i != lj){
			temp[j++] = nodes[i];
		}
	}
	nnodes--; //lower index of last node
	temp[nnodes-1] = top; //place new top node in last position
	For(i,nnodes){ //replace node pointers in nodes array
		nodes[i] = temp[i];
	}
	if(nnodes > 1){ //if more than one node left, do it again
		Join_Nodes_Balanced(nodes,nnodes,index);
	}
}


/*********************************************************/

void printTree(t_node* node){
	int i;
	if(!node->tax){
		For(i,3){
			if(node->b[i]->num != node->anc_edge->num){
				printf("tree %d\t%d\t%lf\t",node->num,node->v[i]->num,node->b[i]->l);
				if(node->v[i]->tax)printf("%s",node->v[i]->name);
				printf("\n");
				printTree(node->v[i]);
			}
		}
	}
}

/********************************************************
 * Attach node a below node top in direction ai from a and
 * direction topi from top using an edge of specified length
 * */
void Attach_Edge(t_node* top, t_node* a, int topi, int ai, phydbl length){
	a->v[ai] = top;
	a->anc = top;
	top->v[topi] = a;
	top->b[topi] = a->b[ai];
	a->anc_edge = a->b[ai];
	a->anc_edge->l = length;
	// set top node on appropriate sides of ancestral edge
	if(a->b[ai]->rght->num == a->num)a->b[ai]->left=top;
	else a->b[ai]->rght = top;
}


/*********************************************************
 * Join two nodes together with a new internal node, return
 * their top node. nindex is the index of this new node.
 * Returns top node
 */
t_node* Join_Nodes(t_node* a, t_node* b, int nindex){
	int i;
	t_node* top = Make_Node_Light(nindex); //make new top node
	top->tax = 0;
	top->b[2] = Make_Edge_Light(NULL,NULL,nindex); //ancestral edge
	top->b[2]->rght = top;
	top->b[2]->left = NULL;
	top->anc_edge = top->b[2];
	top->anc_edge->l = 0.0;
	top->polytomy=1;
	//printf("joining %d\t%d\t%d\t%d\n",a->num,b->num,top->num,nindex);
	For(i,3){
		if(!a->tax || i == 0){ //attach top node to node a via ancestral node
			if(a->b[i]->num == a->anc_edge->num){
				// printf("joining %d and %d along edge %d\n",a->num,top->num,i);
				Attach_Edge(top, a, 0, i, a->anc_edge->l);
			}
		}
		if(!b->tax || i == 0){ //attach top node to node b via ancestral node
			if(b->b[i]->num == b->anc_edge->num){
				//  printf("joining %d and %d along edge %d\n",b->num,top->num,i);
				Attach_Edge(top,b,1,i,b->anc_edge->l);
			}
		}
	}
	return top;
}


/*********************************************************
 * Get the nodes at the bottom of the polytomy
 * Starts at polytomy top, store edge number of descendant nodes
 * node numbers of internal nodes, and array of nodes that
 * are the bottom of the polytomy.
 * Returns index of nums array, equal to the number of bottom
 * nodes in polytomy
 */
int Prune_Polytomy(t_node* d, t_node** nodes, int* scores, int* nums, int* edges, int index, int* numindex, int* edgeindex, int root, phydbl thresh, t_tree* tree){
	int i;
	edges[*edgeindex] = d->anc_edge->num; //store edge number of current node
    *edgeindex = *edgeindex + 1;
    // if at a bottom node of the polytomy
	if((d->tax || d->anc_edge->l > thresh) && !root){
		nodes[index++] = d; //place node in nodes array
		scores[d->pstate]++; //add tally to the score of the state at that node
	}else{ //if not at a bottom of the polytomy
		if(tree->mod->optDebug)printf("Prune_Polytomy, node %d %d %lf\n",d->num,*numindex,d->anc_edge->l);
		d->polytomy = 1; //identify this node as an internal polytomy node
		nums[*numindex] = d->num; //store node number
		*numindex = *numindex+1;
		For(i,3){ //recurse downward, storing the index of the last bottom node
			if(d->v[i]->num != d->anc->num){
				index = Prune_Polytomy(d->v[i],nodes,scores, nums, edges, index, numindex, edgeindex, 0, thresh,tree);
			}
		}
	}
	return index;
}

/*********************************************************
 * Count number of states among the bottom nodes of the polytomy
 * The mark parameter causes the scores to be stored
 */
void Count_Polytomy_States(t_node* d, int* scores, int root, phydbl thresh, int mark, t_tree* tree){
  int i;
  if((d->tax || d->anc_edge->l > thresh) && !root){
    if(mark)scores[d->pstate]++;
  }else{
    d->polytomy = 1;
    if(mark){
      d->polystates = scores;
    }
    For(i,3){
      if(d->v[i]->num != d->anc->num){
        Count_Polytomy_States(d->v[i], scores, 0, thresh, mark, tree);
      }
    }
  }
}

/********************************************************
 * Score a polytomy by how far away the first node of each state
 * found at the polytomy tips is from the top of the polytomy.
 * In max ambiguity state, this distance should be minimized.
 * Relative score is relative to the best possible, should be 0 in
 * max ambiguity state.
 * input is top of polytomy, threshold, maxtrees not necessary
 * */
phydbl Score_Polytomy(t_node* top, phydbl thresh, int maxtrees, int relative, t_tree* tree){
	int debug = 1;
	maxtrees = 1;
	int i,j,k;
	int nedges = (tree->n_otu-1)*2;
	t_node* r = tree->noeud[tree->mod->startnode]; //root of tree

  	int pars = Fill_Sankoff(r,tree,1); //max parsimony score
  	Set_Pars_Counters(r,tree,1); //set counters to leftmost configuration

	t_node* tr2 = tree->noeud[top->num]; //get top num (just set to top?)
	int taxinfo[nedges]; //store taxa labels for all tips
	For(i,nedges)taxinfo[i] = tree->noeud[i]->tax;
	For(i,3){ // Isolate polytomy on the descendant sides of the top node
		if(tr2->anc->num != tr2->v[i]->num) Isolate_Polytomy(tr2->v[i], thresh, tree);
	}

	// Find a max parsimony state at the root of the tree
	int smin = INT_MAX;
	int mins = 0;
	For(i,tree->nstate){
		if(r->sroot[i] < smin){
			smin = r->sroot[i];
			mins = i;
		}
	}
	// Set states to first path given isolated polytomy and starting state
	Get_First_Path(tree->noeud[tree->mod->startnode], mins, tree, 1);

	// Score each state by how far away the first polytomy node with that
	// state is from the polytomy top. Find the maximum of these scores.
	int mscore = 0; //max score
	int nzero = 0; //number of non-zero scores (states at polytomy tips)
	int* scores = mCalloc(tree->nstate,sizeof(int));
	Score_Mono(tree->noeud[top->num], 1, scores, debug, tree);
	For(j,tree->nstate){
		if(tree->mod->optDebug) printf("%d\t",scores[j]);
		if(mscore < scores[j]) mscore = scores[j];
		if(scores[j] > 0) nzero++;
	}
	if(tree->mod->optDebug)printf("\n");
	free(scores);

	//find score relative to the best possible score, which is
	//the minimum necessary number of ambiguous internal nodes in the polytomy
	//divided by two, plus any remainder
	if(relative){
		int ninternal = nzero - 2; //number of internal nodes in max ambig state
		int best = 1;
		if(ninternal >= 0)best = floor(ninternal/2.0) + (ninternal % 2) + 2;
		mscore -= best;
		if(tree->mod->optDebug)
			printf("\nRelative polytomy score! states:%d nodes:%d best:%d relative:%d %f %d\n",
				nzero, ninternal, best, mscore, floor(ninternal/2.0), (ninternal % 2));
	}
	//reset taxa designations for all nodes (undo Isolate_Polytomy)
	For(i,nedges)tree->noeud[i]->tax = taxinfo[i];

	return(mscore);
}

/********************************************************
 * Recurse down polytomy, find the first time a node with a particular
 * state is found and record how far away it is from the polytomy top
 * only counts if the state is later found at one of the polytomy ends
 * */
void Score_Mono(t_node* d, int level, int* scores, int debug, t_tree* tree){
	int i;
	//set scores of this node's state to the negative of the current level
	//if the state hasn't been found in the polytomy yet.
	if(scores[d->pstate] == 0){
		scores[d->pstate] = -level;
	}
	if(d->tax){ //if at a tip or end of polytomy, set the score to a positive
		if(scores[d->pstate] < 0) scores[d->pstate] = -scores[d->pstate];
	}
	if(!d->tax){ // recurse down to the end of the polytomy
		For(i,3){
			if(d->v[i]->num != d->anc->num){
				Score_Mono(d->v[i], level+1, scores, debug, tree);
			}
		}
	}
}

/********************************************************
 * Isolate the polytomy from the rest of the tree by setting
 * it's last nodes as tips (tax=1). Do this by recusively
 * moving down until you reach either a tip or an internal node
 * with an ancestral branch length > thresh
 * */
void Isolate_Polytomy(t_node* d, phydbl thresh, t_tree* tree){
	int i;
	if(tree->mod->optDebug)printf("isolating %d\n",d->num);
	if(d->tax && d->num != tree->mod->startnode)return;
	if(d->anc_edge->l < thresh){
		For(i,3){
			if(d->v[i]->num != d->anc->num){
				Isolate_Polytomy(d->v[i],thresh,tree);
			}
		}
	}else{
		if(tree->mod->optDebug)printf("%d is now taxa\n",d->num);
		d->tax = 1;
	}
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

	// try two different NNIs to get rearranged parsimony score
	int pars1 = NNI_ParsSwaps(a,c,d,e,tree);
	int pars2 = NNI_ParsSwaps(b,c,d,e,tree);

	if(tree->mod->optDebug)printf("%d\t%d\t%d\n",pars0,pars1,pars2);
	if(pars0 <= MIN(pars1,pars2)){ //if neither NNI improves score, do nothing
		if(tree->mod->optDebug)printf("Not swapping!\n");
	}else if(pars1 < MIN(pars2,pars0)){ //if first NN smallest, commit to it
	  if(tree->mod->optDebug)printf("0 Swapping!\n");
	  Swap(a,c,d,e,tree);
	  return 1;
	}else if(pars2 < MIN(pars0,pars1)){ //if second, swap that
	  if(tree->mod->optDebug)printf("1 Swapping!\n");
	  Swap(b,c,d,e,tree);
	  return 1;
	}else if(pars1 == pars2){//if both options are equally better than the original, do the first
		if(tree->mod->optDebug)printf("2 Swapping!\n");
		Swap(a,c,d,e,tree);
		return 1;
	}else{ //shouldn't happen
		Warn_And_Exit("Inconsistency in parsimony swapping detected\n");
	}
	// if no NNI found, recurse along another edge
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

	  //get initial parsimony score
	  pars0 = Fill_Sankoff(tree->noeud[tree->mod->startnode],tree,1);
	  if(tree->mod->optDebug){
	  	  Get_First_Path(tree->noeud[tree->mod->startnode],0,tree,1);
	  	  printTreeState(tree->noeud[tree->mod->startnode],tree,1);printf(" 0\n");
	  }
	  // swap nodes, get new score
	  Swap(a,b,c,d,tree);
	  pars1 = Fill_Sankoff(tree->noeud[tree->mod->startnode],tree,1);
	  if(tree->mod->optDebug){
		  Get_First_Path(tree->noeud[tree->mod->startnode],0,tree,1);
		  printTreeState(tree->noeud[tree->mod->startnode],tree,1);printf(" 1\n");
	  }
	  // swap them back
	  Swap(d,b,c,a,tree); //swap nodes
	  pars2 = Fill_Sankoff(tree->noeud[tree->mod->startnode],tree,1);
	  if(tree->mod->optDebug){
	  	  Get_First_Path(tree->noeud[tree->mod->startnode],0,tree,1);
	  	  printTreeState(tree->noeud[tree->mod->startnode],tree,1);printf(" 2\n");
	  }

	  if(pars0 != pars2){ //oops
		  printf("pars score %d\t%d\t%d\n",pars0,pars1,pars2);
		  printf("tree %d\t%s\n",tree->mod->num,tree->mod->rootname);
		  printf("\n.\tParsimony reconstruction swap inconsistent!\n");
		  exit(EXIT_FAILURE);
	  }
	  return pars1;
}

/*********************************************************
 * Fill in dynamic programming tables of the Sankoff algorithm at a
 * given node.
 * Recursive method that begins at root node
 * */
int Fill_Sankoff(t_node *d, t_tree *tree, int root){
  int i,j,k,dir1,dir2; //Recurse to a tip!
  if(!d->tax || root){ //if root or internal node
	  if(tree->mod->optDebug)printf("Sankoff at node %d\n",d->num);
	  if(!root){ // get indexes of descendant nodes
		  t_node* a = d->anc;
		  dir1=dir2=-1;
		  For(i,3) if(d->v[i]->num != a->num) (dir1<0)?(dir1=i):(dir2=i);
	  }else{ //if at root, only have a single descendant
		  dir1=0;
	  }
	  For(j,tree->nstate){ //set state and minimum subtree scores to their highest values
		  if(!root)d->s[j]=1000;
		  else d->sroot[j]=1000;
		  d->lmin[j]=MAX_PARS;
		  d->rmin[j]=MAX_PARS;
	  }
	  // Recurse to descendant nodes!
	  if(!root && tree->mod->optDebug)printf("%d\t%d\t%d\n",d->num,d->v[dir1]->num,d->v[dir2]->num);
      Fill_Sankoff(d->v[dir1],tree,0); //left = 1, right = 2
      if(!root)Fill_Sankoff(d->v[dir2],tree,0);
      //fill in pointers and minimums of dynamic programming table
      For(i,tree->nstate){
    		For(j,tree->nstate){
    			if(root){
    				//root is same as below, but also add the state score since the root is a tip
    				d->pl[i*tree->nstate+j] = d->s[i] + tree->step_mat[i*tree->nstate+j] + d->v[dir1]->s[j];
    				if(d->pl[i*tree->nstate+j] < d->lmin[i]) d->lmin[i] = d->pl[i*tree->nstate+j];
    			}else{
    				//score from i to j in left table is score of j in left node + steps from i to j
    				//similar for right, record minimums on either side along the way
    				d->pl[i*tree->nstate+j] = d->v[dir1]->s[j] + tree->step_mat[i*tree->nstate+j];
    				d->pr[i*tree->nstate+j] = d->v[dir2]->s[j] + tree->step_mat[i*tree->nstate+j];
    				if(d->pl[i*tree->nstate+j] < d->lmin[i]) d->lmin[i] = d->pl[i*tree->nstate+j];
    				if(d->pr[i*tree->nstate+j] < d->rmin[i]) d->rmin[i] = d->pr[i*tree->nstate+j];
    			}
    			if(tree->mod->optDebug)printf("%d\t%d\t%d\t%d\t%d\t%d\n",i,j,d->pl[i*tree->nstate+j],d->pr[i*tree->nstate+j],d->lmin[i],d->rmin[i]);
    		}
    		// score at state i is that resulting in the minimum combined score on both sides.
    		if(!root)d->s[i] = d->rmin[i] + d->lmin[i];
    		else d->sroot[i] = d->lmin[i];
      }
      if(tree->mod->optDebug)printf("\n");
  }
  // find the minimum parsimony score of your state array
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
 * Get maximally parsimonious labeling of tree according to
 * plc and prc arrays for each node
 * "index" is the state at the current node, passed by previous
 * recursion call. "root" indicates whether the node is the root (1/0)
 * */
void Get_First_Path(t_node *d, int index, t_tree *tree, int root){
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
 * Set counters to "left most" minimum scoring position for a node
 * Recursive method, start at root node
 * */
void Set_Pars_Counters(t_node *d, t_tree *tree, int root){
	int i,j,k,dir1,dir2;
	if(!d->tax || root){ //only set counters at internal/root node
	  if(!root){ //find indexes of descendants
		  t_node* a = d->anc;
		  dir1=dir2=-1;
		  For(i,3) if(d->v[i]->num != a->num) (dir1<0)?(dir1=i):(dir2=i);
	  }else{ //if at root, only have a single descendant
		  dir1=0;
	  }
	  int index;
	  int lfound,rfound;
	  For(index,tree->nstate){
		  lfound = rfound = 0;
		  // Counter for each state index is held in plc/prc array
		  // minimum values stored in lmin/rmin
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
		  // set locks (deprecate?)
		  if(lfound==1)d->llock[index]=1;
		  if(!root)if(rfound==1)d->rlock[index]=1;
		  if(!root)if(lfound==0 || rfound==0){
			  printf("No valid pointer found in Sankoff algorithm!\n");
			  exit(EXIT_FAILURE);
		  }
	  }
	  // Recurse to each descendant node
	  Set_Pars_Counters(d->v[dir1],tree,0);
	  if(!root)Set_Pars_Counters(d->v[dir2],tree,0);
	}
}

/********************************************************
 * Get all possible maximum parsimony labels of internal nodes
 *
int Get_All_Paths(t_node *d, int index, t_tree *tree, t_tree** btrees, int root,int maxtrees,int treeindex, int repindex, int rootstate){
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
