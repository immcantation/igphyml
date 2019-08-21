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

/*********************************************************/
// TODO: Edit to output and then free each tree at a time
void Pars_Reconstructions(option* io){
	int j,i,k;
  printf("doing pars reconstructions\n");
    io->maxparsotu = 1;
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

  		  //array of all possible reconstructions
  		  t_tree** trees = mCalloc(maxtrees,sizeof(t_tree*));
  		  trees[0]=tree;
  		  t_node* r = tree->noeud[tree->mod->startnode]; //root node
  		  Init_Class_Tips(tree,io->precon); //initialize tip states and data structures
  		  int pars = Fill_Sankoff(r,tree,1); //max parsimony score
  		  Set_Pars_Counters(r,tree,1); //set counters to their first minimum position

  		  //Data structures for counting the number of type switches and branch lengths in each state

        tree->polytomy_states = mCalloc((tree->n_otu-1)*2,sizeof(int));
        tree->polytomy_swaps = mCalloc((tree->n_otu-1)*2,sizeof(int*));
  	    phydbl* switches = mCalloc(tree->nstate*tree->nstate,sizeof(phydbl));
  	    phydbl* classl = mCalloc(tree->nstate,sizeof(phydbl));
        tree->charindex = mCalloc(tree->nstate,sizeof(int));

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
         //printf("%s\t%d\t%d\t%d\n",tree->chars[i],r->sroot[i],pars,minroots);
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
                  For(k,(trees[i]->n_otu-1)*2){
                    printf("%d",trees[i]->noeud[k]->num);
                    if(trees[i]->noeud[k]->num != trees[i]->mod->startnode)
                      printf("\t%d\n",trees[i]->noeud[k]->anc->num);
                    if(trees[i]->noeud[k]->tax)printf("\t%s\n",trees[i]->noeud[k]->name);
                  }
                  For(k,(trees[i]->n_otu-1)*2){
                    t_node* d = trees[i]->noeud[k];
                    printf("%d\n",d->num);
                    if(!d->tax && !d->anc->tax && d->anc_edge->l < tree->io->thresh && d->polytomy == 2)
                       Count_Polytomy_Switches(d, switches, tree->io->thresh, trees[i]);
                  }
  		  	  		  Fill_Pars_Stats(trees[i]->noeud[trees[i]->mod->startnode],trees[i],switches,classl,1); //summarize statistics of the reconstruction
                  //exit(1);
                  For(k,(trees[i]->n_otu-1)*2)trees[i]->noeud[k]->polytomy = 0;
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
  	  		  FILE* treeout1 = Openfile(fout, 1 );
  			  printf("\n. Sampling %d trees instead for tree %d %s.",sample,j,tree->mod->rootname);
  			  for(i=1;i<maxtrees;i++){ //free trees from attempt to solve for all
  				  if(tree->n_otu < maxotu){
  					  Clean_Tree(trees[i]);
  					  Free_Tree(trees[i]);
  				  }
  			  }
  			  fprintf(treeout1,"#NEXUS\nBegin taxa;\nDimensions ntax=%d;\nTaxlabels\n",tree->n_otu);
  		  	  For(i,nedges){
  		  	  	  if(tree->noeud[i]->tax){
  		  	  		  fprintf(treeout1,"%s\n",tree->noeud[i]->name);
  		  	  	  }
  		  	  }
  		  	  fprintf(treeout1,";\nEnd\nBegin trees;\n");
            t_tree* tree2 = Make_Tree_From_Scratch(tree->n_otu,tree->data);
            Copy_Tree(tree,tree2);
            //Copy_Sankoff_Tree(tree,tree2);
            tree2->chars = tree->chars;
            tree2->mod = tree->mod;
            tree2->io = tree->io;
            Set_Pars_Counters(r,tree,1);
              Get_First_Path(tree->noeud[tree->mod->startnode],minrootsar[0],tree,1);
              if(tree->mod->polytomyresolve >= 2){
              int iters = 0;
              int nnifound = 1;
              while(nnifound){ //execute if a rearrangement is found in a tree
                nnifound=0;
                For(k,nedges){ //if no taxa on either side of edge and length is below threshold, search for NNIs
                  t_node* node = tree->noeud[k];
                  if(!node->tax && !node->anc->tax && node->anc_edge->l < tree->io->thresh && node->polytomy == 0){
                    Resolve_Polytomy_Mono(node,tree->io->thresh,1,tree);
                    nnifound++;
                    if(nnifound){
                      break;
                    }
                  }
                }
                iters++;
              }
              }
              For(k,(tree->n_otu-1)*2)tree->noeud[k]->polytomy = 0;
  			  For(i,sample){
            //printf("\nsample %d",i);
	  	  		  fprintf(treeout1,"Tree TREE%d = [&R] ",i+1);
	  	  		  int n = rand() % minroots;
	  	  		  int minroot = minrootsar[n];
              For(k,(tree->n_otu-1)*2)tree->noeud[k]->polytomy = 0;
              int poly;
              if(tree->mod->polytomyresolve >= 2){
              For(poly,tree->polytomies){
                int states = tree->polytomy_states[poly];
                if(states > 2){
                  int* swaps = tree->polytomy_swaps[poly];
                  For(k,states){
                    //printf("\n%d\t%d\t%d",k,swaps[k],states);
                    //printf("\n%d\t%d\t%d",k,swaps[k]>anc->num,states);
                  }
                  int f = rand() % states;
                  int s;
                  t_node* first = tree->noeud[swaps[f]];
                  t_node* second = first;
                  while(first->anc->num == second->anc->num){
                    s = rand() % states;
                    second = tree->noeud[swaps[s]];
                  }
                  //printf("\nswapping %d %d %d %d",first->num,second->num,f,s);
                  //printf("\nswapping2 %d %d",first->anc->num,second->anc->num);
                  Swap(first,first->anc,second->anc,second,tree);
                }
              }
              }
              //Init_Class_Tips(tree,io->precon);
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
  			  }
  			  fprintf(treeout1,"END;\n");
  			  fclose(treeout1);
          Clean_Tree(tree2);
          Free_Tree(tree2);
  		  }
  		 int tposi, tposj;
       phydbl total = 0.0;
  		 For(tposi,tree->nstate){
  		 	classl[tposi] = classl[tposi]/(npars*1.0);
  		 	For(tposj,tree->nstate){
  		 		//printf("%lf\t%d\n",switches[tposi*tree->nstate+tposj],npars);
  		 		switches[tposi*tree->nstate+tposj]=switches[tposi*tree->nstate+tposj]/(npars*1.0);
          total += switches[tposi*tree->nstate+tposj];
  		 	}
  		 }
       //printf("\nTOTAL %lf\n",total);
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


/********************************************************
 * Recurse down tree, randomly choosing ambiguous pointers
 * TODO: Make descendant nodes from polytomies consider all possible 
 * states within a polytomy as a possible ancestor
 * if(d->anc->polytomy == 1) xyz...
 * */
void Get_Rand_Path(t_node *d, int index, t_tree *tree, int root){
  int i,j,dir1,dir2;
  int ldraw=0;
  int rdraw=0;
  int lfound=0;
  int rfound=0;
  if(d->polytomy == 0){
    d->pstate=index;
  }
  if(!d->tax&&tree->mod->optDebug)printf("%s\t%d\n",tree->chars[index],index);
  if(!d->tax || root){
    if(!root){
      t_node* a = d->anc;
      dir1=dir2=-1;
      For(i,3) if(d->v[i]->num != a->num) (dir1<0)?(dir1=i):(dir2=i);
      if((d->polytomy == 0) && (d->b[dir1]->l < tree->io->thresh || 
          d->b[dir2]->l < tree->io->thresh) && tree->mod->polytomyresolve >= 2){
        int* scores = mCalloc(tree->nstate,sizeof(int));
        scores[d->pancstate]++;
        Get_Rand_Path_Polytomy(d, scores, index, tree, 1);
        d->polytomy = 1;
      }
    }else{
      dir1=0;
    }
    int* scores;
    if(d->polytomy == 0){
      //free(d->polystates);
      scores = mCalloc(tree->nstate,sizeof(int));
      scores[index] = 1;
    }else{
      scores = d->polystates;
      //printf("\n");
      //For(i,tree->nstate)printf("%d\t",scores[i]);
      //printf("\n");
      //For(i,tree->nstate)printf("%s\t",tree->chars[i]);
    }
    int lmins=0;int rmins=0; //tally up number of minimal paths
    For(j,tree->nstate){ //allow for any state in a polytomy to be used
      if(scores[j] > 0){
        For(i,tree->nstate){
         if(d->pl[j*tree->nstate+i] == d->lmin[d->pstate])lmins++;
         if(!root)if(d->pr[j*tree->nstate+i] == d->rmin[d->pstate])rmins++;
        }
      }
    }
    if(lmins>1)ldraw = rand() % lmins; //not perfectly random but probably okay
    if(!root && rmins>1)rdraw = rand() % rmins;
    For(j,tree->nstate){
      if(scores[j] > 0){
        For(i,tree->nstate){
        /*if(d->num == 13){
          printf("\noptions %s\t%s\t%d\t%d\t%d\t%d",tree->chars[j],tree->chars[i],
          d->pl[j*tree->nstate+i],d->lmin[d->pstate],d->pr[j*tree->nstate+i],d->rmin[d->pstate]);
        }*/
         if(d->pl[j*tree->nstate+i] == d->lmin[d->pstate]){
          //if(d->polytomy==1 && d->v[dir1]->polytomy == 0 && !d->v[dir1]->tax)
          //  printf("\n!!!adjacent tip %d\t%s\t%s",d->num,tree->chars[j],tree->chars[i]);
          if(lfound == ldraw){
            if(d->v[dir1]->polytomy==0)d->v[dir1]->pancstate = j;
            Get_Rand_Path(d->v[dir1],i,tree,0);
          }
            lfound++;
         }
         if(d->pr[j*tree->nstate+i] == d->rmin[d->pstate] && !root ){
          //if(!root)if(d->polytomy==1 && d->v[dir2]->polytomy == 0  && !d->v[dir2]->tax)
          //  printf("\n!!!adjacent tip %d\t%s\t%s",d->num,tree->chars[j],tree->chars[i]);
          if(rfound == rdraw){
            if(d->v[dir2]->polytomy==0)d->v[dir2]->pancstate = j;
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
  /*if(!root)printf("\nTREE %d\t%d\t%s",d->anc->num,d->num,tree->chars[d->pstate]);
  if(d->tax)printf("\t%s",d->name);*/
}

/********************************************************
 * Recurse down polytomy, randomly choosing ambiguous pointers
 * TODO: Make descendant nodes from polytomies consider all possible 
 * states within a polytomy as a possible ancestor
 * if(d->anc->polytomy == 1) xyz...
 * */
void Get_Rand_Path_Polytomy(t_node *d, int* scores, int index, t_tree *tree, int top){
  int i,j,dir1,dir2;
  int ldraw=0;
  int rdraw=0;
  int lfound=0;
  int rfound=0;
  int root = 0;
  if(d->num == tree->mod->startnode)root=1;
  d->pstate=index;
  d->polystates = scores;
  if((!d->tax && d->anc_edge->l < tree->io->thresh) || top){
    d->polytomy = 1;
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
          if(lfound == ldraw){
            d->v[dir1]->pancstate = index;
            Get_Rand_Path_Polytomy(d->v[dir1],scores,i,tree,0);
          }
          lfound++;
        }
        if(d->pr[index*tree->nstate+i] == d->rmin[index] && !root ){
          if(rfound == rdraw){
            d->v[dir2]->pancstate = index;
            Get_Rand_Path_Polytomy(d->v[dir2],scores,i,tree,0);
          }
          rfound++;
        }
    }
  }else{
      scores[d->pstate]++;
      if(tree->mod->optDebug)printf("%s\t%s\n",tree->chars[index],d->name);
  }
}

/********************************************************
 * Recursively add values to switch and length matrices for each type
 * */
void Fill_Pars_Stats(t_node* d,t_tree* tree, phydbl* switches, phydbl* classl, int root){
  int i,j,dir1,dir2;
  dir1=dir2=-1;
  //printf("%d rank %d\n",d->num,d->polytomy);
  int polytip = 0;
  if(!root){
    if(d->anc->polytomy != 0){
      if(d->anc_edge->l >= tree->io->thresh || d->tax){
        polytip = 1;
      }
    }  
  }
  //if(!root)printf("%d\t%d\t%d\t%lf\n",d->pancstate,d->pstate,d->anc->polytomy,d->anc_edge->l);
  if(!root){
    if(d->pstate != d->pancstate){
      //printf("\nswitch! %d\t%s\t%s\n",d->num,tree->chars[d->pancstate],tree->chars[d->pstate]);
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
 * */
void Count_Polytomy_Switches(t_node* d, phydbl* switches, phydbl thresh, t_tree* tree){
  int i,j;
 // printf("counting polytomy switches\n");
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
          nswitches++;
        }
      }
    }
  }
  For(i,tree->nstate){
    For(j,tree->nstate){
      if(scores[i] > 0 && scores[j] > 0 && i != j){
        if(tree->step_mat[i*tree->nstate+j] <= 1 && j != ancstate){
          //printf("switch counts %s\t%s\t%lf\n",tree->chars[i],tree->chars[j],nsteps/nswitches);
          switches[i*tree->nstate+j] += nsteps/nswitches;
        }
      }
    }
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
	//For(i,tree->n_otu)printf("1 %d\t%s\n",i,tree->noeud[i]->name);
	int n = tree->n_otu;
	char temp[T_MAX_NAME];
    For(i,n - 1) {
    		if(strcmp(tree->noeud[i]->name,tree->mod->rootname) != 0){
    			 //   			printf("%s\n",tree->noeud[i]->name);
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
void Permute_All_MetaData(option* io, int pos){
	int i, j;
	int n_otu = 0;
	For(i,io->ntrees){
		n_otu += io->tree_s[i]->n_otu-1;
	}
	int* indexes = mCalloc(n_otu,sizeof(int));
	char** names = mCalloc(n_otu,sizeof(char*));

	int index = 0;
	For(i,io->ntrees){
		For(j,io->tree_s[i]->n_otu){
			if(strcmp(io->tree_s[i]->noeud[j]->name,io->mod_s[i]->rootname) != 0){
				indexes[index] = index;
				names[index] = mCalloc(T_MAX_NAME,sizeof(char));
				strcpy(names[index],io->tree_s[i]->noeud[j]->name);
				//printf("%d %s\n",index,names[index]);
				index++;
			}
		}
	}
	int n = n_otu;
	int temp;
	For(i,n - 1){
    	j = i + rand() / (RAND_MAX / (n - i) + 1);
    	temp = indexes[i];
    	indexes[i] = indexes[j];
    	indexes[j] = temp;
    }

	index= 0;
	For(i,io->ntrees){
		For(j,io->tree_s[i]->n_otu){
			if(strcmp(io->tree_s[i]->noeud[j]->name,io->mod_s[i]->rootname) != 0){
				strcpy(io->tree_s[i]->noeud[j]->name,names[indexes[index++]]);
				//printf("x %s\n",io->tree_s[i]->noeud[j]->name);
			}
		}
	}

    free(indexes);
    For(i,n_otu)free(names[i]);
	free(names);
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
	 /*For(j,statecount){
		 //printf("state %d\t%s\t%d\n",j,ambigstatesfrom[j],ambigstatesto[j]);
	 }*/
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

/*********************************************************/
// Resolve polytomies based on parsimony score
int Resolve_Polytomies_Pars(t_tree* tree, phydbl thresh){
	int i,j;
	int nni=1;
	int nnifound=1; //initial parsimony score
	int pars0 = Fill_Sankoff(tree->noeud[tree->mod->startnode],tree,1);
	int iters = 0;
	int nedges = (tree->n_otu-1)*2;
	while(nnifound){ //execute if a rearrangement is found in a tree
		nnifound=0;
		For(i,nedges){ //if no taxa on either side of edge and length is below threshold, search for NNIs
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
				Resolve_Polytomy_Mono(node,thresh,1,tree);
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


/*********************************************************/

t_node* Resolve_Polytomy_Mono(t_node* b, phydbl thresh, int randomize, t_tree* tree){
	int i,j;
	int nedges = (tree->n_otu-1)*2;
	t_node* top = b;
	while(top->anc->num != tree->mod->startnode && top->anc_edge->l < thresh){
		top = top->anc;
	}
	t_node* anc = top->anc;
	//Init_Class_Tips(tree,tree->io->precon);
	Fill_Sankoff(tree->noeud[tree->mod->startnode],tree,1);
	int pscore = Score_Polytomy(top,thresh,1000,0,tree);
  /*tree->io->precon *= 10;
    printf("%s\n",Write_Tree(tree));
    tree->io->precon /= 10;*/

	t_node** nodes = mCalloc(nedges,sizeof(t_node*));
	t_node** tops = mCalloc(tree->nstate,sizeof(t_node*));
	int* scores = mCalloc(tree->nstate,sizeof(int));
	int* nums = malloc(nedges*sizeof(int));
  int* edges = malloc(nedges*sizeof(int));
	int* numindex = malloc(sizeof(int));
  int* edgeindex = malloc(sizeof(int));
	*numindex = 0;
  *edgeindex = 0;
	int nnodes = Prune_Polytomy(top,nodes,scores,nums,edges,0,numindex,edgeindex,1,thresh,tree);
	int nzero = 0;
	For(i,tree->nstate)if(scores[i] > 0)nzero++;
	if(nnodes < 3 || nzero == 1)return top; //if polytomy only two nodes, ignore
  int n = tree->nstate;
  if(randomize){
    For(i,tree->nstate){
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

	int node;
	int tindex = -1;
	int nindex = nedges+1;
	For(n,tree->nstate){
    j = tree->charindex[n];
		if(scores[j] == 0)continue;
		For(node,nnodes){
			if(nodes[node]->pstate == j){
				nodes[node]->polytomy = 1;
				if(scores[j] > 0){
					tops[++tindex] = nodes[node];
					scores[j] = -1*scores[j];
				}else{
					tops[tindex] = Join_Nodes(tops[tindex],nodes[node],nindex+=3);
				}
			}
		}
	}

  int* swaplist = mCalloc(tindex+1,sizeof(int));
  tree->polytomy_states[tree->polytomies] = tindex+1;
  For(n,tindex+1)swaplist[n] = tops[n]->num;
  tree->polytomy_swaps[tree->polytomies] = swaplist;
  
  /*For(n,tindex+1){
      printf("\nnodes to swap1: %d\t%d",n,tree->polytomy_swaps[tree->polytomies][n]->num);
  }*/

	Join_Nodes_Balanced(tops,tindex+1,nindex);
	t_node* newtop = tops[0];

	For(i,3){
		For(j,3){
			if(anc->v[i]->num == top->num && newtop->b[j]->num == newtop->anc_edge->num){
				Attach_Edge(anc,newtop,i,j,anc->b[i]->l);
				i=3;j=3;
				break;
			}
		}
	}

  //printf("\npolytomy top %d\t%d\t%d",top->num,newtop->num,nums[0]);

	Fix_Node_Numbers(newtop,nums,0,1,thresh,tree);
	Fix_Edge_Numbers(newtop,edges,0,1,thresh,tree);

  /*For(n,tindex+1){
      printf("\nnodes to swap2: %d\t%d",n,tree->polytomy_swaps[tree->polytomies][n]);
  }*/

  //Clean_Tree(tree);
	//Init_Class_Tips(tree,tree->io->precon);
  int firstscore = Score_Polytomy(newtop,thresh,1000,0,tree);
	
  //!!TODO: detect if the underlying characters have changed with this rearrangement
  int redo = 0;
  int* nscores = mCalloc(tree->nstate,sizeof(int));
  Count_Polytomy_States(newtop,nscores,1,thresh,1,tree);
  For(i,tree->nstate)if(nscores[i] != -scores[i])redo=1;
  
  if(redo){
    printf("\nstate change detected at node %d - recursively fixing\n",newtop->num);
    For(i,tree->nstate)printf("%s\t",tree->chars[i]);
    printf("\n");
    For(i,tree->nstate)printf("%d\t",-scores[i]);
    printf("\n");
    For(i,tree->nstate)printf("%d\t",nscores[i]);
    printf("\n");
    For(i,nnodes-1)tree->noeud[i]->polytomy = 0;
    newtop = Resolve_Polytomy_Mono(newtop,thresh,0,tree);
  }

	int score = Score_Polytomy(newtop,thresh,1000,0,tree);
	int rscore = Score_Polytomy(newtop,thresh,1000,1,tree);

  Update_Dirs(tree);
  Update_Ancestors_Edge(tree->noeud[tree->mod->startnode],tree->noeud[tree->mod->startnode]->v[0],tree->noeud[tree->mod->startnode]->b[0],tree);
  tree->polytomies++;

 /* printf("printing tree\n");
  printTree(tree->noeud[tree->mod->startnode]->v[0]);
printf("printing tree\n");*/

	if(rscore > 0){
		printf("\n!! POLYTOMY NOT FULLY RESOLVED. Trying to resursively fix..");
		printf("\nPolytomy at %d, PSCORE: %d, SCORE: %d, RSCORE %d",newtop->num,pscore,score,rscore);
		tree->io->precon *= 10;
    printf("%s\n",Write_Tree(tree));
    tree->io->precon /= 10;
    exit(1);
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
 */
int Fix_Node_Numbers(t_node* d, int* nums, int index, int root, phydbl thresh, t_tree* tree){
	int i;
	if(!(d->tax || d->anc_edge->l > thresh) || root){
		//printf("replacing %d %d\n",d->num, nums[index]);
    //if(!root){
      For(i,tree->polytomy_states[tree->polytomies]){
        if(tree->polytomy_swaps[tree->polytomies][i] == d->num){
          //printf("replacing swap %d %d\n",d->num, nums[index]);
          tree->polytomy_swaps[tree->polytomies][i] = nums[index];
        }
      }
		  d->num = nums[index++];
      Clean_Sankoff_Node(tree->noeud[d->num]);
		  Free_Node(tree->noeud[d->num]);
		  tree->noeud[d->num] = d;
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
    /*}else{
      if(d->num != nums[index]){
        printf("\n POLYTOMY TOP NUMBER HAS CHANGED!\n");
        exit(1);
      }
      index++;
    }*/
		For(i,3){
			if(d->v[i]->num != d->anc->num){
				index = Fix_Node_Numbers(d->v[i], nums, index, 0, thresh, tree);
			}
		}
	}
	return index;
}

/********************************************************
 * Fixes edges numbers to those in the original tree
 */
int Fix_Edge_Numbers(t_node* d, int* edges, int index, int root, phydbl thresh, t_tree* tree){
  int i;
  //if(!root){
    d->anc_edge->num = edges[index++];
    tree->t_edges[d->anc_edge->num] = d->anc_edge;
  /*}else{
    if(d->anc_edge->num != edges[index]){
        printf("\n POLYTOMY TOP EDGE NUMBER HAS CHANGED!\n");
        exit(1);
    }
  }*/
  if(!(d->tax || d->anc_edge->l > thresh) || root){
    //printf("replacing edges %d %d\n",d->anc_edge->num, edges[index]);
  //  d->anc_edge->num = edges[index++];
   // Free_Node(tree->noeud[d->num]);
    //tree->noeud[d->num] = d;
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

	For(i,nnodes){
		for(j=i+1;j<nnodes;j++){
			if(nodes[i]->polytomy + nodes[j]->polytomy < lowest){
				lowest = nodes[i]->polytomy + nodes[j]->polytomy;
				li = i;
				lj = j;
			}
		}
	}
	t_node* top = Join_Nodes(nodes[li],nodes[lj], index+=3);
	top->polytomy = lowest;

	t_node* temp[nnodes-1];
	j = 0;
	For(i,nnodes)if(i != li && i != lj)temp[j++] = nodes[i];

	nnodes--;
	temp[nnodes-1] = top;

	For(i,nnodes){
		nodes[i] = temp[i];
	}

	if(nnodes > 1){
		//printf("NEXT LEVEL!\n");
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

/*********************************************************/

void Attach_Edge(t_node* top, t_node* a, int topi, int ai, phydbl length){
	a->v[ai] = top;
	a->anc = top;
	top->v[topi] = a;

	top->b[topi] = a->b[ai];
	a->anc_edge = a->b[ai];
	a->anc_edge->l = length;
	if(a->b[ai]->rght->num == a->num)a->b[ai]->left=top;
	else a->b[ai]->rght = top;
}


/*********************************************************
 * Join two nodes together, return their top
 */
t_node* Join_Nodes(t_node* a, t_node* b, int nindex){
	int i;
	t_node* top = Make_Node_Light(nindex);
	top->tax = 0;
	top->b[2] = Make_Edge_Light(NULL,NULL,nindex);
	top->b[2]->rght = top;
	top->b[2]->left = NULL;
	top->anc_edge = top->b[2];
	top->anc_edge->l = 0.0;
	top->polytomy=1;
  //printf("joining %d\t%d\t%d\t%d\n",a->num,b->num,top->num,nindex);
	For(i,3){
		if(!a->tax || i == 0){
			if(a->b[i]->num == a->anc_edge->num){
        //printf("%d\t%d\t%d\n",a->b[i]->num,a->anc_edge->num,nindex);
       // printf("joining %d and %d along edge %d\n",a->num,top->num,i);
				Attach_Edge(top,a,0,i,a->anc_edge->l);
			}
		}
		if(!b->tax || i == 0){
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
 */
int Prune_Polytomy(t_node* d, t_node** nodes, int* scores, int* nums, int* edges, int index, int* numindex, int* edgeindex, int root, phydbl thresh, t_tree* tree){
	int i;
  edges[*edgeindex] = d->anc_edge->num;
  *edgeindex = *edgeindex + 1;
	if((d->tax || d->anc_edge->l > thresh) && !root){
		nodes[index++] = d;
		scores[d->pstate]++;
	}else{
		//printf("at node %d %d %lf\n",d->num,*numindex,d->anc_edge->l);
		d->polytomy = 1;
    nums[*numindex] = d->num;
		*numindex = *numindex+1;
		For(i,3){
			if(d->v[i]->num != d->anc->num){
				index = Prune_Polytomy(d->v[i],nodes,scores, nums, edges, index, numindex, edgeindex, 0, thresh,tree);
			}
		}
	}
	return index;
}

/*********************************************************
 * Get the nodes at the bottom of the polytomy
 */
void Count_Polytomy_States(t_node* d, int* scores, int root, phydbl thresh, int mark, t_tree* tree){
  int i;
  if((d->tax || d->anc_edge->l > thresh) && !root){
    if(mark)scores[d->pstate]++;
  }else{
    d->polytomy = 1;
    if(mark){
      d->polytomy = 1;
      d->polystates = scores;
    }
    For(i,3){
      if(d->v[i]->num != d->anc->num){
        Count_Polytomy_States(d->v[i], scores, 0, thresh, mark, tree);
      }
    }
  }
}

/*********************************************************/

phydbl Score_Polytomy(t_node* top,phydbl thresh, int maxtrees, int relative, t_tree* tree){
	int debug = 1;
	maxtrees = 1;
	int i,j,k;
	int nedges = (tree->n_otu-1)*2;

	t_tree** btrees = mCalloc(maxtrees,sizeof(t_tree*));
	t_node* r = tree->noeud[tree->mod->startnode];
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

	int mscore = 0;
	int nzero = 0;
		int* scores = mCalloc(tree->nstate,sizeof(int));
		Score_Mono(tree->noeud[top->num],1,scores,debug,tree);
		For(j,tree->nstate){
      //printf("%d\t",scores[j]);
			if(mscore < scores[j])mscore = scores[j];
			if(scores[j] > 0)nzero++;
		}
    //printf("\n");
		free(scores);


	if(relative){
		int ninternal = nzero - 2;
		int best = 1;
		if(ninternal >= 0)best = floor(ninternal/2.0) + (ninternal % 2) + 2;
		mscore -= best;
		//printf("\n%d\t%d\t%d",ninternal,best,mscore);
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
	if(scores[d->pstate] == 0){
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

/*********************************************************/
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
      printf("tree %d\t%s\n",tree->mod->num,tree->mod->rootname);
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
	  //printf("at node %d\n",d->num);
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
