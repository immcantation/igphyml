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
extern int     stopCodons[64];
extern int   senseCodons[64];
extern char aminoAcidmap[65];
extern int indexSenseCodons[64];




/*********************************************************/
void Init_Class_Tips(t_tree* tree){
	 int i,j;
	 char* mtemp;
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
	 For(i,(tree->n_otu-1)*2){
		// printf("%d\n",i);
	 	tree->noeud[i]->pl = (int *)mCalloc(tree->nstate*tree->nstate,sizeof(int ));
	 	tree->noeud[i]->pr = (int *)mCalloc(tree->nstate*tree->nstate,sizeof(int ));
	 	tree->noeud[i]->s = (int *)mCalloc(tree->nstate,sizeof(int));
	 	tree->noeud[i]->sroot = (int *)mCalloc(tree->nstate,sizeof(int));
	 	tree->noeud[i]->lmin = (int *)mCalloc(tree->nstate,sizeof(int));
	 	tree->noeud[i]->rmin = (int *)mCalloc(tree->nstate,sizeof(int));
	 	tree->noeud[i]->prc = (int *)mCalloc(tree->nstate,sizeof(int));
	 	tree->noeud[i]->plc = (int *)mCalloc(tree->nstate,sizeof(int));
	 	tree->noeud[i]->llock = (int *)mCalloc(tree->nstate,sizeof(int));
	 	tree->noeud[i]->rlock = (int *)mCalloc(tree->nstate,sizeof(int));
	 	For(j,tree->nstate){
	 		tree->noeud[i]->s[j]=1000;
	 		tree->noeud[i]->lmin[j]=MAX_PARS;
	 		tree->noeud[i]->rmin[j]=MAX_PARS;
	 	}
	 }
	 //printf("here\n");
	 free(tree->step_mat);
	 tree->step_mat = mCalloc(tree->nstate*tree->nstate,sizeof(int));
	 For(i,tree->nstate){ //set up step mat
	 	 For(j,tree->nstate){
	 			 if(i < j && tree->io->precon>0)tree->step_mat[j*tree->nstate+i]=10000; //if precon is negative, no costraint
	 			 else if(i == j)tree->step_mat[j*tree->nstate+i]=0;
	 			 else tree->step_mat[j*tree->nstate+i] =1;
	 	 }
	 }
	 //printf("here\n");
	 For(i,tree->n_otu){
		// tree->noeud[i]->b[0]->ui_r[0] = 0;
		 int nelements=0;
		 char* state = mCalloc(T_MAX_OPTION,sizeof(char));
		 char* minfo1 = strdup(tree->noeud[i]->name);
		 char* minfo2 = strdup(tree->noeud[i]->name);
		 while ((mtemp = strsep(&minfo1, "_")) != NULL){nelements++;}
		 For(j,nelements){
			 strcpy(state,strsep(&minfo2, "_"));
		 }
		 if(tree->mod->optDebug)printf("\n%s\t%s",tree->noeud[i]->name,state);
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
		 For(j,tree->nstate){
			 if(tree->mod->optDebug)printf(" %d",tree->noeud[i]->s[j]);
		 }
		 if(tree->mod->optDebug)printf("\n");
	 }
}

/*********************************************************/
//Free extraneous data structures
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
	 //free(tree->step_mat);
}


/*********************************************************/
//Free extraneous data structures
void Copy_Sankoff_Tree(t_tree* tree1,t_tree* tree2){
	int i,j,k;
	tree2->nstate=tree1->nstate;
    For(i,(tree1->n_otu-1)*2){
    	//printf("here\n");
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
    		//printf("here\n");
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
void Fill_Sankoff_Root(t_tree *tree, int root){
	int i,j;
	t_node* d = tree->noeud[tree->mod->startnode];
	Fill_Sankoff(d->v[0],tree,1);
	if(tree->mod->optDebug)printf("root\n");
	For(i,tree->nstate){
	   For(j,tree->nstate){//fill in pointers and minimums
		   d->pl[i*tree->nstate+j]=d->s[i]+tree->step_mat[i*tree->nstate+j]+d->v[0]->s[j];
		   if(d->pl[i*tree->nstate+j] < d->lmin[i])d->lmin[i]=d->pl[i*tree->nstate+j];
		  // printf("%d\t%d\t%d\t%d\n",i,j,d->pl[i*tree->nstate+j],d->lmin[i]);
	   }
	   d->sroot[i]=d->lmin[i];
	   //printf("%d ",d->sroot[i]);
	}
	//printf("\n");
}
/*********************************************************/
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
	  //if(!root)printf("%d\t%d\t%d\n",d->num,d->v[dir1]->num,d->v[dir2]->num);
      Fill_Sankoff(d->v[dir1],tree,0); //left = 1, right = 2
      if(!root)Fill_Sankoff(d->v[dir2],tree,0);
      //Fill in arrays
      For(i,tree->nstate){
    		For(j,tree->nstate){//fill in pointers and minimums
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
/*********************************************************/
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
		/*For(i,tree->nstate){
			if(d->pl[index*tree->nstate+i] == d->lmin[index] && !lfound){
				Get_First_Path(d->v[dir1],i,tree,0);
				lfound=1;
			}
			if(d->pr[index*tree->nstate+i] == d->rmin[index] && !root && !rfound){
				Get_First_Path(d->v[dir2],i,tree,0);
				rfound=1;
			}
		}*/
	}
}
/*********************************************************/
void Get_Rand_Path(t_node *d, int index, t_tree *tree,int root){
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

/*********************************************************/
//set counters
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
		  if(!root)if(lfound==0 || rfound==0)exit(EXIT_FAILURE);
	  }
	  Set_Pars_Counters(d->v[dir1],tree,0);
	  if(!root)Set_Pars_Counters(d->v[dir2],tree,0);
	}
}

/*********************************************************/

int Get_All_Paths(t_node *d, int index, t_tree *tree, t_tree** btrees, int root,int maxtrees,int treeindex){
	int i,j,dir1,dir2;
	int lfound=0;
	int rfound=0;
	d->pstate=index;
	//printf("%s\t%d\t%d\n",tree->chars[index],index,treeindex);
	if(!d->tax || root){
		if(!root){
			t_node* a = d->anc;
			dir1=dir2=-1;
			For(i,3) if(d->v[i]->num != a->num) (dir1<0)?(dir1=i):(dir2=i);
		}else{
			dir1=0;
		}
		int lm = d->lmin[index];
		int rm = d->rmin[index];
		int lc = d->plc[index];
		int rc = d->prc[index];

		//printf("Left treeb! %s\t%d\t%d\t%d\n",tree->chars[lc],treeindex,d->llock[index],d->plc[index]);
		//exit(EXIT_FAILURE);
		if(d->llock[index]){ //left node
			treeindex=Get_All_Paths(d->v[dir1],d->plc[index],tree,btrees,0,maxtrees,treeindex);
		}else{
			For(i,tree->nstate){
				if(d->pl[index*tree->nstate+i] == d->lmin[index]){
					if(lfound>0 && treeindex+1 < maxtrees){ //if alternate left path found
						//printf("Left tree!\n");
						t_tree* tree2;
						tree2 = Read_User_Tree(tree->io->tree_s[0]->data,tree->io->mod_s[0],tree->io);
						tree2->mod=tree->mod;
						tree2->mod->startnode = tree->mod->startnode;
						t_node* r2 = tree2->noeud[tree->mod->startnode];
						tree2->io=tree->io;
						Update_Ancestors_Edge(r2,r2->v[0],r2->b[0],tree);
						d->plc[index]=i;
						d->llock[index]=1;
						Copy_Sankoff_Tree(tree,tree2);
						treeindex++;
						btrees[treeindex]=tree2;
						treeindex = Get_All_Paths(r2,0,tree2,btrees,1,maxtrees,treeindex);
					}
					lfound++;
				}
			}
			d->plc[index]=lc;//continue down original tree
			d->llock[index]=1; //with path through node fixed
			treeindex=Get_All_Paths(d->v[dir1],d->plc[index],tree,btrees,0,maxtrees,treeindex);
		}
		if(!root){
			//printf("Right treeb! %s\t%d\t%d\t%d\n",tree->chars[rc],treeindex,d->rlock[index],d->prc[index]);
			if(d->rlock[index]){ //right node
				treeindex=Get_All_Paths(d->v[dir2],d->prc[index],tree,btrees,0,maxtrees,treeindex);
			}else{
				For(i,tree->nstate){
					if(d->pr[index*tree->nstate+i] == d->rmin[index]){
						if(rfound>0 && treeindex+1 < maxtrees){ //if alternate left path found
							//printf("Right tree! %s\t%s\t%s\n",tree->chars[index],tree->chars[rc],tree->chars[i]);
							t_tree* tree2;
							tree2 = Read_User_Tree(tree->io->tree_s[0]->data,tree->io->mod_s[0],tree->io);
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
							treeindex = Get_All_Paths(r2,0,tree2,btrees,1,maxtrees,treeindex);
						}
						rfound++;
					}
				}
				d->prc[index]=rc;//continue down original tree
				d->rlock[index]=1; //with path through node fixed
				treeindex=Get_All_Paths(d->v[dir2],d->prc[index],tree,btrees,0,maxtrees,treeindex);
			}
		}
	}else{
		if(tree->mod->optDebug)printf("%s\t%s\n",tree->chars[index],d->name);
	}
	return treeindex;
}

/*********************************************************/

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
