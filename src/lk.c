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

#include "lk.h"
extern int senseCodons[64];
extern int indexSenseCodons[64];
extern phydbl ecmK07freq[61];
extern phydbl ecmS05freq[64];
/* int    LIM_SCALE; */
/* phydbl LIM_SCALE_VAL; */
/* phydbl BIG; */
/* phydbl SMALL; */

/*********************************************************/

void Init_Tips_At_One_Site_Nucleotides_Float(char state, int pos, phydbl *p_lk)
{
  switch(state)
    {
    case 'A' : p_lk[pos+0]=1.; p_lk[pos+1]=p_lk[pos+2]=p_lk[pos+3]=.0;
      break;
    case 'C' : p_lk[pos+1]=1.; p_lk[pos+0]=p_lk[pos+2]=p_lk[pos+3]=.0;
      break;
    case 'G' : p_lk[pos+2]=1.; p_lk[pos+1]=p_lk[pos+0]=p_lk[pos+3]=.0;
      break;
    case 'T' : p_lk[pos+3]=1.; p_lk[pos+1]=p_lk[pos+2]=p_lk[pos+0]=.0;
      break;
    case 'U' : p_lk[pos+3]=1.; p_lk[pos+1]=p_lk[pos+2]=p_lk[pos+0]=.0;
      break;
    case 'M' : p_lk[pos+0]=p_lk[pos+1]=1.; p_lk[pos+2]=p_lk[pos+3]=.0;
      break;
    case 'R' : p_lk[pos+0]=p_lk[pos+2]=1.; p_lk[pos+1]=p_lk[pos+3]=.0;
      break;
    case 'W' : p_lk[pos+0]=p_lk[pos+3]=1.; p_lk[pos+1]=p_lk[pos+2]=.0;
      break;
    case 'S' : p_lk[pos+1]=p_lk[pos+2]=1.; p_lk[pos+0]=p_lk[pos+3]=.0;
      break;
    case 'Y' : p_lk[pos+1]=p_lk[pos+3]=1.; p_lk[pos+0]=p_lk[pos+2]=.0;
      break;
    case 'K' : p_lk[pos+2]=p_lk[pos+3]=1.; p_lk[pos+0]=p_lk[pos+1]=.0;
      break;
    case 'B' : p_lk[pos+1]=p_lk[pos+2]=p_lk[pos+3]=1.; p_lk[pos+0]=.0;
      break;
    case 'D' : p_lk[pos+0]=p_lk[pos+2]=p_lk[pos+3]=1.; p_lk[pos+1]=.0;
      break;
    case 'H' : p_lk[pos+0]=p_lk[pos+1]=p_lk[pos+3]=1.; p_lk[pos+2]=.0;
      break;
    case 'V' : p_lk[pos+0]=p_lk[pos+1]=p_lk[pos+2]=1.; p_lk[pos+3]=.0;
      break;
    case 'N' : case 'X' : case '?' : case 'O' : case '-' :
      p_lk[pos+0]=p_lk[pos+1]=p_lk[pos+2]=p_lk[pos+3]=1.;break;
    default :
      {
	PhyML_Printf("\n. Unknown character state : %c\n",state);
	Exit("\n. Init failed (check the data type)\n");
	break;
      }
    }
}

/*********************************************************/

void Init_Tips_At_One_Site_Nucleotides_Int(char state, int pos, short int *p_pars)
{
  switch(state)
    {
    case 'A' : p_pars[pos+0]=1; p_pars[pos+1]=p_pars[pos+2]=p_pars[pos+3]=0;
      break;
    case 'C' : p_pars[pos+1]=1; p_pars[pos+0]=p_pars[pos+2]=p_pars[pos+3]=0;
      break;
    case 'G' : p_pars[pos+2]=1; p_pars[pos+1]=p_pars[pos+0]=p_pars[pos+3]=0;
      break;
    case 'T' : p_pars[pos+3]=1; p_pars[pos+1]=p_pars[pos+2]=p_pars[pos+0]=0;
      break;
    case 'U' : p_pars[pos+3]=1; p_pars[pos+1]=p_pars[pos+2]=p_pars[pos+0]=0;
      break;
    case 'M' : p_pars[pos+0]=p_pars[pos+1]=1; p_pars[pos+2]=p_pars[pos+3]=0;
      break;
    case 'R' : p_pars[pos+0]=p_pars[pos+2]=1; p_pars[pos+1]=p_pars[pos+3]=0;
      break;
    case 'W' : p_pars[pos+0]=p_pars[pos+3]=1; p_pars[pos+1]=p_pars[pos+2]=0;
      break;
    case 'S' : p_pars[pos+1]=p_pars[pos+2]=1; p_pars[pos+0]=p_pars[pos+3]=0;
      break;
    case 'Y' : p_pars[pos+1]=p_pars[pos+3]=1; p_pars[pos+0]=p_pars[pos+2]=0;
      break;
    case 'K' : p_pars[pos+2]=p_pars[pos+3]=1; p_pars[pos+0]=p_pars[pos+1]=0;
      break;
    case 'B' : p_pars[pos+1]=p_pars[pos+2]=p_pars[pos+3]=1; p_pars[pos+0]=0;
      break;
    case 'D' : p_pars[pos+0]=p_pars[pos+2]=p_pars[pos+3]=1; p_pars[pos+1]=0;
      break;
    case 'H' : p_pars[pos+0]=p_pars[pos+1]=p_pars[pos+3]=1; p_pars[pos+2]=0;
      break;
    case 'V' : p_pars[pos+0]=p_pars[pos+1]=p_pars[pos+2]=1; p_pars[pos+3]=0;
      break;
    case 'N' : case 'X' : case '?' : case 'O' : case '-' :
      p_pars[pos+0]=p_pars[pos+1]=p_pars[pos+2]=p_pars[pos+3]=1;break;
    default :
      {
	PhyML_Printf("\n. Unknown character state : %c\n",state);
	Exit("\n. Init failed (check the data type)\n");
	break;
      }
    }
}

/*********************************************************/



void Init_Tips_At_One_Site_AA_Float(char aa, int pos, phydbl *p_lk)
{
  int i;

  For(i,20) p_lk[pos+i] = .0;

  switch(aa){
  case 'A' : p_lk[pos+0]= 1.; break;/* Alanine */
  case 'R' : p_lk[pos+1]= 1.; break;/* Arginine */
  case 'N' : p_lk[pos+2]= 1.; break;/* Asparagine */
  case 'D' : p_lk[pos+3]= 1.; break;/* Aspartic acid */
  case 'C' : p_lk[pos+4]= 1.; break;/* Cysteine */
  case 'Q' : p_lk[pos+5]= 1.; break;/* Glutamine */
  case 'E' : p_lk[pos+6]= 1.; break;/* Glutamic acid */
  case 'G' : p_lk[pos+7]= 1.; break;/* Glycine */
  case 'H' : p_lk[pos+8]= 1.; break;/* Histidine */
  case 'I' : p_lk[pos+9]= 1.; break;/* Isoleucine */
  case 'L' : p_lk[pos+10]=1.; break;/* Leucine */
  case 'K' : p_lk[pos+11]=1.; break;/* Lysine */
  case 'M' : p_lk[pos+12]=1.; break;/* Methionine */
  case 'F' : p_lk[pos+13]=1.; break;/* Phenylalanin */
  case 'P' : p_lk[pos+14]=1.; break;/* Proline */
  case 'S' : p_lk[pos+15]=1.; break;/* Serine */
  case 'T' : p_lk[pos+16]=1.; break;/* Threonine */
  case 'W' : p_lk[pos+17]=1.; break;/* Tryptophan */
  case 'Y' : p_lk[pos+18]=1.; break;/* Tyrosine */
  case 'V' : p_lk[pos+19]=1.; break;/* Valine */

  case 'B' : p_lk[pos+2]= 1.; break;/* Asparagine */
  case 'Z' : p_lk[pos+5]= 1.; break;/* Glutamine */

  case 'X' : case '?' : case '-' : For(i,20) p_lk[pos+i] = 1.; break;
  default :
    {
      PhyML_Printf("\n. Unknown character state : %c\n",aa);
      Exit("\n. Init failed (check the data type)\n");
      break;
    }
  }
}

/*********************************************************/

void Init_Tips_At_One_Site_AA_Int(char aa, int pos, short int *p_pars)
{
  int i;

  For(i,20) p_pars[pos+i] = .0;

  switch(aa){
  case 'A' : p_pars[pos+0]  = 1; break;/* Alanine */
  case 'R' : p_pars[pos+1]  = 1; break;/* Arginine */
  case 'N' : p_pars[pos+2]  = 1; break;/* Asparagine */
  case 'D' : p_pars[pos+3]  = 1; break;/* Aspartic acid */
  case 'C' : p_pars[pos+4]  = 1; break;/* Cysteine */
  case 'Q' : p_pars[pos+5]  = 1; break;/* Glutamine */
  case 'E' : p_pars[pos+6]  = 1; break;/* Glutamic acid */
  case 'G' : p_pars[pos+7]  = 1; break;/* Glycine */
  case 'H' : p_pars[pos+8]  = 1; break;/* Histidine */
  case 'I' : p_pars[pos+9]  = 1; break;/* Isoleucine */
  case 'L' : p_pars[pos+10] = 1; break;/* Leucine */
  case 'K' : p_pars[pos+11] = 1; break;/* Lysine */
  case 'M' : p_pars[pos+12] = 1; break;/* Methionine */
  case 'F' : p_pars[pos+13] = 1; break;/* Phenylalanin */
  case 'P' : p_pars[pos+14] = 1; break;/* Proline */
  case 'S' : p_pars[pos+15] = 1; break;/* Serine */
  case 'T' : p_pars[pos+16] = 1; break;/* Threonine */
  case 'W' : p_pars[pos+17] = 1; break;/* Tryptophan */
  case 'Y' : p_pars[pos+18] = 1; break;/* Tyrosine */
  case 'V' : p_pars[pos+19] = 1; break;/* Valine */

  case 'B' : p_pars[pos+2]  = 1; break;/* Asparagine */
  case 'Z' : p_pars[pos+5]  = 1; break;/* Glutamine */

  case 'X' : case '?' : case '-' : For(i,20) p_pars[pos+i] = 1; break;
  default :			//printf("here\n");

    {
      PhyML_Printf("\n. Unknown character state : %c\n",aa);
      Exit("\n. Init failed (check the data type)\n");
      break;
    }
  }
}

/*********************************************************/

void Init_Tips_At_One_Site_Generic_Float(char *state, int ns, int state_len, int pos, phydbl *p_lk)
{
  int i;
  int state_int;
 
  For(i,ns) p_lk[pos+i] = 0.;

  if(Is_Ambigu(state,GENERIC,state_len)) For(i,ns) p_lk[pos+i] = 1.;
  else
    {
      char format[6];
      sprintf(format,"%%%dd",state_len);
      if(!sscanf(state,format,&state_int))
	{
	  PhyML_Printf("\n. state='%c'",state);
	  PhyML_Printf("\n. Err in file %s at line %d\n\n",__FILE__,__LINE__);
	  Warn_And_Exit("");
	}
      if(state_int > ns)
	{
	  PhyML_Printf("\n. %s %d cstate: %.2s istate: %d state_len: %d.\n",__FILE__,__LINE__,state,state_int,state_len);
	  PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
	  Warn_And_Exit("");	  
	}
      p_lk[pos+state_int] = 1.;
      /*       PhyML_Printf("\n. %s %d cstate: %.2s istate: %d state_len: %d ns: %d pos: %d",__FILE__,__LINE__,state,state_int,state_len,ns,pos); */
    }
}

/*********************************************************/

void Init_Tips_At_One_Site_Generic_Int(char *state, int ns, int state_len, int pos, short int *p_pars)
{
  int i;
  int state_int;

  For(i,ns) p_pars[pos+i] = 0;
  
  if(Is_Ambigu(state,GENERIC,state_len)) For(i,ns) p_pars[pos+i] = 1;
  else 
    {
      char format[6];
      sprintf(format,"%%%dd",state_len);
      if(!sscanf(state,format,&state_int))
	{
	  PhyML_Printf("\n. state='%c'",state);
	  PhyML_Printf("\n. Err in file %s at line %d\n\n",__FILE__,__LINE__);
	  Warn_And_Exit("");
	}
      if(state_int > ns)
	{
	  PhyML_Printf("\n. %s %d cstate: %.2s istate: %d state_len: %d.\n",__FILE__,__LINE__,state,state_int,state_len);
	  PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
	  Warn_And_Exit("");	  
	}
      p_pars[pos+state_int] = 1;
    }
}

/*********************************************************/

void Get_All_Partial_Lk_Scale(t_tree *tree, t_edge *b_fcus, t_node *a, t_node *d)
{
  if(d->tax) return;
  else Update_P_Lk(tree,b_fcus,d);
}

/*********************************************************/

void Post_Order_Lk(t_node *a, t_node *d, t_tree *tree)
{
  int i,dir;

  dir = -1;

  if(d->tax) return;
  else
    {
      For(i,3)
	{
	  if(d->v[i] != a)
	    Post_Order_Lk(d,d->v[i],tree);
	  else dir = i;
	}      
      Get_All_Partial_Lk_Scale(tree,d->b[dir],a,d);
    }
}


/*********************************************************/

void Pre_Order_Lk(t_node *a, t_node *d, t_tree *tree)
{
  int i;

  if(d->tax) return;
  else
    {
      For(i,3)
	{
	  if(d->v[i] != a)
	    {
	      Get_All_Partial_Lk_Scale(tree,d->b[i],d->v[i],d);
	      Pre_Order_Lk(d,d->v[i],tree);
	    }
	}
    }
}
/*********************************************************
 * Optimization code:
 * 0: do not optimize upper model, copy from upper to lower models
 * 1: optimize upper, copy from upper to lower
 * 2: optimize lower, do not copy from upper to lower
 * 3: do not optimize lower model, do not copy from upper to lower
 */
//added by Ken 17/1/2018
phydbl Lk_rep(option *io){
	phydbl replnL = 0.0;
	//copy parameters to submodels
	int i,j;
	For(i,io->ntrees){
		model* mod=io->mod_s[i];
		t_tree* tree=io->tree_s[i];
		tree->both_sides=io->both_sides;
		mod->update_eigen=io->mod->update_eigen;
		For(j,io->mod->nomega_part){ //copy omega values if specified
			if(io->mod->omega_part_opt[j]<2 || io->mod->optIter==0){
				mod->omega_part[j]=io->mod->omega_part[j];
				//printf("omega %lf\n",mod->omega_part[j]);
			}
		}
		if(io->mod->optFreq<2 || io->mod->optIter==0){ //copy frequencies if specified
			For(j,12){//assumes CF3X4 -\_(:/)_/-
				mod->base_freq[j]=io->mod->base_freq[j];
				mod->uns_base_freq[j]=io->mod->uns_base_freq[j];
			}
		}
		if(io->mod->optKappa<2 || io->mod->optIter==0){//copy kappa if specified
			mod->kappa=io->mod->kappa;
		}
		if(io->mod->opthotness || io->mod->optIter==0){ //copy hotness if specified
			For(j,io->mod->nhotness){
				if(io->mod->hoptindex[j]<2 || io->mod->optIter==0){
					mod->hotness[j]=io->mod->hotness[j];
				}
			}
		}
	}
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
		Lk(tree);
	}
	For(i,io->ntrees){
		replnL += io->tree_s[i]->c_lnL;
	}
	io->replnL=replnL;
	//if(io->mod->optDebug)printf("\nreplikelihood %lf",replnL);
	return replnL;
}



/*********************************************************/
phydbl Lk(t_tree *tree)
{
//printf("\nboth sides: %d %d",tree->both_sides,tree->mod->update_eigen);
  int br, n_patterns, n_edges;  
  n_edges=2*tree->n_otu-3;
  tree->old_lnL = tree->c_lnL;
 
  if(tree->rates && tree->rates->lk_approx == NORMAL)
  {
    tree->c_lnL = Dnorm_Multi_Given_InvCov_Det(tree->rates->u_cur_l, 
					       tree->rates->u_ml_l,
					       tree->rates->invcov,
					       tree->rates->covdet,
					       2*tree->n_otu-3,YES);
					       return tree->c_lnL;
  }
  
  n_patterns = tree->n_pattern;
  Set_Model_Parameters(tree->mod);
  #ifdef TIMES
  if((tree->rates) && (tree->rates->bl_from_rt)) RATES_Update_Cur_Bl(tree);
  if(tree->bl_from_node_stamps) TIMES_Bl_From_T(tree);
  #endif
 
  #if defined OMP || defined BLAS_OMP
  #pragma omp parallel for if(tree->mod->datatype==CODON && !tree->mod->splitByTree)
  #endif
  For(br,n_edges) Update_PMat_At_Given_Edge(tree->t_edges[br],tree);

 // For(br,2*tree->n_otu-3) printf("%lf\n",tree->t_edges[br]->l);
  //added and modified by Ken to start likelihood calculation at specified node
  int startnode = 0;
  if(tree->mod->whichrealmodel == HLP17){
	  startnode = tree->mod->startnode;
  }
  Post_Order_Lk(tree->noeud[startnode],tree->noeud[startnode]->v[0],tree);

  if(tree->both_sides)
    Pre_Order_Lk(tree->noeud[startnode], tree->noeud[startnode]->v[0], tree);

  tree->c_lnL             = .0;
  tree->sum_min_sum_scale = .0;

  if(tree->mod->whichrealmodel==HLP17){
	  Fill_UPP_root(tree,tree->noeud[startnode]->b[0]);
  }

	tree->c_lnL             = .0;
	tree->sum_min_sum_scale = .0;

  For(tree->curr_site,n_patterns){
	  if(tree->data->wght[tree->curr_site] > SMALL && tree->mod->whichrealmodel==HLP17) Lk_Core_UPP(tree->noeud[startnode]->b[0],tree,tree->noeud[startnode],tree->noeud[startnode]->b[0]->des_node);
	  else if(tree->data->wght[tree->curr_site] > SMALL) Lk_Core(tree->noeud[startnode]->b[0],tree);
  }
  Adjust_Min_Diff_Lk(tree);
  if(tree->mod->optDebug && tree->mod->whichrealmodel==HLP17){
#if defined OMP || defined BLAS_OMP
#pragma omp critical
#endif
	  printf("\nomega %d %lf %lf %lf %lf %lf %lf",tree->both_sides,tree->mod->omega_part[0],tree->mod-> kappa,tree->c_lnL,tree->mod->qmat_part[0][1],tree->mod->hotness[0],tree->mod->hotness[1]);
	  int i;
#if defined OMP || defined BLAS_OMP
#pragma omp critical
#endif
	  For(i,12)printf(" %lf",tree->mod->uns_base_freq[i]);
  }

  return tree->c_lnL;
}


void upAllPmats(t_tree *tree)
{
  int br, n_patterns, n_edges;
  n_edges=2*tree->n_otu-3;

  For(br,n_edges) Update_PMat_At_Given_Edge(tree->t_edges[br],tree);
}


/*********************************************************/
//Added by Ken 20/10/2016
void Get_UPP(t_node *a, t_node *d,  t_tree *tree)
{
  int i;

  if(d->tax) return;
  else{
	  //recurse
	  int desc[2];
	  int count=0;
      For(i,3){//figure out which nodes are the descendants
    	  if(d->v[i]->num != d->anc->num){
			  desc[count]=i;
			  count++;
    	  }
      }
      //fill in upp matrix
      if(a->num == tree->mod->startnode){//if next to the root
    	  Fill_UPP_root(tree,d->anc_edge);
      }
      Fill_UPP_single(tree,d->b[desc[0]]);
      Fill_UPP_single(tree,d->b[desc[1]]);

      Get_UPP(d,d->v[desc[0]],tree);
      Get_UPP(d,d->v[desc[1]],tree);
    }
}

void Fill_UPP_single(t_tree *tree,  t_edge *target)
{
/*	  anc
       |
	   |<- b
	   |
	   d
      / \
u_e->/   \ <-c_e
  	/     \
n_v1   n_v2
*/
//printf("about to do this\n");
t_node *d = target->anc_node;

  t_node *adjn, *targn;
  t_edge *b, *adje;
  phydbl padj_lk,ptar_lk;
  phydbl *p_lk,*p_lk_adj,*p_lk_tar;
  int *sum_scale, *sum_scale_v1, *sum_scale_v2, *sum_scale_adj, *sum_scale_tar, *sum_scale_upp;
  int sum_scale_adj_val, sum_scale_upp_val, sum_scale_tar_val, sum_scale_val;
  int i,j;
  int catg,site;
  int dir1,dir2,adjdir,tardir;
  int n_patterns;
  short int ambiguity_check_adj,ambiguity_check_tar;
  int state_adj,state_tar;
  int dim1, dim2, dim3;
  phydbl curr_scaler;
  int curr_scaler_pow, piecewise_scaler_pow;
  phydbl p_lk_lim_inf;
  phydbl smallest_p_lk;

  p_lk_lim_inf = (phydbl)P_LK_LIM_INF;

  dim1 = tree->mod->n_catg * tree->mod->ns;
  dim2 = tree->mod->ns;
  dim3 = tree->mod->ns * tree->mod->ns;

  state_adj = state_tar = -1;
  ambiguity_check_adj = ambiguity_check_tar = NO;
  sum_scale_adj_val = sum_scale_tar_val = 0;

  if(d->tax){
      PhyML_Printf("\n. t_node %d is a leaf...",d->num);
      PhyML_Printf("\n. Err in file %s at line %d\n\n",__FILE__,__LINE__);
      Warn_And_Exit("\n");
    }

  n_patterns = tree->n_pattern;

  b=d->anc_edge;
  dir1=dir2=-1;
  For(i,3){
	  if(d->b[i] != b){
		  (dir1<0)?(dir1=i):(dir2=i);
	  }
  }

  if(target->num == d->b[dir1]->num){
	  adjdir=dir2;
	  tardir=dir1;
  }else if(target->num == d->b[dir2]->num){
	  adjdir=dir1;
	  tardir=dir2;
  }else{
	  printf("didn't work\n");
	  printf("Edge %d error as descendant of node d: %d b: %d target_des: %d | %d %d | %d %d %d\n",target->num,d->num,b->num,target->des_node->num,dir1,dir2, d->b[0]->num, d->b[1]->num, d->b[2]->num);
	  exit(EXIT_FAILURE);
  }

  adje = d->b[adjdir];
  adjn = d->v[adjdir];
  targn = d->v[tardir];

  if(d->num == b->left->num){
      p_lk = b->p_lk_left;
    }
  else{
      p_lk = b->p_lk_rght;
    }

  //upp scaling
  sum_scale = target->sum_scale_upp;
  sum_scale_upp = b->sum_scale_upp;

  if(d->num == d->b[adjdir]->left->num){
        p_lk_adj = d->b[adjdir]->p_lk_rght;
        sum_scale_adj = d->b[adjdir]->sum_scale_rght;
      }else{
        p_lk_adj = d->b[adjdir]->p_lk_left;
        sum_scale_adj = d->b[adjdir]->sum_scale_left;
      }

  if(d->num == d->b[tardir]->left->num){
       p_lk_tar = d->b[tardir]->p_lk_rght;
       sum_scale_tar = d->b[tardir]->sum_scale_rght;
     }else{
       p_lk_tar = d->b[tardir]->p_lk_left;
       sum_scale_tar = d->b[tardir]->sum_scale_left;
     }

  phydbl tarsum=0;

      //int interest = 29;

#if defined OMP || BLAS_OMP
//#pragma omp parallel for if(!tree->mod->splitByTree)
#endif

  For(site,n_patterns){
      phydbl tsitesum=0.0;

      state_adj = -1;
      ambiguity_check_adj = NO;
      if(!tree->mod->s_opt->greedy){
	   /* n_v1 and n_v2 are tip nodes */
	      if(adjn->tax){
	      /* Is the state at this tip ambiguous? */
	    	  if(adjn->num == target->left->num){
	    		  printf("left taxa %d %s\n",b->num,d->name);
	    	      	 exit(EXIT_FAILURE);
	    	  }
	       ambiguity_check_adj = tree->data->c_seq[adjn->num]->is_ambigu[site];
	       if(ambiguity_check_adj == NO) state_adj = Get_State_From_P_Pars(adje->p_lk_tip_r,site*dim2,tree);
	      }
      }

      if(tree->mod->use_m4mod){
    	  ambiguity_check_adj = YES;
      }

      /* For all the rate classes */
      For(catg,tree->mod->n_catg){
	 smallest_p_lk  =  BIG;

	  /* For all the state at node d */
	  For(i,tree->mod->ns){
		  target->upp[site][i]=0.0;
	      /* n_v1 is a tip */
	      if((adjn->tax) && (!tree->mod->s_opt->greedy)){
		  if(ambiguity_check_adj == NO){
		      /* For the (non-ambiguous) state at node n_v1 */
			  target->upp[site][i] = d->b[adjdir]->bPmat_part[tree->mod->partIndex[site]][catg*dim3+i*dim2+state_adj]; // Ken 17/8/2016
		    }else{
		      For(j,tree->mod->ns){
		    	phydbl val = d->b[adjdir]->bPmat_part[tree->mod->partIndex[site]][catg*dim3+i*dim2+j]*(phydbl)adje->p_lk_tip_r[site*dim2+j];
		    	target->upp[site][i] += val;
		      }
		    }
		  }
	      /* n_v1 is an internal node */
	   else{
		    For(j,tree->mod->ns){
		    	phydbl val = d->b[adjdir]->bPmat_part[tree->mod->partIndex[site]][catg*dim3+i*dim2+j] * p_lk_adj[site*dim1+catg*dim2+j];
		    	target->upp[site][i] += val;
		    }
	   }
	      //add in upper likelihood of b
	   phydbl sum=0.0;
	   For(j,tree->mod->ns){ /* For all the states at node n_v2 */
          phydbl val = b->bPmat_part[tree->mod->partIndex[site]][catg*dim3+j*dim2+i]*b->upp[site][j];
          sum += val;
        }
	   target->upp[site][i] *= sum;


	     tsitesum += target->upp[site][i];


	   	   if(target->upp[site][i] < smallest_p_lk) smallest_p_lk = target->upp[site][i] ;
	   }//For(i,tree->mod->ns)

	   	  sum_scale_adj_val = (sum_scale_adj)?(sum_scale_adj[catg*n_patterns+site]):(0);
	  	  sum_scale_upp_val = (sum_scale_upp)?(sum_scale_upp[catg*n_patterns+site]):(0);
	  	  sum_scale[catg*n_patterns+site] = sum_scale_adj_val + sum_scale_upp_val;

	  	  if(smallest_p_lk < p_lk_lim_inf){
	  	      curr_scaler_pow = (int)(LOG(p_lk_lim_inf)-LOG(smallest_p_lk))/LOG2;
	  	      curr_scaler     = (phydbl)((unsigned long long)(1) << curr_scaler_pow);
	  	      sum_scale[catg*n_patterns+site] += curr_scaler_pow;
	  	      do{
	  		  piecewise_scaler_pow = MIN(curr_scaler_pow,63);
	  		  curr_scaler = (phydbl)((unsigned long long)(1) << piecewise_scaler_pow);
	  		  For(i,tree->mod->ns){
	  			  target->upp[site][i] *= curr_scaler;
	  		      if(target->upp[site][i] > BIG){
	  			  PhyML_Printf("\n. p_lk_lim_inf = %G smallest_p_lk = %G",p_lk_lim_inf,smallest_p_lk);
	  			  PhyML_Printf("\n. curr_scaler_pow = %d",curr_scaler_pow);
	  			  PhyML_Printf("\n. Err in file %s at line %d.",__FILE__,__LINE__);
	  			  Warn_And_Exit("\n");
	  			}
	  		  }
	  		  curr_scaler_pow -= piecewise_scaler_pow;
	  		}
	  	      while(curr_scaler_pow != 0);
	  	 }

	  }//For(catg,tree->mod->n_catg){
    }//For(n,sites
}

void Fill_UPP_root(t_tree *tree, t_edge *b)
{
/*	  anc
       |
	   |<- b
	   |
	   d
      / \
u_e->/   \ <-c_e
  	/     \
n_v1   n_v2
*/

	t_node *d = b->des_node;
	t_node *anc = b->anc_node;
  t_node *n_v1, *n_v2;
  phydbl p1_lk1,p2_lk2;
  phydbl *p_lk,*p_lk_v1,*p_lk_v2;
  int *sum_scale, *sum_scale_adj, *sum_scale_upp;
  int sum_scale_v1_val, sum_scale_v2_val;
  int i,j;
  int catg,site;
  int dir1,dir2,adjdir,tardir;
  int n_patterns;
  short int ambiguity_check_anc;
  int state_anc;
  int dim1, dim2, dim3;
  phydbl curr_scaler;
  int curr_scaler_pow, piecewise_scaler_pow;
  phydbl p_lk_lim_inf;
  phydbl smallest_p_lk;

  p_lk_lim_inf = (phydbl)P_LK_LIM_INF;

  dim1 = tree->mod->n_catg * tree->mod->ns;
  dim2 = tree->mod->ns;
  dim3 = tree->mod->ns * tree->mod->ns;

  state_anc = -1;
  ambiguity_check_anc = NO;

  sum_scale = b->sum_scale_upp;

  int sum_scale_adj_val;


  if(d->tax)
    {
      PhyML_Printf("\n. t_node %d is a leaf...",d->num);
      PhyML_Printf("\n. Err in file %s at line %d\n\n",__FILE__,__LINE__);
      Warn_And_Exit("\n");
    }

  n_patterns = tree->n_pattern;

  /* Get the partial likelihood vectors on edge b and the two pendant
     edges (i.e., the two other edges connected to d) */
  if(d == b->left){
      p_lk = b->p_lk_left;

    }
  else{
      p_lk = b->p_lk_rght;

    }

#if defined OMP || BLAS_OMP
//#pragma omp parallel for if(!tree->mod->splitByTree)
#endif

  For(site,tree->n_pattern){
      state_anc = -1;
      ambiguity_check_anc= NO;
      if(!tree->mod->s_opt->greedy){
	  /* n_v1 and n_v2 are tip nodes */
	  if(anc->tax){
	      /* Is the state at this tip ambiguous? */
	      ambiguity_check_anc = tree->data->c_seq[anc->num]->is_ambigu[site];
	      if(ambiguity_check_anc == NO) state_anc = Get_State_From_P_Pars(anc->b[0]->p_lk_tip_r,site*dim2,tree);
	    }
      }

      if(tree->mod->use_m4mod){
    	  ambiguity_check_anc = YES;
      }

      /* For all the rate classes */
      For(catg,tree->mod->n_catg){
	  smallest_p_lk  =  BIG;

	  /* For all the state at node d */
	  For(i,tree->mod->ns){
		  b->upp[site][i]=0.0;
	      /* n_v1 is a tip */
	      if((anc->tax) && (!tree->mod->s_opt->greedy)){
		    if(ambiguity_check_anc == NO){
		      /* For the (non-ambiguous) state at node n_v1 */
		    	phydbl state = 0.0;
		    	if(state_anc==i){state=1;}
			     b->upp[site][i] = state; // Ken 17/8/2016
		    }else{
		    	if(!tree->mod->rootpi) b->upp[site][i] = anc->b[0]->p_lk_tip_r[site*tree->mod->ns+i];//Changed 12/Feb/2018 by Ken
		    	else b->upp[site][i]=tree->mod->pi[i];
		    	//char* s1=mCalloc(4,sizeof(char));
		    	//Sprint_codon(s1,senseCodons[i]);
		    	// printf("\nUPP %d\t%d\t%s\t%lf",site,i,s1,b->upp[site][i]);

           //b->upp[site][i]=1/61;
		    }
		  }

	      //if not at root node, add in upper likelihood of b
	      if(anc->num != tree->mod->startnode){
	    	  printf("something weird happened\n");
	    	  exit(EXIT_FAILURE);
	      }
	   	   if(b->upp[site][i] < smallest_p_lk) smallest_p_lk = b->upp[site][i] ;
	    }
	  sum_scale[catg*tree->n_pattern+site]=0;
    }
  }
}
/*********************************************************/

void Get_Lhood(t_node *a, t_node *d,  t_tree *tree)
{
  int i;

  phydbl lk = Lk_At_Given_Edge(d->anc_edge,tree);
  if(d->tax) return;
  else{
	  //recurse
	  int desc[2];
	  int count=0;
      For(i,3){//figure out which nodes are the descendants
    	  if(d->v[i]->num != d->anc->num){
			  desc[count]=i;
			  count++;
    	  }
      }
      Get_Lhood(d,d->v[desc[0]],tree);
      Get_Lhood(d,d->v[desc[1]],tree);
    }
}
/*********************************************************/
void BSort(phydbl* probs, int* codons){
int i, j;
phydbl temp;
int tempc;
for (i = 0; i <61-1; i++){
     for (j = 0; j < 61 - i-1; j++){
          if (probs[j] < probs[j+1]){
               temp = probs[j+1];
               probs[j+1] = probs[j];
               probs[j] = temp;
               tempc = codons[j+1];
               codons[j+1] = codons[j];
               codons[j] = tempc;
          	  }
     	 }
	}
}
/*********************************************************/
char ambigAssign(char** codons,int ncodons, int site){
	char result;
	int i,a,c,g,t;
	a=c=g=t=0;
	For(i,ncodons){
		//printf("here: %c\n",codons[i][site]);
		if(codons[i][site] == 'A')a++;
		if(codons[i][site] == 'C')c++;
		if(codons[i][site] == 'G')g++;
		if(codons[i][site] == 'T')t++;
	}

	//https://www.bioinformatics.org/sms/iupac.html
	if(a && !c && !g && !t)result='A';
	if(!a && c && !g && !t)result='C';
	if(!a && !c && g && !t)result='G';
	if(!a && !c && !g && t)result='T';

	if(a && !c && g && !t)result='R';//R 	A or G
	if(!a && c && !g && t)result='Y';//Y 	C or T
	if(!a && c && g && !t)result='S';//S 	G or C
	if(a && !c && !g && t)result='W';//W 	A or T
	if(!a && !c && g && t)result='K';//K 	G or T
	if(a && c && !g && !t)result='M';//M 	A or C

	if(!a && c && g && t)result='B';//B 	C or G or T
	if(a && !c && g && t)result='D';//D 	A or G or T
	if(a && c && !g && t)result='H';//H 	A or C or T
	if(a && c && g && !t)result='V';//V 	A or C or G

	if(a && c && g && t)result='N';//N 	any base
	if(!a && !c && !g && !t){
		printf("\nSomething went wrong 0_0 in ASR %s",codons[0]);
		exit(EXIT_FAILURE);
	}
	return result;
}



phydbl logAdd(phydbl a, phydbl b){
	return a + log(1+exp(b-a));
}


char* printThreeNum(phydbl G1p){
	char* s1=malloc(4*sizeof(char));
	int per=roundf(G1p*100);
	//printf("%d\n",per);
	if(per == 100)sprintf(s1,"%d",per);
	else if(per >=10)sprintf(s1," %d",per);
	else sprintf(s1,"  %d",per);
	return s1;
}
/*********************************************************/
//Reconstruct most likely set of germline genotypes and alleles
void reconGermline(option *io, char* v, int site,FILE* f){
	int i,j,k,c;

	int matches[io->ntrees];
	int	imgt[io->ntrees];
	int nmatch=0;
	For(i,io->ntrees){//make a list of trees/models with specified V and IMGT site
		if(io->mod->optDebug)printf("\n%d\t%s\t%s\t%d\n",i,io->mod_s[i]->germlineV,v,site);
		if(strcmp(io->mod_s[i]->germlineV,v)==0){
			int found=0;
			For(j,io->mod_s[i]->init_len/3){
				if(io->mod_s[i]->imgt[j] == site){
					imgt[i]=j;
					found++;
				}
			}
			if(found>0)matches[nmatch++]=i;
			//else printf("IMGT site %d not found in %d!\n",site,i);
		}
	}
	if(!nmatch){
		printf("\nNo matches were found to %s at IMGT site %d!\n",v,site);
		return;
		//exit(EXIT_FAILURE);
	}
	phydbl** rprobs=mCalloc(nmatch,sizeof(phydbl*));
	For(i,nmatch)rprobs[i]=mCalloc(61,sizeof(phydbl));
	For(i,nmatch){//calculate relative probabilities of different codons at the root
			t_tree* tree = io->tree_s[matches[i]];
			model* mod = io->mod_s[matches[i]];
			t_node* root = tree->noeud[mod->startnode];
			int isite = imgt[matches[i]];
			int ori[61];
			For(j,61)ori[j]=root->b[0]->p_lk_tip_r[isite*tree->mod->ns+j];
			phydbl flk=Lk(tree); //full ML using starting sequence
			tree->curr_site=isite; //find the LK of this site under starting sequence
			phydbl slk = Lk_Core_UPP(root->b[0],tree,root,root->b[0]->des_node);
			phydbl blk=flk-slk;
			For(j,61)root->b[0]->p_lk_tip_r[isite*tree->mod->ns+j]=0;//set all codons to zero
			For(j,61){
				root->b[0]->p_lk_tip_r[isite*tree->mod->ns+j]=1;
				Fill_UPP_root(tree,root->b[0]);
				rprobs[i][j]=Lk_Core_UPP(root->b[0],tree,root,root->b[0]->des_node)+blk+log(0.5);//rprobs[i][j]
				//rprobs[i][j]=Lk(tree)+log(0.5);
				char s[4];
				Sprint_codon(s,io->senseCodons[j]);
				//printf("%d\t%s\t%lf\t%lf\n",isite,s,flk,rprobs[i][j]);
				root->b[0]->p_lk_tip_r[isite*tree->mod->ns+j]=0;
			}
			For(j,61)root->b[0]->p_lk_tip_r[isite*tree->mod->ns+j]=ori[j];
	}
	/*For(i,nmatch){//calculate relative probabilities of different codons at the root
		t_tree* tree = io->tree_s[matches[i]];
		model* mod = io->mod_s[matches[i]];
		t_node* root = tree->noeud[mod->startnode];
		int isite = imgt[i];
		int ori = Get_State_From_P_Pars(root->b[0]->p_lk_tip_r,isite*61,tree); //store original codon
		phydbl flk=Lk(tree);
		For(j,61)root->b[0]->p_lk_tip_r[isite*tree->mod->ns+j]=0;//set all codons to zero
		For(j,61){
			root->b[0]->p_lk_tip_r[isite*tree->mod->ns+j]=1;
			//rprobs[i][j]=Lk(tree)+log(0.5);
			rprobs[i][j]=Lk(tree);
			printf("%d\t%lf\t%lf\n",isite,flk,rprobs[i][j]);
			root->b[0]->p_lk_tip_r[isite*tree->mod->ns+j]=0;
		}
		root->b[0]->p_lk_tip_r[isite*tree->mod->ns+ori]=1; //restore original codon
		if(io->mod->optDebug)printf("%d\t%d\n",ori,Get_State_From_P_Pars(root->b[0]->p_lk_tip_r,isite*61,tree));
		if(io->mod->optDebug)printf("\n%d\t%d\t%lf\t%lf",i,j,flk,Lk(tree));
	}*/

	phydbl* recons=mCalloc(61*61,sizeof(phydbl));
	phydbl* mprobs=mCalloc(61,sizeof(phydbl));
	phydbl sum=0;

	For(i,61){
		for(j=i;j<61;j++){
			For(k,nmatch)recons[i*61+j] += logAdd(rprobs[k][j],rprobs[k][i]);
			if(i==0&&j==0)sum=recons[i*61+j];
			else sum=logAdd(sum,recons[i*61+j]);
			if(mprobs[i]==0)mprobs[i]=recons[i*61+j];
			else mprobs[i]=logAdd(mprobs[i],recons[i*61+j]);
			if(mprobs[j]==0)mprobs[j]=recons[i*61+j];
			else mprobs[j]=logAdd(mprobs[j],recons[i*61+j]);
		}
	}
	//Find ML genotypes
	int G1=-1;
	phydbl G1p=-DBL_MAX;
	For(i,61){
		for(j=i;j<61;j++){
			if(recons[i*61+j]>G1p){
				G1p=recons[i*61+j];
				G1=i*61+j;
			}
		}
	}
	int G2=-1;
		phydbl G2p=-DBL_MAX;
		For(i,61){
			for(j=i;j<61;j++){
				if(recons[i*61+j]>G2p && (i*61+j) != G1){
					G2p=recons[i*61+j];
					G2=i*61+j;
				}
			}
		}
		int G3=-1;
		phydbl G3p=-DBL_MAX;
		For(i,61){
			for(j=i;j<61;j++){
				if(recons[i*61+j]>G3p && (i*61+j) != G1&& (i*61+j) != G2){
					G3p=recons[i*61+j];
					G3=i*61+j;
				}
			}
		}

	int mi=(int)G1/61;
	int mj=G1%61;
	int m2i=(int)G2/61;
	int m2j=G2%61;
	int m3i=(int)G3/61;
	int m3j=G3%61;
	//Find number of genotypes with LRs close enough to the MLE
	int ci_size=0;
	For(i,61){
		for(j=i;j<61;j++){
			if(recons[i*61+j]>=(G1p-1.96))ci_size++;
		}
	}
	//place in CI intervals
	int* ci1 = mCalloc(ci_size,sizeof(int));
	int* ci2 = mCalloc(ci_size,sizeof(int));
	c=0;
	For(i,61){
		for(j=i;j<61;j++){
			if(recons[i*61+j]>=(G1p-1.96)){
				if(mprobs[i] > mprobs[j]){
					ci1[c]=i;
					ci2[c]=j;
				}else if(mprobs[i] < mprobs[j]){
					ci1[c]=j;
					ci2[c]=i;
				}else{
					if(i!=j && (i*61+j) != G1)printf("Exactly equal marginals for a plausible genotype. Might mess up CI: %d\t%d",i,j);
					ci1[c]=i;
					ci2[c]=j;
				}
				char s1[4];
				char s2[4];
				Sprint_codon(s1,io->senseCodons[ci1[c]]);
				Sprint_codon(s2,io->senseCodons[ci2[c]]);
				//printf("%s\t%s\t%lf\t%lf\t%lf\t%lf\t%d\t%d\t%d\n",s1,s2,mprobs[ci1[c]],mprobs[ci2[c]],recons[i*61+j],G1p-1.96,i,j,ci1[c]);
				c++;
			}
		}
	}
	char** codonset1=malloc(c*sizeof(char*));
	char** codonset2=malloc(c*sizeof(char*));
	For(i,c){
		char* s1=mCalloc(4,sizeof(char));
		char* s2=mCalloc(4,sizeof(char));
		Sprint_codon(s1,io->senseCodons[ci1[i]]);
		Sprint_codon(s2,io->senseCodons[ci2[i]]);
		codonset1[i]=s1;
		codonset2[i]=s2;
	}

	char cigen[8];
	For(i,3)cigen[i]=ambigAssign(codonset1,c,i);
	cigen[3]=':';
	for(i=4;i<7;i++)cigen[i]=ambigAssign(codonset2,c,i-4);
	cigen[7]='\0';


	char s1[4];
	char s2[4];
	char s2_1[4];
	char s2_2[4];
	char s3_1[4];
	char s3_2[4];
	Sprint_codon(s1,io->senseCodons[mi]);
	Sprint_codon(s2,io->senseCodons[mj]);
	Sprint_codon(s2_1,io->senseCodons[m2i]);
	Sprint_codon(s2_2,io->senseCodons[m2j]);
	Sprint_codon(s3_1,io->senseCodons[m3i]);
	Sprint_codon(s3_2,io->senseCodons[m3j]);
	fprintf(f,"%s\t%d\t%d\t%s\t%s:%s\t%s:%s:%.2lf\t%s:%s:%.2lf\n",v,site,nmatch,cigen,s1,s2,s2_1,s2_2,G1p-G2p,s3_1,s3_2,G1p-G3p);
	//printf("Genotype: %s\n",cigen);
}

/*********************************************************/
//Calculate the marginal probabilities of each codon at each site for node d
/*	  anc
       |
	   |<- b
	   |
	   d
      / \
u_e->/   \ <-c_e
  	/     \
n_v1   n_v2
*/

phydbl ASR_Core(t_edge *b, t_tree *tree, t_node *anc, t_node *d)
{
	phydbl log_site_lk;
	  phydbl site_lk_cat, site_lk;
	  int fact_sum_scale;
	  phydbl max_sum_scale;
	  phydbl sum;
	  int ambiguity_check,state;
	  int catg,ns,k,l,site;
	  int dim1,dim2,dim3;
	  int *sum_scale_down_cat,*sum_scale_upp_cat,*sum_scale_down,*sum_scale_upp;
	  phydbl multiplier;
	  int exponent, piecewise_exponent;
	  phydbl tmp;
	  phydbl logbig;
	  phydbl inv_site_lk;

	  logbig = LOG((phydbl)BIG);
	  //printf("logbig: %lf\n",logbig);

	  dim1 = tree->mod->n_catg * tree->mod->ns;
	  dim2 = tree->mod->ns;
	  dim3 = tree->mod->ns * tree->mod->ns;

	  log_site_lk     = .0;
	  site_lk         = .0;
	  site_lk_cat     = .0;
	  ambiguity_check = -1;
	  state           = -1;
	  site            = tree->curr_site;
	  ns              = tree->mod->ns;

	  phydbl *p_lk;

	  if(d->num == b->left->num){
	      p_lk = b->p_lk_left;
	      sum_scale_down = b->sum_scale_left;
	      sum_scale_down_cat = b->sum_scale_left_cat;
	     if((d->tax) && (!tree->mod->s_opt->greedy)){
	    	 if(site==16){printf("left taxa %d %s\n",b->num,d->name);}
	    	 exit(EXIT_FAILURE);
	        ambiguity_check = tree->data->c_seq[d->num]->is_ambigu[site];
	        if(!ambiguity_check)state = Get_State_From_P_Pars(b->p_lk_tip_l,site*dim2,tree);
	      }
	    }
	  else{
	      p_lk = b->p_lk_rght;
	      sum_scale_down = b->sum_scale_rght;
	      sum_scale_down_cat = b->sum_scale_rght_cat;
	      if((d->tax) && (!tree->mod->s_opt->greedy)){
	        ambiguity_check = tree->data->c_seq[d->num]->is_ambigu[site];
	        if(!ambiguity_check)state = Get_State_From_P_Pars(b->p_lk_tip_r,site*dim2,tree);
	      }
	    }
	  sum_scale_upp_cat = b->sum_scale_upp_cat;
	  sum_scale_upp = b->sum_scale_upp;

	  if(tree->mod->use_m4mod) ambiguity_check = 1;


	  phydbl* probc = mCalloc(ns,sizeof(phydbl));
	  For(catg,tree->mod->n_catg){
	      site_lk_cat = .0;
	      /* b is an external edge */
	      if((d->tax) && (!tree->mod->s_opt->greedy)){
	          /* If the character observed at the tip is NOT ambiguous: ns x 1 terms to consider */
	          if(!ambiguity_check){
	        	 // printf("a\n");
	              sum = .0;
	              For(l,ns){
	                phydbl val = b->bPmat_part[tree->mod->partIndex[site]][catg*dim3+l*dim2+state]*b->upp[site][l];
	                sum+=val;
	              }
	              probc[state]=sum;
	              site_lk_cat = sum;
	            }
	          /* If the character observed at the tip is ambiguous: ns x ns terms to consider */
	          else{
	        	//  printf("b\n");
	              For(l,ns){
	                  sum = .0;
	                      For(k,ns){
	                        phydbl val = b->bPmat_part[tree->mod->partIndex[site]][catg*dim3+k*dim2+l]*b->p_lk_tip_r[site*dim2+l]*b->upp[site][k];
	                        sum+=val;
	                      }
	                      probc[l]=sum;
	                      site_lk_cat += sum;
	                    }
	            }
	        }
	      /* b is an internal edge: ns x ns terms to consider */
	      else{
	    	  //printf("c\n");
	          For(l,ns){
	        	  	sum = .0;
	                For(k,ns){
	                  	phydbl val = b->bPmat_part[tree->mod->partIndex[site]][catg*dim3+k*dim2+l]*p_lk[site*dim1+catg*dim2+l]*b->upp[site][k];//Modified by Ken 17/8/2016
	                    sum+=val;
	                }
	                probc[l]=sum;
	                site_lk_cat+=sum;
	           }
	        }
	      tree->site_lk_cat[catg] = site_lk_cat;
	    }

	  phydbl max=-DBL_MAX;
	  int maxc=-1;
	  int i;
	  For(i,ns){
		  char s1[4];
		  Sprint_codon(s1,tree->io->senseCodons[i]);
		  /*if(site==0 && (i==0 || i==8)){
			  printf("\n%d\t%lf\t%s\t%lf\t%lf\t%lf\t%lf\t%lf",b->num,b->l,s1,probc[i]/tree->site_lk_cat[0],b->bPmat_part[0][0],b->bPmat_part[0][8],b->bPmat_part[0][8*61],b->bPmat_part[0][8*61+8]);
		  }*/
		  if(probc[i]>max){
			  max=probc[i];
			  maxc=i;
		  }
	  }
	 /* int sense[61];
	  For(i,61){
		  probc[i]=probc[i]/tree->site_lk_cat[0];
		  sense[i]=tree->io->senseCodons[i];
	  }
	  BSort(probc,sense);
	  int ncodons=0;
	  phydbl probsum=0;
	  For(i,61){
		  char s1[4];
		  Sprint_codon(s1,sense[i]);
		 // printf("\npossible: %d\t%d\t%s\t%lf",b->num,site,s1,probc[i]);
		  probsum+=probc[i];
		  ncodons++;
		  if(probsum>=tree->mod->ASRcut)break;
	  }*/

	  int ncodons=0;
	  For(i,61){
		 // printf("\n%d\t%lf\t%lf",i,probc[i],probc[i]/tree->site_lk_cat[0]);
		  if(probc[i]==0)probc[i]=-DBL_MAX;
		  else{
			  probc[i]=log(probc[i]);
			  //printf("\n%d\t%d\t%lf\t%lf",site,i,probc[i],log(max)-1.96);
		  }
		  if(probc[i]>=(log(max)-1.96))ncodons++;
	  }

	  char** codonset=malloc(ncodons*sizeof(char*));
	  int c=0;
	  //For(i,ncodons){
	  For(i,61){
		  if(probc[i]>=(log(max)-1.96)){
			  codonset[c]=mCalloc(4,sizeof(char));
		  	  char* s1=mCalloc(4,sizeof(char));
		  	  Sprint_codon(s1,tree->io->senseCodons[i]);
		  	  strcpy(codonset[c],s1);
		  	//printf("\n%c\t%c\t%c\t%s",codonset[c][0],codonset[c][1],codonset[c][2],codonset[c]);
		  	c++;
		  }
	  }
	  For(i,3){
		  //printf("\n%s\t%d\n",codonset[0],ncodons);
		  tree->mod->mlCodon[b->num][site*3+i]=ambigAssign(codonset,ncodons,i);
	  }


	  //tree->mod->mlASR[b->num][site]=maxc;
	  //tree->mod->probASR[b->num][site]=max/tree->site_lk_cat[0];
	 // printf("%d\t%lf\t%lf\n",maxc,max,tree->site_lk_cat[0]);
	  if(site==0){
		//  printf("labels: %d\n",b->n_labels);
	  	  if(!b->n_labels ){
		  	  b->labels=mCalloc(1,sizeof(char*));
		  	  b->labels[0]=mCalloc(T_MAX_LABEL,sizeof(char));
		  	  b->n_labels++;
	  	  }
	  	  char s1[4];
	  	  sprintf(s1,"%d",b->num);
	  	  strcpy(b->labels[0],s1);
	  	  //printf("edge and des: %d\t%d\n",b->num,b->des_node->num);
	  }
	  free(probc);
	 // exit(EXIT_FAILURE);

	  max_sum_scale =  (phydbl)BIG;
	  For(catg,tree->mod->n_catg){

	      sum_scale_down_cat[catg] = //evaluates as false when sum_scale_down is NULL
	        (sum_scale_down)?
	        (sum_scale_down[catg*tree->n_pattern+site]):
	        (0.0);

	      sum_scale_upp_cat[catg] =
	        (sum_scale_upp)?
	        (sum_scale_upp[catg*tree->n_pattern+site]):
	        (0.0);

	      sum = sum_scale_down_cat[catg] + sum_scale_upp_cat[catg];

	      if(sum < .0){
	    	  PhyML_Printf("\n. sum = %G",sum);
	    	  PhyML_Printf("\n. Err in file %s at line %d\n\n",__FILE__,__LINE__);
	    	  Warn_And_Exit("\n");
	      }
	      tmp = sum + (logbig - LOG(tree->site_lk_cat[catg]))/(phydbl)LOG2;
	      if(tmp < max_sum_scale) max_sum_scale = tmp; /* min of the maxs */
	    }
	  fact_sum_scale = (int)(max_sum_scale / 2);

	  // Apply scaling factors
	  For(catg,tree->mod->n_catg){
	      exponent = -(sum_scale_down_cat[catg]+sum_scale_upp_cat[catg])+fact_sum_scale;
	      site_lk_cat = tree->site_lk_cat[catg];
	      if(exponent > 0){
	    	  do{
	    		  piecewise_exponent = MIN(exponent,63);
	    		  multiplier = (phydbl)((unsigned long long)(1) << piecewise_exponent);
	    		  site_lk_cat *= multiplier;
	    		  exponent -= piecewise_exponent;
	    	  }
	    	  while(exponent != 0);
	      }
	      else{
	    	  do{
	    		  piecewise_exponent = MAX(exponent,-63);
	    		  multiplier = 1. / (phydbl)((unsigned long long)(1) << -piecewise_exponent);
	    		  site_lk_cat *= multiplier;
	    		  exponent -= piecewise_exponent;
	    	  }
	    	  while(exponent != 0);
	      }

	      if(isinf(site_lk_cat)){
	    	  PhyML_Printf("\n+ site=%4d site_lk_cat=%G sum_scale=%d max=%G fact=%d expo=%d dbl=%G",
			       tree->curr_site,
			       tree->site_lk_cat[catg],
			       sum_scale_down_cat[catg]+sum_scale_upp_cat[catg],
			       max_sum_scale,
			       fact_sum_scale,
			       -(sum_scale_down_cat[catg]+sum_scale_upp_cat[catg])+fact_sum_scale,
			       (double)tree->site_lk_cat[catg] * pow(2.,-(sum_scale_down_cat[catg]+sum_scale_upp_cat[catg])+fact_sum_scale));
	    	  site_lk_cat = BIG / 10;
	      }

	      if(site_lk_cat < SMALL){
	    	  site_lk_cat = SMALL;
	    	  printf("site lk = 0 %d %d %lf %d %d %d.\nSetting underflow to DBL_MIN. See Manual.\n",anc->num,d->num,site_lk_cat,site,sum_scale_upp_cat[catg],sum_scale_down_cat[catg]);
	    	  printf("%lf %lf",tree->mod->omega_part[0],tree->mod->kappa);
	       	  printf("\n");
	      }
	      tree->site_lk_cat[catg] = site_lk_cat;
	    }


	  site_lk = .0;
	  if(tree->mod->datatype!=CODON ) For(catg,tree->mod->n_catg) site_lk += tree->site_lk_cat[catg]*tree->mod->gamma_r_proba[catg];
	  else //!< Added by Marcelo.
	  {
	    if(tree->mod->n_catg>1 && tree->mod->n_w_catg==1) For(catg,tree->mod->n_catg) site_lk += tree->site_lk_cat[catg]*tree->mod->gamma_r_proba[catg];
	    //else For(catg,tree->mod->n_catg) site_lk += tree->site_lk_cat[catg]*tree->mod->prob_omegas[catg];
	    else site_lk = tree->site_lk_cat[0]; //changed by Ken 7/2/2018
	  }

	  /* The substitution model does include invariable sites */
	   if(tree->mod->invar)
	     {
	       /* The site is invariant */
	       if(tree->data->invar[site] > -0.5)
	         {
	 	  /* Multiply P(D|r=0) by 2^(fact_sum_scale) */

	 	  inv_site_lk = tree->mod->pi[tree->data->invar[site]];
	 	  exponent = fact_sum_scale;
	 	  do
	 	    {
	 	      piecewise_exponent = MIN(exponent,63);
	 	      multiplier = (phydbl)((unsigned long long)(1) << piecewise_exponent);
	 	      inv_site_lk *= multiplier;
	 	      exponent -= piecewise_exponent;
	 	    }
	 	  while(exponent != 0);
	 	  inv_site_lk *= tree->mod->pinvar;
	 	  /* Update the value of site_lk */
	 	  site_lk *= (1. - tree->mod->pinvar);
	 	  site_lk += inv_site_lk;
	         }
	       else
	         {
	           /* Same formula as above with P(D | subs, rate = 0) = 0 */
	 	  site_lk *= (1. - tree->mod->pinvar);
	         }
	     }

	   //printf("site_lk: %lf\n",site_lk);
	   //printf("tree site lk: %lf\n",tree->site_lk_cat[0]);

	  log_site_lk = LOG(site_lk) - (phydbl)LOG2 * fact_sum_scale;
	  For(catg,tree->mod->n_catg) tree->log_site_lk_cat[catg][site] = LOG(tree->site_lk_cat[catg]) - (phydbl)LOG2 * fact_sum_scale;
	  if(isinf(log_site_lk) || isnan(log_site_lk))
	    {
	      PhyML_Printf("\n. site = %d",site);
	      PhyML_Printf("\n. invar = %f",tree->data->invar[site]);
	      PhyML_Printf("\n. scale_down = %d scale_upp = %d %lf",sum_scale_down[0],sum_scale_upp[0],fact_sum_scale);
	      PhyML_Printf("\n. Lk = %G LOG(Lk) = %f < %G",site_lk,log_site_lk,-BIG);
	      PhyML_Printf("\n. Err in file %s at line %d\n\n",__FILE__,__LINE__);
	      PhyML_Printf("\n. %d %d %d\n\n",b->num,d->num,anc->num);
	      Warn_And_Exit("\n");
	    }
	  tree->cur_site_lk[site] = log_site_lk;
	  /* Multiply log likelihood by the number of times this site pattern is found in the data */
	  tree->c_lnL_sorted[site] = tree->data->wght[site]*log_site_lk;
	  tree->c_lnL += tree->data->wght[site]*log_site_lk;
	  return log_site_lk;
}

phydbl ASR_Core_root(t_edge *b, t_tree *tree, t_node *anc, t_node *d)
{
	phydbl log_site_lk;
	  phydbl site_lk_cat, site_lk;
	  int fact_sum_scale;
	  phydbl max_sum_scale;
	  phydbl sum;
	  int ambiguity_check,state;
	  int catg,ns,k,l,site;
	  int dim1,dim2,dim3;
	  int *sum_scale_down_cat,*sum_scale_upp_cat,*sum_scale_down,*sum_scale_upp;
	  phydbl multiplier;
	  int exponent, piecewise_exponent;
	  phydbl tmp;
	  phydbl logbig;
	  phydbl inv_site_lk;

	  logbig = LOG((phydbl)BIG);
	  //printf("logbig: %lf\n",logbig);

	  dim1 = tree->mod->n_catg * tree->mod->ns;
	  dim2 = tree->mod->ns;
	  dim3 = tree->mod->ns * tree->mod->ns;

	  log_site_lk     = .0;
	  site_lk         = .0;
	  site_lk_cat     = .0;
	  ambiguity_check = -1;
	  state           = -1;
	  site            = tree->curr_site;
	  ns              = tree->mod->ns;

	  phydbl *p_lk;

	  if(d->num == b->left->num){
	      p_lk = b->p_lk_left;
	      sum_scale_down = b->sum_scale_left;
	      sum_scale_down_cat = b->sum_scale_left_cat;
	     if((d->tax) && (!tree->mod->s_opt->greedy)){
	    	 if(site==16){printf("left taxa %d %s\n",b->num,d->name);}
	    	 exit(EXIT_FAILURE);
	        ambiguity_check = tree->data->c_seq[d->num]->is_ambigu[site];
	        if(!ambiguity_check)state = Get_State_From_P_Pars(b->p_lk_tip_l,site*dim2,tree);
	      }
	    }
	  else{
	      p_lk = b->p_lk_rght;
	      sum_scale_down = b->sum_scale_rght;
	      sum_scale_down_cat = b->sum_scale_rght_cat;
	      if((d->tax) && (!tree->mod->s_opt->greedy)){
	        ambiguity_check = tree->data->c_seq[d->num]->is_ambigu[site];
	        if(!ambiguity_check)state = Get_State_From_P_Pars(b->p_lk_tip_r,site*dim2,tree);
	      }
	    }
	  sum_scale_upp_cat = b->sum_scale_upp_cat;
	  sum_scale_upp = b->sum_scale_upp;

	  if(tree->mod->use_m4mod) ambiguity_check = 1;


	  phydbl* probc = mCalloc(ns,sizeof(phydbl));
	  For(catg,tree->mod->n_catg){
	      site_lk_cat = .0;
	      /* b is an external edge */
	      if((d->tax) && (!tree->mod->s_opt->greedy)){
	          /* If the character observed at the tip is NOT ambiguous: ns x 1 terms to consider */
	          if(!ambiguity_check){
	        	 // printf("a\n");
	              sum = .0;
	              For(l,ns){//re-arrange direction of change
	                phydbl val = b->bPmat_part[tree->mod->partIndex[site]][catg*dim3+state*dim2+l]*b->upp[site][state];
	                sum+=val;
	              }
	              probc[state]=sum;
	              site_lk_cat = sum;
	            }
	          /* If the character observed at the tip is ambiguous: ns x ns terms to consider */
	          else{
	        	 // printf("b\n");
	              For(l,ns){//direction is changed from l to k
	                  sum = .0;
	                      For(k,ns){
	                        phydbl val = b->bPmat_part[tree->mod->partIndex[site]][catg*dim3+l*dim2+k]*b->p_lk_tip_r[site*dim2+k]*b->upp[site][l];
	                        sum+=val;
	                      }
	                      probc[l]=sum;
	                      site_lk_cat += sum;
	                    }
	            }
	        }
	      /* b is an internal edge: ns x ns terms to consider */
	      else{
	    	//  printf("c\n");
	          For(l,ns){
	              sum = .0;
	                  For(k,ns){
	                    	phydbl val = b->bPmat_part[tree->mod->partIndex[site]][catg*dim3+l*dim2+k]*p_lk[site*dim1+catg*dim2+k]*b->upp[site][l];//Modified by Ken 17/8/2016
	                      sum+=val;
	                    }
	                  probc[l]=sum;
	                  site_lk_cat+=sum;
	                }
	        }
	      tree->site_lk_cat[catg] = site_lk_cat;
	    }

	  phydbl max=-DBL_MAX;
	  int maxc=-1;
	  int i;
	  For(i,ns){
		  char s1[4];
		  		  Sprint_codon(s1,tree->io->senseCodons[i]);
		  		  //printf("%s\t%lf\n",s1,probc[i]/tree->site_lk_cat[0]);
		  		  if(site==0 && (i==0 || i==8)){
		  			  //printf("%d\t%d\t%s\t%lf\n",-1,i,s1,probc[i]/tree->site_lk_cat[0]);
					  //printf("\n%d\t%s\t%lf\t%lf\t%lf\t%lf\t%lf",-1,s1,probc[i]/tree->site_lk_cat[0],b->bPmat_part[0][0],b->bPmat_part[0][8],b->bPmat_part[0][8*61],b->bPmat_part[0][8*61+8]);

		  		  }
		  //printf("%lf\t%lf\n",probc[i],max);
		  if(probc[i]>max){
			  max=probc[i];
			  maxc=i;
		  }
	  }
	  tree->mod->mlASR[tree->mod->nedges][site]=maxc;
	  tree->mod->probASR[tree->mod->nedges][site]=max/tree->site_lk_cat[0];
	  /*int sense[61];
	  	  For(i,61){
	  		  probc[i]=probc[i]/tree->site_lk_cat[0];
	  		  sense[i]=tree->io->senseCodons[i];
	  	  }
	  	  BSort(probc,sense);
	  	  int ncodons=0;
	  	  phydbl probsum=0;
	  	  For(i,61){
	  		  char s1[4];
	  		  Sprint_codon(s1,sense[i]);
	  		  printf("\npossible: %d\t%d\t%s\t%lf",b->num,site,s1,probc[i]);
	  		  probsum+=probc[i];
	  		  ncodons++;
	  		  if(probsum>=tree->mod->ASRcut)break;
	  	  }
	  	  char** codonset=malloc(ncodons*sizeof(char*));
	  	  For(i,ncodons){
	  		  char* s1=mCalloc(4,sizeof(char));
	  		  Sprint_codon(s1,sense[i]);
	  		  codonset[i]=s1;
	  		  //printf("\n%c\t%c\t%c",codonset[i][0],codonset[i][1],codonset[i][2]);

	  	  }
	  	  //printf("\n%c",ambigAssign(codonset,ncodons,1));
	  	  //printf("\n%c\t%c\t%c",ambigAssign(codonset,ncodons,0),ambigAssign(codonset,ncodons,1),ambigAssign(codonset,ncodons,2));
	  	  For(i,3)tree->mod->mlCodon[tree->mod->nedges][site*3+i]=ambigAssign(codonset,ncodons,i);*/


		  int ncodons=0;
		  For(i,61){
			 // printf("\n%d\t%lf\t%lf",i,probc[i],probc[i]/tree->site_lk_cat[0]);
			  if(probc[i]==0)probc[i]=-DBL_MAX;
			  else{
				  probc[i]=log(probc[i]);
				  //printf("\n%d\t%d\t%lf\t%lf",site,i,probc[i],log(max)-1.96);
			  }
			  if(probc[i]>=(log(max)-1.96))ncodons++;
		  }

		  char** codonset=malloc(ncodons*sizeof(char*));
		  int c=0;
		  //For(i,ncodons){
		  For(i,61){
			  if(probc[i]>=(log(max)-1.96)){
				  codonset[c]=mCalloc(4,sizeof(char));
			  	  char* s1=mCalloc(4,sizeof(char));
			  	  Sprint_codon(s1,tree->io->senseCodons[i]);
			  	  strcpy(codonset[c],s1);
			  	//printf("\n%c\t%c\t%c\t%s",codonset[c][0],codonset[c][1],codonset[c][2],codonset[c]);
			  	c++;
			  }
		  }
		  For(i,3){
			 // printf("\n%s\t%d\n",codonset[0],ncodons);
			  For(i,3)tree->mod->mlCodon[tree->mod->nedges][site*3+i]=ambigAssign(codonset,ncodons,i);
		  }
	  free(probc);

	  max_sum_scale =  (phydbl)BIG;
	  For(catg,tree->mod->n_catg){

	      sum_scale_down_cat[catg] = //evaluates as false when sum_scale_down is NULL
	        (sum_scale_down)?
	        (sum_scale_down[catg*tree->n_pattern+site]):
	        (0.0);

	      sum_scale_upp_cat[catg] =
	        (sum_scale_upp)?
	        (sum_scale_upp[catg*tree->n_pattern+site]):
	        (0.0);

	      sum = sum_scale_down_cat[catg] + sum_scale_upp_cat[catg];

	      if(sum < .0){
	    	  PhyML_Printf("\n. sum = %G",sum);
	    	  PhyML_Printf("\n. Err in file %s at line %d\n\n",__FILE__,__LINE__);
	    	  Warn_And_Exit("\n");
	      }
	      tmp = sum + (logbig - LOG(tree->site_lk_cat[catg]))/(phydbl)LOG2;
	      if(tmp < max_sum_scale) max_sum_scale = tmp; /* min of the maxs */
	    }
	  fact_sum_scale = (int)(max_sum_scale / 2);

	  // Apply scaling factors
	  For(catg,tree->mod->n_catg){
	      exponent = -(sum_scale_down_cat[catg]+sum_scale_upp_cat[catg])+fact_sum_scale;
	      site_lk_cat = tree->site_lk_cat[catg];
	      if(exponent > 0){
	    	  do{
	    		  piecewise_exponent = MIN(exponent,63);
	    		  multiplier = (phydbl)((unsigned long long)(1) << piecewise_exponent);
	    		  site_lk_cat *= multiplier;
	    		  exponent -= piecewise_exponent;
	    	  }
	    	  while(exponent != 0);
	      }
	      else{
	    	  do{
	    		  piecewise_exponent = MAX(exponent,-63);
	    		  multiplier = 1. / (phydbl)((unsigned long long)(1) << -piecewise_exponent);
	    		  site_lk_cat *= multiplier;
	    		  exponent -= piecewise_exponent;
	    	  }
	    	  while(exponent != 0);
	      }

	      if(isinf(site_lk_cat)){
	    	  PhyML_Printf("\n+ site=%4d site_lk_cat=%G sum_scale=%d max=%G fact=%d expo=%d dbl=%G",
			       tree->curr_site,
			       tree->site_lk_cat[catg],
			       sum_scale_down_cat[catg]+sum_scale_upp_cat[catg],
			       max_sum_scale,
			       fact_sum_scale,
			       -(sum_scale_down_cat[catg]+sum_scale_upp_cat[catg])+fact_sum_scale,
			       (double)tree->site_lk_cat[catg] * pow(2.,-(sum_scale_down_cat[catg]+sum_scale_upp_cat[catg])+fact_sum_scale));
	    	  site_lk_cat = BIG / 10;
	      }

	      if(site_lk_cat < SMALL){
	    	  site_lk_cat = SMALL;
	    	  printf("site lk = 0 %d %d %lf %d %d %d.\nSetting underflow to DBL_MIN. See Manual.\n",anc->num,d->num,site_lk_cat,site,sum_scale_upp_cat[catg],sum_scale_down_cat[catg]);
	    	  printf("%lf %lf",tree->mod->omega_part[0],tree->mod->kappa);
	       	  printf("\n");
	      }
	      tree->site_lk_cat[catg] = site_lk_cat;
	    }


	  site_lk = .0;
	  if(tree->mod->datatype!=CODON ) For(catg,tree->mod->n_catg) site_lk += tree->site_lk_cat[catg]*tree->mod->gamma_r_proba[catg];
	  else //!< Added by Marcelo.
	  {
	    if(tree->mod->n_catg>1 && tree->mod->n_w_catg==1) For(catg,tree->mod->n_catg) site_lk += tree->site_lk_cat[catg]*tree->mod->gamma_r_proba[catg];
	    //else For(catg,tree->mod->n_catg) site_lk += tree->site_lk_cat[catg]*tree->mod->prob_omegas[catg];
	    else site_lk = tree->site_lk_cat[0]; //changed by Ken 7/2/2018
	  }

	  /* The substitution model does include invariable sites */
	   if(tree->mod->invar)
	     {
	       /* The site is invariant */
	       if(tree->data->invar[site] > -0.5)
	         {
	 	  /* Multiply P(D|r=0) by 2^(fact_sum_scale) */

	 	  inv_site_lk = tree->mod->pi[tree->data->invar[site]];
	 	  exponent = fact_sum_scale;
	 	  do
	 	    {
	 	      piecewise_exponent = MIN(exponent,63);
	 	      multiplier = (phydbl)((unsigned long long)(1) << piecewise_exponent);
	 	      inv_site_lk *= multiplier;
	 	      exponent -= piecewise_exponent;
	 	    }
	 	  while(exponent != 0);
	 	  inv_site_lk *= tree->mod->pinvar;
	 	  /* Update the value of site_lk */
	 	  site_lk *= (1. - tree->mod->pinvar);
	 	  site_lk += inv_site_lk;
	         }
	       else
	         {
	           /* Same formula as above with P(D | subs, rate = 0) = 0 */
	 	  site_lk *= (1. - tree->mod->pinvar);
	         }
	     }

	   //printf("site_lk: %lf\n",site_lk);
	   //printf("tree site lk: %lf\n",tree->site_lk_cat[0]);

	  log_site_lk = LOG(site_lk) - (phydbl)LOG2 * fact_sum_scale;
	  For(catg,tree->mod->n_catg) tree->log_site_lk_cat[catg][site] = LOG(tree->site_lk_cat[catg]) - (phydbl)LOG2 * fact_sum_scale;
	  if(isinf(log_site_lk) || isnan(log_site_lk))
	    {
	      PhyML_Printf("\n. site = %d",site);
	      PhyML_Printf("\n. invar = %f",tree->data->invar[site]);
	      PhyML_Printf("\n. scale_down = %d scale_upp = %d %lf",sum_scale_down[0],sum_scale_upp[0],fact_sum_scale);
	      PhyML_Printf("\n. Lk = %G LOG(Lk) = %f < %G",site_lk,log_site_lk,-BIG);
	      PhyML_Printf("\n. Err in file %s at line %d\n\n",__FILE__,__LINE__);
	      PhyML_Printf("\n. %d %d %d\n\n",b->num,d->num,anc->num);
	      Warn_And_Exit("\n");
	    }
	  tree->cur_site_lk[site] = log_site_lk;
	  /* Multiply log likelihood by the number of times this site pattern is found in the data */
	  tree->c_lnL_sorted[site] = tree->data->wght[site]*log_site_lk;
	  tree->c_lnL += tree->data->wght[site]*log_site_lk;
	  return log_site_lk;
}

/*********************************************************/
//Updated by Ken 3/8/2016
phydbl ASR_At_Given_Edge(t_edge *b_fcus, t_tree *tree, int root)
{
	if(b_fcus->left->tax){
	   PhyML_Printf("\n. Err in file %s at line %d\n\n",__FILE__,__LINE__);
	   Warn_And_Exit("");
	}

	if(tree->mod->whichrealmodel == HLP17){//Added by Ken 3/8/2016
	  Update_PMat_At_Given_Edge(b_fcus,tree);

	  if(tree->mod->s_opt->opt_topo){
		  if(b_fcus->anc_node->num != tree->mod->startnode){
			  Fill_UPP_single(tree,b_fcus);
		  }else{
			  Fill_UPP_root(tree,b_fcus);
		  }
	  }
	  tree->c_lnL             = .0;
	   tree->sum_min_sum_scale = .0;
	   For(tree->curr_site,tree->n_pattern){
	   	 if(!root)ASR_Core(b_fcus,tree, b_fcus->anc_node, b_fcus->des_node);
	   	 else ASR_Core_root(b_fcus,tree, b_fcus->anc_node, b_fcus->des_node);
		   //Lk_Core_UPP(b_fcus,tree, b_fcus->anc_node, b_fcus->des_node);
	   }
	  return tree->c_lnL;
  }else{

  int n_patterns;
  n_patterns = tree->n_pattern;

  #ifdef TIMES
  if((tree->rates) && (tree->rates->bl_from_rt)) RATES_Update_Cur_Bl(tree);
  if(tree->bl_from_node_stamps) TIMES_Bl_From_T(tree);
  #endif
  Update_PMat_At_Given_Edge(b_fcus,tree);

  tree->c_lnL             = .0;
  tree->sum_min_sum_scale = .0;
  For(tree->curr_site,n_patterns)
  {
    if(tree->data->wght[tree->curr_site] > SMALL) Lk_Core(b_fcus,tree);

  }
  //printf("%d tree LK3: %lf \n",b_fcus->num,tree->c_lnL);

  return tree->c_lnL;
  }
}



/*********************************************************/
//Updated by Ken 3/8/2016
phydbl Lk_At_Given_Edge(t_edge *b_fcus, t_tree *tree)
{
	if(b_fcus->left->tax){
	   PhyML_Printf("\n. Err in file %s at line %d\n\n",__FILE__,__LINE__);
	   Warn_And_Exit("");
	}

	if(tree->mod->whichrealmodel == HLP17){//Added by Ken 3/8/2016
	  Update_PMat_At_Given_Edge(b_fcus,tree);

	  if(tree->mod->s_opt->opt_topo){
		  if(b_fcus->anc_node->num != tree->mod->startnode){
			  Fill_UPP_single(tree,b_fcus);
		  }else{
			  Fill_UPP_root(tree,b_fcus);
		  }
	  }
	  tree->c_lnL             = .0;
	   tree->sum_min_sum_scale = .0;
	   For(tree->curr_site,tree->n_pattern){
	   	 Lk_Core_UPP(b_fcus,tree, b_fcus->anc_node, b_fcus->des_node);
	   }
	  return tree->c_lnL;
  }else{

  int n_patterns;
  n_patterns = tree->n_pattern;

  #ifdef TIMES
  if((tree->rates) && (tree->rates->bl_from_rt)) RATES_Update_Cur_Bl(tree);
  if(tree->bl_from_node_stamps) TIMES_Bl_From_T(tree);
  #endif
  Update_PMat_At_Given_Edge(b_fcus,tree);

  tree->c_lnL             = .0;
  tree->sum_min_sum_scale = .0;
  For(tree->curr_site,n_patterns)
  {
    if(tree->data->wght[tree->curr_site] > SMALL) Lk_Core(b_fcus,tree);

  }
  //printf("%d tree LK3: %lf \n",b_fcus->num,tree->c_lnL);

  return tree->c_lnL;
  }
}

/*********************************************************/
/* Core of the likelihood calculcation. Assume that the partial likelihoods on both
   sides of t_edge *b are up-to-date. Calculate the log-likelihood at one site.
*/


phydbl Lk_Core_UPP(t_edge *b, t_tree *tree, t_node *anc, t_node *d)
{
	phydbl log_site_lk;
	  phydbl site_lk_cat, site_lk;
	  int fact_sum_scale;
	  phydbl max_sum_scale;
	  phydbl sum;
	  int ambiguity_check,state;
	  int catg,ns,k,l,site;
	  int dim1,dim2,dim3;
	  int *sum_scale_down_cat,*sum_scale_upp_cat,*sum_scale_down,*sum_scale_upp;
	  phydbl multiplier;
	  int exponent, piecewise_exponent;
	  phydbl tmp;
	  phydbl logbig;
	  phydbl inv_site_lk;

	  logbig = LOG((phydbl)BIG);

	  dim1 = tree->mod->n_catg * tree->mod->ns;
	  dim2 = tree->mod->ns;
	  dim3 = tree->mod->ns * tree->mod->ns;

	  log_site_lk     = .0;
	  site_lk         = .0;
	  site_lk_cat     = .0;
	  ambiguity_check = -1;
	  state           = -1;
	  site            = tree->curr_site;
	  ns              = tree->mod->ns;

	  phydbl *p_lk;

	  if(d->num == b->left->num){
	      p_lk = b->p_lk_left;
	      sum_scale_down = b->sum_scale_left;
	      sum_scale_down_cat = b->sum_scale_left_cat;
	     if((d->tax) && (!tree->mod->s_opt->greedy)){
	    	 if(site==16){printf("left taxa %d %s\n",b->num,d->name);}
	    	 exit(EXIT_FAILURE);
	        ambiguity_check = tree->data->c_seq[d->num]->is_ambigu[site];
	        if(!ambiguity_check)state = Get_State_From_P_Pars(b->p_lk_tip_l,site*dim2,tree);
	      }
	    }
	  else{
	      p_lk = b->p_lk_rght;
	      sum_scale_down = b->sum_scale_rght;
	      sum_scale_down_cat = b->sum_scale_rght_cat;
	      if((d->tax) && (!tree->mod->s_opt->greedy)){
	        ambiguity_check = tree->data->c_seq[d->num]->is_ambigu[site];
	        if(!ambiguity_check)state = Get_State_From_P_Pars(b->p_lk_tip_r,site*dim2,tree);
	      }
	    }
	  sum_scale_upp_cat = b->sum_scale_upp_cat;
	  sum_scale_upp = b->sum_scale_upp;

	  if(tree->mod->use_m4mod) ambiguity_check = 1;


	  #if defined OMP || defined BLAS_OMP
	  #pragma omp parallel for if(tree->mod->n_w_catg>1) private(site_lk_cat,sum,k,l)
	  #endif


	  For(catg,tree->mod->n_catg){
	      site_lk_cat = .0;
	      /* b is an external edge */
	      if((d->tax) && (!tree->mod->s_opt->greedy)){
	          /* If the character observed at the tip is NOT ambiguous: ns x 1 terms to consider */
	          if(!ambiguity_check){
	              sum = .0;
	              For(l,ns){
	                phydbl val = b->bPmat_part[tree->mod->partIndex[site]][catg*dim3+l*dim2+state]*b->upp[site][l];
	                sum+=val;
	              }
	              site_lk_cat = sum;
	            }
	          /* If the character observed at the tip is ambiguous: ns x ns terms to consider */
	          else{
	              For(k,ns){
	                  sum = .0;
	                //  if(b->p_lk_tip_r[site*dim2+l] > .0){
	                      For(l,ns){
	                        phydbl val = b->bPmat_part[tree->mod->partIndex[site]][catg*dim3+k*dim2+l]*b->p_lk_tip_r[site*dim2+l]*b->upp[site][k];
	                        sum+=val;
	                      }
	                      site_lk_cat += sum;
	                    }
	            //    }
	            }
	        }
	      /* b is an internal edge: ns x ns terms to consider */
	      else{
	          For(k,ns){
	              sum = .0;
	         //     if(b->p_lk_rght[site*dim1+catg*dim2+k] > .0){
	                  For(l,ns){
	                    	phydbl val = b->bPmat_part[tree->mod->partIndex[site]][catg*dim3+k*dim2+l]*p_lk[site*dim1+catg*dim2+l]*b->upp[site][k];//Modified by Ken 17/8/2016
	                      sum+=val;
	                    }
	                  site_lk_cat+=sum;
	                }
	        //    }
	        }
	      tree->site_lk_cat[catg] = site_lk_cat;
	    }

	  max_sum_scale =  (phydbl)BIG;
	  For(catg,tree->mod->n_catg){

	      sum_scale_down_cat[catg] =
	        (sum_scale_down)?
	        (sum_scale_down[catg*tree->n_pattern+site]):
	        (0.0);

	      sum_scale_upp_cat[catg] =
	        (sum_scale_upp)?
	        (sum_scale_upp[catg*tree->n_pattern+site]):
	        (0.0);

	      sum = sum_scale_down_cat[catg] + sum_scale_upp_cat[catg];
	      if(sum < .0){
	    	  PhyML_Printf("\n. sum = %G",sum);
	    	  PhyML_Printf("\n. Err in file %s at line %d\n\n",__FILE__,__LINE__);
	    	  Warn_And_Exit("\n");
	      }
	      tmp = sum + (logbig - LOG(tree->site_lk_cat[catg]))/(phydbl)LOG2;
	      if(tmp < max_sum_scale) max_sum_scale = tmp; /* min of the maxs */
	    }

	  fact_sum_scale = (int)(max_sum_scale / 2);

	  // Apply scaling factors
	  For(catg,tree->mod->n_catg){
	      exponent = -(sum_scale_down_cat[catg]+sum_scale_upp_cat[catg])+fact_sum_scale;
	      site_lk_cat = tree->site_lk_cat[catg];
	      if(exponent > 0){
	    	  do{
	    		  piecewise_exponent = MIN(exponent,63);
	    		  multiplier = (phydbl)((unsigned long long)(1) << piecewise_exponent);
	    		  site_lk_cat *= multiplier;
	    		  exponent -= piecewise_exponent;
	    	  }
	    	  while(exponent != 0);
	      }
	      else{
	    	  do{
	    		  piecewise_exponent = MAX(exponent,-63);
	    		  multiplier = 1. / (phydbl)((unsigned long long)(1) << -piecewise_exponent);
	    		  site_lk_cat *= multiplier;
	    		  exponent -= piecewise_exponent;
	    	  }
	    	  while(exponent != 0);
	      }

	      if(isinf(site_lk_cat)){
	    	  PhyML_Printf("\n+ site=%4d site_lk_cat=%G sum_scale=%d max=%G fact=%d expo=%d dbl=%G",
			       tree->curr_site,
			       tree->site_lk_cat[catg],
			       sum_scale_down_cat[catg]+sum_scale_upp_cat[catg],
			       max_sum_scale,
			       fact_sum_scale,
			       -(sum_scale_down_cat[catg]+sum_scale_upp_cat[catg])+fact_sum_scale,
			       (double)tree->site_lk_cat[catg] * pow(2.,-(sum_scale_down_cat[catg]+sum_scale_upp_cat[catg])+fact_sum_scale));
	    	  site_lk_cat = BIG / 10;
	      }

	      if(site_lk_cat < SMALL){
	    	  site_lk_cat = SMALL;
	    	  printf("site lk = 0 %d %d %lf %d %d %d.\nSetting underflow to DBL_MIN. See Manual.\n",anc->num,d->num,site_lk_cat,site,sum_scale_upp_cat[catg],sum_scale_down_cat[catg]);
	    	  printf("%lf %lf",tree->mod->omega_part[0],tree->mod->kappa);
	       	  printf("\n");
	      }
	      tree->site_lk_cat[catg] = site_lk_cat;
	    }


	  site_lk = .0;
	  if(tree->io->datatype!=CODON ) For(catg,tree->mod->n_catg) site_lk += tree->site_lk_cat[catg]*tree->mod->gamma_r_proba[catg];
	  else //!< Added by Marcelo.
	  {
	    if(tree->mod->n_catg>1 && tree->mod->n_w_catg==1) For(catg,tree->mod->n_catg) site_lk += tree->site_lk_cat[catg]*tree->mod->gamma_r_proba[catg];
	    else For(catg,tree->mod->n_catg) site_lk += tree->site_lk_cat[catg]*tree->mod->prob_omegas[catg];
	  }

	  /* The substitution model does include invariable sites */
	   if(tree->mod->invar)
	     {
	       /* The site is invariant */
	       if(tree->data->invar[site] > -0.5)
	         {
	 	  /* Multiply P(D|r=0) by 2^(fact_sum_scale) */

	 	  inv_site_lk = tree->mod->pi[tree->data->invar[site]];
	 	  exponent = fact_sum_scale;
	 	  do
	 	    {
	 	      piecewise_exponent = MIN(exponent,63);
	 	      multiplier = (phydbl)((unsigned long long)(1) << piecewise_exponent);
	 	      inv_site_lk *= multiplier;
	 	      exponent -= piecewise_exponent;
	 	    }
	 	  while(exponent != 0);
	 	  inv_site_lk *= tree->mod->pinvar;
	 	  /* Update the value of site_lk */
	 	  site_lk *= (1. - tree->mod->pinvar);
	 	  site_lk += inv_site_lk;
	         }
	       else
	         {
	           /* Same formula as above with P(D | subs, rate = 0) = 0 */
	 	  site_lk *= (1. - tree->mod->pinvar);
	         }
	     }

	  log_site_lk = LOG(site_lk) - (phydbl)LOG2 * fact_sum_scale;
	  For(catg,tree->mod->n_catg) tree->log_site_lk_cat[catg][site] = LOG(tree->site_lk_cat[catg]) - (phydbl)LOG2 * fact_sum_scale;
	  if(isinf(log_site_lk) || isnan(log_site_lk))
	    {
	      PhyML_Printf("\n. site = %d",site);
	      PhyML_Printf("\n. invar = %f",tree->data->invar[site]);
	      PhyML_Printf("\n. scale_down = %d scale_upp = %d %lf",sum_scale_down[0],sum_scale_upp[0],fact_sum_scale);
	      PhyML_Printf("\n. Lk = %G LOG(Lk) = %f < %G",site_lk,log_site_lk,-BIG);
	      PhyML_Printf("\n. Err in file %s at line %d\n\n",__FILE__,__LINE__);
	      PhyML_Printf("\n. %d %d %d\n\n",b->num,d->num,anc->num);
	      Warn_And_Exit("\n");
	    }
	  tree->cur_site_lk[site] = log_site_lk;
	  /* Multiply log likelihood by the number of times this site pattern is found in the data */
	  tree->c_lnL_sorted[site] = tree->data->wght[site]*log_site_lk;
	  tree->c_lnL += tree->data->wght[site]*log_site_lk;
	  return log_site_lk;
}

#ifndef USE_OLD_LK
phydbl Lk_Core(t_edge *b, t_tree *tree)
{
  phydbl log_site_lk;
  phydbl site_lk_cat, site_lk;
  int sum_scale_left, sum_scale_rght;
  int fact_sum_scale;
  phydbl max_sum_scale;
  phydbl sum;
  int ambiguity_check,state;
  int catg,ns,k,l,site;
  int dim1,dim2,dim3;
  int *sum_scale_left_cat,*sum_scale_rght_cat;
  phydbl multiplier;
  int exponent, piecewise_exponent;
  phydbl tmp;
  phydbl logbig;
  phydbl inv_site_lk;

  logbig = LOG((phydbl)BIG);


  dim1 = tree->mod->n_catg * tree->mod->ns;
  dim2 = tree->mod->ns;
  dim3 = tree->mod->ns * tree->mod->ns;

  log_site_lk     = .0;
  site_lk         = .0;
  site_lk_cat     = .0;
  ambiguity_check = -1;
  state           = -1;
  site            = tree->curr_site;
  ns              = tree->mod->ns;



  if((b->rght->tax) && (!tree->mod->s_opt->greedy))
    {
      ambiguity_check = tree->data->c_seq[b->rght->num]->is_ambigu[site];
      if(!ambiguity_check) state = Get_State_From_P_Pars(b->p_lk_tip_r,site*dim2,tree);
    }

  sum_scale_left = .0;
  sum_scale_rght = .0;
  
  if(tree->mod->use_m4mod) ambiguity_check = 1;
  
  sum_scale_left_cat = b->sum_scale_left_cat;
  sum_scale_rght_cat = b->sum_scale_rght_cat;

  /* Actual likelihood calculation */
  /* For all classes of rates */
 
  #if defined OMP || defined BLAS_OMP
  
  #pragma omp parallel for if(tree->mod->n_w_catg>1) private(site_lk_cat,sum,k,l) 
     
  #endif

 
  For(catg,tree->mod->n_catg)
     {

      site_lk_cat = .0;

      /* b is an external edge */
      if((b->rght->tax) && (!tree->mod->s_opt->greedy))
        {
          /* If the character observed at the tip is NOT ambiguous: ns x 1 terms to consider */
          if(!ambiguity_check)
            {
              sum = .0;
              For(l,ns)
                {
                  sum +=
                    b->bPmat_part[tree->mod->partIndex[site]][catg*dim3+state*dim2+l] * //Modified by Ken 17/8/2016
                    b->p_lk_left[site*dim1+catg*dim2+l];
                }
              site_lk_cat += sum * tree->mod->pi[state];
            }
          /* If the character observed at the tip is ambiguous: ns x ns terms to consider */
          else
            {
              For(k,ns)
                {
                  sum = .0;
                  if(b->p_lk_tip_r[site*dim2+k] > .0)
                    {
                      For(l,ns)
                        {
                          sum +=
                        	b->bPmat_part[tree->mod->partIndex[site]][catg*dim3+k*dim2+l] *			//Modified by Ken 17/8/2016
                            b->p_lk_left[site*dim1+catg*dim2+l];
                        }
                      site_lk_cat +=
                        sum *
                        tree->mod->pi[k] *
                        b->p_lk_tip_r[site*dim2+k];
                    }
                }
            }
        }
      /* b is an internal edge: ns x ns terms to consider */
      else
        {
          For(k,ns)
            {
              sum = .0;
              if(b->p_lk_rght[site*dim1+catg*dim2+k] > .0)
                {
                  For(l,ns)
                    {
                      sum +=
                    	b->bPmat_part[tree->mod->partIndex[site]][catg*dim3+k*dim2+l] * //Modified by Ken 17/8/2016
                        b->p_lk_left[site*dim1+catg*dim2+l];
                      //printf("lk internal!\n");
                    }
                  site_lk_cat +=
                    sum *
                    tree->mod->pi[k] *
                    b->p_lk_rght[site*dim1+catg*dim2+k];
                }
            }
        }
      tree->site_lk_cat[catg] = site_lk_cat;
    }
    
  max_sum_scale =  (phydbl)BIG;
  For(catg,tree->mod->n_catg)
    {
      sum_scale_left_cat[catg] =
        (b->sum_scale_left)?
        (b->sum_scale_left[catg*tree->n_pattern+site]):
        (0.0);

      sum_scale_rght_cat[catg] =
        (b->sum_scale_rght)?
        (b->sum_scale_rght[catg*tree->n_pattern+site]):
        (0.0);

      sum = sum_scale_left_cat[catg] + sum_scale_rght_cat[catg];

      if(sum < .0)
	{
	  PhyML_Printf("\n. sum = %G",sum);
	  PhyML_Printf("\n. Err in file %s at line %d\n\n",__FILE__,__LINE__);
	  Warn_And_Exit("\n");
	}

      tmp = sum + (logbig - LOG(tree->site_lk_cat[catg]))/(phydbl)LOG2;
      if(tmp < max_sum_scale) max_sum_scale = tmp; /* min of the maxs */
    }

/*   fact_sum_scale = (int)((max_sum_scale + min_sum_scale) / 2); */


  fact_sum_scale = (int)(max_sum_scale / 2);


  /* Apply scaling factors */
  For(catg,tree->mod->n_catg)
    {
      exponent = -(sum_scale_left_cat[catg]+sum_scale_rght_cat[catg])+fact_sum_scale;

      site_lk_cat = tree->site_lk_cat[catg];

      if(exponent > 0)
	{
	  do
	    {
	      piecewise_exponent = MIN(exponent,63);
	      multiplier = (phydbl)((unsigned long long)(1) << piecewise_exponent);
	      site_lk_cat *= multiplier;
	      exponent -= piecewise_exponent;
	    }
	  while(exponent != 0);
	}
      else
	{
	  do
	    {
	      piecewise_exponent = MAX(exponent,-63);
	      multiplier = 1. / (phydbl)((unsigned long long)(1) << -piecewise_exponent);
	      site_lk_cat *= multiplier;
	      exponent -= piecewise_exponent;
	    }
	  while(exponent != 0);
	}

      if(isinf(site_lk_cat))
	{
	  PhyML_Printf("\n+ site=%4d site_lk_cat=%G sum_scale=%d max=%G fact=%d expo=%d dbl=%G",
		       tree->curr_site,
		       tree->site_lk_cat[catg],
		       sum_scale_left_cat[catg]+sum_scale_rght_cat[catg],
		       max_sum_scale,
		       fact_sum_scale,
		       -(sum_scale_left_cat[catg]+sum_scale_rght_cat[catg])+fact_sum_scale,
		       (double)tree->site_lk_cat[catg] * pow(2.,-(sum_scale_left_cat[catg]+sum_scale_rght_cat[catg])+fact_sum_scale));
	  site_lk_cat = BIG / 10;
	}

      //if(site_lk_cat < SMALL) site_lk_cat = 0;
      if(site_lk_cat < SMALL){
    	  site_lk_cat = SMALL;
    	  printf("site lk = 0 %lf %d %d %d.\nSetting underflow to DBL_MIN. Probably not a big problem during initial parameter searching..\n",site_lk_cat,site,sum_scale_left_cat[catg],sum_scale_rght_cat[catg]);
    	  printf("%lf %lf",tree->mod->omega_part[0],tree->mod->kappa);
       	  printf("\n");
      }
      tree->site_lk_cat[catg] = site_lk_cat;
    }


  site_lk = .0;
  if(tree->mod->datatype!=CODON ) For(catg,tree->mod->n_catg) site_lk += tree->site_lk_cat[catg]*tree->mod->gamma_r_proba[catg];
  else //!< Added by Marcelo.
  {
    if(tree->mod->n_catg>1 && tree->mod->n_w_catg==1) For(catg,tree->mod->n_catg) site_lk += tree->site_lk_cat[catg]*tree->mod->gamma_r_proba[catg];
    else For(catg,tree->mod->n_catg) site_lk += tree->site_lk_cat[catg]*tree->mod->prob_omegas[catg];
  }

 /* The substitution model does include invariable sites */
  if(tree->mod->invar)
    {
      /* The site is invariant */
      if(tree->data->invar[site] > -0.5)
        {
	  /* Multiply P(D|r=0) by 2^(fact_sum_scale) */
	  inv_site_lk = tree->mod->pi[tree->data->invar[site]];
	  exponent = fact_sum_scale;
	  do
	    {
	      piecewise_exponent = MIN(exponent,63);
	      multiplier = (phydbl)((unsigned long long)(1) << piecewise_exponent);
	      inv_site_lk *= multiplier;
	      exponent -= piecewise_exponent;
	    }
	  while(exponent != 0);
	  inv_site_lk *= tree->mod->pinvar;
	  /* Update the value of site_lk */
	  site_lk *= (1. - tree->mod->pinvar);
	  site_lk += inv_site_lk;
        }
      else
        {
          /* Same formula as above with P(D | subs, rate = 0) = 0 */
	  site_lk *= (1. - tree->mod->pinvar);
        }
    }

  log_site_lk = LOG(site_lk) - (phydbl)LOG2 * fact_sum_scale;

  For(catg,tree->mod->n_catg) tree->log_site_lk_cat[catg][site] = LOG(tree->site_lk_cat[catg]) - (phydbl)LOG2 * fact_sum_scale;
  //HERE!!!
  if(isinf(log_site_lk) || isnan(log_site_lk))
    {
      PhyML_Printf("\n. site = %d",site);
      PhyML_Printf("\n. invar = %f",tree->data->invar[site]);
      PhyML_Printf("\n. scale_left = %d scale_rght = %d",sum_scale_left,sum_scale_rght);
      PhyML_Printf("\n. Lk = %G LOG(Lk) = %f < %G",site_lk,log_site_lk,-BIG);
      PhyML_Printf("\n. Err in file %s at line %d\n\n",__FILE__,__LINE__);
      Warn_And_Exit("\n");
    }

  tree->cur_site_lk[site] = log_site_lk;

  /* Multiply log likelihood by the number of times this site pattern is found in the data */
  tree->c_lnL_sorted[site] = tree->data->wght[site]*log_site_lk;

  tree->c_lnL += tree->data->wght[site]*log_site_lk;

  return log_site_lk;
}


/*********************************************************/

/* Update partial likelihood on edge b on the side of b where
   node d lies.
*/
void Update_P_Lk(t_tree *tree, t_edge *b, t_node *d)
{
/*
       |
	   |<- b
	   |
	   d
      / \
   	 /   \
  	/     \
n_v1   n_v2
*/
  t_node *n_v1, *n_v2;
  phydbl p1_lk1,p2_lk2;
  phydbl *p_lk,*p_lk_v1,*p_lk_v2;
 // phydbl *Pij1,*Pij2;
  int *sum_scale, *sum_scale_v1, *sum_scale_v2;
  int sum_scale_v1_val, sum_scale_v2_val;
  int i,j;
  int catg,site;
  int dir1,dir2;
  int n_patterns;
  short int ambiguity_check_v1,ambiguity_check_v2;
  int state_v1,state_v2;
  int dim1, dim2, dim3;
  phydbl curr_scaler;
  int curr_scaler_pow, piecewise_scaler_pow;
  phydbl p_lk_lim_inf;
  phydbl smallest_p_lk;

  p_lk_lim_inf = (phydbl)P_LK_LIM_INF;
 
  dim1 = tree->mod->n_catg * tree->mod->ns;
  dim2 = tree->mod->ns;
  dim3 = tree->mod->ns * tree->mod->ns;

  state_v1 = state_v2 = -1;
  ambiguity_check_v1 = ambiguity_check_v2 = NO;
  sum_scale_v1_val = sum_scale_v2_val = 0;
  p1_lk1 = p2_lk2 = .0;

  if(d->tax)
    {
      PhyML_Printf("\n. t_node %d is a leaf...",d->num);
      PhyML_Printf("\n. Err in file %s at line %d\n\n",__FILE__,__LINE__);
      Warn_And_Exit("\n");
    }

  n_patterns = tree->n_pattern;
  
  /* TO DO: Might be worth keeping these directions in memory instead of
     calculating them every time... */
  dir1=dir2=-1;
  For(i,3) if(d->b[i] != b) (dir1<0)?(dir1=i):(dir2=i);
  
  n_v1 = d->v[dir1];
  n_v2 = d->v[dir2];

  /* Get the partial likelihood vectors on edge b and the two pendant
     edges (i.e., the two other edges connected to d) */
  if(d == b->left)
    {
      p_lk = b->p_lk_left;
      sum_scale = b->sum_scale_left;
    }
  else
    {
      p_lk = b->p_lk_rght;
      sum_scale = b->sum_scale_rght;
    }
      
  if(d == d->b[dir1]->left)
    {
      p_lk_v1 = d->b[dir1]->p_lk_rght;
      sum_scale_v1 = d->b[dir1]->sum_scale_rght;
    }
  else
    {
      p_lk_v1 = d->b[dir1]->p_lk_left;
      sum_scale_v1 = d->b[dir1]->sum_scale_left;
    }
  
  if(d == d->b[dir2]->left)
    {
      p_lk_v2 = d->b[dir2]->p_lk_rght;
      sum_scale_v2 = d->b[dir2]->sum_scale_rght;
    }
  else
    {
      p_lk_v2 = d->b[dir2]->p_lk_left;
      sum_scale_v2 = d->b[dir2]->sum_scale_left;
    }
  

  phydbl **Ppart1 = d->b[dir1]->bPmat_part;
  phydbl **Ppart2 = d->b[dir2]->bPmat_part;
  
  /* For every site in the alignment */
  
  #if defined OMP || defined BLAS_OMP
  #pragma omp parallel for if(tree->mod->datatype==CODON  && !tree->mod->splitByTree) private(state_v1, state_v2, ambiguity_check_v1, ambiguity_check_v2, smallest_p_lk, p1_lk1, catg, i, j, p2_lk2, sum_scale_v1_val, sum_scale_v2_val, curr_scaler_pow, curr_scaler, piecewise_scaler_pow )
  #endif
  
  For(site,n_patterns)
    {
	  //printf("%d\n",omp_get_thread_num());
      state_v1 = state_v2 = -1;
      ambiguity_check_v1 = ambiguity_check_v2 = NO;
      if(!tree->mod->s_opt->greedy)
	{
	  /* n_v1 and n_v2 are tip nodes */
	  if(n_v1->tax)
	    {
	      /* Is the state at this tip ambiguous? */
	      ambiguity_check_v1 = tree->data->c_seq[n_v1->num]->is_ambigu[site];
	      if(ambiguity_check_v1 == NO) state_v1 = Get_State_From_P_Pars(n_v1->b[0]->p_lk_tip_r,site*dim2,tree);
	    }
	      
	  if(n_v2->tax)
	    {
	      /* Is the state at this tip ambiguous? */
	      ambiguity_check_v2 = tree->data->c_seq[n_v2->num]->is_ambigu[site];
	      if(ambiguity_check_v2 == NO) state_v2 = Get_State_From_P_Pars(n_v2->b[0]->p_lk_tip_r,site*dim2,tree);

	    }
	}
      
      if(tree->mod->use_m4mod)
	{
	  ambiguity_check_v1 = YES;
	  ambiguity_check_v2 = YES;
	}

      /* For all the rate classes */
      For(catg,tree->mod->n_catg)
	{
	  smallest_p_lk  =  BIG;

	  /* For all the state at node d */
	  For(i,tree->mod->ns)
	    {
	      p1_lk1 = .0;
	      
	      /* n_v1 is a tip */
	      if((n_v1->tax) && (!tree->mod->s_opt->greedy))
		{
		  if(ambiguity_check_v1 == NO)
		    {
		      /* For the (non-ambiguous) state at node n_v1 */
		      p1_lk1 = d->b[dir1]->bPmat_part[tree->mod->partIndex[site]][catg*dim3+i*dim2+state_v1]; // Ken 17/8/2016
		    }
		  else
		    {
		      /* For all the states at node n_v1 */
		      For(j,tree->mod->ns)
			{
			  p1_lk1 +=d->b[dir1]->bPmat_part[tree->mod->partIndex[site]][catg*dim3+i*dim2+j] * (phydbl)n_v1->b[0]->p_lk_tip_r[site*dim2+j]; //!!!!!
			  phydbl temp = (phydbl)n_v1->b[0]->p_lk_tip_r[site*dim2+j];
			}
		    }
		}
	      /* n_v1 is an internal node */
	      else
		{
		  /* For the states at node n_v1 */
		  For(j,tree->mod->ns)
		    {
		      p1_lk1 += d->b[dir1]->bPmat_part[tree->mod->partIndex[site]][catg*dim3+i*dim2+j] * p_lk_v1[site*dim1+catg*dim2+j];//!!!!!!
		    }
		}
	      
	      p2_lk2 = .0;
	      
	      /* We do exactly the same as for node n_v1 but for node n_v2 this time.*/
	      /* n_v2 is a tip */
	      if((n_v2->tax) && (!tree->mod->s_opt->greedy))
		{
		  if(ambiguity_check_v2 == NO)
		    {
		      /* For the (non-ambiguous) state at node n_v2 */
		      p2_lk2 =  d->b[dir2]->bPmat_part[tree->mod->partIndex[site]][catg*dim3+i*dim2+state_v2];

		    }
		  else
		    {
		      /* For all the states at node n_v2 */
		      For(j,tree->mod->ns)
			{
			  p2_lk2 += d->b[dir2]->bPmat_part[tree->mod->partIndex[site]][catg*dim3+i*dim2+j] * (phydbl)n_v2->b[0]->p_lk_tip_r[site*dim2+j]; //!!!!!
			  phydbl temp = (phydbl)n_v2->b[0]->p_lk_tip_r[site*dim2+j];
			}
		    }
		}
	      /* n_v2 is an internal node */
	      else
		{
		  /* For all the states at node n_v2 */
		  For(j,tree->mod->ns)
		    {
		      p2_lk2 += d->b[dir2]->bPmat_part[tree->mod->partIndex[site]][catg*dim3+i*dim2+j] * p_lk_v2[site*dim1+catg*dim2+j];//!!!!!
		    }
		}
	      
	      p_lk[site*dim1+catg*dim2+i] = p1_lk1 * p2_lk2;	    

	      if(p_lk[site*dim1+catg*dim2+i] < smallest_p_lk) smallest_p_lk = p_lk[site*dim1+catg*dim2+i] ; 
	    }
	      
	  /* Current scaling values at that site */
 	  sum_scale_v1_val = (sum_scale_v1)?(sum_scale_v1[catg*n_patterns+site]):(0);
	  sum_scale_v2_val = (sum_scale_v2)?(sum_scale_v2[catg*n_patterns+site]):(0);
	  sum_scale[catg*n_patterns+site] = sum_scale_v1_val + sum_scale_v2_val;
	  
	  /* Scaling */
	  if(smallest_p_lk < p_lk_lim_inf)
	    {
	      curr_scaler_pow = (int)(LOG(p_lk_lim_inf)-LOG(smallest_p_lk))/LOG2;
	      curr_scaler     = (phydbl)((unsigned long long)(1) << curr_scaler_pow);

	      sum_scale[catg*n_patterns+site] += curr_scaler_pow;

	      do
		{
		  piecewise_scaler_pow = MIN(curr_scaler_pow,63);
		  curr_scaler = (phydbl)((unsigned long long)(1) << piecewise_scaler_pow);
		  For(i,tree->mod->ns)
		    {
		      p_lk[site*dim1+catg*dim2+i] *= curr_scaler;
		      
		      if(p_lk[site*dim1+catg*dim2+i] > BIG)
			{
			  PhyML_Printf("\n. p_lk_lim_inf = %G smallest_p_lk = %G",p_lk_lim_inf,smallest_p_lk);
			  PhyML_Printf("\n. curr_scaler_pow = %d",curr_scaler_pow);
			  PhyML_Printf("\n. Err in file %s at line %d.",__FILE__,__LINE__);
			  Warn_And_Exit("\n");
			}
		    }
		  curr_scaler_pow -= piecewise_scaler_pow;
		}
	      while(curr_scaler_pow != 0);
	    }
	}
    }
  if(tree->mod->s_opt->opt_topo && tree->mod->whichrealmodel == HLP17){
	  if(b->anc_node->num != tree->mod->startnode){
		  Fill_UPP_single(tree,b);
	  }else{
		  Fill_UPP_root(tree,b);
	  }
  }
}
#endif

/*********************************************************/

matrix *ML_Dist(calign *data, model *mod)
{
  int i,j,k,l;
  phydbl init;
  int n_catg;
  phydbl d_max,sum;
  matrix *mat;
  calign *twodata,*tmpdata;
  int state0, state1,len;
  phydbl *F;
  eigen *eigen_struct;
  time_t tbegin, tend; //! Added by Marcelo
  
  tmpdata         = (calign *)mCalloc(1,sizeof(calign));
  tmpdata->c_seq  = (align **)mCalloc(2,sizeof(align *));
  tmpdata->b_frq  = (phydbl *)mCalloc(mod->ns,sizeof(phydbl));
  tmpdata->ambigu = (short int *)mCalloc(data->crunch_len,sizeof(short int));
  F               = (phydbl *)mCalloc(mod->ns*mod->ns,sizeof(phydbl ));
  eigen_struct    = (eigen *)Make_Eigen_Struct(mod);
  
  tmpdata->n_otu  = 2;
  
  tmpdata->crunch_len = data->crunch_len;
  tmpdata->init_len   = data->init_len;
  
  
  mat = NULL;

  For(i,mod->n_catg) /* Don't use the discrete gamma distribution */
  {
    mod->gamma_rr[i]      = 1.0;
    mod->gamma_r_proba[i] = 1.0;
  }
  
  n_catg = mod->n_catg;
  mod->n_catg = 1;
  
  time(&tbegin); //! Added by Marcelo.
  
  int progress = 0; //!Added by Marcelo.
  double loops = ((data->n_otu)*(data->n_otu)-(data->n_otu))/2.00;//!Added by Marcelo.
  For(j,data->n_otu-1)
  {
    tmpdata->c_seq[0]       = data->c_seq[j];
    tmpdata->c_seq[0]->name = data->c_seq[j]->name;
    tmpdata->wght           = data->wght;
    
    for(k=j+1;k<data->n_otu;k++)
    {
      
      tmpdata->c_seq[1]       = data->c_seq[k];
      tmpdata->c_seq[1]->name = data->c_seq[k]->name;
      twodata = Compact_Cdata(tmpdata,mod->io,mod);
      For(l,mod->ns) twodata->b_frq[l] = data->b_frq[l];
      Check_Ambiguities(twodata,mod->datatype,mod->state_len);
      
      Hide_Ambiguities(twodata);
      
      init = mat->dist[j][k];
      
      if((init > DIST_MAX-SMALL) || (init < .0)) init = 0.1;
      
      d_max = init;
      
      For(i,mod->ns*mod->ns) F[i]=.0;
      len = 0;
      For(l,twodata->c_seq[0]->len)
      {
	state0 = Assign_State(twodata->c_seq[0]->state+l*mod->state_len,mod->datatype,mod->state_len); //was originall mod-> io-> mod->state_len Ken 9/1/2018
	state1 = Assign_State(twodata->c_seq[1]->state+l*mod->state_len,mod->datatype,mod->state_len);
	
	if((state0 > -1) && (state1 > -1))
	{
	  if(mod->datatype == CODON) //! Added by Marcelo.
	  { 
	    state0=indexSenseCodons[state0];
	    state1=indexSenseCodons[state1];
	  }
	  F[mod->ns*state0+state1] += (int)twodata->wght[l];
	  len += (int)twodata->wght[l];
	}
      }
      sum=0;
      For(i,mod->ns*mod->ns) sum += F[i];
           
      if(len > .0) {For(i,mod->ns*mod->ns) F[i] /= (phydbl)len;}
      
      sum = 0.;
      For(i,mod->ns*mod->ns) sum += F[i];
      
      
      if(sum < .001) d_max = -1.;
      else if((sum > 1. - .001) && (sum < 1. + .001)) Opt_Dist_F(&(d_max),F,mod);
      else
      {
	PhyML_Printf("\n. sum = %f\n",sum);
	PhyML_Printf("\n. Err in file %s at line %d\n\n",__FILE__,__LINE__);
	Exit("");
      }
      
      if(d_max >= DIST_MAX) d_max = DIST_MAX;
      
      
      /* Do not correct for dist < BL_MIN, otherwise Fill_Missing_Dist
      *  will not be called
      */
      
      mat->dist[j][k] = d_max;
      mat->dist[k][j] = mat->dist[j][k];
      Free_Cseq(twodata);
      progress++;
      if(!mod->quiet) PhyML_Printf("\r. Computing pairwise distances...%3.2f%c concluded.", progress*100.00/loops, '%');
    }
  }
  if(!mod->quiet) PhyML_Printf("\n");
  time(&tend); //! Added by Marcelo.
  
  if((tend-tbegin)>=1) //! Added by Marcelo.
  {
    PhyML_Printf("\n. Calculation of pairwise distances completed in %ld sec.\n",(long int)(tend-tbegin));
  }
  else
  {
    PhyML_Printf("\n. Calculation of pairwise distances completed in less than 1 sec.\n");
  }
  
  mod->n_catg = n_catg;
  
  free(tmpdata->ambigu);
  free(tmpdata->b_frq);
  free(tmpdata->c_seq);
  free(tmpdata);
  Free_Eigen(eigen_struct);
  free(F);
  
  return mat;
}
/*********************************************************/

void Unconstraint_Lk(t_tree *tree)
{
  int i;

  tree->unconstraint_lk = .0;

  For(i,tree->data->crunch_len)
    {
      tree->unconstraint_lk +=
	tree->data->wght[i]*(phydbl)LOG(tree->data->wght[i]);
    }
  tree->unconstraint_lk -=
    tree->data->init_len*(phydbl)LOG(tree->data->init_len);
}

/*********************************************************/

void Make_Tree_4_Lk(t_tree *tree, calign *cdata, int n_site)
{
  int i;
  
  tree->c_lnL_sorted = (phydbl *)mCalloc(tree->n_pattern, sizeof(phydbl));
  tree->cur_site_lk  = (phydbl *)mCalloc(cdata->crunch_len,sizeof(phydbl));
  tree->old_site_lk  = (phydbl *)mCalloc(cdata->crunch_len,sizeof(phydbl));
    
  tree->site_lk_cat          = (phydbl *)mCalloc(tree->mod->n_catg,sizeof(phydbl));  
  tree->log_site_lk_cat      = (phydbl **)mCalloc(tree->mod->n_catg,sizeof(phydbl *));
  For(i,tree->mod->n_catg) tree->log_site_lk_cat[i] = (phydbl *)mCalloc(cdata->crunch_len,sizeof(phydbl));
      
  //printf("%d %d catg cdata crunch_len %d\n",tree->n_pattern,tree->mod->n_catg,cdata->crunch_len);

  tree->log_lks_aLRT = (phydbl **)mCalloc(3,sizeof(phydbl *));
  For(i,3) tree->log_lks_aLRT[i] = (phydbl *)mCalloc(tree->data->init_len,sizeof(phydbl));

  For(i,2*tree->n_otu-3)
  {
    Make_Edge_Lk(tree->t_edges[i],tree);
    Make_Edge_NNI(tree->t_edges[i]);
  }
    
  //Does nothing? Ken 4/1/2017
  For(i,2*tree->n_otu-2) Make_Node_Lk(tree->noeud[i]);

  if(tree->mod->s_opt->greedy) 
  {
    Init_P_Lk_Tips_Double(tree);
  }
  else 
  {
    Init_P_Lk_Tips_Int(tree);
  }
}

/*********************************************************/

void Init_P_Lk_Tips_Double(t_tree *tree)
{
  int curr_site,i,j,k,dim1,dim2;
  
  dim1 = tree->mod->n_catg * tree->mod->ns;
  dim2 = tree->mod->ns;

  For(curr_site,tree->data->crunch_len)
    {
      For(i,tree->n_otu)
	{
	  if (tree->mod->datatype == NT)
	    Init_Tips_At_One_Site_Nucleotides_Float(tree->data->c_seq[i]->state[curr_site],
						    curr_site*dim1+0*dim2,
						    tree->noeud[i]->b[0]->p_lk_rght);

	  else if(tree->mod->datatype == AA)
	    Init_Tips_At_One_Site_AA_Float(tree->data->c_seq[i]->state[curr_site],
					   curr_site*dim1+0*dim2,
					   tree->noeud[i]->b[0]->p_lk_rght);

	  else if(tree->mod->datatype == GENERIC)
	    Init_Tips_At_One_Site_Generic_Float(tree->data->c_seq[i]->state+curr_site*tree->mod->state_len,
						tree->mod->ns,
						tree->mod->state_len,
						curr_site*dim1+0*dim2,
						tree->noeud[i]->b[0]->p_lk_rght);
	  else if (tree->mod->datatype == CODON)//!< Added by Marcelo.
	    Init_Tips_At_One_Site_Codons_Float(tree->data->c_seq[i]->state[curr_site],
						    curr_site*dim1+0*dim2,//!< Changed by Marcelo.
						    tree->noeud[i]->b[0]->p_lk_rght,
					             tree->data->c_seq[i]->alternativeCodons[curr_site]);

	 for(j=1;j<tree->mod->n_catg;j++)
	    {
	      For(k,tree->mod->ns)
		{
		  tree->noeud[i]->b[0]->p_lk_rght[curr_site*dim1+j*dim2+k] = 
		    tree->noeud[i]->b[0]->p_lk_rght[curr_site*dim1+0*dim2+k];
		}
	    }
	}
    }

  #ifdef M4
  if(tree->mod->m4mod) M4_Init_P_Lk_Tips_Double(tree);
  #endif
}

/*********************************************************/

void Init_P_Lk_Tips_Int(t_tree *tree)
{
  int curr_site,i,dim1;

  dim1 = tree->mod->ns;

  For(curr_site,tree->data->crunch_len)
    {
      For(i,tree->n_otu)
	{
	  if(tree->mod->datatype == NT)
	    Init_Tips_At_One_Site_Nucleotides_Int(tree->data->c_seq[i]->state[curr_site],
						  curr_site*dim1,
						  tree->noeud[i]->b[0]->p_lk_tip_r);

	  else if(tree->mod->datatype == AA)
	    Init_Tips_At_One_Site_AA_Int(tree->data->c_seq[i]->state[curr_site],
					 curr_site*dim1,					   
					 tree->noeud[i]->b[0]->p_lk_tip_r);

	  else if(tree->mod->datatype == GENERIC)
	    {
	      Init_Tips_At_One_Site_Generic_Int(tree->data->c_seq[i]->state+curr_site*tree->mod->state_len,
						tree->mod->ns,
						tree->mod->state_len,
						curr_site*dim1,
						tree->noeud[i]->b[0]->p_lk_tip_r);
	    }
	    else if(tree->mod->datatype == CODON)//!< Added by Marcelo.
	    {
	      Init_Tips_At_One_Site_Codons_Int(tree->data->c_seq[i]->state[curr_site],
						  curr_site*dim1,
						  tree->noeud[i]->b[0]->p_lk_tip_r,
					          tree->data->c_seq[i]->alternativeCodons[curr_site]);
	    }
	}
    }
  #ifdef M4
  if(tree->mod->m4mod) M4_Init_P_Lk_Tips_Int(tree);
  #endif
}

/*********************************************************/

void Update_PMat_At_Given_Edge(t_edge *b_fcus, t_tree *tree)
{
  phydbl len;  
  int i,k;
  
  if(b_fcus->l < BL_MIN)      b_fcus->l = BL_MIN;
  else if(b_fcus->l > BL_MAX) b_fcus->l = BL_MAX;
   
  if(tree->mod->datatype==CODON) //!< Added by Marcelo.
  {
    if(b_fcus->has_zero_br_len) 
    {
    	printf("zero branch length\nKen should have fixed this by now.\nPlease complain to him.\n");
    	exit(EXIT_FAILURE);
      For(k,tree->mod->n_w_catg) For(i,tree->mod->ns*tree->mod->ns) b_fcus->Pij_rr[i+k*tree->mod->ns*tree->mod->ns]=tree->mod->A0_part[0][i];//Ken 22/8 Need to fix
    } 
    else
    {
      if(tree->mod->n_catg>1 && tree->mod->n_w_catg==1)
      {
	For(i,tree->mod->n_catg)
	{
	  len = b_fcus->l*tree->mod->gamma_rr[i];	  
	  if(len < BL_MIN)      len = BL_MIN;
	  else if(len > BL_MAX) len = BL_MAX;
	  PMat_CODON(len,tree->mod,0,b_fcus->Pij_rr+tree->mod->ns*tree->mod->ns*i);
	}
      }
      else
      {
    	  len = b_fcus->l;
    	  if(tree->mod->testcondition){len=BL_MIN;}
    	  int modeli;
    	  for(modeli=0;modeli<tree->mod->nparts;modeli++){
    		  PMat_CODON_part(len,tree->mod,0,tree->mod->qmat_part[modeli],b_fcus->bPmat_part[modeli],modeli);
    	  }
      }
    }
  }
  else
  {
    
    len = -1.0;
    
    For(i,tree->mod->n_catg)
    {
      if(b_fcus->has_zero_br_len) len = -1.0;
      else
      {
	len = b_fcus->l*tree->mod->gamma_rr[i];	  
	if(len < BL_MIN)      len = BL_MIN;
	else if(len > BL_MAX) len = BL_MAX;
      }
      PMat(len,tree->mod,tree->mod->ns*tree->mod->ns*i,b_fcus->Pij_rr);
    }
  }
}

/*********************************************************/

void Update_P_Lk_Along_A_Path(t_node **path, int path_length, t_tree *tree)
{
  int i,j;

  For(i,path_length-1)
    {
      For(j,3)
	if(path[i]->v[j] == path[i+1])
	  {
	    if(path[i] == path[i]->b[j]->left)
	      {
		Update_P_Lk(tree,path[i]->b[j],path[i]->b[j]->left);		    
	      }

	    else if(path[i] == path[i]->b[j]->rght)
	      {
		Update_P_Lk(tree,path[i]->b[j],path[i]->b[j]->rght);
	      }
	    else
	      {
		PhyML_Printf("\n. Err in file %s at line %d\n\n",__FILE__,__LINE__);
		Exit("");
	      }
	    break;
	  }      
#ifdef DEBUG
      if(j == 3)
	{
	  PhyML_Printf("\n. Err in file %s at line %d\n\n",__FILE__,__LINE__);
	  Exit("");
	}
#endif
    }
}

/*********************************************************/

phydbl Lk_Dist(phydbl *F, phydbl dist, model *mod)
{
  int i,j,k;
  phydbl lnL,len;
  int dim1,dim2;

  For(k,mod->n_catg)
    {
      len = dist*mod->gamma_rr[k];
      if(len < BL_MIN)      len = BL_MIN;
      else if(len > BL_MAX) len = BL_MAX;
      PMat(len,mod,mod->ns*mod->ns*k,mod->Pij_rr);
    }
  
  dim1 = mod->ns*mod->ns;
  dim2 = mod->ns;
  lnL = .0;
  For(i,mod->ns)
    {
      For(j,mod->ns)
	{
	  For(k,mod->n_catg)
	    {
 	      lnL += F[dim1*k+dim2*i+j] * LOG(mod->Pij_rr[dim1*k+dim2*i+j]);
	    }
	}
    }

  return lnL;
}

/*********************************************************/
//changed by Ken
phydbl Update_Lk_At_Given_Edge(t_edge *b_fcus, t_tree *tree)
{
  if(!b_fcus->left->tax) Update_P_Lk(tree,b_fcus,b_fcus->left);
  if(!b_fcus->rght->tax) Update_P_Lk(tree,b_fcus,b_fcus->rght);
  tree->c_lnL = Lk_At_Given_Edge(b_fcus,tree);
  return tree->c_lnL;
}

/*********************************************************/

/*       root
           \
           /
          a
	  |
	  |
	  d
	 / \
        /   \
       w     x    	

       d->t has changed and we need to compute
       the likelihood.
       (1) update the three branch lengths l(ad), l(dw) and l(dx)
       (2) update the change proba matrices along these branches
       (3) update the likelihood of subtree (w,x) (WARNING: (a,x) and (a,w) are not updated)
*/
phydbl Lk_Triplet(t_node *a, t_node *d, t_tree *tree)
{
  int i;
  phydbl max_height;
  phydbl up_bound, low_bound;

  if(d->tax)
    {
      PhyML_Printf("\n. Err in file %s at line %d\n\n",__FILE__,__LINE__);
      Warn_And_Exit("");
    }

  up_bound = low_bound = -1.0;
  max_height = -1.0;
  For(i,3)
    {
      if((d->v[i] != a) && (d->b[i] != tree->e_root))
	{
	  if(tree->rates->nd_t[d->v[i]->num] > max_height)
	    {
	      max_height = tree->rates->nd_t[d->v[i]->num];
	    }
	}
      else
	{
	  up_bound = 
	    (a == tree->n_root)?
	    (tree->rates->nd_t[a->num]):
	    (tree->rates->nd_t[d->v[i]->num]);
	}
    }

  low_bound = max_height;

  if(up_bound < low_bound - 1.E-10)
    {
      PhyML_Printf("\n. a->num=%d d->num=%d",a->num,d->num);
      PhyML_Printf("\n. up_bound = %f, low_bound = %f",up_bound,low_bound);
      Warn_And_Exit("\n");
    }

  if(tree->rates->nd_t[d->num] < low_bound) tree->rates->nd_t[d->num] = low_bound;
  else if(tree->rates->nd_t[d->num] > up_bound) tree->rates->nd_t[d->num] = up_bound;

  /* Step (1) */
  For(i,3)
    {
      if((d->v[i] != a) && (d->b[i] != tree->e_root)) 
	{
	  d->b[i]->l = 
	    (tree->rates->nd_t[d->num] - tree->rates->nd_t[d->v[i]->num]) * 
	    tree->rates->clock_r * 
	    tree->rates->br_r[d->b[i]->num];
	}
      else
	{
	  if(a == tree->n_root)
	    {
	      d->b[i]->l = 
		(tree->rates->nd_t[tree->n_root->num] - tree->rates->nd_t[tree->n_root->v[0]->num] + 
		 tree->rates->nd_t[tree->n_root->num] - tree->rates->nd_t[tree->n_root->v[1]->num]) * tree->rates->clock_r;
	    }
	  else
	    {
	      d->b[i]->l = (tree->rates->nd_t[a->num] - tree->rates->nd_t[d->num]) * tree->rates->clock_r * tree->rates->br_r[d->b[i]->num];	    
	    }
	}
    }
  
  /* Step (2) */
  For(i,3) Update_PMat_At_Given_Edge(d->b[i],tree);
  
  For(i,3) 
    if((d->v[i] == a) || (d->b[i] == tree->e_root))
      {
	Update_P_Lk(tree,d->b[i],d); 
	Lk_At_Given_Edge(d->b[i],tree);
	break;
      }

  return tree->c_lnL;
}

/*********************************************************/

void Print_Lk_Given_Edge_Recurr(t_node *a, t_node *d, t_edge *b, t_tree *tree)
{
  PhyML_Printf("\n___ Edge %3d (left=%3d rght=%3d) lnL=%f",
	 b->num,
	 b->left->num,
	 b->rght->num,
	 Lk_At_Given_Edge(b,tree));

  if(d->tax) return;
  else
    {
      int i;
      For(i,3)
	if(d->v[i] != a)
	  Print_Lk_Given_Edge_Recurr(d,d->v[i],d->b[i],tree);
    }
}

/*********************************************************/
void Init_Tips_At_One_Site_Codons_Float(char state, int pos, phydbl *p_lk, char * alternatives) //!< Added by Marcelo.
{
  int i,j,numSenseCodons;
  numSenseCodons=Genetic_code_ns();
  
  For(i,numSenseCodons) p_lk[pos+i]=0.0;
  
  if((state>=(char)0) && (state<(char)64)) p_lk[pos+indexSenseCodons[(int)state]]=1.0; 
  else
    if(state==(char)88)
    {
      j=-1;
      while(alternatives[++j]<(char)64) p_lk[pos+indexSenseCodons[(int)alternatives[j]]]=1.0;
    }else 
    {
      PhyML_Printf("\n. Unknown character state : %c\n",state);
      Exit("\n. Init failed (check the data type)\n");
    }
}
/*********************************************************/

void Init_Tips_At_One_Site_Codons_Int(char state, int pos, short int *p_pars, char * alternatives) //!< Added by Marcelo.
{
  int i,j,numSenseCodons;
  numSenseCodons=Genetic_code_ns();
  
  For(i,numSenseCodons) p_pars[pos+i]=0;
  
  if((state>=(char)0) && (state<(char)64)) p_pars[pos+indexSenseCodons[(int)state]]=1; 
  else
    if(state==(char)88)
    {
      j=-1;
      while(alternatives[++j]<(char)64)
      {
	p_pars[pos+indexSenseCodons[(int)alternatives[j]]]=1;
      }    
    }else 
    {
      PhyML_Printf("\n. Unknown character state : %c\n",state);
      Exit("\n. Init failed (check the data type)\n");
    }
}
/*********************************************************/
matrix *ML_CODONDist_Pairwise(calign *data, option *io, model *mod) //!<Added by Marcelo.
{
  int j,k,tmpkappa;
  model *mod_tmp;
  optimiz *s_opt;
  matrix *mat;
  calign *twodata,*tmpdata;
  phydbl lnL, *expt, *uexpt, brLen, mr;
  time_t tbegin, tend;
 
  //printf("\nml codon pairwise1\n");
  mod_tmp = Make_Model_Basic();
  s_opt = Make_Optimiz();
  Set_Defaults_Model(mod_tmp);
  mod_tmp->genetic_code = io->genCode;
  Genetic_code_index_stopCodons(mod_tmp->genetic_code);
  Set_Defaults_Optimiz(s_opt);
  mod_tmp->io = io;
  mod_tmp->s_opt = s_opt;
  mod_tmp->ns = mod->ns;
  mod_tmp->whichmodel = GY;
  mod_tmp->whichrealmodel = GY;
  mod_tmp->initqrates = io->init_DistanceTreeCD;
  mod_tmp->freq_model = FMODEL;
  mod_tmp->omegaSiteVar = NOOMEGA;
  tmpkappa = io->kappaECM;
  io->kappaECM = kap1;

  mod_tmp->constB = mod->constB;
  mod_tmp->nparts=mod->nparts;
  mod_tmp->n_otu = mod->n_otu;
  mod_tmp->init_len = mod->init_len;
  mod_tmp->state_len = mod->state_len;
  mod_tmp->whichrealmodel = mod->whichrealmodel;

  copyIOtoMod(mod_tmp->io,mod_tmp);


   mod_tmp->expm=mod->expm; /*!< 0 Eigenvalue; 1 scaling and squaring Pade. approx. COPIED FROM OPTION*/ //!< Added by Marcelo.
   mod_tmp->datatype=mod->datatype; //copied from option
   mod_tmp->init_DistanceTreeCD=mod->init_DistanceTreeCD; /*!< 0: ML JC69; 1: ML M0; 2: Schneider 2005 CodonPam; 3: Kosiol 2007 Empirical. COPIED FROM OPTION*/ //!<Added by Marcelo.
   mod_tmp->kappaECM=mod->kappaECM; /*!< According to Kosiol 2007.COPIED FROM OPTION*///!<Added by Marcelo.
   mod_tmp->n_termsTaylor=mod->n_termsTaylor; //! Added by Marcelo. COPIED FROM OPTION
   mod_tmp->heuristicExpm=mod->heuristicExpm; /*!< TAYLOR in parameter opt?.COPIED FROM OPTION*/
   mod_tmp->modeltypeOpt=mod->modeltypeOpt; //COPIED FROM OPTION
   mod_tmp->freqmodelOpt=mod->freqmodelOpt; //COPIED FROM OPTION
   mod_tmp->omegaOpt=mod->omegaOpt; //COPIED FROM OPTION
   mod_tmp->eq_freq_handling=mod->eq_freq_handling; //COPIED FROM OPTION
   mod_tmp->quiet=mod->quiet;
   mod_tmp->optParam=mod->optParam;
	  mod_tmp->nomega_part=mod_tmp->nparts;

   /*mod_tmp->opt_heuristic_manuel=mod->io->opt_heuristic_manuel;
   mod_tmp->opt_heuristic_manuel=mod->io->opt_heuristic_manuel;
   mod_tmp->roundMax_start=mod->io->roundMax_start;
   mod_tmp->roundMax_end=mod->io->roundMax_end;
   mod_tmp->logtree=mod->io->logtree;
   mod_tmp->print_site_lnl = mod->io->print_site_lnl;
   mod_tmp->out_stats_format = mod->io->out_stats_format;
   mod_tmp->random_boot_seq_order = mod->io->random_boot_seq_order;
   mod_tmp->minParam = mod->io->minParam;
   mod_tmp->maxParam = mod->io->maxParam;
   mod_tmp->lkExpStepSize = mod->io->lkExpStepSize;*/

  //mod_tmp->initqrates = mod->initqrates;
   //printf("expm %d\n",mod_tmp->expm);
  Make_Model_Complete(mod_tmp);
  //printf("expm %d\n",mod_tmp->expm);
  mod_tmp->s_opt->opt_kappa = NO;
  //printf("\nml codon pairwise2\n");
  mod_tmp->structTs_and_Tv  = (ts_and_tv *)mCalloc(64*64,sizeof(ts_and_tv));//!< Added by Marcelo. 64 possible codons ... will be initialized together with the model parameters in set model default.
  if(mod_tmp->datatype==CODON) Make_Ts_and_Tv_Matrix(mod_tmp->io,mod_tmp);



  expt            = (phydbl *)mCalloc(mod->ns,sizeof(phydbl));
  uexpt           = (phydbl *)mCalloc(mod->ns*mod->ns,sizeof(phydbl));
  
  tmpdata         = (calign *)mCalloc(1,sizeof(calign));
  tmpdata->c_seq  = (align **)mCalloc(2,sizeof(align *));
  tmpdata->b_frq  = (phydbl *)mCalloc(mod->ns,sizeof(phydbl));
  tmpdata->ambigu = (short int *)mCalloc(data->crunch_len,sizeof(short int));
  
  tmpdata->n_otu  = 2;
  tmpdata->crunch_len = data->crunch_len;
  tmpdata->init_len   = data->init_len;
  
  int i=0;
  For(i,mod_tmp->ns){
	  mod_tmp->pi[i]  = mod->pi[i];
	//  printf("first pi %lf\n",mod_tmp->pi[i]);
  }

  mat = Make_Mat(data->n_otu);
  Init_Mat(mat,data);
  if(mod->optDebug)printf("inited mat\n");
  //Added by Ken
  //Make Bmat array equivalent to GY94
  if(mod->whichrealmodel == HLP17){ //was io-> mod
	  /*int combinations = 3721;
	  phydbl hot[combinations];
	  //phydbl *tfr = ecmK07freq;
	  int i1, i2;
	  for(i1=0;i1<61;i1++){
		  for(i2=0;i2<61;i2++){
			  hot[i1*61+i2]=0;
		  }
	  }
	  mod_tmp->Bmat = hot;*/
	  mod_tmp->Bmat = (phydbl *)mCalloc(3721,sizeof(phydbl));
	  mod_tmp->nmotifs = 1;
	  mod_tmp->nhotness = 1;
	  //mod_tmp->hotspotcmps = malloc(sizeof(int *) * mod_tmp->nmotifs);
	  mod_tmp->hotness = malloc(sizeof(double ) * mod_tmp->nhotness);
	  mod_tmp->motif_hotness = malloc(sizeof(int ) * mod_tmp->nmotifs);
	  mod_tmp->hotness[0] = 0.0;
	  mod_tmp->motif_hotness[0]=0;
	  //mod_tmp->hotspotcmps[0] = mod->hotspotcmps[0]; //was io-> mod
  }

  //Add in stuff for partitioned model
  int modeli;
  for(modeli=0;modeli<mod->nparts;modeli++){ //was io-> mod
	  mod_tmp->nparts=1;
	  mod_tmp->qmat_part = (phydbl**)mCalloc(1,sizeof(phydbl*));
	  mod_tmp->partNames = (char**)mCalloc(1,sizeof(char*));
	  mod_tmp->omega_part = (phydbl*)mCalloc(1,sizeof(phydbl));
	  mod_tmp->qmat_part[0]=(phydbl *)mCalloc(mod->n_w_catg*mod->ns*mod->ns,sizeof(phydbl)); //was io-> mod
	  mod_tmp->partNames[0]="SINGLE";
	  mod_tmp->omega_part[0]=0.4;
  }
  //printf("first pi %lf %lf\n",mod_tmp->pi[0],mod->pi[0]);
  if(io->mod->optDebug)printf("initqs %d %d %d\n",mod_tmp->initqrates,mod->initqrates,io->mod->initqrates);
  mr = Update_Qmat_Codons(mod_tmp, 0,0);
  if(io->mod->optDebug)printf("first pi %lf %lf\n",mod_tmp->pi[0],mod->pi[0]);

  //printf("\n. FYI: Using single partitioned GY94 for initial tree.\n\n");
  EigenQREV(mod_tmp->qmat_part[0], mod_tmp->pi, mod->ns, mod_tmp->eigen->e_val, mod_tmp->eigen->r_e_vect, mod_tmp->eigen->l_e_vect, mod_tmp->eigen->space);
 //printf("first pi %lf %lf\n",mod_tmp->pi[0],mod->pi[0]);

  For(j, mod->ns) mod_tmp->eigen->e_val[j] /= mr;
  
  brLen=1.0;
  
  #if defined OMP || defined BLAS_OMP

  tbegin= omp_get_wtime();
  
  #else

  time(&tbegin);
 
  #endif
  int progress = 0; //!Added by Marcelo.
  double loops = ((data->n_otu)*(data->n_otu)-(data->n_otu))/2.00;//!Added by Marcelo.
  if(io->mod->optDebug)PhyML_Printf("\n. Computing pairwise distances...");
  For( j, data->n_otu-1 ) {
    tmpdata->c_seq[0]       = data->c_seq[j];
    tmpdata->wght           = data->wght;
    
    for(k=j+1;k<data->n_otu;k++) {
      tmpdata->c_seq[1]       = data->c_seq[k];
      twodata = Compact_Cdata(tmpdata,io,mod);

      //printf("first pi %lf %lf\n",mod_tmp->pi[0],mod->pi[0]);
      lnL=Br_Len_Brent_Codon_Pairwise(BL_MIN, brLen, BL_MAX, mod->s_opt->min_diff_lk_local, &brLen, mod_tmp->Pij_rr, mod_tmp->pi, mod_tmp->eigen, twodata, mod->ns, mod->s_opt->brent_it_max, mod->s_opt->quickdirty, uexpt, expt);
            //was io-> mod
      mat->dist[j][k] = brLen;
      mat->dist[k][j] = mat->dist[j][k];
      Free_Cseq(twodata);
      progress++;
      //if(!io->quiet) PhyML_Printf("\r. Computing pairwise distances...%3.2f%c concluded.", progress*100.00/loops, '%');
      //PhyML_Printf("\n. Computing pairwise distances...");
    }
  }
  if(!io->quiet) PhyML_Printf("\n");
  #if defined OMP || defined BLAS_OMP

  tend= omp_get_wtime();
  
  #else

  time(&tend);
 
  #endif
  
  if(io->ntrees == 1 || io->mod->optDebug){
	  if((tend-tbegin)>=1) {
	  	  PhyML_Printf("\n. Calculation of pairwise distances completed in %ld sec.\n",(long int)(tend-tbegin));
  	  } else {
  		  PhyML_Printf("\n. Calculation of pairwise distances completed in less than 1 sec.\n");
  	  }
  }
  
  tmpdata->c_seq[0] = NULL;
  tmpdata->c_seq[1] = NULL;
  tmpdata->wght     = NULL;
  
  free(tmpdata->ambigu);
  free(tmpdata->b_frq);
  free(tmpdata->c_seq);
  free(tmpdata);
  free(expt); 
  free(uexpt);
  Free_Optimiz(s_opt);
  Free_Model(mod_tmp);
  
  io->kappaECM = tmpkappa;

  return mat;
}

/*********************************************************/
phydbl LK_Codon_Pairwise(calign *data, phydbl *Pij, phydbl *pi, int ns, phydbl len, eigen *eigenStruct, phydbl *uexpt, phydbl *expt) //!< Added by Marcelo.
{
  phydbl cn_lk;
  // phydbl p_sum0;
  int i;
  // int j,k,m,n;
  char *seq0, *seq1;
  // char state0, state1;
  
  PMat_CODON_Pairwise(len, Pij, eigenStruct->r_e_vect, eigenStruct->l_e_vect, eigenStruct->e_val, ns, uexpt, expt);

  cn_lk=0.0;
  seq0=data->c_seq[0]->state;
  seq1=data->c_seq[1]->state;
  
  #if defined OMP || defined BLAS_OMP
  #pragma omp parallel for reduction(+:cn_lk)
  #endif
  
  For(i,data->crunch_len)
  {
    if((seq0[i]!=(char)88)&&(seq1[i]!=(char)88))
    {
      cn_lk+=data->wght[i]*LOG(pi[indexSenseCodons[(int)seq0[i]]]*Pij[indexSenseCodons[(int)seq0[i]]*ns+indexSenseCodons[(int)seq1[i]]]);
    }
  } 
  

  return cn_lk;
}


/*********************************************************/
phydbl LK_BFGS_from_CODEML(option* io, phydbl *x, int n){
  int  i, numParams=0;
  //if(io->mod->optDebug)printf("lk bfgs\n");
  if (io->mod->datatype==CODON){
    if(io->mod->s_opt->opt_kappa && (io->mod->optKappa==1||(io->mod->optKappa==2 && io->mod->optIter==0))){
      if(io->mod->whichmodel!=GYECMK07WK && io->mod->whichmodel!=GYECMK07WKF &&
		io->mod->whichmodel!=GYECMS05WK && io->mod->whichmodel!=GYECMS05WKF
		&& io->mod->whichmodel!=MGECMK07WK && io->mod->whichmodel!=MGECMK07WKF &&
		io->mod->whichmodel!=MGECMS05WK && io->mod->whichmodel!=MGECMS05WKF
		&& io->mod->whichmodel!=YAPECMK07WK && io->mod->whichmodel!=YAPECMK07WKF &&
		io->mod->whichmodel!=YAPECMS05WK && io->mod->whichmodel!=YAPECMS05WKF &&
		io->mod->whichmodel!=GYECMUSRWK && io->mod->whichmodel!=GYECMUSRWKF &&
		io->mod->whichmodel!=MGECMUSRWKF && io->mod->whichmodel!=YAPECMUSRWKF &&
		io->mod->whichmodel!=MGECMUSRWK &&io->mod->whichmodel!=YAPECMUSRWK ){
			io->mod->kappa = x[numParams++];
      }else{
		switch(io->mod->kappaECM){
		  case kap6:
		  case kap2:
		  case kap3:          io->mod->pkappa[0] = x[numParams++];break;
		  case kap4:          io->mod->pkappa[0] = x[numParams++]; io->mod->pkappa[1] = x[numParams++]; break;
		  case kap5: For(i,io->mod->nkappa-1) io->mod->unspkappa[i] = x[numParams++]; break;
		  default: break;
		}
      }
    }
    
    if(io->mod->s_opt->opt_omega){
      if(io->mod->omegaSiteVar==DM0){
    	  int omegai; //added by Ken 18/8/2016
    	  for(omegai=0;omegai<io->mod->nomega_part;omegai++){
    		  if(io->mod->omega_part_opt[omegai] == 1 || (io->mod->optIter==0 && io->mod->omega_part_opt[omegai] == 2)){
    			  io->mod->omega_part[omegai] = x[numParams++];
    		  }
    		  //printf("opting o %lf\n",io->mod->omega_part[omegai]);
    	  }
      }else if(io->mod->omegaSiteVar==DMODELK){
		For(i,io->mod->n_w_catg) io->mod->omegas[i] = x[numParams++];
		For(i,io->mod->n_w_catg-1) io->mod->prob_omegas_uns[i] = x[numParams++];
      }else if(io->mod->omegaSiteVar==DGAMMAK){
		io->mod->alpha                = x[numParams++];
		io->mod->beta                 = x[numParams++];
      }
    }

    if(io->mod->opthotness){ //added by Kenneth Hoehn 3/6/2016
    	int c;
    	for(c=0;c<io->mod->nhotness;c++){
    		if(io->mod->hoptindex[c] == 1 || (io->mod->optIter==0 && io->mod->hoptindex[c] == 2)){
    			io->mod->hotness[c] = x[numParams++];
    		}
    	}
     }

    if(io->mod->s_opt->opt_state_freq){
      switch(io->mod->freq_model){
		case F1XSENSECODONS: For(i,io->mod->num_base_freq-1) io->mod->pi_unscaled[i] = x[numParams++]; break;
		case F1X4: For(i,io->mod->num_base_freq-1) io->mod->uns_base_freq[i] = x[numParams++]; break;
		case F3X4:
		case CF3X4:{
	  		for(i=0;i<3;i++) io->mod->uns_base_freq[i] = x[numParams++];
	  		for(i=4;i<7;i++) io->mod->uns_base_freq[i] = x[numParams++];
	  		for(i=8;i<11;i++)io->mod->uns_base_freq[i] = x[numParams++];
	  		break;
		}
		default:
		  break;
      }
    }

    if(io->mod->s_opt->opt_pinvar) io->mod->pinvar = x[numParams++];
    
    if(io->mod->n_catg>1 && io->mod->n_w_catg==1 && io->mod->s_opt->opt_alphaCD) io->mod->alpha = x[numParams++];
    
    if(io->mod->pcaModel==1) For(i,io->mod->npcs) io->mod->pcsC[i] = x[numParams++];
    
    if(n!=numParams) Warn_And_Exit("Number of parameters optimized does not match.");
    
  }else if(io->mod->datatype==AA){
    
    if(io->mod->s_opt->opt_state_freq) For(i,io->mod->ns-1) io->mod->pi_unscaled[i] = x[numParams++];
    
    if(n!=numParams) Warn_And_Exit("Number of parameters optimized does not match.");
  }
  //if(io->mod->optDebug)printf("calling lk_rep\n");
  return -Lk_rep(io);
}
/*********************************************************/
