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

#include "optimiz.h"
#include "io.h"

/*********************************************************/

void CI_Wrapper(option* io){
	int i;
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


/*********************************************************/

void findCIs(model* mod, option *io, FILE* CI){
	int i;
	phydbl ori_min = io->mod->s_opt->min_diff_lk_global;
	io->mod->s_opt->min_diff_lk_global *= io->roughCI;

	For(i,mod->nomega_part){
	  		  if(mod->omega_part_ci[i]==1){
	  			  printf("\n\n. Estimating 95%% confidence intervals for omega %d model %d",i,mod->num);
	  			  printf("\n. Target likelihood: %.4f",io->replnL-1.92);
	  			  int opt=mod->omega_part_opt[i];
	  	  	  	  if(mod->primary)mod->omega_part_opt[i]=0;
	  	  	  	  else mod->omega_part_opt[i]=3;

	  	  	  	  char* buf= mCalloc(T_MAX_OPTION,sizeof(char));
	  	  	  	  int temp = asprintf(&buf, "omega %d", i);
	  	  	  	  phydbl upper=binarySearchCI(&mod->omega_part[i],io,0.01,0.5,0.0,100.0,CI,buf);
	  	  	  	  phydbl lower=binarySearchCI(&mod->omega_part[i],io,0.01,-0.2,0.0,100.0,CI,buf);
	  	  	  	  printf("\n. Estimated CI: %.4f (%.4f, %.4f)\n",mod->omega_part[i],lower,upper);
	  	  	  	  mod->omega_part_uci[i]=upper;
	  	  	  	  mod->omega_part_lci[i]=lower;
	  	  	  	  mod->omega_part_opt[i]=opt;
	  		  }
	  }
	  For(i,mod->nhotness){
	  	if(mod->hoptci[i]==1){
	  			printf("\nEstimating 95%% confidence intervals for h %d model %d\n",i,mod->num);
	  			printf("\n. Target likelihood: %.4f",io->replnL-1.92);
	  			 int opt=mod->hoptindex[i];
	  			 if(mod->primary)mod->hoptindex[i]=0;
	  			 else mod->hoptindex[i]=3;

	  			 char* buf= mCalloc(T_MAX_OPTION,sizeof(char));
	  			 int temp = asprintf(&buf, "hotness %d", i);
	  			 phydbl upper=binarySearchCI(&mod->hotness[i],io,0.05,2.0,-1.0,100.0,CI,buf);
	  			 phydbl lower=binarySearchCI(&mod->hotness[i],io,0.05,-0.5,-1.0,100.0,CI,buf);
	  			 printf("\n. Estimated CI: %.4f (%.4f, %.4f)\n",mod->hotness[i],lower,upper);
	  			 mod->hoptuci[i]=upper;
	  			 mod->hoptlci[i]=lower;
	  			 mod->hoptindex[i]=opt;
	  		}
	  	  }
	  	  if(mod->kappaci==1){
	  		printf("\nEstimating 95%% confidence intervals for kappa model %d\n",mod->num);
	  		printf("\n. Target likelihood: %.4f",io->replnL-1.92);
	  		  	int opt=mod->s_opt->opt_kappa;
	  		  	if(mod->primary)mod->optKappa=0;
	  		  	else mod->optKappa=3;
	  		  	phydbl upper=binarySearchCI(&mod->kappa,io,0.05,0.5,-1,100.0,CI,"kappa");
	  		  	phydbl lower=binarySearchCI(&mod->kappa,io,0.05,-1.0,-1,100.0,CI,"kappa");
	  		  	printf("\n. Estimated CI: %.4f (%.4f, %.4f)\n",mod->kappa,lower,upper);
	  		  	mod->kappauci=upper;
	  		  	mod->kappalci=lower;
	  		  	mod->s_opt->opt_kappa=opt;
	  	  }
	  	io->mod->s_opt->min_diff_lk_global=ori_min;
}

/*********************************************************/

/* Binary search for upper/lower CI
 * Before this can be run:
 * 1. Must have the initial value of the parameter as its MLE
 * 2. Optimization for this parameter must be turned off
 */
phydbl binarySearchCI(phydbl* param,option* io,phydbl tol,phydbl delta,phydbl lowerb,phydbl upperb,FILE* CI, char* ID){
	io->both_sides=1;
	io->mod->update_eigen=1;
	Lk_rep(io);
	phydbl mle=*param;
	phydbl target=io->replnL-1.92;//95% CI value
	phydbl ol = io->replnL;
	phydbl nl=0;
	fprintf(CI,"%s %lf %lf %lf %lf %lf %lf\n",ID,target,io->replnL,*param,0.0,0.0,0.0);
	phydbl* paramStore=mCalloc(io->nparams,sizeof(phydbl));

	//!!!Record params and branch lengths
	storeParams(io,1,paramStore);
	io->mod->quiet=1;
	//find a point on the other side of CI
	phydbl d1 = delta;
	phydbl bound;
	int breaknext=0;
	int iter=1;
	resetSubstParams(io);
	do{//need to include lower boundary condition!
		*param=mle+d1;
		if(*param <= lowerb){
			*param=lowerb+SMALL;
			breaknext=1;
		}
		if(*param >= upperb){
			*param=upperb;
			breaknext=2;
		}
		bound=*param;
		Round_Optimize(io,ROUND_MAX*2);
		nl=Lk_rep(io);
		if(nl > ol){
			printf("Optimization failed! Interval is higher lhood than MLE!\n");
			printf("%lf\t%lf\n",nl,ol);
			exit(EXIT_FAILURE);
		}
		Print_Lk_rep(io,"[parameter value    ]");
		PhyML_Printf("[%.2f ]",*param);
		resetSubstParams(io);
		phydbl nl2=Lk_rep(io);
		if(breaknext>0){
			if(breaknext==1){
				resetSubstParams(io);
				//restoreParams(io,1,paramStore);
				if(nl>=target)return lowerb;
				else break;
			}else{
				//restoreParams(io,1,paramStore);
				resetSubstParams(io);
				if(nl>=target)return upperb;
				else break;
			}
		}
		//!Restore parameter values and branch lengths if changed!
		d1+=delta*iter; //make iterations progressively larger
		iter++;
	}while(nl >= target);
	if(io->mod->optDebug)printf("\nFound boundary %lf",bound);

	phydbl b1=mle;
	phydbl b2=bound;
	do{
		*param=(b2+b1)/2;
		//else *param=(b1-b2)/2;
		Round_Optimize(io,ROUND_MAX*2);
		nl=Lk_rep(io); //needs to be Round_Optimize
		Print_Lk_rep(io,"[parameter value    ]");
		PhyML_Printf("[%.2f ]",*param);
		fprintf(CI,"%s %lf %lf %lf %lf %lf %lf\n",ID,target,nl,*param,b1,b2,fabs(b1-b2));
		if(nl>=target)b1=*param;//want most conservative estimate of CI
		else b2=*param;
		//Reset parameter values and branch lengths if changed!
		resetSubstParams(io);
		phydbl nl2=Lk_rep(io);
		fabs(b1-b2);
	}while(fabs(b1-b2) >= tol);

	phydbl cb = (b1+b2)/2;
	*param=cb;
	Round_Optimize(io,ROUND_MAX*2);
	//return midpoint between the two intervals
	if(io->mod->optDebug)printf("\n%lf %lf %lf %lf %lf %lf",target,nl,*param,b1,b2,fabs(b1-b2));
	//reset params
	restoreParams(io,1,paramStore);
	Lk_rep(io);
	return (b1+b2)/2;
}
/*********************************************************/

/* Store parameters for optimization*/
int storeParams(option* io, int reseto, phydbl* ar){
	//phydbl* ar=io->paramStore;
	int c=0;
	int i,j,k;
	//store repertoire parameters
	ar[c++]=io->replnL;
	ar[c++]=io->mod->kappa;
	For(i,io->mod->nomega_part)ar[c++]=io->mod->omega_part[i];
	For(i,12)ar[c++]=io->mod->uns_base_freq[i];
	if(io->mod->whichrealmodel<=HLP17)For(i,io->mod->nhotness)ar[c++]=io->mod->hotness[i];
	//store subtree model parameters and branch lengths
	For(j,io->ntrees){
		model* mod=io->mod_s[j];
		t_tree* tree=io->tree_s[j];
		ar[c++]=tree->c_lnL;
		ar[c++]=mod->kappa;
		For(i,io->mod->nomega_part)ar[c++]=mod->omega_part[i];
		For(i,12)ar[c++]=mod->uns_base_freq[i];
		if(io->mod->whichrealmodel<=HLP17)For(i,io->mod->nhotness)ar[c++]=mod->hotness[i];
		For(k,2*tree->n_otu-3){
			  ar[c++]=tree->t_edges[k]->l;
		 }
	}
	//assumes you want to reset optimization params
	if(reseto){
		io->SIZEp=0;
		  io->noisy=0;
		  io->Iround=0;
		  io->NFunCall=0;
		  io->AlwaysCenter=0;
		  io->gemin=1e-6;
		  io->Small_Diff=.5e-6;
		  io->both_sides=1;
	}
	return c;
}
/*********************************************************/

/* Set substitution parameters to defaults in HLP19*/
int resetSubstParams(option* io){
	//phydbl* ar=io->paramStore;
	int c=0;
	int i,j,k;
	//store repertoire parameters
	io->mod->kappa=1.0;
	For(i,io->mod->nomega_part)io->mod->omega_part[i]=0.4;
	if(io->mod->whichrealmodel<=HLP17)For(i,io->mod->nhotness)io->mod->hotness[i]=0.0;

	//resotre original branch lengths
	For(i,io->ntrees){
		For(j,io->mod_s[i]->nedges){
			io->tree_s[i]->t_edges[j]->l = io->tree_s[i]->t_edges[j]->ol;
		}
	}

	if(io->repwidefreqs){ //calculate repertoire-wide codon frequencies
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
		  For(j,12)io->mod->base_freq[j]=roundf(io->mod->base_freq[j]*10000.0f)/10000.0f; //round base freqs to 5 decimal places
		  if(io->mod->optDebug)For(j,12){printf("%lf\t%lf\n",io->mod->baseCounts[j],io->mod->base_freq[j]);}
		  CF3x4(io->mod->base_freq, io->mod->genetic_code);
  	  }else{
	  	  For(i,12)io->mod->base_freq[i]=io->mod->user_b_freq[i];
	  	  CF3x4(io->mod->base_freq, io->mod->genetic_code);
  	  }
  	  if(io->mod->optDebug)printf("cf3x4\n");
  	  if(io->mod->optDebug)For(j,12){printf("%lf\t%lf\n",io->mod->baseCounts[j],io->mod->base_freq[j]);}
  	  Freq_to_UnsFreq(io->mod->base_freq,   io->mod->uns_base_freq,   4, 1);
  	  Freq_to_UnsFreq(io->mod->base_freq+4, io->mod->uns_base_freq+4, 4, 1);
  	  Freq_to_UnsFreq(io->mod->base_freq+8, io->mod->uns_base_freq+8, 4, 1);
  	  if(io->mod->optDebug)For(j,12)printf("%lf\t%lf\t%lf\t%lf\n",io->mod->baseCounts[j],io->mod->base_freq[j],io->mod->uns_base_freq[j],io->mod_s[0]->base_freq[j]);
	}else if(io->mod->whichrealmodel == HLP19){ //calculate midpoint divergence along tree
		For(i,io->ntrees){ //calculate midpoint divergence for each tree
			t_tree* tree = io->tree_s[i];
			phydbl tl = Get_Total_Divergence(tree->noeud[tree->mod->startnode], tree->noeud[tree->mod->startnode]->v[0],
					tree->noeud[tree->mod->startnode]->b[0],0.0, tree);
			tree->mod->midpoint_div = tl/(tree->n_otu*2.0);
			if(io->mod->optDebug)printf("Total divergence: %lf %lf\n",tl,tree->mod->midpoint_div);
		}
		if(io->mod->freq_model==MROOT)Setup_Midpoint_Flux(io);
		else Warn_And_Exit("HLP19 must be run with freq_model == MROOT");
	}else{ //set up empirical frequencies for individual lineages (GY, or if -f empirical)
		 For(i,io->ntrees){
			model* mod = io->tree_s[i]->mod;
			CF3x4(mod->base_freq, mod->genetic_code); //cut these segments?
			EqFrequencies(mod->freq_model, mod->pi, mod->base_freq, mod->ns);
		}
	}
	return c;
}
/*********************************************************/

/*restore parameters after optimization*/
int restoreParams(option* io, int reseto, phydbl* ar){
	int c=0;
	int i,j,k;
	//restore repertoire parameters
	io->replnL=ar[c++];
	io->mod->kappa=ar[c++];
	For(i,io->mod->nomega_part)io->mod->omega_part[i]=ar[c++];
	For(i,12)io->mod->uns_base_freq[i]=ar[c++];
	if(io->mod->whichrealmodel<=HLP17)For(i,io->mod->nhotness)io->mod->hotness[i]=ar[c++];
	//restore subtree model parameters and branch lengths
	For(j,io->ntrees){
		model* mod=io->mod_s[j];
		t_tree* tree=io->tree_s[j];
		tree->c_lnL=ar[c++];
		mod->kappa=ar[c++];
		For(i,io->mod->nomega_part)mod->omega_part[i]=ar[c++];
		For(i,12)mod->uns_base_freq[i]=ar[c++];
		if(io->mod->whichrealmodel<=HLP17)For(i,io->mod->nhotness)mod->hotness[i]=ar[c++];
		For(k,2*tree->n_otu-3){
			  tree->t_edges[k]->l=ar[c++];
		 }
	}
	io->mod->update_eigen=1;
		//assumes you want to reset optimization params
	if(reseto){
	  io->SIZEp=0;
	  io->noisy=0;
	  io->Iround=0;
	  io->NFunCall=0;
	  io->AlwaysCenter=0;
	  io->gemin=1e-6;
	  io->Small_Diff=.5e-6;
	  io->both_sides=1;
	}
	return c;
}



/*********************************************************/

phydbl Br_Len_Brent(phydbl ax, phydbl bx, phydbl cx, phydbl tol,
		    t_edge *b_fcus, t_tree *tree, int n_iter_max, int quickdirty){
  int iter;
  phydbl a,b,d,etemp,fu,fv,fw,fx,p,q,r,tol1,tol2,u,v,w,x,xm;
  phydbl e=0.0;
  phydbl old_lnL, init_lnL;

  d=0.0;
  a=((ax < cx) ? ax : cx);
  b=((ax > cx) ? ax : cx);
  x=w=v=bx;
  old_lnL = UNLIKELY;
  b_fcus->l = FABS(bx);
  fw=fv=fx=fu=-Lk_At_Given_Edge(b_fcus,tree);
  init_lnL = -fw;

  for(iter=1;iter<=BRENT_ITMAX;iter++){
      xm=0.5*(a+b);
      tol2=2.0*(tol1=tol*FABS(x)+BRENT_ZEPS);

      if((tree->c_lnL > init_lnL + tol) && (quickdirty)){
    	  b_fcus->l = x;
    	  Lk_At_Given_Edge(b_fcus,tree);
    	  return tree->c_lnL;
      }
      if((FABS(tree->c_lnL-old_lnL) < tol) || (iter > n_iter_max - 1)){
    	  b_fcus->l=x;
    	  Lk_At_Given_Edge(b_fcus,tree);
    	  return tree->c_lnL;
      }
      if(FABS(e) > tol1){
    	  r=(x-w)*(fx-fv);
    	  q=(x-v)*(fx-fw);
    	  p=(x-v)*q-(x-w)*r;
    	  q=2.0*(q-r);
    	  if(q > 0.0) p = -p;
    	  q=FABS(q);
    	  etemp=e;
    	  e=d;
    	  if(FABS(p) >= FABS(0.5*q*etemp) || p <= q*(a-x) || p >= q*(b-x))
    		  d=BRENT_CGOLD*(e=(x >= xm ? a-x : b-x));
    	  else{
    		  d=p/q;
    		  u=x+d;
    		  if (u-a < tol2 || b-u < tol2)
    			  d=SIGN(tol1,xm-x);
    	  }
	}else{
	  d=BRENT_CGOLD*(e=(x >= xm ? a-x : b-x));
	}
    u=(FABS(d) >= tol1 ? x+d : x+SIGN(tol1,d));
    b_fcus->l=FABS(u);
    old_lnL = tree->c_lnL;
    fu=-Lk_At_Given_Edge(b_fcus,tree);

    if(fu < fx){
	  if(u > x) a=x; else b=x;
	  SHFT(v,w,x,u)
	  SHFT(fv,fw,fx,fu)
	}else{
	  if (u < x) a=u; else b=u;
	  if (fu < fw || FABS(w-x) < SMALL){
	      v=w;
	      w=u;
	      fv=fw;
	      fw=fu;
	  }else if (fu < fv || FABS(v-x) < SMALL || FABS(v-w) < SMALL){
	      v=u;
	      fv=fu;
	    }
	}
  }
  if(iter > BRENT_ITMAX) PhyML_Printf("\n. Too many iterations in BRENT (%d) (%f)",iter,b_fcus->l);
  return(-1);
}

/*********************************************************/
//Modified by Ken
void Round_Optimize(option *io, int n_round_max){
  int n_round,each;
  phydbl lk_old, lk_new, tol;

  lk_new = io->replnL;
  lk_old = UNLIKELY;
  n_round = 0;
  each = 0;
  tol = 1.e-2;
  int i,j;

  while(n_round < n_round_max){
      if(io->tree_s[0]->has_branch_lengths){//!< Added by Marcelo ... in this case opt parameters first
      	if(!each){
	  		each = 1;
	  		if(io->mod->optDebug)printf("optimizing all free params\n");
	  		Optimiz_All_Free_Param(io,(io->mod->quiet)?(0):(io->mod->s_opt->print),0);
		}
		if(io->mod->s_opt->opt_bl){
#if defined OMP || defined BLAS_OMP
#pragma omp parallel for if(io->splitByTree)
#endif
			For(i,io->ntrees){
				int n_roundt=n_round;
				t_tree* tree = io->tree_s[i];
				t_node* root = tree->noeud[tree->mod->startnode];
				if(tree->mod->whichrealmodel > HLP17)(!((n_roundt+2)%2))?(root=tree->noeud[0]):(root=tree->noeud[tree->n_otu-1]);
				Lk(tree);
				if(tree->mod->whichrealmodel <= HLP17){Get_UPP(root, root->v[0], tree);}
				//printf("startnode %d %d\n",root->num, root->v[0]->num);
				Optimize_Br_Len_Serie(root,root->v[0],root->b[0],tree,tree->data);
			}
			Lk_rep(io);
			if(!io->mod->quiet)Print_Lk_rep(io,"[Branch lengths     ]");
		}
		if(io->ntrees>1)Lk_rep(io);
		lk_new = io->replnL;
      }else{
		if(io->mod->s_opt->opt_bl){
#if defined OMP || defined BLAS_OMP
#pragma omp parallel for if(io->splitByTree)
#endif
			For(i,io->ntrees){
				int n_roundt=n_round;
				t_tree* tree = io->tree_s[i];
				t_node* root = tree->noeud[tree->mod->startnode];
				if(tree->mod->whichrealmodel > HLP17)(!((n_roundt+2)%2))?(root=tree->noeud[0]):(root=tree->noeud[tree->n_otu-1]);
				tree->both_sides = 1;
				Lk(tree);
				if(tree->mod->whichrealmodel <= HLP17){Get_UPP(root, root->v[0], tree);}
				Optimize_Br_Len_Serie(root,root->v[0],root->b[0],tree,tree->data);
			}
			Lk_rep(io);
			if(!io->mod->quiet)Print_Lk_rep(io,"[Branch lengths     ]");
		}
		if(!each){
	  		each = 1;
	  		Optimiz_All_Free_Param(io,(io->mod->quiet)?(0):(io->mod->s_opt->print),0);
		}
      }
      io->both_sides = 1;
      Lk_rep(io);
      lk_new = io->replnL;
      if(lk_new < lk_old - io->mod->s_opt->min_diff_lk_global){
    	  printf("Old: %lf, New: %lf\n",lk_old,lk_new);
    	  Exit("\n. Optimisation failed ! (Round_Optimize)\n");
      }
      if(FABS(lk_new - lk_old) < io->mod->s_opt->min_diff_lk_global)  break;
      else lk_old  = lk_new;
      n_round++;
      each--;
   }
  Optimiz_All_Free_Param(io,(io->mod->quiet)?(0):(io->mod->s_opt->print),0);
}

/*********************************************************/
//Edited by Ken
void Optimize_Br_Len_Serie(t_node *a, t_node *d, t_edge *b_fcus, t_tree *tree, calign *cdata){
  int i;
  phydbl l_infa,l_max,l_infb;
  phydbl lk_init;
  lk_init = tree->c_lnL;
  l_infa = l_max  = l_infb = BL_MIN;
  l_infa = BL_MAX;
  l_max  = b_fcus->l;
  l_infb = BL_MIN;

  if(a->num == tree->mod->startnode){
	  l_infa = tree->mod->maxtrunkl;
	  if(b_fcus->l > tree->mod->maxtrunkl){
		  // printf("\n%d\t%d\t%lf",a->num,tree->mod->startnode,b_fcus->l);
		  b_fcus->l = tree->mod->maxtrunkl;
		  l_max = b_fcus->l;
		  Lk(tree);
		  lk_init = tree->c_lnL;
	  }
  }
  //printf("\np %d\t%d\t%lf",b_fcus->anc_node->num,tree->mod->startnode,b_fcus->l);

  if(tree->mod->optDebug)printf("\n%d\t%lf\t%lf",b_fcus->num,b_fcus->l,tree->c_lnL);
  Br_Len_Brent(l_infa,l_max,l_infb,
	       tree->mod->s_opt->min_diff_lk_local,
	       b_fcus,tree,
	       tree->mod->s_opt->brent_it_max,
	       tree->mod->s_opt->quickdirty);
  //Added to catch potential issues with branch optimization
  if(tree->mod->optDebug)printf("\n%d\t%lf\t%lf",b_fcus->num,b_fcus->l,tree->c_lnL);
  if(tree->mod->whichrealmodel <= HLP17){
	  if(tree->c_lnL < lk_init - tree->mod->s_opt->min_diff_lk_local){
		  PhyML_Printf("\n. %f %f %f %f %d %d %d\n",l_infa,l_max,l_infb,b_fcus->l,b_fcus->num,a->num,d->num);
		  PhyML_Printf("\n. %f -- %f\n",lk_init,tree->c_lnL);
		  Warn_And_Exit("\n. Err. in Optimize_Br_Len_Serie\n");
	  }
  }
  //Lk(tree);
    /*if(b_fcus->anc_node->num == tree->mod->startnode){
  	  if(b_fcus->l > tree->mod->maxtrunkl){
  		  b_fcus->l = tree->mod->maxtrunkl;
  	  }
  }*/
  //printf("\na %d\t%d\t%lf",b_fcus->anc_node->num,tree->mod->startnode,b_fcus->l);

  if(d->tax) return;

  else For(i,3) if(d->v[i] != a){
	   if(tree->mod->freq_model!=ROOT){
		   Update_P_Lk(tree,d->b[i],d);
		   if(tree->mod->whichrealmodel <= HLP17){
			   Fill_UPP_single(tree,d->b[i]);
		   }
	   }
	   Optimize_Br_Len_Serie(d,d->v[i],d->b[i],tree,cdata);
   }
  if(tree->mod->freq_model!=ROOT)For(i,3) if((d->v[i] == a) && !(d->v[i]->tax)) Update_P_Lk(tree,d->b[i],d);
}

/*********************************************************/

void Optimiz_Ext_Br(t_tree *tree){
  int i;
  t_edge *b;
  phydbl l_infa,l_max,l_infb;
  phydbl lk, lk_init,l_init;
  lk_init = tree->c_lnL;
  For(i,2*tree->n_otu-3){
      b = tree->t_edges[i];
      if((b->left->tax) || (b->rght->tax)){
    	  l_init = b->l;
    	  l_infa = 10.*b->l;
    	  l_max  = b->l;
    	  l_infb = BL_MIN;
	  	  lk = Br_Len_Brent(l_infa,l_max,l_infb,
			    tree->mod->s_opt->min_diff_lk_local,
			    b,tree,
			    tree->mod->s_opt->brent_it_max,
			    tree->mod->s_opt->quickdirty);

	  	  b->nni->best_l    = b->l;
	  	  b->nni->l0        = b->l;
	  	  b->nni->best_conf = 0;
	  	  b->l              = l_init;
	  	  if(tree->mod->whichrealmodel <= HLP17)
	  		  Update_PMat_At_Given_Edge(b,tree);
	}
  }
  tree->c_lnL = lk_init; 
}

/*********************************************************/
//Optimize submodel parameters using Brent
void Optimiz_Submodel_Params(option* io,int verbose){
	int i,n,c;

	io->threads=0;
#if defined OMP || defined BLAS_OMP
#pragma omp parallel for private(i,n,c)
#endif

	For(i,io->ntrees){
		t_tree* tree;

#pragma omp critical
		{
			tree=io->tree_s[io->threads++];
		}

		tree = io->tree_s[i];
		tree->mod->update_eigen=1;
		For(n, tree->mod->s_opt->nBrentCycles){
			if(tree->mod->optKappa==2){
				Generic_Brent_Lk(&(tree->mod->kappa),
						     TREEBOUNDLOW,TREEBOUNDHIGH,
						     tree->mod->s_opt->min_diff_lk_global,
						     tree->mod->s_opt->brent_it_max,
						     tree->mod->s_opt->quickdirty,
						     Optwrap_Lk,NULL,tree,NULL);
			}
		if(tree->mod->s_opt->opt_omega){
			  if(tree->mod->omegaSiteVar==DM0){
				  For(c,tree->mod->nomega_part){
					  if(tree->mod->omega_part_opt[c]==2){
					  	  Generic_Brent_Lk(&(tree->mod->omega_part[c]),
					  			  TREEBOUNDLOW,TREEBOUNDHIGH,
								  tree->mod->s_opt->min_diff_lk_global,
								  tree->mod->s_opt->brent_it_max,
								  tree->mod->s_opt->quickdirty,
								  Optwrap_Lk,NULL,tree,NULL);
					  }
				  }
			 }
		}
		if(tree->mod->whichrealmodel<=HLP17){ //added by Kenneth Hoehn
		  	 For(c,tree->mod->nhotness){
		  		if(tree->mod->hoptindex[c] == 2){
		  			Generic_Brent_Lk(&(tree->mod->hotness[c]),
		  					TREEBOUNDLOW,TREEBOUNDHIGH,
		  					tree->mod->s_opt->min_diff_lk_global,
		  					tree->mod->s_opt->brent_it_max,
		  					tree->mod->s_opt->quickdirty,
		  					Optwrap_Lk,NULL,tree,NULL);
		  		}
		  	 }
		 }
	  }
	}
	if(verbose){
	    For(i,io->ntrees){
	    	t_tree* tree = io->tree_s[i];
			if(tree->mod->optKappa==2){
		  			Print_Lk(tree,"[ts/tv ratio        ]");
		  			PhyML_Printf("[%.2f ]",tree->mod->kappa);
			}
			if(tree->mod->s_opt->opt_omega){
				For(c,tree->mod->nomega_part){
					if(tree->mod->omega_part_opt[c]==2){
						Print_Lk(tree,"[dn/ds ratio        ]");
						PhyML_Printf("[%.2f ]",tree->mod->omega_part[c]);
					}
				}
			}
			if(tree->mod->whichrealmodel<=HLP17){ //added by Kenneth Hoehn
				    int d;
				    For(c,tree->mod->nmotifs){
				    	if(tree->mod->hoptindex[tree->mod->motif_hotness[c]]==2){
				    	  char *info = malloc(22);
				    	  char motifh[10];
				    	  sprintf(motifh,"%d",tree->mod->motif_hotness[c]);
				    	  strcpy(info, "[");
				    	  strcat(info, tree->mod->motifs[c]);
				    	  strcat(info, " h");
				    	  strcat(info, motifh);
				    	  for(d=strlen(info);d<20;d++){
				    		  strcat(info, " ");
				    	  }
				    	  strcat(info,"]");
				    	  Print_Lk(tree,info);
				    	  PhyML_Printf("[%.3f]",tree->mod->hotness[tree->mod->motif_hotness[c]]);
				    	}
				    }
			}
	    }
	}
}




/*********************************************************/

//Global variable made threadprivate
//phydbl gemin=1e-6; //!< Added by Marcelo.

/*#if defined OMP || defined BLAS_OMP
#pragma omp threadprivate(gemin)
#endif*/

//Modified by Ken to update h
void Optimiz_All_Free_Param(option* io, int verbose, int recurse){

  if(io->mod->datatype==CODON){ //!< Added by Marcelo.
    char s[100],r[100];
    int  init_both_sides, numParams = 0, i, n;
    if(io->mod->s_opt->opt_method==optPAML){
      phydbl x2min[120], x2minbound[120][2], fx, intf,newf, *space; //!< 120 is a very pessimistic number of parameters adopted to simplify the selection for optimiazion.
      
      init_both_sides  = io->both_sides;
      io->both_sides = 0;

      //! Those models are updated only once ... need to clean up and make this more robust to consider the case where no parameters are optimized
      if((io->mod->whichmodel==GYECMS05  ||
		io->mod->whichmodel==YAPECMS05 ||
		io->mod->whichmodel==GYECMK07  ||
		io->mod->whichmodel==YAPECMK07 ||
		io->mod->whichmodel==GYECMS05F   ||
		io->mod->whichmodel==GYECMK07F   ||
		io->mod->whichmodel==YAPECMS05F  ||
		io->mod->whichmodel==YAPECMK07F ||
		io->mod->whichmodel==GYECMUSR ||
		io->mod->whichmodel==GYECMUSRF ||
		io->mod->whichmodel==MGECMUSRF ||
		io->mod->whichmodel==MGECMUSR ||
		io->mod->whichmodel==YAPECMUSR ||
		io->mod->whichmodel==YAPECMUSRF ) &&
		io->mod->s_opt->opt_state_freq==NO &&
		(io->mod->n_catg>1 && io->mod->n_w_catg==1)&&
		(io->mod->pcaModel==0)){
			io->mod->update_eigen = 0;
      }else{
		io->mod->update_eigen = 1;
      }
      
      if(io->mod->heuristicExpm){io->mod->expm=TAYLOR; io->mod->optParam=1;}
      if(io->mod->optDebug)printf("\nopt iter %d",io->mod->optIter);
      if(io->mod->s_opt->opt_kappa && (io->mod->optKappa==1||(io->mod->optKappa==2 && io->mod->optIter==0))){
		if(io->mod->whichmodel!=GYECMK07WK  && io->mod->whichmodel!=GYECMK07WKF &&
	  	io->mod->whichmodel!=GYECMS05WK  && io->mod->whichmodel!=GYECMS05WKF &&
	  	io->mod->whichmodel!=MGECMK07WK  && io->mod->whichmodel!=MGECMK07WKF &&
	  	io->mod->whichmodel!=MGECMS05WK  && io->mod->whichmodel!=MGECMS05WKF &&
	  	io->mod->whichmodel!=YAPECMK07WK && io->mod->whichmodel!=YAPECMK07WKF&&
	  	io->mod->whichmodel!=YAPECMS05WK && io->mod->whichmodel!=YAPECMS05WKF &&
	  	io->mod->whichmodel!=GYECMUSRWK && io->mod->whichmodel!=GYECMUSRWKF &&
	  	io->mod->whichmodel!=MGECMUSRWKF && io->mod->whichmodel!=YAPECMUSRWKF &&
	  	io->mod->whichmodel!=MGECMUSRWK && io->mod->whichmodel!=YAPECMUSRWK){
			if(io->mod->optDebug)printf("optimizing kappa\n");
	  		x2min[numParams++]               = io->mod->kappa;
	  		x2minbound[numParams-1][0]       = TREEBOUNDLOW;
	  		x2minbound[numParams-1][1]       = TREEBOUNDHIGH;
	  	}else{
	  		switch(io->mod->kappaECM){
	    		case kap1:
	    		  break;
	    		case kap2:
	    		case kap3:
	    		case kap6:{
	    		  x2min[numParams++]           = io->mod->pkappa[0];
	    		  x2minbound[numParams-1][0]   = TREEBOUNDLOW;
	    		  x2minbound[numParams-1][1]   = TREEBOUNDHIGH;
	    		  break;
	    		}
	    		case kap4:{
	    		  x2min[numParams++]           = io->mod->pkappa[0];
	    		  x2minbound[numParams-1][0]   = TREEBOUNDLOW;
	    		  x2minbound[numParams-1][1]   = TREEBOUNDHIGH;
	    		  x2min[numParams++]           = io->mod->pkappa[1];
	    		  x2minbound[numParams-1][0]   = TREEBOUNDLOW;
	    		  x2minbound[numParams-1][1]   = TREEBOUNDHIGH;
	    		  break;
	    		}
	    		case kap5:{
	    		  For(i,io->mod->nkappa-1){
					x2min[numParams++]         = io->mod->unspkappa[i];
					x2minbound[numParams-1][0] = -99.00;
					x2minbound[numParams-1][1] = 99.00;
	    		  }
	    		  break;
	    		}
	    		default: Warn_And_Exit("Error in Kappa assignment.");
	    		break;
	  		}
		}
    }
      
    if(io->mod->s_opt->opt_omega){
		if(io->mod->omegaSiteVar==DM0){
			int omegai; //added by Ken 18/8/2016
		 	for(omegai=0;omegai<io->mod->nomega_part;omegai++){
		    	if(io->mod->omega_part_opt[omegai]==1 || (io->mod->omega_part_opt[omegai]==2 && io->mod->optIter==0)){
		 			x2min[numParams++]           = io->mod->omega_part[omegai];
		    		x2minbound[numParams-1][0]   = TREEBOUNDLOW*100; //changed by Ken 9/2/2017 due to underflow issues with highly polymorphic lineages
		    		x2minbound[numParams-1][1]   = TREEBOUNDHIGH;
		    	}
		  	}
		}else if(io->mod->omegaSiteVar==DMODELK){
	  		For(i,io->mod->n_w_catg){
	  		  x2min[numParams++]         = io->mod->omegas[i];
	  		  x2minbound[numParams-1][0] = TREEBOUNDLOW;
	  		  x2minbound[numParams-1][1] = TREEBOUNDHIGH;
	  		}
	  		For(i,io->mod->n_w_catg-1){
	  		  x2min[numParams++]         = io->mod->prob_omegas_uns[i];
	  		  x2minbound[numParams-1][0] = -99.0;
	  		  x2minbound[numParams-1][1] = 99.0;
	  		}
		}else if(io->mod->omegaSiteVar==DGAMMAK){
	  		x2min[numParams++]           = io->mod->alpha;
	  		x2minbound[numParams-1][0]   = 2e-3;
	  		x2minbound[numParams-1][1]   = TREEBOUNDHIGH;
	  		x2min[numParams++]           = io->mod->beta;
	  		x2minbound[numParams-1][0]   = 2e-3;
	  		x2minbound[numParams-1][1]   = TREEBOUNDHIGH;
		}
    }
    if(io->mod->whichrealmodel<=HLP17){ //added by Kenneth Hoehn //changed from opthotness==1 on 31/1/2018
  		int c;
  	  	for(c=0;c<io->mod->nhotness;c++){
  		if(io->mod->hoptindex[c] == 1  || (io->mod->optIter==0 && io->mod->hoptindex[c] == 2)){
  		  x2min[numParams++]           = io->mod->hotness[c];
  	  	  x2minbound[numParams-1][0]   = -0.99; //new minimum value as of 21/Nov/2018.
  	  	  x2minbound[numParams-1][1]   = TREEBOUNDHIGH;
  	 	}
  	   }
     }
      
    if(io->mod->s_opt->opt_state_freq){
    	//printf("optimizing state freq\n");
		switch(io->mod->freq_model){
		  case F1XSENSECODONS:{
		    For(i,io->mod->num_base_freq-1){
		      x2min[numParams++]         = io->mod->pi_unscaled[i];
		      x2minbound[numParams-1][0] = -99.0;
		      x2minbound[numParams-1][1] = 99.0;
		    }
		    break;
		  }
		  case F1X4:{
		    For(i,io->mod->num_base_freq-1){
		      x2min[numParams++]         = io->mod->uns_base_freq[i];
		      x2minbound[numParams-1][0] = -99.0;
		      x2minbound[numParams-1][1] = 99.0;
		    }
		    break;
		  }
		  case F3X4:
		  case CF3X4:{
		    for(i=0;i<3;i++){
		      x2min[numParams++]         = io->mod->uns_base_freq[i];
		      x2minbound[numParams-1][0] = -99.0;
		      x2minbound[numParams-1][1] = 99.0;
		    }
		    for(i=4;i<7;i++){
		      x2min[numParams++]         = io->mod->uns_base_freq[i];
		      x2minbound[numParams-1][0] = -99.0;
		      x2minbound[numParams-1][1] = 99.0;
		    }
		    for(i=8;i<11;i++){
		      x2min[numParams++]         = io->mod->uns_base_freq[i];
		      x2minbound[numParams-1][0] = -99.0;
		      x2minbound[numParams-1][1] = 99.0;
		    }
		    break;
		  }
		  default:
		    break;
		}
    }
      
    if(io->mod->s_opt->opt_pinvar){
		x2min[numParams++]               = io->mod->pinvar;
		x2minbound[numParams-1][0]       = 0.0;
		x2minbound[numParams-1][1]       = 1.0;
    }
      
    if(io->mod->n_catg>1 && io->mod->n_w_catg==1 && io->mod->s_opt->opt_alphaCD){
		x2min[numParams++]           = io->mod->alpha;
		x2minbound[numParams-1][0]   = 2e-3;
		x2minbound[numParams-1][1]   = TREEBOUNDHIGH;
    }
      
    if(io->mod->pcaModel==1){
		For(i,io->mod->npcs){
		  x2min[numParams++]         = io->mod->pcsC[i];
		  x2minbound[numParams-1][0] = (-1) * TREEBOUNDHIGH;
		  x2minbound[numParams-1][1] = TREEBOUNDHIGH ;
		}
      }
      
    intf=fx=io->replnL;
      
    if(numParams>0){
    	storeParams(io,0,io->paramStore);
		space = (phydbl *) mCalloc((20+5*numParams)*numParams, sizeof(phydbl));
		if(io->mod->optDebug)printf("\ngemin in: %lf, numparams: %d",io->gemin,numParams);
		int result = BFGS_from_CODEML(&fx, io, x2min, x2minbound, space, io->gemin, numParams);
		if(io->mod->optDebug)printf("\n\nRESULT: %d\n\n",result);
		if(io->mod->optDebug)printf("\ngemin out: %lf",io->gemin);
		//Recursive function to save parameter estimates if they fail
		if(isnan(io->mod->omega_part[0])){
			restoreParams(io,0,io->paramStore);
			For(i,io->mod->nomega_part){
				phydbl r = (rand()*1.0)/RAND_MAX;
				io->mod->omega_part[i]+=r*0.5;
			}
			printf("\nOptimization failed - jiggling omega and trying again %lf\n",
					io->mod->omega_part[0]);
			if(recurse > 5)Warn_And_Exit("\\n\n\nCouldn't get BFGS to work - sorry :-(\n\n\n");
			Optimiz_All_Free_Param(io, verbose,++recurse);
		}
		newf=io->replnL;
		io->gemin/=2;  if(FABS(newf-intf)<1) io->gemin/=2;
		if(FABS(newf-intf)<0.5)     io->gemin = min2(io->gemin,1e-3);
		else if(FABS(newf-intf)>10) io->gemin = max2(io->gemin,0.1);
		io->gemin = max2(io->gemin,1e-6);
		if(io->mod->optDebug)printf("\ngemin mod: %lf %lf %lf",io->gemin,newf,intf);
		free(space);
    }
      
    io->both_sides = init_both_sides;
    if(io->mod->heuristicExpm){
		io->mod->update_eigen = 1;
		io->mod->optParam=0;
		io->mod->expm=EIGEN;
		Lk_rep(io);
		io->mod->update_eigen = 0;
    }else{
		io->mod->update_eigen = 0;
		if(io->both_sides){
			Lk_rep(io); /* Needed to update all partial likelihoods.*/
		}
      }
    }/*
    else
    {
      //!< use the old phyml method of parameter by parameter with brent.
      
      init_both_sides  = tree->both_sides;
      tree->both_sides = 0;
      For(n, tree->mod->s_opt->nBrentCycles)
      {
	if ((tree->mod->whichmodel==GYECMS05  || //! Those models are updated only once ... need to clean up and make this more robust to consider the case where no parameters are optimized
	  tree->mod->whichmodel==YAPECMS05 || 
	  tree->mod->whichmodel==GYECMK07  || 
	  tree->mod->whichmodel==YAPECMK07 ||
	  tree->mod->whichmodel==GYECMS05F   ||
	  tree->mod->whichmodel==GYECMK07F   ||
	  tree->mod->whichmodel==YAPECMS05F  ||
	  tree->mod->whichmodel==YAPECMK07F ||
	  tree->mod->whichmodel==GYECMUSR ||
	  tree->mod->whichmodel==GYECMUSRF ||
	  tree->mod->whichmodel==MGECMUSRF ||
	  tree->mod->whichmodel==MGECMUSR ||
	  tree->mod->whichmodel==YAPECMUSR ||
	  tree->mod->whichmodel==YAPECMUSRF ) && 
	  tree->mod->s_opt->opt_state_freq==NO &&
	  (tree->mod->n_catg>1 && tree->mod->n_w_catg==1)&&
	  (tree->mod->pcaModel==0))
	{
	  tree->mod->update_eigen = 0;
	}
	else
	{
	  tree->mod->update_eigen = 1;
	}
	
	if(tree->mod->heuristicExpm){ tree->mod->expm=TAYLOR; tree->mod->optParam=1; }
	
	if(tree->mod->s_opt->opt_kappa)
	{
	  if(tree->mod->whichmodel!=GYECMK07WK  && tree->mod->whichmodel!=GYECMK07WKF && 
	    tree->mod->whichmodel!=GYECMS05WK  && tree->mod->whichmodel!=GYECMS05WKF && 
	    tree->mod->whichmodel!=MGECMK07WK  && tree->mod->whichmodel!=MGECMK07WKF && 
	    tree->mod->whichmodel!=MGECMS05WK  && tree->mod->whichmodel!=MGECMS05WKF && 
	    tree->mod->whichmodel!=YAPECMK07WK && tree->mod->whichmodel!=YAPECMK07WKF&& 
	    tree->mod->whichmodel!=YAPECMS05WK && tree->mod->whichmodel!=YAPECMS05WKF &&
	    tree->mod->whichmodel!=GYECMUSRWK && tree->mod->whichmodel!=GYECMUSRWKF &&
	    tree->mod->whichmodel!=MGECMUSRWKF && tree->mod->whichmodel!=YAPECMUSRWKF &&
	    tree->mod->whichmodel!=MGECMUSRWK && tree->mod->whichmodel!=YAPECMUSRWK)
	  {
	    Generic_Brent_Lk(&(tree->mod->kappa),
			     TREEBOUNDLOW,TREEBOUNDHIGH,
			     tree->mod->s_opt->min_diff_lk_global,
			     tree->mod->s_opt->brent_it_max,
			     tree->mod->s_opt->quickdirty,
			     Optwrap_Lk,NULL,tree,NULL);
	  }
	  else
	  {
	    switch(tree->mod->kappaECM)
	    {
	      case kap1: 
		break;
	      case kap2:
	      case kap3: 
	      case kap6:  
	      {
		Generic_Brent_Lk(&(tree->mod->pkappa[0]),
				 TREEBOUNDLOW,TREEBOUNDHIGH,
				 tree->mod->s_opt->min_diff_lk_global,
				 tree->mod->s_opt->brent_it_max,
				 tree->mod->s_opt->quickdirty,
				 Optwrap_Lk,NULL,tree,NULL);
				 break;
	      }
	      case kap4:
	      {
		Generic_Brent_Lk(&(tree->mod->pkappa[0]),
				 TREEBOUNDLOW,TREEBOUNDHIGH,
				 tree->mod->s_opt->min_diff_lk_global,
				 tree->mod->s_opt->brent_it_max,
				 tree->mod->s_opt->quickdirty,
				 Optwrap_Lk,NULL,tree,NULL);
				 Generic_Brent_Lk(&(tree->mod->pkappa[1]),
						  TREEBOUNDLOW,TREEBOUNDHIGH,
						  tree->mod->s_opt->min_diff_lk_global,
						  tree->mod->s_opt->brent_it_max,
						  tree->mod->s_opt->quickdirty,
						  Optwrap_Lk,NULL,tree,NULL);
						  break;
	      }
	      case kap5:
	      {
		For(i,tree->mod->nkappa-1)
		{
		  Generic_Brent_Lk(&(tree->mod->unspkappa[i]),
				   -99.0,99.0,
				   tree->mod->s_opt->min_diff_lk_global,
				   tree->mod->s_opt->brent_it_max,
				   tree->mod->s_opt->quickdirty,
				   Optwrap_Lk,NULL,tree,NULL);
		}
		break;
	      }
	      default: Warn_And_Exit("Error in Kappa assignment.");
	      break;
	    }
	  }
	}
	
	if(tree->mod->s_opt->opt_omega)
	{
	  if(tree->mod->omegaSiteVar==DM0)
	  {
		  int omegai; //added by Ken 18/8/2016
		  for(omegai=0;omegai<tree->mod->nomega_part;omegai++){ //Ken 18/8

		  Generic_Brent_Lk(&(tree->mod->omega_part[omegai]),
			     TREEBOUNDLOW,TREEBOUNDHIGH,
			     tree->mod->s_opt->min_diff_lk_global,
			     tree->mod->s_opt->brent_it_max,
			     tree->mod->s_opt->quickdirty,
			     Optwrap_Lk,NULL,tree,NULL);
		  }
	  }
	  else if(tree->mod->omegaSiteVar==DMODELK)
	  {
	    For(i,tree->mod->n_w_catg)
	    {
	      Generic_Brent_Lk(&(tree->mod->omegas[i]),
			       TREEBOUNDLOW,TREEBOUNDHIGH,
			       tree->mod->s_opt->min_diff_lk_global,
			       tree->mod->s_opt->brent_it_max,
			       tree->mod->s_opt->quickdirty,
			       Optwrap_Lk,NULL,tree,NULL);
			       
	      if(i<tree->mod->n_w_catg-1)
	      {
		
		Generic_Brent_Lk(&(tree->mod->prob_omegas_uns[i]),
				-99.0,99.0,
				tree->mod->s_opt->min_diff_lk_global,
				tree->mod->s_opt->brent_it_max,
				tree->mod->s_opt->quickdirty,
				Optwrap_Lk,NULL,tree,NULL);
	      }
	    }
	  }
	  else if(tree->mod->omegaSiteVar==DGAMMAK)
	  {
	    Generic_Brent_Lk(&(tree->mod->alpha),
			     2e-3,TREEBOUNDHIGH,
			     tree->mod->s_opt->min_diff_lk_global,
			     tree->mod->s_opt->brent_it_max,
			     tree->mod->s_opt->quickdirty,
			     Optwrap_Lk,NULL,tree,NULL);
			     Generic_Brent_Lk(&(tree->mod->beta),
					      2e-3,999.0,
					      tree->mod->s_opt->min_diff_lk_global,
					      tree->mod->s_opt->brent_it_max,
					      tree->mod->s_opt->quickdirty,
					      Optwrap_Lk,NULL,tree,NULL);
	  }
	}
	printf("opt state freq: %d\n",tree->mod->s_opt->opt_state_freq);
	if(tree->mod->s_opt->opt_state_freq)
	{
	  switch(tree->mod->freq_model)
	  {
	    case F1XSENSECODONS:
	    {
	      For(i,tree->mod->num_base_freq-1)
	      {
		Generic_Brent_Lk(&(tree->mod->pi_unscaled[i]),
				 -99.0,99.0,
				 tree->mod->s_opt->min_diff_lk_global,
				 tree->mod->s_opt->brent_it_max,
				 tree->mod->s_opt->quickdirty,
				 Optwrap_Lk,NULL,tree,NULL);
	      }
	      break;
	    }
	    case F1X4:
	    {
	      For(i,tree->mod->num_base_freq-1)
	      {
		Generic_Brent_Lk(&(tree->mod->uns_base_freq[i]),
				 -99.0,99.0,
				 tree->mod->s_opt->min_diff_lk_global,
				 tree->mod->s_opt->brent_it_max,
				 tree->mod->s_opt->quickdirty,
				 Optwrap_Lk,NULL,tree,NULL);
	      }
	     PhyML_Printf("%lf %lf %lf %lf\n",tree->mod->uns_base_freq[0],tree->mod->uns_base_freq[1],tree->mod->uns_base_freq[2],tree->mod->uns_base_freq[3]);
	      break;
	    }
	    case F3X4:
	    case CF3X4:
	    {
	      for(i=0;i<3;i++)
	      {
	    	  printf("doing optimizations\n");

		Generic_Brent_Lk(&(tree->mod->uns_base_freq[i]),
				 -99.0,99.0,
				 tree->mod->s_opt->min_diff_lk_global,
				 tree->mod->s_opt->brent_it_max,
				 tree->mod->s_opt->quickdirty,
				 Optwrap_Lk,NULL,tree,NULL);
	      }
	      for(i=4;i<7;i++)
	      {
		Generic_Brent_Lk(&(tree->mod->uns_base_freq[i]),
				 -99.0,99.0,
				 tree->mod->s_opt->min_diff_lk_global,
				 tree->mod->s_opt->brent_it_max,
				 tree->mod->s_opt->quickdirty,
				 Optwrap_Lk,NULL,tree,NULL);
	      }
	      for(i=8;i<11;i++)
	      {
		Generic_Brent_Lk(&(tree->mod->uns_base_freq[i]),
				 -99.0,99.0,
				 tree->mod->s_opt->min_diff_lk_global,
				 tree->mod->s_opt->brent_it_max,
				 tree->mod->s_opt->quickdirty,
				 Optwrap_Lk,NULL,tree,NULL);
	      }
	      break;
	    }
	    default:
	      break;
	  }
	}
	
	if(tree->mod->s_opt->opt_pinvar)
	{
	  Generic_Brent_Lk(&(tree->mod->pinvar),
			   0.0,1.0,
			   tree->mod->s_opt->min_diff_lk_global,
			   tree->mod->s_opt->brent_it_max,
			   tree->mod->s_opt->quickdirty,
			   Optwrap_Lk,NULL,tree,NULL);
	}
	
	if(tree->mod->n_catg>1 && tree->mod->n_w_catg==1 && tree->mod->s_opt->opt_alphaCD)
	{ 
	  Generic_Brent_Lk(&(tree->mod->alpha),
			   2e-3,TREEBOUNDHIGH,
			   tree->mod->s_opt->min_diff_lk_global,
			   tree->mod->s_opt->brent_it_max,
			   tree->mod->s_opt->quickdirty,
			   Optwrap_Lk,NULL,tree,NULL);
	}
	
	if(tree->mod->pcaModel==1)
	{
            printf("hier\n");
	  For(i,tree->mod->npcs)
	  {
	    Generic_Brent_Lk(&(tree->mod->pcsC[i]),
			     ((-1) * TREEBOUNDHIGH),TREEBOUNDHIGH,
			     tree->mod->s_opt->min_diff_lk_global,
			     tree->mod->s_opt->brent_it_max,
			     tree->mod->s_opt->quickdirty,
			     Optwrap_Lk,NULL,tree,NULL);
	  }	
	}
	
	tree->both_sides = init_both_sides;
	
	if(tree->mod->heuristicExpm)
	{
	  tree->mod->update_eigen = 1; 
	  tree->mod->optParam=0;
	  tree->mod->expm=EIGEN;
	  Lk(tree);
	  tree->mod->update_eigen = 0;  
	}
	else
	{
	  tree->mod->update_eigen = 0;
	  if(tree->both_sides) Lk(tree); // Needed to update all partial likelihoods.
	}
      }
    }*/
    if(verbose){
      if(io->mod->s_opt->opt_kappa && (io->mod->optKappa==1||io->mod->optIter==0)){
		if(io->mod->whichmodel!=GYECMK07WK  && io->mod->whichmodel!=GYECMK07WKF &&
	  		io->mod->whichmodel!=GYECMS05WK  && io->mod->whichmodel!=GYECMS05WKF &&
	  		io->mod->whichmodel!=MGECMK07WK  && io->mod->whichmodel!=MGECMK07WKF &&
	  		io->mod->whichmodel!=MGECMS05WK  && io->mod->whichmodel!=MGECMS05WKF &&
	  		io->mod->whichmodel!=YAPECMK07WK && io->mod->whichmodel!=YAPECMK07WKF &&
	  		io->mod->whichmodel!=YAPECMS05WK && io->mod->whichmodel!=YAPECMS05WKF &&
	  		io->mod->whichmodel!=GYECMUSRWK && io->mod->whichmodel!=GYECMUSRWKF &&
	  		io->mod->whichmodel!=MGECMUSRWKF && io->mod->whichmodel!=YAPECMUSRWKF &&
	  		io->mod->whichmodel!=MGECMUSRWK && io->mod->whichmodel!=YAPECMUSRWK ){
	  			Print_Lk_rep(io,"[ts/tv ratio        ]");
	  			PhyML_Printf("[%.2f ]",io->mod->kappa);
		}else{
	  		switch(io->mod->kappaECM){
			    case kap1:
			      break;
			    case kap2:
			    case kap3:{
			      Print_Lk_rep(io,"[Emp. ts/tv ratio   ]");
			      PhyML_Printf("[%.2f ]",io->mod->pkappa[0]);
			      break;
			    }
			    case kap4:{
			      Print_Lk_rep(io,"[Emp. ts/tv ratio   ]");
			      PhyML_Printf("[ %.2f %.2f ]",io->mod->pkappa[0], io->mod->pkappa[1]);
			      break;
			    }
			    case kap5:{
			      Print_Lk_rep(io,"[Emp. ts/tv ratio   ]");
			      break;
			    }
			    case kap6:{
			      Print_Lk_rep(io,"[Multi-NT parameter ]");
			      PhyML_Printf("[%.2f ]",io->mod->pkappa[0]);
			      break;
			    }
			    default: Warn_And_Exit("Error in Kappa assignment.");
			    break;
			}
		}
      }
      if(io->mod->s_opt->opt_omega){
		if(io->mod->omegaSiteVar==DM0){
	  		if(io->mod->whichmodel!=GYECMK07WK  && io->mod->whichmodel!=GYECMK07WKF &&
	    		io->mod->whichmodel!=GYECMS05WK  && io->mod->whichmodel!=GYECMS05WKF &&
	    		io->mod->whichmodel!=MGECMK07WK  && io->mod->whichmodel!=MGECMK07WKF &&
	    		io->mod->whichmodel!=MGECMS05WK  && io->mod->whichmodel!=MGECMS05WKF &&
	    		io->mod->whichmodel!=YAPECMK07WK && io->mod->whichmodel!=YAPECMK07WKF&&
	    		io->mod->whichmodel!=YAPECMS05WK && io->mod->whichmodel!=YAPECMS05WKF &&
	    		io->mod->whichmodel!=GYECMUSRWK && io->mod->whichmodel!=GYECMUSRWKF &&
	    		io->mod->whichmodel!=MGECMUSRWKF && io->mod->whichmodel!=YAPECMUSRWKF &&
	    		io->mod->whichmodel!=MGECMUSRWK &&io->mod->whichmodel!=YAPECMUSRWK){
		 			int omegai; //added by Ken 18/8/2016
		  			for(omegai=0;omegai<io->mod->nomega_part;omegai++){
			  			if(io->mod->omega_part_opt[omegai]==1 || io->mod->optIter==0){
		  					Print_Lk_rep(io,"[dn/ds ratio        ]");
			  				PhyML_Printf("[%.2f ]",io->mod->omega_part[omegai]);
		  				}
		  			}
	  			}else{
	    			Print_Lk_rep(io,"[Emp. dn/ds ratio   ]");
    				if(io->mod->nomega_part > 1){printf("options not compatible with partitioned model error 5\n");exit(EXIT_FAILURE);}
	  			}
		}else if(io->mod->omegaSiteVar==DMODELK){
			  if(io->mod->n_w_catg<5){
			    if((io->mod->whichmodel!=MGECMUSRWK)&&(io->mod->whichmodel!=GYECMK07WK)&&(io->mod->whichmodel!=GYECMK07WKF)&&(io->mod->whichmodel!=GYECMS05WK)&&(io->mod->whichmodel!=GYECMS05WKF) && (io->mod->whichmodel!=MGECMK07WK)&&(io->mod->whichmodel!=MGECMK07WKF)&&(io->mod->whichmodel!=MGECMS05WK)&&(io->mod->whichmodel!=MGECMS05WKF) && (io->mod->whichmodel!=YAPECMK07WK)&&(io->mod->whichmodel!=YAPECMK07WKF)&&(io->mod->whichmodel!=YAPECMS05WK)&&(io->mod->whichmodel!=YAPECMS05WKF) &&(io->mod->whichmodel!=GYECMUSRWK && io->mod->whichmodel!=GYECMUSRWKF && io->mod->whichmodel!=MGECMUSRWKF && io->mod->whichmodel!=YAPECMUSRWKF && io->mod->whichmodel!=YAPECMUSRWK)){
			      int m;
			      Print_Lk_rep(io,"[dn/ds ratio        ]");
			      printf("[");
			      For(m,io->mod->n_w_catg) printf("%.2f ",io->mod->omegas[m]);
			      printf("]");
			      Print_Lk_rep(io,"[dn/ds frequencies  ]");
			      printf("[");
			      For(m,io->mod->n_w_catg) printf("%.2f ",io->mod->prob_omegas[m]);
			      printf("]");
		          fflush(NULL);
			    }else{
			      int m;
			      Print_Lk_rep(io,"[Emp. dn/ds ratio   ]");
			      printf("[");
			      if(io->mod->nomega_part > 1){printf("options not compatible with partitioned model error 6\n");exit(EXIT_FAILURE);}
			      printf("]");
			      Print_Lk_rep(io,"[dn/ds frequencies  ]");
			      printf("[");
			      For(m,io->mod->n_w_catg) printf("%.2f ",io->mod->prob_omegas[m]);
			      printf("]");
			    }
			  }else{
			    Print_Lk_rep(io,"[dn/ds ratio        ]");
			    Print_Lk_rep(io,"[dn/ds frequencies  ]");
			  }
		}
		else if(io->mod->omegaSiteVar==DGAMMAK){
		  if(io->mod->n_w_catg<5){
		    int m;
		    if((io->mod->whichmodel!=MGECMK07WK)&&(io->mod->whichmodel!=GYECMK07WK)&&(io->mod->whichmodel!=GYECMK07WKF)&&(io->mod->whichmodel!=GYECMS05WK)&&(io->mod->whichmodel!=GYECMS05WKF)&&(io->mod->whichmodel!=MGECMK07WK)&&(io->mod->whichmodel!=MGECMK07WKF)&&(io->mod->whichmodel!=MGECMS05WK)&&(io->mod->whichmodel!=MGECMS05WKF)&&(io->mod->whichmodel!=YAPECMK07WK)&&(io->mod->whichmodel!=YAPECMK07WKF)&&(io->mod->whichmodel!=YAPECMS05WK)&&(io->mod->whichmodel!=YAPECMS05WKF)&&(io->mod->whichmodel!=GYECMUSRWK && io->mod->whichmodel!=GYECMUSRWKF && io->mod->whichmodel!=MGECMUSRWKF && io->mod->whichmodel!=YAPECMUSRWKF && io->mod->whichmodel!=YAPECMUSRWK)){
		      Print_Lk_rep(io,"[dn/ds ratio        ]");
		      printf("[");
		      For(m,io->mod->n_w_catg) printf("%.2f ",io->mod->omegas[m]);
		      printf("]");
		      Print_Lk_rep(io,"[Alpha and Beta     ]");
		      PhyML_Printf("[%.2f %.2f ]",io->mod->alpha, io->mod->beta);
		    }else{
		      Print_Lk_rep(io,"[Emp. dn/ds ratio   ]");
		      printf("[");
	      	if(io->mod->nomega_part > 1){printf("options not compatible with partitioned model error 7\n");exit(EXIT_FAILURE);}
		      printf("]");
		      Print_Lk_rep(io,"[Alpha and Beta     ]");
		      PhyML_Printf("[%.2f %.2f ]",io->mod->alpha, io->mod->beta);
		    }
		  }else{
		    Print_Lk_rep(io,"[dn/ds ratio        ]");
		    Print_Lk_rep(io,"[Alpha and Beta     ]");
		    PhyML_Printf("[%.2f %.2f ]",io->mod->alpha, io->mod->beta);
		  }
		}
      }

      if(io->mod->opthotness==1){
    	 int c,d;
    	 for(c=0;c<io->mod->nmotifs;c++){
    		 if(io->mod->hoptindex[io->mod->motif_hotness[c]]==1 || io->mod->optIter==0){
    			 char *info = malloc(22);
    			 char motifh[10];
    			 sprintf(motifh,"%d",io->mod->motif_hotness[c]);
    			 strcpy(info, "[");
    			 strcat(info, io->mod->motifs[c]);
    			 strcat(info, " h");
    			 strcat(info, motifh);
    			 for(d=strlen(info);d<20;d++){
    				 strcat(info, " ");
    			 }
    			 strcat(info,"]");

    			 Print_Lk_rep(io,info);
    			 PhyML_Printf("[%.3f]",io->mod->hotness[io->mod->motif_hotness[c]]);
    		 }
    	 }
      }

      if(io->mod->s_opt->opt_state_freq){
		switch(io->mod->freq_model){
		  case F1XSENSECODONS:{
			    strcpy(s,"[F1x");
			    sprintf(r,"%d",io->mod->ns);
			    strcat(s,r);
			    strcat(s,"        freqs.]");
			    i=-1;
			    while(s[++i]!=']');
			    s[++i]=0;
			    Print_Lk_rep(io,s);
			    break;
		  }
		  case F1X4: Print_Lk_rep(io,"[F1x4 freqs.        ]"); break;
		  case F3X4: Print_Lk_rep(io,"[F3x4 freqs.        ]"); break;
		  case CF3X4: Print_Lk_rep(io,"[CF3x4 freqs.       ]"); break;
		  default:
	    break;
		}
      }
      
      if(io->mod->s_opt->opt_pinvar){
		Print_Lk_rep(io,"[P-inv              ]");
		PhyML_Printf("[%.2f]",io->mod->pinvar);
      }
      
      if(io->mod->n_catg>1 && io->mod->n_w_catg==1 && io->mod->s_opt->opt_alphaCD){
		Print_Lk_rep(io,"[Alpha              ]");
		PhyML_Printf("[%.2f]",io->mod->alpha);
      }
      if(io->mod->pcaModel==1){
		Print_Lk_rep(io,"[PCA linear coeff.  ]");
      }
    }          
  } ///////!<End of CODON MODELS///////////

  //Optimize remaining submodel parameters!
  Optimiz_Submodel_Params(io,verbose);
  Lk_rep(io);
  io->mod->optIter++;
}
  
//Global variable made private in OMP
static phydbl sqrarg;

#if defined OMP || defined BLAS_OMP
#pragma omp threadprivate(sqrarg)
#endif

/*********************************************************/

void BFGS(t_tree *tree, 
	  phydbl *p, 
	  int n, 
	  phydbl gtol, 
	  phydbl step_size,
	  phydbl(*func)(t_tree *tree), 
	  int(*dfunc)(t_tree *tree,phydbl *param,int n_param,phydbl stepsize,phydbl(*func)(t_tree *tree),phydbl *derivatives), 
	  int(*lnsrch)(t_tree *tree, int n, phydbl *xold, phydbl fold,phydbl *g, phydbl *p, phydbl *x,phydbl *f, phydbl stpmax, int *check),
	  int *failed)
{
	Warn_And_Exit("BFGS\n");
  int check,i,its,j;
  phydbl den,fac,fad,fae,fp,stpmax,sum=0.0,sumdg,sumxi,temp,test,fret;
  phydbl *dg,*g,*hdg,**hessin,*pnew,*xi;
  
  hessin = (phydbl **)mCalloc(n,sizeof(phydbl *));
  For(i,n) hessin[i] = (phydbl *)mCalloc(n,sizeof(phydbl));
  dg   = (phydbl *)mCalloc(n,sizeof(phydbl ));
  g    = (phydbl *)mCalloc(n,sizeof(phydbl ));
  pnew = (phydbl *)mCalloc(n,sizeof(phydbl ));
  hdg  = (phydbl *)mCalloc(n,sizeof(phydbl ));
  xi   = (phydbl *)mCalloc(n,sizeof(phydbl ));
  
  fp=(*func)(tree);
  (*dfunc)(tree,p,n,step_size,func,g);

  for (i=0;i<n;i++){
      for (j=0;j<n;j++) hessin[i][j]=0.0;
      hessin[i][i]=1.0;
      xi[i] = -g[i];
      sum += p[i]*p[i];
  }

  stpmax=STPMX*MAX(SQRT(sum),(phydbl)n);

  for(its=1;its<=ITMAX;its++){
      lnsrch(tree,n,p,fp,g,xi,pnew,&fret,stpmax,&check);
      fp = fret;
      for (i=0;i<n;i++){
    	  xi[i]=pnew[i]-p[i];
	  	  p[i]=pnew[i];
      }
      test=0.0;
      for (i=0;i<n;i++){
    	  temp=FABS(xi[i])/MAX(FABS(p[i]),1.0);
    	  if(temp > test) test=temp;
      }
      if (test < TOLX){
    	  (*func)(tree);
	  	  For(i,n) free(hessin[i]);
	  	  free(hessin);
	  	  free(xi);
	  	  free(pnew);
	  	  free(hdg);
	  	  free(g);
	  	  free(dg);

	  	  if(its == 1)*failed = 1;
	  	 return;
      }

      for (i=0;i<n;i++) dg[i]=g[i];
      (*dfunc)(tree,p,n,step_size,func,g);
      test=0.0;
      den=MAX(fret,1.0);
      for (i=0;i<n;i++){
    	  temp=FABS(g[i])*MAX(FABS(p[i]),1.0)/den;
    	  if (temp > test) test=temp;
      }
      if (test < gtol){
    	  (*func)(tree);
	  	  For(i,n) free(hessin[i]);
	  	  free(hessin);
	  	  free(xi);
	  	  free(pnew);
	  	  free(hdg);
	  	  free(g);
	  	  free(dg);
	  	  return;
      }

      for (i=0;i<n;i++) dg[i]=g[i]-dg[i];

      for (i=0;i<n;i++){
    	  hdg[i]=0.0;
    	  for (j=0;j<n;j++) hdg[i] += hessin[i][j]*dg[j];
      }

      fac=fae=sumdg=sumxi=0.0;
      for (i=0;i<n;i++){
    	  fac += dg[i]*xi[i];
    	  fae += dg[i]*hdg[i];
    	  sumdg += SQR(dg[i]);
    	  sumxi += SQR(xi[i]);
      }
    
      if(fac*fac > EPS*sumdg*sumxi){
    	  fac=1.0/fac;
    	  fad=1.0/fae;
    	  for (i=0;i<n;i++) dg[i]=fac*xi[i]-fad*hdg[i];
    	  for (i=0;i<n;i++){
    		  for (j=0;j<n;j++){
    			  hessin[i][j] += fac*xi[i]*xi[j]
					-fad*hdg[i]*hdg[j]+fae*dg[i]*dg[j];
    		  }
    	  }
      }
      for (i=0;i<n;i++){
    	  xi[i]=0.0;
    	  for(j=0;j<n;j++) xi[i] -= hessin[i][j]*g[j];
      }
    }
  *failed = 1;
  For(i,n) free(hessin[i]);
  free(hessin);
  free(xi);
  free(pnew);
  free(hdg);
  free(g);
  free(dg);   
}

#undef ITMAX
#undef EPS
#undef TOLX
#undef STPMX


#define ALF 1.0e-4
#define TOLX 1.0e-7
#undef ALF
#undef TOLX
#undef NRANSI

#define ALF 1.0e-4
#define TOLX 1.0e-7

#undef ALF
#undef TOLX
#undef NRANSI

/*********************************************************/

phydbl Dist_F_Brent(phydbl ax, phydbl bx, phydbl cx, phydbl tol, int n_iter_max, 
		    phydbl *param, phydbl *F, model *mod){
  int iter;
  phydbl a,b,d,etemp,fu,fv,fw,fx,p,q,r,tol1,tol2,u,v,w,x,xm;
  phydbl e=0.0;
  phydbl old_lnL,init_lnL, curr_lnL;

  d=0.0;
  a=((ax < cx) ? ax : cx);
  b=((ax > cx) ? ax : cx);
  x = w = v = bx;
  old_lnL = UNLIKELY;
  fw = fv = fx = -Lk_Dist(F,FABS(bx),mod);
  curr_lnL = init_lnL = -fw;

  for(iter=1;iter<=BRENT_ITMAX;iter++){
      xm=0.5*(a+b);

      tol2=2.0*(tol1=tol*FABS(x)+BRENT_ZEPS);

      if(((FABS(curr_lnL-old_lnL) < mod->s_opt->min_diff_lk_local) &&
	  (curr_lnL > init_lnL - mod->s_opt->min_diff_lk_local)) ||	 
	  (iter > n_iter_max - 1)){
    	  *param = x;
    	  curr_lnL = Lk_Dist(F,*param,mod);
    	  return -curr_lnL;
      }
      
      if(FABS(e) > tol1){
    	  r=(x-w)*(fx-fv);
    	  q=(x-v)*(fx-fw);
    	  p=(x-v)*q-(x-w)*r;
    	  q=2.0*(q-r);
    	  if(q > 0.0) p = -p;
    	  q=FABS(q);
    	  etemp=e;
    	  e=d;
    	  if(FABS(p) >= FABS(0.5*q*etemp) || p <= q*(a-x) || p >= q*(b-x)){
    		  d=BRENT_CGOLD*(e=(x >= xm ? a-x : b-x));
    	  }else{
    		  d=p/q;
    		  u=x+d;
    		  if (u-a < tol2 || b-u < tol2)
    		  d=SIGN(tol1,xm-x);
    	  }
      }else{
    	  d=BRENT_CGOLD*(e=(x >= xm ? a-x : b-x));
      }
      u=(FABS(d) >= tol1 ? x+d : x+SIGN(tol1,d));
      if(u<BL_MIN) u = BL_MIN;
      (*param) = FABS(u);
      old_lnL = curr_lnL;
      fu = -Lk_Dist(F,FABS(u),mod);
      curr_lnL = -fu;      
      if(fu < fx){
    	  if(iter > n_iter_max) return -fu;
    	  if(u >= x) a=x; else b=x;
    	  SHFT(v,w,x,u)
    	  SHFT(fv,fw,fx,fu)
      }else{
    	  if (u < x) a=u; else b=u;
    	  if (fu < fw || FABS(w-x) < SMALL){
    		  v=w;
    		  w=u;
    		  fv=fw;
    		  fw=fu;
	    }else if (fu < fv || FABS(v-x) < SMALL
	    		|| FABS(v-w) < SMALL){
	      v=u;
	      fv=fu;
	    }
	}
  }
  Exit("\n. Too many iterations in BRENT !");
  return(-1);
}

/*********************************************************/

void Opt_Dist_F(phydbl *dist, phydbl *F, model *mod)
{
  phydbl ax,bx,cx;

  if(*dist < BL_MIN) *dist = BL_MIN;

  ax = BL_MIN;
  bx =  (*dist);
  cx = BL_MAX;

/*   Dist_F_Brak(&ax,&bx,&cx,F,dist,mod); */
  Dist_F_Brent(ax,bx,cx,1.E-10,1000,dist,F,mod);
}

/*********************************************************/

phydbl Optwrap_Lk(t_edge *b, t_tree *tree, supert_tree *stree)
{
  Lk(tree);
  return tree->c_lnL;
}

/*********************************************************/

phydbl Optwrap_Lk_At_Given_Edge(t_edge *b, t_tree *tree, supert_tree *stree)
{
  Lk_At_Given_Edge(b,tree);
  return tree->c_lnL;
}

/*********************************************************/

phydbl Generic_Brent_Lk(phydbl *param, phydbl ax, phydbl cx, phydbl tol, 
			int n_iter_max, int quickdirty,
			phydbl (*obj_func)(t_edge *,t_tree *,supert_tree *), 
			t_edge *branch, t_tree *tree, supert_tree *stree){
  int iter;
  phydbl a,b,d,etemp,fu,fv,fw,fx,p,q,r,tol1,tol2,u,v,w,x,xm;
  phydbl e=0.0;
  phydbl old_lnL,init_lnL,cur_lnL;
  phydbl bx = *param;

  d=0.0;
  a=((ax < cx) ? ax : cx);
  b=((ax > cx) ? ax : cx);
  x=w=v=bx;
  old_lnL = UNLIKELY;
  (*param) = FABS(bx);
  fw=fv=fx=fu=-(*obj_func)(branch,tree,stree);
  init_lnL = -fw;

  for(iter=1;iter<=BRENT_ITMAX;iter++){
      xm=0.5*(a+b);
      tol2=2.0*(tol1=tol*FABS(x)+BRENT_ZEPS);

      cur_lnL = (stree)?(stree->tree->c_lnL):(tree->c_lnL);

      if((cur_lnL > init_lnL + tol) && (quickdirty)){
    	  (*param) = x;
    	  (*obj_func)(branch,tree,stree);
    	  return (stree)?(stree->tree->c_lnL):(tree->c_lnL);
      }
      if((FABS(cur_lnL-old_lnL) < tol) || (iter > n_iter_max - 1)){
    	  (*param) = x;
    	  (*obj_func)(branch,tree,stree);
    	  return (stree)?(stree->tree->c_lnL):(tree->c_lnL);
      }
      if(FABS(e) > tol1){
    	  r=(x-w)*(fx-fv);
    	  q=(x-v)*(fx-fw);
    	  p=(x-v)*q-(x-w)*r;
    	  q=2.0*(q-r);
    	  if(q > 0.0) p = -p;
    	  q=FABS(q);
    	  etemp=e;
    	  e=d;
    	  if(FABS(p) >= FABS(0.5*q*etemp) || p <= q*(a-x) || p >= q*(b-x)){
    		  d=BRENT_CGOLD*(e=(x >= xm ? a-x : b-x));
    	 }else{
	      d=p/q;
	      u=x+d;
	      if (u-a < tol2 || b-u < tol2) d=SIGN(tol1,xm-x);
	    }
     }else{
	  d=BRENT_CGOLD*(e=(x >= xm ? a-x : b-x));
     }
      u=(FABS(d) >= tol1 ? x+d : x+SIGN(tol1,d));
      (*param) = FABS(u);
      old_lnL = (stree)?(stree->tree->c_lnL):(tree->c_lnL);
      fu = -(*obj_func)(branch,tree,stree);

      if(fu <= fx){
    	  if(u >= x) a=x; else b=x;
    	  SHFT(v,w,x,u)
    	  SHFT(fv,fw,fx,fu)
      }else{
    	  if (u < x) a=u; else b=u;
    	  if (fu < fw || FABS(w-x) < SMALL){
    		  v=w;
    		  w=u;
    		  fv=fw;
    		  fw=fu;
    	  }
	  else if (fu < fv || FABS(v-x) < SMALL || FABS(v-w) < SMALL){
	      v=u;
	      fv=fu;
	  }
	}
  }
  Exit("\n. Too many iterations in BRENT !");
  return(-1);
}

/*********************************************************/
phydbl Br_Len_Brent_Codon_Pairwise(phydbl ax, phydbl bx, phydbl cx, phydbl tol, phydbl *b_fcus,
		phydbl *Pij, phydbl *pi, eigen *eigenStruct, calign *data,
		int ns, int n_iter_max, int quickdirty, phydbl *uexpt, phydbl *expt){
  int iter;
  phydbl a,b,d,etemp,fu,fv,fw,fx,p,q,r,tol1,tol2,u,v,w,x,xm;
  phydbl e=0.0;
  phydbl old_lnL, init_lnL, curr_lnL;

  d=0.0;
  a=((ax < cx) ? ax : cx);
  b=((ax > cx) ? ax : cx);
  x=w=v=bx;
  old_lnL = UNLIKELY;
  (*b_fcus) = FABS(bx);
  fw=fv=fx=fu=-LK_Codon_Pairwise(data, Pij, pi, ns, (*b_fcus), eigenStruct, uexpt, expt);
  curr_lnL = init_lnL = -fw;

  for(iter=1;iter<=BRENT_ITMAX;iter++){
      xm=0.5*(a+b);
      tol2=2.0*(tol1=tol*FABS(x)+BRENT_ZEPS);

      if((curr_lnL > init_lnL + tol) && (quickdirty)){
    	  (*b_fcus) = x;
	  	  curr_lnL=LK_Codon_Pairwise(data, Pij, pi, ns, (*b_fcus), eigenStruct, uexpt, expt);
	  	  return curr_lnL;
      }

      if((FABS(curr_lnL-old_lnL) < tol) || (iter > n_iter_max - 1)){
    	  (*b_fcus)=x;
    	  curr_lnL=LK_Codon_Pairwise(data, Pij, pi, ns, (*b_fcus), eigenStruct, uexpt, expt);
    	  return curr_lnL;
      }
      
      if(FABS(e) > tol1){
    	  r=(x-w)*(fx-fv);
    	  q=(x-v)*(fx-fw);
    	  p=(x-v)*q-(x-w)*r;
    	  q=2.0*(q-r);
    	  if(q > 0.0) p = -p;
    	  q=FABS(q);
    	  etemp=e;
    	  e=d;
    	  if(FABS(p) >= FABS(0.5*q*etemp) || p <= q*(a-x) || p >= q*(b-x))
    		  d=BRENT_CGOLD*(e=(x >= xm ? a-x : b-x));
    	  else{
    		  d=p/q;
    		  u=x+d;
    		  if (u-a < tol2 || b-u < tol2)
    			  d=SIGN(tol1,xm-x);
    	  }
      }else{
    	  d=BRENT_CGOLD*(e=(x >= xm ? a-x : b-x));
      }
      u=(FABS(d) >= tol1 ? x+d : x+SIGN(tol1,d));
      (*b_fcus)=FABS(u);
      old_lnL = curr_lnL;
      fu=-LK_Codon_Pairwise(data, Pij, pi, ns, (*b_fcus), eigenStruct, uexpt, expt);
      curr_lnL=-fu;
      if(fu < fx){
    	  if(u > x) a=x; else b=x;
    	  SHFT(v,w,x,u)
    	  SHFT(fv,fw,fx,fu)
      }else{
    	  if (u < x) a=u; else b=u;
    	  if (fu < fw || FABS(w-x) < SMALL){
    		  v=w;
    		  w=u;
    		  fv=fw;
    		  fw=fu;
    	  }else if (fu < fv || FABS(v-x) < SMALL || FABS(v-w) < SMALL){
    		  v=u;
    		  fv=fu;
    	  }
      }
  }
  if(iter > BRENT_ITMAX) PhyML_Printf("\n. Too many iterations in BRENT (%d) (%f)",iter,(*b_fcus));
  return(-1);
}

/*********************************************************/
#define BFGS

int BFGS_from_CODEML (phydbl *f, option* io, phydbl *x, phydbl xb[120][2], phydbl space[], phydbl e, int n)
{
/* n-variate minimization with bounds using the BFGS algorithm
     g0[n] g[n] p[n] x0[n] y[n] s[n] z[n] H[n*n] C[n*n] tv[2*n]
     xmark[n],ix[n]
   Size of space should be (check carefully?)
      #define spaceming2(n) ((n)*((n)*2+9+2)*sizeof(phydbl))
   nfree: # free variables
   xmark[i]=0 for inside space; -1 for lower boundary; 1 for upper boundary.
   x[] has initial values at input and returns the estimates in return.
   ix[i] specifies the i-th free parameter
  
   ALL CREDIT GOES TO PROF YANG and CODEML

   Ken: returns -100 if optimization fails due to parameter issues
*/
   int i,j, i1,i2,it, maxround=10000, fail=0, *xmark, *ix, nfree;
   int Ngoodtimes=2, goodtimes=0;
   phydbl small=1.e-30, sizep0=0;     /* small value for checking |w|=0 */
   phydbl f0, *g0, *g, *p, *x0, *y, *s, *z, *H, *C, *tv;
   phydbl w,v, alpha, am, h, maxstep=8;

   if(n==0) return(0);
   g0=space;   g=g0+n;  p=g+n;   x0=p+n;
   y=x0+n;     s=y+n;   z=s+n;   H=z+n;  C=H+n*n, tv=C+n*n;
   xmark=(int*)(tv+2*n);  ix=xmark+n;

   for(i=0; i<n; i++)  { xmark[i]=0; ix[i]=i; }
   
   for(i=0,nfree=0;i<n;i++) {
      if(x[i]<=xb[i][0]) { x[i]=xb[i][0]; xmark[i]=-1; continue; }
      if(x[i]>=xb[i][1]) { x[i]=xb[i][1]; xmark[i]= 1; continue; }
      ix[nfree++]=i;
   }
   f0=*f=LK_BFGS_from_CODEML(io,x,n);
   if(f0 != f0)return -100; //failure!

   xtoy(x,x0,n);
   io->SIZEp=99;

   gradientB (n, x0, f0, g0, io, tv, xmark);

   identity (H,nfree);
   for(io->Iround=0; io->Iround<maxround; io->Iround++) {
     

      for (i=0,zero(p,n); i<nfree; i++)  For (j,nfree)
         p[ix[i]] -= H[i*nfree+j]*g0[ix[j]];
      sizep0 = io->SIZEp;
      io->SIZEp  = norm(p,n);      /* check this */

      for (i=0,am=maxstep; i<n; i++) {  /* max step length */
         if (p[i]>0 && (xb[i][1]-x0[i])/p[i]<am) am=(xb[i][1]-x0[i])/p[i];
         else if (p[i]<0 && (xb[i][0]-x0[i])/p[i]<am) am=(xb[i][0]-x0[i])/p[i];
      }

      if (io->Iround==0) {
         h=fabs(2*f0*.01/innerp(g0,p,n));  /* check this?? */
         h=min2(h,am/2000);

      }
      else {
         h=norm(s,nfree)/io->SIZEp;
         h=max2(h,am/500);
      }
      h = max2(h,1e-5);   h = min2(h,am/5);
      *f = f0;
      alpha = LineSearch2(io,f,x0,p,h,am, min2(1e-3,e), tv,n); /* n or nfree? */
      if(alpha != alpha)return -100; //failure!

      if (alpha<=0) {
         if (fail) {
            if (io->AlwaysCenter) { io->Iround=maxround;  break; }
            else { io->AlwaysCenter=1; identity(H,n); fail=1; }
         }
         else   
            { if(io->noisy>2) printf(".. ");  identity(H,nfree); fail=1; }
      }
      else  {
         fail=0;
         For(i,n)  x[i]=x0[i]+alpha*p[i];
         w=min2(2,e*1000); if(e<1e-4 && e>1e-6) w=0.01;

         if(io->Iround==0 || io->SIZEp<sizep0 || (io->SIZEp<.001 && sizep0<.001)) goodtimes++;
         else  goodtimes=0;
         if((n==1||goodtimes>=Ngoodtimes) && io->SIZEp<(e>1e-5?1:.001)
            && H_end(x0,x,f0,*f,e,e,n))
            break;
      }
     
      gradientB (n, x, *f, g, io, tv, xmark);

      /* modify the working set */
      for(i=0; i<n; i++) {         /* add constraints, reduce H */
         if (xmark[i]) continue;
         if (fabs(x[i]-xb[i][0])<1e-6 && -g[i]<0)  xmark[i]=-1;
         else if (fabs(x[i]-xb[i][1])<1e-6 && -g[i]>0)  xmark[i]=1;
         if (xmark[i]==0) continue;
         xtoy (H, C, nfree*nfree);
         for(it=0; it<nfree; it++) if (ix[it]==i) break;
         for (i1=it; i1<nfree-1; i1++) ix[i1]=ix[i1+1];
         for (i1=0,nfree--; i1<nfree; i1++) For (i2,nfree)
            H[i1*nfree+i2]=C[(i1+(i1>=it))*(nfree+1) + i2+(i2>=it)];
      }
      for (i=0,it=0,w=0; i<n; i++) {  /* delete a constraint, enlarge H */
         if (xmark[i]==-1 && -g[i]>w)     { it=i; w=-g[i]; }
         else if (xmark[i]==1 && -g[i]<-w) { it=i; w=g[i]; }
      }
      if (w>10*io->SIZEp/nfree) {          /* *** */
         xtoy (H, C, nfree*nfree);
         For (i1,nfree) For (i2,nfree) H[i1*(nfree+1)+i2]=C[i1*nfree+i2];
         For (i1,nfree+1) H[i1*(nfree+1)+nfree]=H[nfree*(nfree+1)+i1]=0;
         H[(nfree+1)*(nfree+1)-1]=1;
         xmark[it]=0;   ix[nfree++]=it;
      }

     
      for (i=0,f0=*f; i<nfree; i++)
        {  y[i]=g[ix[i]]-g0[ix[i]];  s[i]=x[ix[i]]-x0[ix[i]]; }
      For (i,n) { g0[i]=g[i]; x0[i]=x[i]; }


      /* renewal of H varies with different algorithms   */
#if (defined SR1)
      /*   Symmetrical Rank One (Broyden, C. G., 1967) */
      for (i=0,w=.0; i<nfree; i++) {
         for (j=0,v=.0; j<nfree; j++) v += H[i*nfree+j] * y[j];
         z[i]=s[i] - v;
         w += y[i]*z[i];
      }
      if (fabs(w)<small)   { identity(H,nfree); fail=1; continue; }
      For (i,nfree)  For (j,nfree)  H[i*nfree+j] += z[i]*z[j]/w;
#elif (defined DFP)
      /* Davidon (1959), Fletcher and Powell (1963). */
      for (i=0,w=v=0.; i<nfree; i++) {
         for (j=0,z[i]=0; j<nfree; j++) z[i] += H[i*nfree+j] * y[j];
         w += y[i]*z[i];  v += y[i]*s[i];
      }
      if (fabs(w)<small || fabs(v)<small)  { identity(H,nfree); fail=1; continue;}
      For (i,nfree)  For (j,nfree)  
         H[i*nfree+j] += s[i]*s[j]/v - z[i]*z[j]/w;
#else /* BFGS */
      for (i=0,w=v=0.; i<nfree; i++) {
         for (j=0,z[i]=0.; j<nfree; j++) z[i]+=H[i*nfree+j]*y[j];
         w+=y[i]*z[i];    v+=y[i]*s[i];
      }
      if (fabs(v)<small)   { identity(H,nfree); fail=1; continue; }
      For (i,nfree)  For (j,nfree)
         H[i*nfree+j] += ((1+w/v)*s[i]*s[j]-z[i]*s[j]-s[i]*z[j])/v;
#endif

   }    /* for (Iround,maxround)  */

   /* try to remove this after updating LineSearch2() */
   *f=LK_BFGS_from_CODEML(io,x,n);
   if(*f != *f)return -100; //failure!

   if (io->Iround==maxround) {
      return(-1);
   }
   if(nfree==n) { 
      xtoy(H, space, n*n);  /* H has variance matrix, or inverse of Hessian */
      return(1);
   }
   return(0);
}
//Original declarations for AlwaysCenter and Small_Diff
//Global variables made threadprivate for OMP


/*********************************************************/

phydbl gradientB (int n, phydbl x[], phydbl f0, phydbl g[],
    option *io, phydbl space[], int xmark[])
{
/* f0=fun(x) is always provided.
   xmark=0: central; 1: upper; -1: down
*/
   int i,j;
   phydbl *x0=space, *x1=space+n, eh0=io->Small_Diff, eh;  /* eh0=1e-6 || 1e-7 */

   For(i,n) {
      eh=eh0*(fabs(x[i])+1);
      if (xmark[i]==0 && (io->AlwaysCenter || io->SIZEp<1)) {   /* central */
         For (j, n)  x0[j]=x1[j]=x[j];
         eh=pow(eh,.67);  x0[i]-=eh; x1[i]+=eh;
         g[i]=(LK_BFGS_from_CODEML(io,x1,n)-LK_BFGS_from_CODEML(io,x0,n))/(eh*2.0);
         if(g[i] != g[i])return NAN; //failure!
      }
      else  {/* forward or backward */
         For (j, n)  x1[j]=x[j];
         if (xmark[i]) eh*=-xmark[i];
         x1[i] += eh;
         g[i]=(LK_BFGS_from_CODEML(io,x1,n)-f0)/eh;
         if(g[i] != g[i])return NAN; //failure!
      }
   }
   return(0);
}

/*********************************************************/

phydbl innerp (phydbl x[], phydbl y[], int n)
{ int i; phydbl t=0;  for(i=0; i<n; i++)  t += x[i]*y[i];  return(t); }

/*********************************************************/

int xtoy (phydbl x[], phydbl y[], int n)
{ int i; for (i=0; i<n; i++) {y[i]=x[i];} return(0); }

/*********************************************************/

phydbl LineSearch2 (option* io, phydbl *f, phydbl x0[],
       phydbl p[], phydbl step, phydbl limit, phydbl e, phydbl space[], int n)
{
/* linear search using quadratic interpolation 
   from x0[] in the direction of p[],
                x = x0 + a*p        a ~(0,limit)
   returns (a).    *f: f(x0) for input and f(x) for output

   x0[n] x[n] p[n] space[n]

   adapted from Wolfe M. A.  1978.  Numerical methods for unconstrained
   optimization: An introduction.  Van Nostrand Reinhold Company, New York.
   pp. 62-73.
   step is used to find the bracket and is increased or reduced as necessary, 
   and is not terribly important.
*/
   int ii=0, maxround=10, status;
   // int i, nsymb=0;
   phydbl *x=space, factor=4, small=1e-10, smallgapa=0.2;
   phydbl a0,a1,a2,a3,a4=-1,a5,a6, f0,f1,f2,f3,f4=-1,f5,f6;

/* look for bracket (a1, a2, a3) with function values (f1, f2, f3)
   step length step given, and only in the direction a>=0 */
   if (step<=0 || limit<small || step>=limit) {
      return (0);
   }
   a0=a1=0; f1=f0=*f;
   a2=a0+step; f2=fun_LineSearch(a2, io,x0,p,x,n);
   if(a2 != a2)return NAN;  //failure!
   if (f2>f1) {  /* reduce step length so the algorithm is decreasing */
      for (; ;) {
         step/=factor;
         if (step<small) return (0);
         a3=a2;    f3=f2;
         a2=a0+step;  f2=fun_LineSearch(a2, io,x0,p,x,n);
         if(f2 != f2)return NAN;  //failure!
         if (f2<=f1) break;
        
      }
   }
   else {       /* step length is too small? */
      for (; ;) {
         step*=factor;
         if (step>limit) step=limit;
         a3=a0+step;  f3=fun_LineSearch(a3, io,x0,p,x,n);
         if(f3 != f3)return NAN;  //failure!
         if (f3>=f2) break;

        
         a1=a2; f1=f2;    a2=a3; f2=f3;
         if (step>=limit) {
            
            *f=f3; return(a3);
         }
      }
   }

   /* iteration by quadratic interpolation, fig 2.2.9-10 (pp 71-71) */
   for (ii=0; ii<maxround; ii++) {
      /* a4 is the minimum from the parabola over (a1,a2,a3)  */
      a4 = (a2-a3)*f1+(a3-a1)*f2+(a1-a2)*f3;
      if(fabs(a4)>1e-100) 
         a4 = ((a2*a2-a3*a3)*f1+(a3*a3-a1*a1)*f2+(a1*a1-a2*a2)*f3)/(2*a4);
      if (a4>a3 || a4<a1) {   /* out of range */
         a4=(a1+a2)/2;
         status='N';
      }
      else {
         if((a4<=a2 && a2-a4>smallgapa*(a2-a1)) || (a4>a2 && a4-a2>smallgapa*(a3-a2)))
            status='Y';
         else 
            status='C';
      }
      f4 = fun_LineSearch(a4, io,x0,p,x,n);
      if(f4 != f4)return NAN;  //failure!
      
      if (fabs(f2-f4)<e*(1+fabs(f2))) {
         
         break;
      }

      /* possible multiple local optima during line search */
      if (a4<=a2) {    /* fig 2.2.10 */
         if (a2-a4>smallgapa*(a2-a1)) {
            if (f4<=f2) { a3=a2; a2=a4;  f3=f2; f2=f4; }
            else        { a1=a4; f1=f4; }
         }
         else {
            if (f4>f2) {
               a5=(a2+a3)/2; f5=fun_LineSearch(a5, io,x0,p,x,n);
               if(f5 != f5)return NAN;  //failure!
               if (f5>f2) { a1=a4; a3=a5;  f1=f4; f3=f5; }
               else       { a1=a2; a2=a5;  f1=f2; f2=f5; }
            }
            else {
               a5=(a1+a4)/2; f5=fun_LineSearch(a5, io,x0,p,x,n);
               if(f5 != f5)return NAN;  //failure!
               if (f5>=f4)
                  { a3=a2; a2=a4; a1=a5;  f3=f2; f2=f4; f1=f5; }
               else {
                  a6=(a1+a5)/2; f6=fun_LineSearch(a6, io,x0,p,x,n);
                  if(f6 != f6)return NAN;  //failure!
                  if (f6>f5)
                       { a1=a6; a2=a5; a3=a4;  f1=f6; f2=f5; f3=f4; }
                  else { a2=a6; a3=a5; f2=f6; f3=f5; }
               }
            }
         }
      }
      else {                     /* fig 2.2.9 */
         if (a4-a2>smallgapa*(a3-a2)) {
            if (f2>=f4) { a1=a2; a2=a4;  f1=f2; f2=f4; }
            else        { a3=a4; f3=f4; }
         }
         else {
            if (f4>f2) {
               a5=(a1+a2)/2; f5=fun_LineSearch(a5, io,x0,p,x,n);
               if(f5 != f5)return NAN;  //failure!
               if (f5>f2) { a1=a5; a3=a4;  f1=f5; f3=f4; }
               else       { a3=a2; a2=a5;  f3=f2; f2=f5; }
            }
            else {
               a5=(a3+a4)/2; f5=fun_LineSearch(a5, io,x0,p,x,n);
               if(f5 != f5)return NAN;  //failure!
               if (f5>=f4)
                  { a1=a2; a2=a4; a3=a5;  f1=f2; f2=f4; f3=f5; }
               else {
                  a6=(a3+a5)/2; f6=fun_LineSearch(a6, io,x0,p,x,n);
                  if(f6 != f6)return NAN;  //failure!
                  if (f6>f5)
                      { a1=a4; a2=a5; a3=a6;  f1=f4; f2=f5; f3=f6; }
                  else { a1=a5; a2=a6;  f1=f5; f2=f6; }
               }
            }
         }
      }
   }

   if (f2>f0 && f4>f0)  a4=0;
   if (f2<=f4)  { *f=f2; a4=a2; }
   else         *f=f4;
   

   return (a4);
}

/*********************************************************/
phydbl fun_LineSearch (phydbl t, option *io, phydbl x0[], phydbl p[], phydbl x[], int n)
{  
  int i;   
  For (i,n) x[i]=x0[i] + t*p[i];   
  return( LK_BFGS_from_CODEML(io,x,n) );

}

/*********************************************************/

int zero (phydbl x[], int n)
{ int i; for(i=0; i<n; i++) x[i]=0; return (0);}

/*********************************************************/

int identity (phydbl x[], int n)
{ int i,j;  for(i=0; i<n; i++)  { for(j=0; j<n; j++)   x[i*n+j]=0;  x[i*n+i]=1; }  return (0); }

/*********************************************************/

phydbl norm (phydbl x[], int n)
{ int i; phydbl t=0;  for(i=0; i<n; i++)  t += x[i]*x[i];  return sqrt(t); }

/*********************************************************/

int H_end (phydbl x0[], phydbl x1[], phydbl f0, phydbl f1,
    phydbl e1, phydbl e2, int n)
/*   Himmelblau termination rule.   return 1 for stop, 0 otherwise.
*/
{
   phydbl r;
   if((r=norm(x0,n))<e2)
      r=1;
   r*=e1;
   if(distance(x1,x0,n)>=r)
      return(0);
   r=fabs(f0);  if(r<e2) r=1;     
   r*=e1;
   if(fabs(f1-f0)>=r) 
      return(0);
   return (1);
}

/*********************************************************/

phydbl distance (phydbl x[], phydbl y[], int n)
{  int i; phydbl t=0;
   for (i=0; i<n; i++) t += square(x[i]-y[i]);
   return(sqrt(t));
}

/*********************************************************/
phydbl Sum (phydbl x[], int n)
{ 
  int i; 
  phydbl t=0;  
  for(i=0; i<n; i++) t += x[i];    
  return(t); 
}
