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
#include "models.h"
#include <stdio.h>


//!< Start added by Marcelo.
#include "pcm.h"
#include "ecm.h"

extern int     stopCodons[64]; 
extern int    senseCodons[64]; 
extern char  aminoAcidmap[65]; 


//!< End added by Marcelo.

/*********************************************************/
/* Handle any number of states (>1) */
/*! Jukes Cantor */
void PMat_JC69(phydbl l, int pos, phydbl *Pij, model *mod){
    int ns;
    int i,j;
    ns = mod->ns;
    For(i,ns) Pij[pos+ ns*i+i] = 1. - ((ns - 1.)/ns)*(1. - EXP(-ns*l/(ns - 1.)));
    For(i,ns-1)
    for(j=i+1;j<ns;j++) {
        Pij[pos+ ns*i+j] = (1./ns)*(1. - EXP(-ns*l/(ns - 1.)));
        if(Pij[pos+ns*i+j] < SMALL_PIJ) Pij[pos+ns*i+j] = SMALL_PIJ;
        Pij[pos+ ns*j+i] = Pij[pos+ ns*i+j];
    }
}


/*********************************************************/

void PMat_Zero_Br_Len(model  *mod, int pos, phydbl *Pij){
    int n = mod->ns;
    int i;
    For (i,n*n) Pij[pos+i] = 0.0;  //!< Changed by Marcelo.
    For(i,n) Pij[pos+n*i+i] = 1.0;//!< Changed by Marcelo.
}

/*********************************************************/
void PMat(phydbl l, model *mod, int pos, phydbl *Pij){
    if(l < BL_MIN-POW(2,-27)){
        PMat_Zero_Br_Len(mod,pos,Pij);
    }else{
        switch(mod->datatype){
            case NT :{
            	Warn_And_Exit("Only codon models supported\n");
                break;
            }
            case AA :{
                Warn_And_Exit("Only codon models supported\n");
                break;
                }

            case CODON :{
                    if(mod->calculate_init_tree && mod->init_DistanceTreeCD==NUCLEO) PMat_JC69(l,pos,Pij,mod);
                    else PMat_CODON(l,mod,0,Pij);
                    break;
                }

            default:{
                    PMat_JC69(l,pos,Pij,mod);
                    break;
                }
        }
    }
}

/*********************************************************/
void Init_Model(calign *data, model *mod, option *io){
    int i,j;
    phydbl sum, aux, mr;
    int result;
    if(io->datatype == CODON) { //!<Added by Marcelo.
        if(!mod->invar) {
            For(i, data->crunch_len) {
                data->invar[i] = 0;
            }
        }
        if(mod->s_opt->opt_pinvar)mod->pinvar = 0.2;
        switch(mod->omegaSiteVar) {
            case DM0: {
                mod->alpha = 1.0;
                mod->omegas[0] = 1.0;
                mod->prob_omegas[0] = 1.0;
                //DiscreteGamma(mod->gamma_r_proba, mod->gamma_rr, mod->alpha, mod->alpha, mod->n_catg, mod->gamma_median);
                break;
            }
            case DGAMMAK: {
                DiscreteGamma(mod->prob_omegas, mod->omegas, mod->alpha, mod->beta, mod->n_w_catg, mod->gamma_median);
                For(i, mod->n_w_catg) {
                    if(mod->omegas[i] < MODELPAREPS) mod->omegas[i] = MODELPAREPS;
                }
                break;
            }
            case DMODELK: {
                Scale_freqs_tol(mod->prob_omegas, mod->n_w_catg, MODELPAREPS, 0.99);
                Freq_to_UnsFreq(mod->prob_omegas, mod->prob_omegas_uns, mod->n_w_catg, 1);
                break;
            }
            default:
                break; 
        } 
        if((mod->initqrates != NOINITMAT) && (mod->omegaSiteVar != NOOMEGA) && (io->kappaECM == kap5)) {
        	Warn_And_Exit("Init_Model options not supported.");
        	/*Scale_freqs(mod->pkappa, mod->nkappa);
            Freq_to_UnsFreq(mod->pkappa, mod->unspkappa, mod->nkappa, 1);*/
        }
        if((mod->freq_model != FUNDEFINED) && (mod->freq_model != FMODEL)) {
            if(mod->s_opt->user_state_freq) {
                mod->s_opt->opt_state_freq = NO;
                if(mod->freq_model<ROOT)EqFrequencies(mod->freq_model, mod->pi, mod->user_b_freq, mod->ns);
            }else{
                if(mod->freq_model<ROOT)EqFrequencies(mod->freq_model, mod->pi, mod->base_freq, mod->ns);
                switch(mod->freq_model) {
                    case F1X4: {
                        Scale_freqs(mod->base_freq, mod->num_base_freq);
                        Freq_to_UnsFreq(mod->base_freq, mod->uns_base_freq, mod->num_base_freq, 1);
                        break;
                    }
                    case F3X4:
                    case CF3X4: {
                        Scale_freqs(mod->base_freq, 4);
                        Freq_to_UnsFreq(mod->base_freq, mod->uns_base_freq, 4, 1);

                        Scale_freqs(mod->base_freq+4, 4);
                        Freq_to_UnsFreq(mod->base_freq+4, mod->uns_base_freq+4, 4, 1);
                        
                        Scale_freqs(mod->base_freq+8, 4);
                        Freq_to_UnsFreq(mod->base_freq+8, mod->uns_base_freq+8, 4, 1);
                        
                        break;
                    }
                    case ROOT: case MROOT:
						break;
                    default:
                        break;
                } 
            }
            if(mod->freq_model<ROOT){
            	Scale_freqs(mod->pi, mod->ns);
            	Freq_to_UnsFreq(mod->pi, mod->pi_unscaled, mod->ns, 1);
            }
        } else { 
            if((mod->initqrates == ECMUSR) && (mod->freq_model == FUNDEFINED || mod->freq_model == FMODEL)) {
                For(i, mod->ns) mod->pi[i] = mod->userfreq[senseCodons[i]]; //was io-> mod->userfreq[senseCodons[i]]; Ken 9/1/2018
                if(mod->whichrealmodel == MG) {
                    For(i, mod->num_base_freq) mod->base_freq[i] = mod->userbfreq[i]; //was io-> mod->userbfreq[i];
                }
                mod->s_opt->user_state_freq = NO;
            }
            
            if(mod->whichrealmodel == PCM) {
                mod->s_opt->user_state_freq = NO;
            } else if(mod->initqrates == KOSI07) {
                For(i, mod->ns) mod->pi[i] = ecmK07freq[i];
                mod->s_opt->user_state_freq = NO;
            } else if(mod->initqrates == SCHN05) {
                For(i, mod->ns) mod->pi[i] = ecmS05freq[senseCodons[i]];
                mod->s_opt->user_state_freq = NO;
            }
            mod->s_opt->opt_state_freq = NO;
        }
        if((mod->s_opt->opt_omega == NO) && 
           (mod->s_opt->opt_kappa == NO) && 
           (mod->s_opt->opt_state_freq == NO)) {
        	if(mod->nomega_part > 1){Warn_And_Exit("options not compatible with partitioned model error 8\n");}
            mr = Update_Qmat_Codons(mod, 0, 0, mod->pi); //modified by Ken 19/8
            EigenQREV(mod->qmat, mod->pi, mod->ns, mod->eigen->e_val, mod->eigen->r_e_vect, mod->eigen->l_e_vect, mod->eigen->space);
            For(i, mod->ns) mod->eigen->e_val[i] /= mr; //was io-> mod->ns 9/1 Ken
            mod->update_eigen = NO;
        } else {
        	if(mod->optDebug)printf("here3.3\n");
            mod->update_eigen = YES;
            Set_Model_Parameters(mod);
            mod->update_eigen = NO;
            if(mod->optDebug)printf("here3.4\n");
        }
        if(mod->optDebug)printf("here4\n");
        mod->beta_old = mod->beta;
    }else{
    	Warn_And_Exit("Currently only codon models supported\n");
    }
    mod->alpha_old   = mod->alpha;
    mod->kappa_old   = mod->kappa;
    mod->lambda_old  = mod->lambda;
    mod->pinvar_old  = mod->pinvar;
}

/*********************************************************/
void Set_Model_Parameters(model *mod) {
    if(mod->datatype == CODON) { //!Added by Marcelo.
        phydbl mr;
        int modeli,i, j, k, n, nn, n_termsTaylor;
        switch(mod->omegaSiteVar) {
            case DGAMMAK: { 
                DiscreteGamma(mod->prob_omegas, mod->omegas, mod->alpha, mod->beta, mod->n_w_catg, mod->gamma_median);
                For(i, mod->n_w_catg) if(mod->omegas[i] < MODELPAREPS) mod->omegas[i] = MODELPAREPS;
                break;
            }
            case DMODELK: { 
                Freq_to_UnsFreq(mod->prob_omegas, mod->prob_omegas_uns, mod->n_w_catg, 0);
                Scale_freqs_tol(mod->prob_omegas, mod->n_w_catg, MODELPAREPS, 0.99);
                break;
            }
            case DM0: {
                //DiscreteGamma(mod->gamma_r_proba, mod->gamma_rr, mod->alpha, mod->alpha, mod->n_catg, mod->gamma_median);
                break;
            }
            default:
                break;
        }
        if((mod->s_opt->opt_state_freq) &&
           (mod->s_opt->opt_omega == NO))Lazy_Exit("Optimizing frequencies without omega",__FILE__,__LINE__);

        //udpate model parameters?
        if(mod->update_eigen) {
            if(mod->n_w_catg == 1) {
               For(modeli,mod->nomega_part){ //Ken 19/8
            	   if(mod->freq_model != MROOT)mod->mr_w[0] = Update_Qmat_Codons(mod, 0, modeli,mod->pi); //Ken 19/8
            	   else mod->mr_w[0] = Update_Qmat_Codons(mod, 0, modeli,mod->mid_pi[modeli]); //Ken 19/8

            	   if(mod->expm == EIGEN) {
            		   if(mod->nomega_part > 1){Lazy_Exit("Eigenvalue exponentiation with partitioned models",__FILE__,__LINE__);}
            		   EigenQREV(mod->qmat_part[0], mod->pi, mod->ns, mod->eigen->e_val, mod->eigen->r_e_vect, mod->eigen->l_e_vect, mod->eigen->space);
            		   For(i,mod->ns) mod->eigen->e_val[i]/=mod->mr_w[0];
            	   }else if(mod->expm == TAYLOR) {
            		   Lazy_Exit("Taylor series expansion",__FILE__,__LINE__);
            	   }else if(mod->expm == SSPADE) { //Modified by Ken 17/8/2016
            		   n=mod->ns;
            		   nn=n*n;
            		   For(i,nn) mod->qmat_part[modeli][i]/=mod->mr_w[0];

#if defined BLAS || defined BLAS_OMP
            		//multiply these matrices together
                    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, n, n, n, 1.0, mod->qmat_part[modeli], n, mod->qmat_part[modeli],  n, 0.0, mod->A2_part[modeli], n);
                    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, n, n, n, 1.0, mod->A2_part[modeli],   n, mod->A2_part[modeli],    n, 0.0, mod->A4_part[modeli], n);
                    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, n, n, n, 1.0, mod->A2_part[modeli],   n, mod->A4_part[modeli],    n, 0.0, mod->A6_part[modeli], n);
                    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, n, n, n, 1.0, mod->A2_part[modeli],   n, mod->A6_part[modeli],    n, 0.0, mod->A8_part[modeli], n);
                    //copy vectors
                    cblas_dcopy(nn,mod->A2_part[modeli], 1, mod->Apowers_part[modeli]+  nn,1);
                    cblas_dcopy(nn,mod->A4_part[modeli], 1, mod->Apowers_part[modeli]+2*nn,1);
                    cblas_dcopy(nn,mod->A6_part[modeli], 1, mod->Apowers_part[modeli]+3*nn,1);
                    cblas_dcopy(nn,mod->A8_part[modeli], 1, mod->Apowers_part[modeli]+4*nn,1);
#else
                    For(i,nn) mod->A2_part[modeli][i]=0.0;
                    For(i,nn) mod->A4_part[modeli][i]=0.0;
                    For(i,nn) mod->A6_part[modeli][i]=0.0;
                    For(i,nn) mod->A8_part[modeli][i]=0.0;
                    //multiply these matrices together
                    For(i,n) For(j,n) For(k,n) mod->A2_part[modeli][n*i+j] += mod->qmat_part[modeli][i*n+k] * mod->qmat_part[modeli][k*n+j];
                    For(i,n) For(j,n) For(k,n) mod->A4_part[modeli][n*i+j] += mod->A2_part[modeli][i*n+k]   * mod->A2_part[modeli][k*n+j];
                    For(i,n) For(j,n) For(k,n) mod->A6_part[modeli][n*i+j] += mod->A2_part[modeli][i*n+k]   * mod->A4_part[modeli][k*n+j];
                    For(i,n) For(j,n) For(k,n) mod->A8_part[modeli][n*i+j] += mod->A2_part[modeli][i*n+k]   * mod->A6_part[modeli][k*n+j];
                    //copy vectors
                    For(i,nn) mod->Apowers_part[modeli][nn+i]  =mod->A2_part[modeli][i];
                    For(i,nn) mod->Apowers_part[modeli][2*nn+i]=mod->A4_part[modeli][i];
                    For(i,nn) mod->Apowers_part[modeli][3*nn+i]=mod->A6_part[modeli][i];
                    For(i,nn) mod->Apowers_part[modeli][4*nn+i]=mod->A8_part[modeli][i];
#endif
            	}
            }//for(modeli)
          }else if(mod->n_w_catg > 1)Lazy_Exit("Multiple omega categories",__FILE__,__LINE__);
        }
    }else{
    	Warn_And_Exit("\nOnly codon data currently supported\n");
    }
}

/*********************************************************/

phydbl F1x4(int codon, phydbl *freq){ //!<Added by marcelo.
    phydbl val=1.0;
    int i, remainder=0;
    For(i,3){
        remainder=codon-((codon>>2)<<2); //!<  n<<k= n*2^k, n>>k=n/2^k.
        codon=codon>>2;
        val*=freq[remainder];
    }
    return val;
}

/*********************************************************/

phydbl F3x4(int codon, phydbl *freq) { //!<Added by marcelo.
    phydbl val=1.0;
    int i, remainder=0;
    for(i=2;i>=0;i--){
        remainder=codon-((codon>>2)<<2); //!<  n<<k= n*2^k, n>>k=n/2^k.
        codon=codon>>2;
        val*=freq[ remainder+i*(2<<1) ];
    }
    return val;
}

/*********************************************************/

void EqFrequencies(int modfreq, phydbl *pi, phydbl *freq, int numSensecodons){ //!<Added by marcelo.
    int i;
    phydbl freqStopcodons;
    
    switch(modfreq) {
        case F1XSENSECODONS: {
            For(i,numSensecodons) pi[i]=freq[i];//!< Initialize the values before optimization
            break;
        }
        case F1X4: {
            freqStopcodons=0.0;//! Calculate the total frequency of stop codons to correct the frequency of the sense codons.
            For(i,64) if(stopCodons[i]) freqStopcodons+=F1x4(i,freq);
            For(i,numSensecodons) pi[i] = F1x4(senseCodons[i], freq)/(1-freqStopcodons); 
            break;
        }
        case F3X4:
        case CF3X4: {
            freqStopcodons=0.0;//! Calculate the total frequency of stop codons to correct the frequency of the sense codons.
            For(i,64) if(stopCodons[i]) freqStopcodons+=F3x4(i,freq);
            For(i,numSensecodons) pi[i] = F3x4(senseCodons[i], freq)/(1-freqStopcodons); 
            break;
        }
        default :{
            Warn_And_Exit("Frequency model not implemented.\n"); break;
        }
    }
}
/**************************************************************/

phydbl Update_Qmat_Codons(model *mod, int cat, int modeli, phydbl* freqs) {
    int numSensecodons, i, j, allocFreqs;
    phydbl sum, mu;
    phydbl *qmat, *mat;
    
    //incorporate equilibrium frequencies in case they have been updated
    switch(mod->freq_model) {
        case F1XSENSECODONS: {
          	Freq_to_UnsFreq(mod->pi, mod->pi_unscaled, mod->ns, 0);
          	break;
        }
        case F1X4: {
            Freq_to_UnsFreq(mod->base_freq, mod->uns_base_freq, mod->num_base_freq, 0);
            EqFrequencies(mod->freq_model, mod->pi, mod->base_freq, mod->ns);
            break;
        }
        case F3X4:
        case CF3X4: {
            Freq_to_UnsFreq(mod->base_freq,   mod->uns_base_freq,   4, 0); //convert to decimal
            Freq_to_UnsFreq(mod->base_freq+4, mod->uns_base_freq+4, 4, 0);
            Freq_to_UnsFreq(mod->base_freq+8, mod->uns_base_freq+8, 4, 0);
            EqFrequencies(mod->freq_model, mod->pi, mod->base_freq, mod->ns);
            break;
        }
        case ROOT:

        	if(mod->tree_loaded) freqs=mod->tree->noeud[mod->freq_node]->partfreqs[modeli];
        	else freqs = mod->root_pi[modeli];
        	break;
        case MROOT:
        	if(mod->tree_loaded>=1){
        		if(mod->tree_loaded==1)Update_Midpoint_Freqs(mod->tree);
        		freqs=mod->mid_pi[modeli];
        	}else{
        		freqs=mod->root_pi[modeli];
        	}
        	break;
        default:
            break;
     }


    qmat = mod->qmat_part[modeli] + (cat * mod->ns * mod->ns);
    numSensecodons = mod->ns;
    allocFreqs = NO;
    
    // Deal with initial rate matrices
    // KOSI07 used when an initial tree topology not provided
    if(mod->initqrates != NOINITMAT) {
        switch(mod->initqrates) {
            case KOSI07:
                mat = (phydbl *) ecmK07;
                if(mod->optDebug)printf("USING KOSI07\n");
                freqs = ecmK07freq;
                break;
                
           case SCHN05:
                mat = (phydbl *) mCalloc(numSensecodons * numSensecodons, sizeof(phydbl));
                if(mod->optDebug)printf("USING SHN05n");
                For(i, numSensecodons) {
                    For(j, numSensecodons) {
                        mat[ i*numSensecodons + j ] = ecmS05[senseCodons[i]][senseCodons[j]];
                    }
                }
                freqs = (phydbl *) mCalloc(numSensecodons, sizeof(phydbl));
                For(i, numSensecodons) {
                    freqs[i] = ecmS05freq[senseCodons[i]];
                }
                allocFreqs = YES;
                break;
                
            case ECMUSR:
                if(mod->calculate_init_tree == 1) {
                    mat = (phydbl *) mCalloc(numSensecodons * numSensecodons, sizeof(phydbl));
                    For(i, numSensecodons) {
                        For(j, numSensecodons) mat[i*numSensecodons+j] = mod->userRatesT[senseCodons[i]][senseCodons[j]];
                    }
                    freqs = (phydbl *) mCalloc(numSensecodons, sizeof(phydbl));
                    For(i, numSensecodons) freqs[i] = mod->userfreqT[senseCodons[i]];
                    allocFreqs = YES;
                } else {
                    mat = (phydbl *) mCalloc(numSensecodons * numSensecodons, sizeof(phydbl));
                    For(i, numSensecodons) {
                        For(j, numSensecodons) mat[i*numSensecodons+j] = mod->userRates[senseCodons[i]][senseCodons[j]];
                    }
                    freqs = (phydbl *) mCalloc(numSensecodons, sizeof(phydbl));
                    For(i, numSensecodons) freqs[i] = mod->userfreq[senseCodons[i]];
                    allocFreqs = YES;
                }
                break;
            default:
                break;
        }
    } else {
    	//if not an initial topology search, start off with a blank Qmat
        mat = (phydbl *) mCalloc(numSensecodons * numSensecodons, sizeof(phydbl));
        For(i, numSensecodons * numSensecodons) mat[i] = 1.0;
        if(mod->freq_model < ROOT){
        	freqs = (phydbl *) mCalloc(numSensecodons, sizeof(phydbl));
        	For(i, numSensecodons) freqs[i] = 1.0;
        }
    }

    // deal with frequency models. Usefd in GY and HLP17
    if((mod->freq_model != FMODEL && mod->freq_model != FUNDEFINED) && mod->freq_model < ROOT) {
        if((mod->initqrates == NOINITMAT && mod->pcaModel == NO) || mod->initqrates == SCHN05) {
            free(freqs);
            allocFreqs = NO;
        }
        freqs = mod->pi;
    }else{ //used in HLP19
        for(i=0; i<numSensecodons; i++) {
        	if(mod->freq_model < ROOT) mod->pi[i] = freqs[i];
        }
    }

    For(i, numSensecodons*numSensecodons) qmat[i] = 0.0;
    
    // calculate the actual Q matrix
    switch(mod->whichrealmodel) {
    	case GY: case PCM:
            Update_Qmat_GY(mat, qmat, freqs, cat, mod);
            break;
        case HLP17:
        case HLP19:
            Update_Qmat_HLP17(mat, qmat, freqs, cat, mod,mod->omega_part[modeli],modeli);
            break;
        default:
            break;
    }

    /*! Calculate the diagonal element which is the negative of the sum of the other elements in a row.*/
    For(i, numSensecodons) {
        sum = 0.0;
        For(j, numSensecodons) sum += qmat[ numSensecodons*i+j ];
        qmat[ numSensecodons*i+i ]= (-1)*sum;    
    }
    
    /*! Normalize Q matrix.*/
    mu = 0;
    For(i, numSensecodons) {
        mu += freqs[i] * (-qmat[numSensecodons*i+i]);
    }

    //Added by Ken - normalizes matrix here.
    if(mod->whichrealmodel <= HLP17){
      For(i, numSensecodons) {
        For(j, numSensecodons) {
    		qmat[numSensecodons*i+j] = qmat[numSensecodons*i+j]/mu;
    	}
       }
       mu = 0;
       For(i, numSensecodons) {
           mu += freqs[i] * (-qmat[numSensecodons*i+i]);
       }
    }

    if(mod->initqrates == NOINITMAT) {
       free(mat);
    }
    
    if(mod->initqrates == SCHN05) {
        if(allocFreqs) {
            free(freqs);
        }
        free(mat);
    }
    return mu;
}

/**************************************************************/
void Update_Qmat_GY(phydbl *mat, phydbl *qmat, phydbl * freqs, int cat, model *mod) {
    int i, j, numSensecodons;
    phydbl value;
    numSensecodons = mod->ns;
    For(i, numSensecodons) {
        For(j, i) {
            value = mat[i*numSensecodons+j] * Kappa_Omega_Factor(i, j, mod, cat,mod->omega_part[0]);
              qmat[ i*numSensecodons+j ] = value * (freqs[j]);
              qmat[ j*numSensecodons+i ] = value * (freqs[i]);
        }
    }
}

/**************************************************************/
//Update parameters of HLP-type models
//Modified to be omega*kappa*pi*(1+b*h) on 12/Jun/2016
//tally up expected number of each type of hotspot mutation
//Cut down on RAM usage by using copy of hotspot tables at upper level
void Update_Qmat_HLP17(phydbl *mat, phydbl *qmat, phydbl * freqs, int cat, model *mod,phydbl omega,int modeli) {
	 int i,j,fi,ti,li,ri,hot,c;
	 double htotal[mod->nmotifs];

	 //skip B matrix construction if constant B matrix used
	 if(!mod->constB){
		 For(fi,61){ //Fill in B matrix
			 For(ti,61){
				For(c,mod->nmotifs)htotal[c]=0; //set htotal array to zero
#if defined OMP || BLAS_OMP
    			omp_lock_t writelock;
    			omp_init_lock(&writelock);
    			omp_set_lock(&writelock);
#endif
    				For(li,61)
    					For(ri,61)
    						For(c,mod->nmotifs)
    							htotal[c] += freqs[li]*freqs[ri]*
    								mod->io->mod->hotspotcmps[c][fi*61*61*61+ti*61*61+li*61+ri];

#if defined OMP || BLAS_OMP
    				omp_unset_lock(&writelock);
    				omp_destroy_lock(&writelock);
#endif
    				//add up hotspot values
    				mod->Bmat[fi*61+ti] = 0.0;
    			    For(c,mod->nmotifs)mod->Bmat[fi*61+ti] += htotal[c]*mod->hotness[mod->motif_hotness[c]];
    			    if(mod->Bmat[fi*61+ti] < -1)mod->Bmat[fi*61+ti]=-1; //constrain total modification to never go below -1
    		}
    	}
    }

    int numSensecodons = mod->ns;
    For(i,numSensecodons){
        For(j, i) {
        	if(mod->constB){
        		mod->Bmat[i*numSensecodons+j]=0.0;
        		mod->Bmat[j*numSensecodons+i]=0.0;
        		For(c,mod->nmotifs){
        			mod->Bmat[i*numSensecodons+j] += mod->cBmat[modeli][i*61+j][c]*mod->hotness[mod->motif_hotness[c]];
        			mod->Bmat[j*numSensecodons+i] += mod->cBmat[modeli][j*61+i][c]*mod->hotness[mod->motif_hotness[c]];
        		}
        		if(mod->Bmat[i*numSensecodons+j] < -1)mod->Bmat[i*numSensecodons+j]=-1;
        		if(mod->Bmat[j*numSensecodons+i] < -1)mod->Bmat[j*numSensecodons+i]=-1;
        	}
            phydbl value = mat[i*numSensecodons+j] * Kappa_Omega_Factor(i, j, mod, cat,omega);
            if(mod->modeltypeOpt==HLP17 && mod->freqsTo){ //HLP17 freqs
            	  qmat[ i*numSensecodons+j ] = value * freqs[j]*(1+mod->Bmat[ i*numSensecodons+j ]);
            	  qmat[ j*numSensecodons+i ] = value * freqs[i]*(1+mod->Bmat[ j*numSensecodons+i ]);
            }else{ //HLP19 freqs
            	qmat[ i*numSensecodons+j ] = value * (1.0/61)*(1+mod->Bmat[ i*numSensecodons+j ]);
            	qmat[ j*numSensecodons+i ] = value * (1.0/61)*(1+mod->Bmat[ j*numSensecodons+i ]);
            }
        }
    }
}

/**************************************************************/
//set up constant Bmat if desired
void Setup_CBmat(model* mod, int uniform, phydbl* pis){
  int modeli,fi,ti,li,ri,hot,c;
  mod->cBmat=(phydbl ***)mCalloc(mod->nomega_part,sizeof(phydbl**));
  For(modeli,mod->nomega_part){
		mod->cBmat[modeli]=(phydbl **)mCalloc(3721,sizeof(phydbl*));
		double htotal[mod->nmotifs];
		  For(fi,61){ //Fill in B matrix
			For(ti,61){
				mod->cBmat[modeli][fi*61+ti]=mCalloc(mod->nmotifs,sizeof(phydbl));
				if(!mod->io->precon){
				For(li,61){
					For(ri,61){
						For(c,mod->nmotifs){
							if(uniform==1) mod->cBmat[modeli][fi*61+ti][c] += (1.0/61)*(1.0/61)*
									mod->io->mod->hotspotcmps[c][fi*61*61*61+ti*61*61+li*61+ri];
							else if(uniform==-1){
								mod->cBmat[modeli][fi*61+ti][c] += mod->root_pi[modeli][li]*
								mod->root_pi[modeli][ri]*mod->io->mod->hotspotcmps[c][fi*61*61*61+ti*61*61+li*61+ri];
							}else{ mod->cBmat[modeli][fi*61+ti][c] += pis[li]*pis[ri]*
									mod->io->mod->hotspotcmps[c][fi*61*61*61+ti*61*61+li*61+ri];
							}
						}
					}
				}
			 }
		  }
		}
  	}
}

/*******************************************************/
// Returns a single factor for a given site that includes both values for kappa
// and/or omega. If neither of those is needed, return value will simply be 1.0 .
phydbl Kappa_Omega_Factor(int senseCodoni, int senseCodonj, model* mod, int cat,phydbl omega) {//modified by Ken 18/8
    phydbl val = 1.0;
    phydbl *kappa = mod->pkappa;
    int numSensecodons = mod->ns;
    ts_and_tv *mat = mod->structTs_and_Tv;
   
    switch(mod->omegaSiteVar) {
        case DM0:
            //omega = omega;
            break;
        case NOOMEGA:
            omega = 1.0;
            break;
        default:
            omega = mod->omegas[ cat ];
            break;
    }
    // First deal with kappa...
        if(mod->initqrates != NOINITMAT) {
          /*  switch(mod->kappaECM) {
                case kap1: { 
                    val = 1.0;
                    break;
                }
                case kap2: {
                    switch(mat[numSensecodons * senseCodoni + senseCodonj].nts) {
                        case 0: val = 1.0; break;
                        case 1: val = kappa[0]; break;
                        case 2: val = kappa[0]*kappa[0]; break;
                        case 3: val = kappa[0]*kappa[0]*kappa[0]; break;
                        default: break;
                    }
                    break;
                }
                case kap3: {
                    switch(mat[numSensecodons * senseCodoni + senseCodonj].ntv) {
                        case 0: val = 1.0; break;
                        case 1: val = kappa[0]; break;
                        case 2: val = kappa[0] * kappa[0]; break;
                        case 3: val = kappa[0] * kappa[0] * kappa[0]; break;
                        default: break;
                    }
                    break;
                }
                case kap4:  {
                    switch(mat[numSensecodons * senseCodoni + senseCodonj].nts) {
                        case 0: val = 1.0; break;
                        case 1: val = kappa[0]; break;
                        case 2: val = kappa[0] * kappa[0]; break;
                        case 3: val = kappa[0] * kappa[0] * kappa[0]; break;
                        default: break;
                    }
                    switch(mat[numSensecodons * senseCodoni + senseCodonj].ntv) {
                        case 0: val *= 1.0; break;
                        case 1: val *= kappa[1]; break;
                        case 2: val *= kappa[1] * kappa[1]; break;
                        case 3: val *= kappa[1] * kappa[1] * kappa[1]; break;
                        default: break;
                    }
                    break;
                }
                case kap5: {
                    val = kappa[ mat[numSensecodons * senseCodoni + senseCodonj].caseK07 ];
                    break;
                }
                case kap6: {
                    if(mat[numSensecodons * senseCodoni + senseCodonj].ndiff > 1) {
                        val = kappa[0];
                    }
                    break;
                }
                default:
                    break;
            }*/
        } else {
            if(mat[numSensecodons * senseCodoni + senseCodonj].ndiff == 1) {
                if(mat[numSensecodons * senseCodoni + senseCodonj].nts) {
                    val = mod->kappa;
                } else {
                    val = 1.0;
                }
            } else {
                val = 0.0;
            }
        }
    // ... Then with omega:
    if(aminoAcidmap[ senseCodons[senseCodoni] ] != aminoAcidmap[ senseCodons[senseCodonj] ]) {
        val *= omega;
    }
    return val;
}

/*********************************************************/
void PMat_CODON(phydbl l, model *mod, int pos, phydbl *Pij){ //!< Added by Marcelo.
    int n = mod->ns, nn = n*n, nw=mod->n_w_catg, i, g, j, k, catg;
    
    if(mod->expm==EIGEN){
        phydbl *U, *V, *R, *P;
        phydbl *expt, expt_n, *pP;
        phydbl *uexpt, uexpt_n;
        
        P=Pij;
        
#if defined OMP || defined BLAS_OMP
        expt  = (phydbl *)malloc(n*nw*sizeof(phydbl));
        uexpt = (phydbl *)malloc(nn*nw*sizeof(phydbl));
#else
        expt  = mod->eigen->e_val_im;
        uexpt = mod->eigen->r_e_vect_im;
#endif
        
        U     = mod->eigen->r_e_vect;
        V     = mod->eigen->l_e_vect;
        R     = mod->eigen->e_val;
        
        For(catg,nw){
#if defined BLAS || defined BLAS_OMP
            /* compute EXP(D/mr * l) into mat_eDmrl */
            For(k,n) expt[k+catg*n] = (phydbl)exp(R[k+catg*n]*l);
            
            /* multiply Vr*EXP(D/mr*l)*Vi into Pij */
            For (i,n) For (k,n) uexpt[i*n+k+catg*nn] = U[i*n+k+catg*nn] * expt[k+catg*n];

            cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, n, n, n, 1.0, uexpt+catg*nn, n, V+catg*nn, n, 0.0, P+catg*nn, n);
#else
            for (k=0,zero(P+catg*nn,nn); k<n; k++) for (i=0,pP=P+catg*nn,expt_n=exp(R[k+catg*n]*l); i<n; i++) for (j=0,uexpt_n=U[i*n+k+catg*nn]*expt_n; j<n; j++) *pP++ += uexpt_n*V[k*n+j+catg*nn];
#endif 
        }
        
        For(i,nn*nw) if(Pij[i]<SMALL_PIJ) Pij[i]=SMALL_PIJ;
        
#if defined OMP || defined BLAS_OMP
        free(expt); 
        free(uexpt);
#endif
    }else if(mod->expm==TAYLOR){
    	if(mod->nomega_part > 1){
    		printf("TAYLOR APPROXIMATION DOESN'T WORK WITH PARITTIONED MODELS IN IGPHYML\n");
    		exit(EXIT_FAILURE);
    	}
        phydbl *Q2;
        int m;
        Q2=mod->A2_part[0]; //won't work with partitioned models
        
        For(i,nn*nw) Pij[i]=0.0;
        
        For(catg,nw){
            For(m,mod->n_termsTaylor){
                For(i,nn) Pij[i+catg*nn]+=pow(l,m)*Q2[i+catg*nn*15+m*nn]/(phydbl)myFactorial(m);
            }
        } 
        For(i,nn*nw) if(Pij[i]<SMALL_PIJ) Pij[i]=SMALL_PIJ;
    } else  if(mod->expm==SSPADE){
        phydbl norm_1_Q, *Q, *A, *B, *F, theta[5] = {1.495585217958292e-002, 2.539398330063230e-001, 9.504178996162932e-001, 2.097847961257068e+000, 5.371920351148152e+000}; 
        int m_vals[5]={3, 5, 7, 9, 13}, z=5;
        
#if defined OMP || defined BLAS_OMP
        
        B=(phydbl *)malloc(nn*sizeof(phydbl));
        A=(phydbl *)malloc(nn*sizeof(phydbl));
        
#endif
        
        int modeli;
		for(modeli=0;modeli<mod->nomega_part;modeli++){ //Added by Ken 22/8
        For(catg,nw){

            pos = nn*catg;
            Q   = mod->qmat_part[modeli]+pos;
            F   = Pij+pos;
            
#if defined OMP || defined BLAS_OMP
            
#else
            
            B = mod->matAux_part[modeli]+pos;
            A = mod->qmat_buff_part[modeli]+pos; //? 22/8 Ken
            
#endif
            
            For(i,nn) A[i]=Q[i]*l;
            
#if defined BLAS || defined BLAS_OMP
            
            norm_1_Q=dlange_("1", &n, &n, A, &n, NULL);
            
#else
            
            phydbl *max;
            max=(phydbl*)mCalloc(n,sizeof(phydbl));
            For(i,n) For(j,n) max[j]+=fabs(A[i*n+j]);
            norm_1_Q=DBL_MIN;
            For(i,n) if(max[i]>norm_1_Q) norm_1_Q=max[i];
            free(max);
            
#endif
            
            if(norm_1_Q<=theta[z-1]){
                For(i,z){
                    if(norm_1_Q<=theta[i]){
                        PadeApprox(n, nn, A, mod, F, pos, l, m_vals[i],modeli);
                        break;
                    }
                }
            }else{
                phydbl Mantissa, sFactor, num;
                int Exponent;
                num=norm_1_Q/theta[z-1];
                Mantissa=frexp(num,&Exponent);
                Exponent-=(Mantissa==0.5);
                sFactor=1.0/pow(2,Exponent);
                
                For(i,nn) A[i]*=sFactor;
                
                PadeApprox(n, nn, A, mod, F, pos, l, m_vals[z-1],modeli);
                
#if defined BLAS || defined BLAS_OMP
                
                Fors(i,Exponent,2){
                    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, n, n, n, 1.0, F, n, F, n, 0.0, B, n);
                    
                    if(i+1<Exponent){
                        cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, n, n, n, 1.0, B, n, B, n, 0.0, F, n);
                    }else{
                        cblas_dcopy(nn,B,1,F,1);
                        break;
                    }
                }
                
#else
                
                Fors(g,Exponent,2){
                    For(i,nn) B[i]=0.0;
                    For(i,n) For(j,n) For(k,n) B[n*i+j] += F[i*n+k] * F[k*n+j];
                    
                    if(g+1<Exponent){
                        For(i,nn) F[i]=0.0;
                        For(i,n) For(j,n) For(k,n) F[n*i+j] += B[i*n+k] * B[k*n+j];
                    }else{
                        For(i,nn) F[i]=B[i];
                        break;
                    }
                }
                
#endif
            }
        }
        For(i,nn*nw) if(Pij[i]<SMALL_PIJ) Pij[i]=SMALL_PIJ; //!< Correct for too small entries.
        
#if defined OMP || defined BLAS_OMP 
        
        free(A); 
        free(B);
        
#endif
    }// for(modeli=0..
    }
}


/*********************************************************/
void PMat_CODON_part(phydbl l, model *mod, int pos, phydbl *Qmat, phydbl *Pij, int modeli){ //!< Added by Marcelo.
    int n = mod->ns, nn = n*n, nw=mod->n_w_catg, i, g, j, k, catg;

    if(mod->expm==EIGEN){
        phydbl *U, *V, *R, *P;
        phydbl *expt, expt_n, *pP;
        phydbl *uexpt, uexpt_n;
        P=Pij;

#if defined OMP || defined BLAS_OMP
        expt  = (phydbl *)malloc(n*nw*sizeof(phydbl));
        uexpt = (phydbl *)malloc(nn*nw*sizeof(phydbl));
#else
        expt  = mod->eigen->e_val_im;
        uexpt = mod->eigen->r_e_vect_im;
#endif

        U     = mod->eigen->r_e_vect;
        V     = mod->eigen->l_e_vect;
        R     = mod->eigen->e_val;

        For(catg,nw){
#if defined BLAS || defined BLAS_OMP

            /* compute EXP(D/mr * l) into mat_eDmrl */
            For(k,n) expt[k+catg*n] = (phydbl)exp(R[k+catg*n]*l);

            /* multiply Vr*EXP(D/mr*l)*Vi into Pij */
            For (i,n) For (k,n) uexpt[i*n+k+catg*nn] = U[i*n+k+catg*nn] * expt[k+catg*n];

            cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, n, n, n, 1.0, uexpt+catg*nn, n, V+catg*nn, n, 0.0, P+catg*nn, n);

#else
            for (k=0,zero(P+catg*nn,nn); k<n; k++) for (i=0,pP=P+catg*nn,expt_n=exp(R[k+catg*n]*l); i<n; i++)
            	for (j=0,uexpt_n=U[i*n+k+catg*nn]*expt_n; j<n; j++) *pP++ += uexpt_n*V[k*n+j+catg*nn];
#endif
        }

        For(i,nn*nw) if(Pij[i]<SMALL_PIJ) Pij[i]=SMALL_PIJ;

#if defined OMP || defined BLAS_OMP

        free(expt);
        free(uexpt);

#endif
    }else if(mod->expm==TAYLOR){
    	if(mod->nomega_part > 1){
    		printf("TAYLOR APPROX DOESNT WORK WITH PARTITIONED MODELS IN IGPHYML YET\n");
    		exit(EXIT_FAILURE);
    	}
        phydbl *Q2;
        int m;
        Q2=mod->A2_part[0];

        For(i,nn*nw) Pij[i]=0.0;

        For(catg,nw){
            For(m,mod->n_termsTaylor){
                For(i,nn) Pij[i+catg*nn]+=pow(l,m)*Q2[i+catg*nn*15+m*nn]/(phydbl)myFactorial(m);
            }
        }
        For(i,nn*nw) if(Pij[i]<SMALL_PIJ) Pij[i]=SMALL_PIJ;
    }else  if(mod->expm==SSPADE){
        phydbl norm_1_Q, *Q, *A, *B, *F, theta[5] = {1.495585217958292e-002, 2.539398330063230e-001, 9.504178996162932e-001, 2.097847961257068e+000, 5.371920351148152e+000};
        int m_vals[5]={3, 5, 7, 9, 13}, z=5;

#if defined OMP || defined BLAS_OMP

        B=(phydbl *)malloc(nn*sizeof(phydbl));
        A=(phydbl *)malloc(nn*sizeof(phydbl));

#endif

        For(catg,nw){
            pos = nn*catg;
            Q   = Qmat+pos;
            F   = Pij+pos;

#if defined OMP || defined BLAS_OMP

#else

            B = mod->matAux_part[modeli]+pos;
            A = mod->qmat_buff_part[modeli]+pos;

#endif

            For(i,nn) A[i]=Q[i]*l;

#if defined BLAS || defined BLAS_OMP

            norm_1_Q=dlange_("1", &n, &n, A, &n, NULL);

#else

            phydbl *max;
            max=(phydbl*)mCalloc(n,sizeof(phydbl));
            For(i,n) For(j,n) max[j]+=fabs(A[i*n+j]);
            norm_1_Q=DBL_MIN;
            For(i,n) if(max[i]>norm_1_Q) norm_1_Q=max[i];
            free(max);

#endif

            if(norm_1_Q<=theta[z-1]){
                For(i,z){
                    if(norm_1_Q<=theta[i]){
                        PadeApprox(n, nn, A, mod, F, pos, l, m_vals[i],modeli);
                        break;
                    }
                }
            }else{
                phydbl Mantissa, sFactor, num;
                int Exponent;
                num=norm_1_Q/theta[z-1];
                Mantissa=frexp(num,&Exponent);
                Exponent-=(Mantissa==0.5);
                sFactor=1.0/pow(2,Exponent);

                For(i,nn) A[i]*=sFactor;

                PadeApprox(n, nn, A, mod, F, pos, l, m_vals[z-1],modeli);

#if defined BLAS || defined BLAS_OMP

                Fors(i,Exponent,2){
                    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, n, n, n, 1.0, F, n, F, n, 0.0, B, n);

                    if(i+1<Exponent){
                        cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, n, n, n, 1.0, B, n, B, n, 0.0, F, n);
                    }else{
                        cblas_dcopy(nn,B,1,F,1);
                        break;
                    }
                }

#else

                Fors(g,Exponent,2){
                    For(i,nn) B[i]=0.0;
                    For(i,n) For(j,n) For(k,n) B[n*i+j] += F[i*n+k] * F[k*n+j];

                    if(g+1<Exponent){
                        For(i,nn) F[i]=0.0;
                        For(i,n) For(j,n) For(k,n) F[n*i+j] += B[i*n+k] * B[k*n+j];
                    }else{
                        For(i,nn) F[i]=B[i];
                        break;
                    }
                }

#endif
            }
        }
        For(i,nn*nw) if(Pij[i]<SMALL_PIJ) Pij[i]=SMALL_PIJ; //!< Correct for too small entries.

#if defined OMP || defined BLAS_OMP

        free(A);
        free(B);

#endif
    }
}

/***********************************************************/
void PadeApprox(int n, int nn, phydbl *A, model *mod, phydbl *F, int pos, phydbl len, int m, int modeli) //!<  Added by Marcelo.
{
    phydbl c3[4]={120.0, 60.0, 12.0, 1.0}, c5[6]={30240.0, 15120.0, 3360.0, 420.0, 30.0, 1.0}, c7[8]={17297280.0, 8648640.0, 1995840.0, 277200.0, 25200.0, 1512.0, 56.0, 1.0}, c9[10]={17643225600.0, 8821612800.0, 2075673600.0, 302702400.0, 30270240.0, 2162160.0, 110880.0, 3960.0, 90.0, 1.0}, c13[14]={64764752532480000.0, 32382376266240000.0, 7771770303897600.0, 1187353796428800.0,  129060195264000.0, 10559470521600.0, 670442572800.0, 33522128640.0, 1323241920.0, 40840800.0, 960960.0, 16380.0, 182.0, 1.0};
    phydbl *ptCoeffPade, *U, *V, *A0, *A2, *A4, *A6, *A8, *Apowers, *Aux, len2, len4, len6, len8;
    phydbl lens[5];
    int j, o, p, q, *Ipiv;
    int info=0;
    
    switch(m)
    {	
        case 3:  ptCoeffPade=c3;  break;
        case 5:  ptCoeffPade=c5;  break;
        case 7:  ptCoeffPade=c7;  break;
        case 9:  ptCoeffPade=c9;  break;
        case 13: ptCoeffPade=c13; break;
        default: Warn_And_Exit("Invalid Pade Coefficient.");break;
    }
    
#if defined OMP || defined BLAS_OMP
    
    Aux=(phydbl *)malloc(nn*sizeof(phydbl));
    U=(phydbl *)malloc(nn*sizeof(phydbl));
    V=(phydbl *)malloc(nn*sizeof(phydbl));
    Ipiv=(int *)malloc(n*sizeof(int));
    
#else
    
    Aux =mod->matAux_part[modeli]+pos;
    U   =mod->U_part[modeli]+pos;
    V   =mod->V_part[modeli]+pos;
    Ipiv=mod->ipiv_part[modeli]+(pos/n);
    
#endif
    
    A0      = mod->A0_part[modeli];
    A2      = mod->A2_part[modeli]+pos;
    A4      = mod->A4_part[modeli]+pos;
    A6      = mod->A6_part[modeli]+pos;
    A8      = mod->A8_part[modeli]+pos;
    Apowers = mod->Apowers_part[modeli]+5*pos;
    
    len2=len*len;
    len4=len2*len2;
    len6=len2*len4;
    len8=len2*len6;
    lens[0]=1.0; lens[1]=len2; lens[2]=len4; lens[3]=len6; lens[4]=len8;
    
    switch(m)
    {
        case 3:  
        case 5:  
        case 7:  
        case 9:  
        {
#if defined BLAS || defined BLAS_OMP
            
            //For(j,nn) U[j]=0;
            For(j,nn) U[j]=lens[m/2]*ptCoeffPade[m]*Apowers[j+nn*(m/2)];      
            //for(j=m;j>0;j-=2) cblas_daxpy(nn, lens[j/2]*ptCoeffPade[j], Apowers+nn*(j/2), 1, U, 1);
            for(j=m-2;j>0;j-=2) cblas_daxpy(nn, lens[j/2]*ptCoeffPade[j], Apowers+nn*(j/2), 1, U, 1);
            
            cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, n, n, n, 1.0, A, n, U, n, 0.0, Aux, n);
            
            //For(j,nn) V[j]=0;
            For(j,nn) V[j]=lens[(m-1)/2]*ptCoeffPade[m-1]*Apowers[j+nn*((m-1)/2)];
            //for(j=m-1;j>-1;j-=2) cblas_daxpy(nn, lens[j/2]*ptCoeffPade[j], Apowers+nn*(j/2), 1, V, 1);
            for(j=m-3;j>-1;j-=2) cblas_daxpy(nn, lens[j/2]*ptCoeffPade[j], Apowers+nn*(j/2), 1, V, 1);
            
            cblas_dcopy(nn, Aux, 1, F, 1);
            cblas_daxpy(nn,  1, V, 1, F, 1);
            
            cblas_daxpy(nn, -1, Aux, 1, V, 1);
            
            dgesv_(&n, &n, V, &n, Ipiv, F, &n, &info);
            if(info!=0) Warn_And_Exit("Matrix cannot be inverted."); 
            
#else
            
            For(j,nn) U[j]=0;
            
            for(j=m;j>0;j-=2) For(o,nn) U[o]+=lens[j/2]*ptCoeffPade[j]*Apowers[nn*(j/2)+o];
            
            For(o,nn) Aux[o]=0.0;
            For(o,n) For(p,n) For(q,n) Aux[n*o+p] += A[o*n+q] * U[q*n+p];
            
            For(j,nn) V[j]=0;
            
            for(j=m-1;j>-1;j-=2) For(o,nn) V[o]+=lens[j/2]*ptCoeffPade[j]*Apowers[nn*(j/2)+o];
            
            For(o,nn) F[o]=Aux[o];
            For(o,nn) F[o]+=V[o];
            
            For(o,nn) V[o]-=Aux[o];
            
            if(MyGaussElimination_gNxRHS(V,F,n,n)!=0) Warn_And_Exit("Matrix cannot be inverted."); 
            
#endif
            
            break;
        }
        case 13: 
        {
#if defined BLAS || defined BLAS_OMP
            
            //For(j,nn) U[j]=0;
            For(j,nn) U[j]=len6*ptCoeffPade[13]*A6[j];
            //!< MATLAB (index -1 to C std) U = A * (A6*(c(14)*A6 + c(12)*A4 + c(10)*A2) + c(8)*A6 + c(6)*A4 + c(4)*A2 + c(2)*eye(n,classA));
            //cblas_daxpy(nn, len6*ptCoeffPade[13], A6, 1, U, 1);
            cblas_daxpy(nn, len4*ptCoeffPade[11], A4, 1, U, 1);
            cblas_daxpy(nn, len2*ptCoeffPade[ 9], A2, 1, U, 1);
            
            cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, n, n, n, len6, A6, n, U, n, 0.0, Aux, n);
            
            cblas_daxpy(nn, len6*ptCoeffPade[ 7], A6, 1, Aux, 1);
            cblas_daxpy(nn, len4*ptCoeffPade[ 5], A4, 1, Aux, 1);
            cblas_daxpy(nn, len2*ptCoeffPade[ 3], A2, 1, Aux, 1);
            cblas_daxpy(nn, ptCoeffPade[ 1], A0, 1, Aux, 1);
            cblas_dcopy(nn,Aux,1,U,1);
            
            cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, n, n, n, 1.0, A, n, U, n, 0.0, Aux, n);
            
            //For(j,nn) V[j]=0;
            For(j,nn) V[j]=len6*ptCoeffPade[12]*A6[j];
            //!< MATLAB (index -1 to C std) V = A6*(c(13)*A6 + c(11)*A4 + c(9)*A2) + c(7)*A6 + c(5)*A4 + c(3)*A2 + c(1)*eye(n,classA);
            //cblas_daxpy(nn, len6*ptCoeffPade[12], A6, 1, V, 1);
            cblas_daxpy(nn, len4*ptCoeffPade[10], A4, 1, V, 1);
            cblas_daxpy(nn, len2*ptCoeffPade[ 8], A2, 1, V, 1);
            
            cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, n, n, n, len6, A6, n, V, n, 0.0, U, n);
            
            cblas_daxpy(nn, len6*ptCoeffPade[ 6], A6, 1, U, 1);
            cblas_daxpy(nn, len4*ptCoeffPade[ 4], A4, 1, U, 1);
            cblas_daxpy(nn, len2*ptCoeffPade[ 2], A2, 1, U, 1);
            cblas_daxpy(nn, ptCoeffPade[ 0], A0, 1, U, 1);
            
            cblas_dcopy(nn, Aux, 1, F, 1);
            cblas_daxpy(nn,  1, U, 1, F, 1);
            
            cblas_daxpy(nn, -1, Aux, 1, U, 1);      
            
            dgesv_(&n, &n, U, &n, Ipiv, F, &n, &info);
            if(info!=0) Warn_And_Exit("Matrix cannot be inverted."); 
            
#else
            
            For(j,nn) U[j]=0;
            
            For(j,nn) U[j]+=len6*ptCoeffPade[13]*A6[j];
            For(j,nn) U[j]+=len4*ptCoeffPade[11]*A4[j];
            For(j,nn) U[j]+=len2*ptCoeffPade[ 9]*A2[j];
            
            For(o,nn) Aux[o]=0.0;
            For(o,n) For(p,n) For(q,n) Aux[n*o+p] +=  len6*A6[o*n+q] *  U[q*n+p];
            
            For(j,nn) Aux[j]+=len6*ptCoeffPade[ 7]*A6[j];
            For(j,nn) Aux[j]+=len4*ptCoeffPade[ 5]*A4[j];
            For(j,nn) Aux[j]+=len2*ptCoeffPade[ 4]*A2[j];
            For(j,nn) Aux[j]+=ptCoeffPade[ 1]*A0[j];
            
            For(o,nn) U[o]=0.0;
            For(o,n) For(p,n) For(q,n) U[n*o+p] +=  A[o*n+q] *  Aux[q*n+p];
            
            For(j,nn) V[j]=0;
            
            For(j,nn) V[j]+=len6*ptCoeffPade[12]*A6[j];
            For(j,nn) V[j]+=len4*ptCoeffPade[10]*A4[j];
            For(j,nn) V[j]+=len2*ptCoeffPade[ 8]*A2[j];
            
            For(o,nn) Aux[o]=0.0;
            For(o,n) For(p,n) For(q,n) Aux[n*o+p] +=  len6*A6[o*n+q] *  V[q*n+p];
            
            For(j,nn) Aux[j]+=len6*ptCoeffPade[ 6]*A6[j];
            For(j,nn) Aux[j]+=len4*ptCoeffPade[ 4]*A4[j];
            For(j,nn) Aux[j]+=len2*ptCoeffPade[ 2]*A2[j];
            For(j,nn) Aux[j]+=ptCoeffPade[ 0]*A0[j];
            
            For(o,nn) F[o]=U[o];
            For(o,nn) F[o]+=Aux[o];
            
            For(o,nn) Aux[o]-=U[o]; 
            
            if(MyGaussElimination_gNxRHS(Aux,F,n,n)!=0) Warn_And_Exit("Matrix cannot be inverted.");     
            
#endif
            
            break;
        }
        default:break;
    }
    
#if defined OMP || defined BLAS_OMP
    
    free(Aux);
    free(U);
    free(V);
    free(Ipiv);
    
#endif
}

/***********************************************************/
void PMat_CODON_Pairwise(phydbl l, phydbl *Pij,  phydbl *U, phydbl *V, phydbl *R, int n, phydbl *uexpt, phydbl *expt) //!< Added by Marcelo.
{
    int nn=n*n;
    int i, j, k;
    // int g;
    phydbl *pP, expt_n, uexpt_n;
    
#if defined BLAS || defined BLAS_OMP
    
    /* compute POW(EXP(D/mr),l) into mat_eDmrl */
    For(k,n) expt[k] = (phydbl)exp(R[k]*l);
    
    /* multiply Vr*POW(EXP(D/mr),l)*Vi into Pij */
    For (i,n) For (k,n) uexpt[i*n+k] = U[i*n+k] * expt[k];
    
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, n, n, n, 1.0, uexpt, n, V, n, 0.0, Pij, n);
    
#else
    
    For(i,nn) Pij[i] = 0.0;
    for (k=0; k<n; k++) for (i=0,pP=Pij,expt_n=exp(R[k]*l); i<n; i++) for (j=0,uexpt_n=U[i*n+k]*expt_n; j<n; j++) *pP++ += uexpt_n*V[k*n+j];
    
#endif
    For(i,nn) if(Pij[i]<SMALL_PIJ) Pij[i]=SMALL_PIJ;
}


/**********************************************************
 * Added by Ken - sets up repertoire-wide model settings
 * declares necessary data structures
 * */
void Setup_Repertoire_Models(option* io){
	int i,j;
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
					tree->noeud[tree->mod->startnode]->b[0], 0.0, tree);
			tree->mod->midpoint_div = tl/(tree->n_otu*2.0);
			//if(io->mod->optDebug)
				//printf("Total divergence: %lf %lf\n",tl,tree->mod->midpoint_div);
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

	//set up data structures to store parameters
	int nparams=0;
	int nmodparams=2+io->mod->nomega_part+io->mod->nhotness+12;
	nparams += nmodparams*(io->ntrees+1);
	For(j,io->ntrees)nparams += 2*io->tree_s[j]->n_otu-3;
	if(io->mod->optDebug)printf("\nNeed to store %d parameters",nparams);
	io->paramStore=mCalloc(nparams,sizeof(phydbl));
	io->nparams=nparams;

	//need to copy this parameter over
	io->mod->s_opt->min_diff_lk_global = io->min_diff_lk_global;

	//set starting time for computations
	#if defined OMP || defined BLAS_OMP
	  io->t_beg=omp_get_wtime();
	#else
	  time(&io->t_beg);
	#endif

	//Set partial upward likelihoods for HLP models
	if(io->mod->whichrealmodel<=HLP17){
		For(i,io->ntrees){
		  if(io->mod->whichrealmodel <= HLP17)
			  Get_UPP(io->tree_s[i]->noeud[io->tree_s[i]->mod->startnode],
				  io->tree_s[i]->noeud[io->tree_s[i]->mod->startnode]->v[0],
				  io->tree_s[i]);
		  }
	}
}




