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


/* Help that needs to be checked: 
 --ts/tv
 --frequencies
 --model
 --nclasses
 
 */


#include "help.h"


/* int  T_MAX_FILE; */
/* phydbl SMALL; */
/* phydbl UNLIKELY; */

/*********************************************************/

void Usage()
{
    
    char *BOLD=(char *)mCalloc(10,sizeof(char));
    char *FLAT=(char *)mCalloc(10,sizeof(char));
    char *LINE=(char *)mCalloc(10,sizeof(char));
    char *cha;
    
    
    cha =getenv("OS");
    
    if(cha!=NULL) 
    {
        strcpy(BOLD, "");
        strcpy(FLAT, "");
        strcpy(LINE, "");
    } 
    else 
    {
        strcpy(BOLD, "\033[00;01m");
        strcpy(FLAT, "\033[00;00m");
        strcpy(LINE, "\033[00;04m");
    }
    
    //! Modified by Louis and Marcelo 25.10.2012
    PhyML_Printf("\n%sIgPhyML %s\n"
    		     "\t%sKB Hoehn, JA Vander Heiden, JQ Zhou, G Lunter, OG Pybus, SH Kleinstein.\n"
    			 "%s\tPlease cite: https://doi.org/10.1073/pnas.1906020116 and 10.1534/genetics.116.196303\n",BOLD,VERSION,BOLD,FLAT);
    PhyML_Printf("%s\n\tFor detail and usage information: https://igphyml.readthedocs.io\n",BOLD);
    

    PhyML_Printf("%s\nUsage:",BOLD);
    PhyML_Printf("%s\n\tigphyml --repfile [lineages file] -m [model] [other options]\n",FLAT);
    PhyML_Printf("%s\n\tigphyml -i [sequence file] -m [model] [other options]\n",FLAT);
    PhyML_Printf("%s\nInput/output options:\n%s",BOLD,FLAT);
    
    PhyML_Printf("\n\t%s--repfile %slineage_file_name%s (required if no -i specified)\n",BOLD,LINE,BOLD);
    PhyML_Printf("\t\t%slineage_file_name%s: File specifying input files for repertoire analysis (see docs).\n",LINE,FLAT);  //! Modified by Louis
    PhyML_Printf("\t\tUse BuildTrees.py (https://changeo.readthedocs.io) to generate lineages file.\n");  //! Modified by Louis
    PhyML_Printf("\n\t%s-i %sseq_file_name%s (required if no --repfile specified)\n",BOLD,LINE,BOLD);
    PhyML_Printf("\t\t%sseq_file_name%s: Codon-aligned sequence file in FASTA or PHYLIP format.\n",LINE,FLAT);  //! Modified by Louis

    PhyML_Printf("%s\n\t--run_id %sID_string%s\n",BOLD,LINE,FLAT);
    PhyML_Printf("\t\t%sAppend the string %sID_string%s at the end of each output file.\n",FLAT,LINE,FLAT);

    PhyML_Printf("%s\n\t-u %suser_tree_file%s (only if -i used)\n",BOLD,LINE,FLAT);
    PhyML_Printf("\t\t%suser_tree_file%s : starting tree filename. The tree must be in Newick format.\n",LINE,FLAT);

    PhyML_Printf("%s\n\t--partfile %spartition_file%s (only if -i used)\n",BOLD,LINE,FLAT);
    PhyML_Printf("\t\t%spartition_file%s : Partition file specifying CDRs/FWRs for sequence file.\n",LINE,FLAT);
    ///////////////////////////////
    PhyML_Printf("%s\nModel options:\n%s",BOLD,FLAT);

    PhyML_Printf("%s\n\t-m %smodel%s %s(required)\n",BOLD,LINE,FLAT,BOLD);
    PhyML_Printf("\t\t%smodel%s : substitution model name.\n",LINE,FLAT);
    PhyML_Printf("\t\t%sCodon%s based models: %sHLP%s (HLP19) | %sGY%s | %sHLP17%s\n",FLAT,FLAT,LINE,FLAT,LINE,FLAT,LINE,FLAT,LINE,FLAT);
    PhyML_Printf("\t\t%sUse GY for quick tree construction.\n\t\tHLP for B cell specific features (see docs).\n",FLAT);
    
    PhyML_Printf("%s\n\tModel parameterization arguments:%s\n",BOLD,LINE,FLAT);
    PhyML_Printf("%s\n\t (!) -t, --omega, and --hotness have the following options.%s\n",BOLD,LINE,FLAT);
    PhyML_Printf("\t\t%sAny combination of these may be specified in a comma separated list\n"
    		"\t\twhen multiple parameters are estimated e.g. --omega ce,1\n",FLAT);
    PhyML_Printf("\t\t%sparameter%s = e: Estimate single value by ML across all lineages (default).\n",LINE,FLAT);
    PhyML_Printf("\t\t%sparameter%s = ce: Same as 'e' but also find 95%% confidence intervals.\n",LINE,FLAT);
    PhyML_Printf("\t\t%sparameter%s = i: Estimate by ML for each lineage individually (experimental).\n",LINE,FLAT);
    PhyML_Printf("\t\t%sparameter%s = ci: Same as 'i' but also find 95%% confidence intervals (experimental).\n",LINE,FLAT);
    PhyML_Printf("\t\t%sparameter%s > -1: Fix parameter to specified value for all lineages (see below).\n",LINE,FLAT);


    PhyML_Printf("%s\n\t-t %sts/tv_ratio%s = [e|ce|i|ci|>0]\n",BOLD,LINE,FLAT);
    PhyML_Printf("\t\tSet the transition/transversion ratio.\n");

    PhyML_Printf("%s\n\t--omega %somega%s = [e|ce|i|ci|>0]\n",BOLD,LINE,FLAT);
    PhyML_Printf("\t\tSet number/value of omegas to estimate.\n");
    PhyML_Printf("\t\tFirst value (0) corresponds to FWRs, second (1) to CDRs.\n");
    PhyML_Printf("\t\tMay specify multiple omegas if partition file(s) specified.\n");

    PhyML_Printf("%s\n\t--hotness %shotness%s = [e|ce|i|ci|>-1]\n",BOLD,LINE,FLAT);
    PhyML_Printf("\t\tSet number hot- and coldspot rates to estimate.\n");
    PhyML_Printf("\t\tMay specify multiple values according to --motifs option.\n");
    PhyML_Printf("\t\t'e,e,e,e,e,e' is default.\n");

    PhyML_Printf("%s\n\t--motifs %smotifs%s\n",BOLD,LINE,FLAT);
    PhyML_Printf("\t\tSpecify hot- and coldspot motifs to be modeled.\n");
    PhyML_Printf("\t\t%smotifs%s = FCH (default) : Free coldspots and hotspots. Estimate separate rates for six canonical motifs.\n",LINE,FLAT);
    PhyML_Printf("\t\tOtherwise, motifs specified by <motif>_<mutable position>:<index_in_hotness>.\n");
    PhyML_Printf("\t\t%smotifs%s = WRC_2:0 | GYW_0:0 | WA_1:0 | TW_0:0 | SYC_2:0 | GRS_0:0 : Model rate specific motif(s).\n",LINE,FLAT,LINE,FLAT,LINE,FLAT,LINE,FLAT,LINE,FLAT,LINE,FLAT,LINE,FLAT);
    PhyML_Printf("\t\te.g. %smotifs%s = WRC_2:0,GYW_0:0 symmetric WR%sC%s/%sG%sYW motifs.\n",LINE,FLAT,LINE,FLAT,LINE,FLAT);
    PhyML_Printf("\t\te.g. %smotifs%s = WRC_2:0,GYW_0:1 asymmetric WR%sC%s/%sG%sYW motifs. Must specify two values in --hotness.\n",LINE,FLAT,LINE,FLAT,LINE,FLAT);
    
    PhyML_Printf("%s\n\t-f (or --frequencies) %sempirical%s, %smodel%s, %soptimized%s, %sfT,fC,fA,fG%s,\n\t\t%sfT1,fC1,fA1,fG1,fT2,fC2,fA2,fG2,fT3,fC3,fA3,fG3%s\n",
                 BOLD,LINE,BOLD,LINE,FLAT,LINE,FLAT,LINE,FLAT,LINE,FLAT,LINE,FLAT);
    PhyML_Printf("\t                       %sor%s %sfC1,fC2, ... ,fC64%s\n",BOLD,FLAT,LINE,FLAT);
    PhyML_Printf("\t\t%sempirical%s: (GY default) the equilibrium codon frequencies are estimated by counting\n"
                 "\t\t the occurence of bases or codons in the alignment.\n",LINE,FLAT);
    PhyML_Printf("\t\t%soptimize%s : (HLP17 default) codon frequencies are estimated using maximum likelihood\n",LINE,FLAT);
    PhyML_Printf("%s\n\t--fmodel %sfrequency model%s\n",BOLD,LINE,FLAT);
    PhyML_Printf("\t\tWhich frequency model to use.\n");
    PhyML_Printf("\t\t%sfrequency model%s = %sF1XCODONS%s | %sF1X4%s | %sF3X4%s | %sCF3X4%s (default)\n",LINE,FLAT,LINE,FLAT,LINE,FLAT,LINE,FLAT,LINE,FLAT);


    PhyML_Printf("%s\nOptimization options:\n%s",BOLD,FLAT);

    PhyML_Printf("%s\n\t-o (or --optimize) %sparams%s\n",BOLD,LINE,FLAT);
    PhyML_Printf("\t\tThis option focuses on specific parameter optimisation.\n");
    PhyML_Printf("\t\t%sparams%s = tlr : (default) tree topology (t), branch length (l) and rate parameters (r) are optimised.\n",LINE,FLAT);
    PhyML_Printf("\t\t%sparams%s = tl  : tree topology and branch length are optimised.\n",LINE,FLAT);
    PhyML_Printf("\t\t%sparams%s = lr  : branch length and rate parameters are optimised.\n",LINE,FLAT);
    PhyML_Printf("\t\t%sparams%s = l   : branch length are optimised.\n",LINE,FLAT);
    PhyML_Printf("\t\t%sparams%s = r   : rate parameters are optimised.\n",LINE,FLAT);
    PhyML_Printf("\t\t%sparams%s = n   : no parameter is optimised.\n",LINE,FLAT);
    
    PhyML_Printf("%s\n\t-s (or --search) %smove%s\n",BOLD,LINE,FLAT);
    PhyML_Printf("\t\tTree topology search operation option.\n");
    PhyML_Printf("\t\tCan be either %sNNI%s (default, fast) or %sSPR%s (thorough, slow).\n",LINE,FLAT,LINE,FLAT,LINE,FLAT);

    PhyML_Printf("%s\n\t--threads %snum_threads%s\n",BOLD,LINE,FLAT);
    PhyML_Printf("\t\tNumber of threads to use for parallelization. Default is 1.\n");

    PhyML_Printf("%s\n\t--minSeq %sminimum_sequences%s\n",BOLD,LINE,FLAT);
    PhyML_Printf("\t\tMinimum number of sequences (including germline) per lineage for inclusion in analysis.\n");

    /*PhyML_Printf("%s\n\t--print_trace%s\n",BOLD,FLAT);
    PhyML_Printf("\t\t%sPrint each phylogeny explored during the tree search process\n",FLAT);
    PhyML_Printf("\t\t%sin file %s*_phyml_trace.txt%s.\n",FLAT,LINE,FLAT);*/

    PhyML_Printf("%s\n\t-h (or --help)%s\n",BOLD,FLAT);
    PhyML_Printf("\t\t%sShow this help message and exit.\n",FLAT);


    exit(0);
}

/*********************************************************/


