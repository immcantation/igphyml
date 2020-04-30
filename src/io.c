/*
 * io.c
 *
 *  Created on: Jan 9, 2018
 *      Author: kenneth
 */

#include "utilities.h"
#include "io.h"

FILE *Openfile(char *filename, int mode)
{
  /* mode = 0 -> read */
  /* mode = 1 -> write */
  /* mode = 2 -> append */

  FILE *fp;
  char *s;
  int open_test=0;

  s = filename;

  fp = NULL;

  switch(mode){
    case 0 :{
      while(!(fp = (FILE *)fopen(s,"r")) && ++open_test<10){
        PhyML_Printf("\n. Can't open file '%s'",s);
        exit(1);
        //Getstring_Stdin(s);
      }
      break;
    }
    case 1 :{
      fp = (FILE *)fopen(s,"w");
      break;
    }
    case 2 :{
      fp = (FILE *)fopen(s,"a");
      break;
    }
    default : break;
  }
  return fp;
}

/*********************************************************/


int Filexists(char *filename)
{
  FILE *fp;
  fp =fopen(filename,"r");
  if (fp) {
    fclose(fp);
    return 1;
  } else
    return 0;
}

/*********************************************************/

void Print_Settings(option *io, model* mod)
{
  int answer;
  char *s,*r;

  s = (char *)mCalloc(100,sizeof(char));
  r = (char *)mCalloc(100,sizeof(char));
  PhyML_Printf("\n\n\n");
  PhyML_Printf("\n\n");
  PhyML_Printf("                               ..................                               \n");
  PhyML_Printf("oooooooooooooooooooooooo       CURRENT   SETTINGS       oooooooooooooooooooooooo\n");
  PhyML_Printf("                               ..................                               \n");

  PhyML_Printf("\n. Sequence filename:\t\t\t\t %s", Basename(mod->in_align_file));


  if(io->datatype == CODON)  strcpy(s,"Codon"); //!< Added by Marcelo.
  else if(io->datatype == NT)       strcpy(s,"DNA");
  else if(io->datatype == AA)        strcpy(s,"Amino Acids");
  else strcpy(s,"generic");

  PhyML_Printf("\n. Data type:\t\t\t\t\t %s",s);
  if(io->datatype == CODON){
    switch(mod->genetic_code){  //was io-> mod Ken 9/1/2018
      case   STANDARD: strcpy(s,"Standard\0");break;
      case       TVMC: strcpy(s,"Vertebrate Mitochondrial\0"); break;
      case       TYMC: strcpy(s,"Yeast Mitochondrial\0"); break;
      case THMPCMCMSC: strcpy(s,"Mold, Protozoan, Coelenterate ...\0"); break;
      case      THIMC: strcpy(s,"Invertebrate Mitochondrial\0"); break;
      case    THCDHNC: strcpy(s,"Ciliate, Dasycladacean ...\0"); break;
      case     THEFMC: strcpy(s,"Echinoderm and Flatworm Mit. ...\0"); break;
      case      THENC: strcpy(s,"Euplotid Nuclear\0"); break;
      case    THBAPPC: strcpy(s,"Bacterial, Archaeal and Plant ...\0"); break;
      case     THAYNC: strcpy(s,"Alternative Yeast Nuclear\0"); break;
      case      THAMC: strcpy(s,"Ascidian Mitochondrial\0"); break;
      case     THAFMC: strcpy(s,"Alternative Flatworm Mit. ...\0"); break;
      case       BLNC: strcpy(s,"Blepharisma Nuclear\0"); break;
      case       CHMC: strcpy(s,"Chlorophycean Mitochondrial\0"); break;
      case       TRMC: strcpy(s,"Trematode Mitochondrial\0"); break;
      case      SCOMC: strcpy(s,"Scenedesmus Obliquus Mit. ...\0"); break;
      case       THMC: strcpy(s,"Thraustochytrium Mitochondrial\0"); break;
      default : Warn_And_Exit("Genetic code not implemented."); break;
        break;
    }
    PhyML_Printf("\n. Alphabet size:\t\t\t\t %d",mod->ns);  //was io-> mod Ken 9/1/2018
    PhyML_Printf("\n. Genetic code:\t\t\t\t\t %s",s);

  }
  else
    PhyML_Printf("\n. Alphabet size:\t\t\t\t %d",mod->ns); //was io-> mod Ken 9/1/2018

  PhyML_Printf("\n. Sequence format:\t\t\t\t %s", io->interleaved ? "Interleaved": "Sequential");
  PhyML_Printf("\n. Number of data sets:\t\t\t\t %d", io->n_data_sets);

  PhyML_Printf("\n. Number of bootstrapped data sets:\t\t %d", mod->bootstrap);  //was io-> mod Ken 9/1/2018

  if (mod->bootstrap > 0)
    PhyML_Printf("\n. Compute aLRT:\t\t\t\t\t No");
  else {
    if(io->ratio_test == 1)
      PhyML_Printf("\n. Compute aLRT:\t\t\t\t\t Yes (aLRT statistics)");
    else if(io->ratio_test == 2)
      PhyML_Printf("\n. Compute aLRT:\t\t\t\t\t Yes (Chi^2-based br. supports)");
    else if(io->ratio_test == 3)
      PhyML_Printf("\n. Compute aLRT:\t\t\t\t\t Yes (Min[SH-like, Chi^2-based]");
    else if(io->ratio_test == 4)
      PhyML_Printf("\n. Compute aLRT:\t\t\t\t\t Yes (SH-like branch supports)");
    else if(io->ratio_test == 5)
      PhyML_Printf("\n. Compute aLRT:\t\t\t\t\t Yes (aBayes branch supports)");

  }
  if(io->datatype == CODON){
    switch(mod->freq_model){ //was io-> mod Ken 9/1/2018
      case F1XSENSECODONS:
        strcpy(s,"F1x");
        sprintf(r,"%d",mod->ns);  //was io-> mod Ken 9/1/2018
        strcat(s,r);
        break;
      case F1X4: strcpy(s,"F1x4");
        break;
      case F3X4: strcpy(s,"F3x4");
        break;
      case CF3X4: strcpy(s,"CF3x4");
        break;
      default: strcpy(s,"\0");break;
    }
    if(mod->whichmodel==MGECMUSR||mod->whichmodel==YAPECMUSR||
    		mod->whichmodel==GYECMUSR ||mod->whichmodel==GYECMK07 ||
			mod->whichmodel==GYECMS05 || mod->whichmodel==YAPECMK07 ||
			mod->whichmodel==YAPECMS05 || mod->whichmodel==MGECMUSRWK||
			mod->whichmodel==YAPECMUSRWK||mod->whichmodel==GYECMUSRWK ||
			mod->whichmodel==GYECMK07WK || mod->whichmodel==GYECMS05WK ||
			mod->whichmodel==YAPECMK07WK || mod->whichmodel==YAPECMS05WK ||
			io->eq_freq_handling == MODEL){
      strcpy(r,"Model-Defined");
    }else{
      if(mod->s_opt->opt_state_freq) strcpy(r,"Optimized");
      else if(mod->s_opt->user_state_freq) strcpy(r,"User-Defined");
      else  strcpy(r,"Empirical");
    }
    PhyML_Printf("\n. Model type:\t\t\t\t\t %s", mod->modelname);
    if(mod->pcaModel) PhyML_Printf("\n. Number of PCs:\t\t\t\t %d", mod->npcs);
    PhyML_Printf("\n. Equilibrium frequencies model:\t\t %s %s", r,s);
    if(io->modeltypeOpt <= HLP17){
        PhyML_Printf("\n. Motifs:\t\t\t\t\t %s", mod->motifstring);
        PhyML_Printf("\n. h parameter(s):\t\t\t\t %s", mod->hotnessstring);
        PhyML_Printf("\n. Root ID:\t\t\t\t\t %s", mod->rootname);
    }
  }else PhyML_Printf("\n. Model name:\t\t\t\t\t %s", mod->modelname);


  if (io->datatype == NT)
	  Warn_And_Exit("Only codon data supported!\n");


  if (mod->s_opt->opt_pinvar)
    PhyML_Printf("\n. Prop. of invariable sites:\t\t\t Estimated");
  else
    PhyML_Printf("\n. Prop. of invariable sites:\t\t\t %.2f", mod->pinvar);

  if(mod->n_w_catg==1 && mod->n_catg>1){ //!< Added By Marcelo.
    PhyML_Printf("\n. Number of subst. rate categs:\t\t\t %d", mod->n_catg);
    if(mod->s_opt->opt_alpha) PhyML_Printf("\n. Gamma distribution parameter:\t\t\t Estimated");
    else PhyML_Printf("\n. Gamma distribution parameter:\t\t\t %f", mod->alpha);
    PhyML_Printf("\n. 'Middle' of each rate class:\t\t\t %s",(mod->gamma_median)?("Median"):("Mean"));
  }else{
    PhyML_Printf("\n. Number of subst. rate categs:\t\t\t 1");
  }

  if(io->datatype==CODON){ //!< Added By Marcelo.
    if(mod->n_w_catg==1) PhyML_Printf("\n. Number of dn/ds rate categs:\t\t\t %d - (M0)", mod->n_w_catg);
    else if(mod->omegaSiteVar==DGAMMAK){
      PhyML_Printf("\n. Number of dn/ds rate ratio categories:\t %d - Discrete Gamma model (M5)", mod->n_w_catg);
      if(mod->s_opt->opt_alpha) PhyML_Printf("\n. Gamma distribution parameter:\t\t\t Estimated");
      else PhyML_Printf("\n. Gamma distribution parameter:\t\t\t %f and %f", mod->alpha, mod->beta);
      PhyML_Printf("\n. 'Middle' of each rate class:\t\t\t %s",(mod->gamma_median)?("Median"):("Mean"));
    }
    else PhyML_Printf("\n. Number of dn/ds rate ratio categories:\t %d - Discrete model (M3)", mod->n_w_catg);
  }

  if(io->datatype == AA)Warn_And_Exit("Only codon data supported!\n");

  PhyML_Printf("\n. Optimise tree topology:\t\t\t %s", (mod->s_opt->opt_topo) ? "Yes": "No");

  if(io->datatype==CODON){
    PhyML_Printf("\n. Heuristic during tree search:\t\t\t %s",
    		(io->opt_heuristic_manuel) ? "Modified": "Original");
    PhyML_Printf("\n. Parameter optimization strategy:\t\t %s",
    		(mod->s_opt->opt_method) ? "Single variate (BRENT)": "Multi-variate (BFGS)");
  }

   switch(io->in_tree){
	case 0: { strcpy(s,"BioNJ");     break; }
	case 1: { strcpy(s,"parsimony"); break; }
	case 2: { strcpy(s,"user tree"); break; }
	default: break;
   }

   if(io->datatype==CODON){
      if((io->in_tree==2 && io->nobl) || io->in_tree!=2){
    	  switch(io->init_DistanceTreeCD){
    	  case NUCLEO: strcat(s," + JC69\0");break;
    	  case KOSI07: strcat(s," + GYECMK07\0");break;
    	  case SCHN05: strcat(s," + GYECMS05\0");break;
    	  case ECMUSR: strcat(s," + ECMUSR\0");break;
    	  default: break;
    	  }
      }
  }

  if(mod->s_opt->opt_topo){
    if(mod->s_opt->topo_search == NNI_MOVE) PhyML_Printf("\n. Tree topology search:\t\t\t\t NNIs");
    else if(mod->s_opt->topo_search == SPR_MOVE) PhyML_Printf("\n. Tree topology search:\t\t\t\t SPRs");
    else if(mod->s_opt->topo_search == BEST_OF_NNI_AND_SPR) PhyML_Printf("\n. Tree topology search:\t\t\t\t Best of NNIs and SPRs");

    PhyML_Printf("\n. Starting tree:\t\t\t\t %s",s);

    PhyML_Printf("\n. Add random input tree:\t\t\t %s", (mod->s_opt->random_input_tree) ? "Yes": "No");
    if(mod->s_opt->random_input_tree)
      PhyML_Printf("\n. Number of random starting trees:\t %d", mod->s_opt->n_rand_starts);
  }else{
    if(!mod->s_opt->random_input_tree)
      PhyML_Printf("\n. Evaluated tree:\t\t\t\t %s",s);
  }
  PhyML_Printf("\n. Optimise branch lengths:\t\t\t %s", (mod->s_opt->opt_bl) ? "Yes": "No");

  answer = 0;
  if(mod->s_opt->opt_alpha      ||
     mod->s_opt->opt_kappa      ||
     mod->s_opt->opt_lambda     ||
     mod->s_opt->opt_pinvar     ||
     mod->s_opt->opt_rr         ||
     mod->s_opt->opt_omega      || //!< Added by Marcelo.
     mod->s_opt->opt_beta       || //!< Added by Marcelo.
     mod->s_opt->opt_prob_omega || //!< Added by Marcelo.
     mod->s_opt->opt_alphaCD) answer = 1;

  PhyML_Printf("\n. Optimise subst. model params:\t\t\t %s", (answer) ? "Yes": "No");
  if(io->expm==SSPADE){
    PhyML_Printf("\n. Likelihood o.f. for subst. mod. param. opt.:\t %s", "Pade approximation\0");
  }else{
    if(answer) PhyML_Printf("\n. Matrix exponential:\t\t\t\t %s", io->heuristicExpm ? "Taylor approximation": "Eigenvalue problem");
  }
  PhyML_Printf("\n. Run ID:\t\t\t\t\t %s", (io->append_run_ID) ? (io->run_id_string): ("None"));
  PhyML_Printf("\n. Random seed:\t\t\t\t\t %d", io->r_seed);

  if(io->datatype==CODON){
#ifdef BLAS_OMP
    PhyML_Printf("\n. Code optimization:\t\t\t\t BLAS, LAPACK and OpenMP");
#elif defined BLAS
    PhyML_Printf("\n. Code optimization:\t\t\t\t BLAS and LAPACK");
#elif defined OMP
    PhyML_Printf("\n. Code optimization:\t\t\t\t OpenMP");
#else
    PhyML_Printf("\n. Code optimization:\t\t\t\t None");
#endif
  }
  PhyML_Printf("\n. Version:\t\t\t\t\t %s", VERSION);
  PhyML_Printf("\n\noooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo\n");
  PhyML_Printf("\n\n");
  fflush(NULL);
  free(s);free(r);
}

/*********************************************************/


align **Get_Seq(option *io, model *mod)
{
  mod->data = NULL;

  Detect_Align_File_Format(io,mod);

  switch(io->data_file_format){
    case PHYLIP:{
      mod->data = Get_Seq_Phylip(io,mod);
      break;
    }
    case NEXUS:{
        Warn_And_Exit("\n. Err: NEXUS cannot be used. Please switch to Phylip format.\n");
        break;
    }
    case FASTA:{
    	scanFasta(mod,mod->fp_in_align);
    	if(mod->n_otu == 0 || mod->init_len == 0){
    		char* tmp = (char *) mCalloc (T_MAX_FILE, sizeof(char));
    		strcpy(tmp, "\n. The file '");
    		strcat(tmp, mod->in_align_file);
    		strcat(tmp, "' is not properly formatted. Has ");
    		sprintf(tmp,"%d",mod->n_otu);
    		strcat(tmp, "sequences of length ");
    		sprintf(tmp,"%d",mod->init_len);
    		strcat(tmp, "\n");
    		Warn_And_Exit(tmp);
    	}
    	mod->data=Read_Seq_Fasta(io,mod);
    	break;
    }
    default:{
      PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
      Warn_And_Exit("");
      break;
    }
  }
  if(mod->n_otu < 3) {
	  printf("Need at least 3 sequences to analyze\n");
	  return NULL;
  }
  if(!mod->data){
    PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
    Warn_And_Exit("");
  }else{
    int i,j;
    char **buff;
    int *remove;
    int n_unkn,pos;
    buff = (char **)mCalloc(mod->n_otu,sizeof(char *));
    For(i,mod->n_otu) buff[i] = (char *)mCalloc(mod->data[0]->len,sizeof(char));
    For(i,mod->data[0]->len){
      For(j,mod->n_otu){
        if((mod->data[j]->state[i] == '?') || (mod->data[j]->state[i] == '-')) mod->data[j]->state[i] = 'X';
        if((io->datatype == NT) && (mod->data[j]->state[i] == 'N')) mod->data[j]->state[i] = 'X';
        if(mod->data[j]->state[i] == 'U') mod->data[j]->state[i] = 'T';
      }
      n_unkn = 0;
      For(j,mod->n_otu) if(mod->data[j]->state[i] == 'X') n_unkn++;

      For(j,mod->n_otu) buff[j][i] = mod->data[j]->state[i];
    }
    pos = 0;
    For(i,mod->data[0]->len){
      For(j,mod->n_otu) mod->data[j]->state[pos] = buff[j][i];
      pos++;
    }
    For(i,mod->n_otu) mod->data[i]->len = pos;
    For(i,mod->n_otu) free(buff[i]);
    free(buff);
    free(remove);
  }
  //Declare number of edges for tree
  mod->nedges=2*mod->n_otu-3;
  return mod->data;
}

/**********************************************************/
// Scan across FASTA file to determine number and length of sequences
void scanFasta(model *mod, FILE* f){
	rewind(f);
	int nseq=0;
	int seql=0;
	char c;
	int maxnl=-1;
	int firstseql=-1;
	int seq=0;
	int id=0;
	int otu=-1;
	int idpos=0;
	do{
	  c=(char)fgetc(mod->fp_in_align);
	  if(c == '>' || c == EOF){
		  id=1;
		  otu++;
		  if(otu==1)firstseql=seql;
		  if(firstseql != seql && otu > 0){
			  printf("\n\nSome sequences in fasta file are not the same length!"
					  " %d vs %d.\nFile: %s\n",firstseql,seql,mod->in_align_file);
			  exit(EXIT_FAILURE);
		  }
		  if(idpos > T_MAX_NAME){
			  printf("\n\nAn ID is too long. Maximum sequence name is %d."
					  "\nFile: %s\n",T_MAX_NAME,mod->in_align_file);
			  exit(EXIT_FAILURE);
		  }
		  idpos=0;
		  seql=0;
	  }else if(c == '\n' && id){
		  id=0;
	  }
	  if(c == EOF)break;
	  if(id){
		  if(c != '>') idpos++;
	  }else if(c!='\n')seql++;
	}while(1);
	mod->n_otu=otu;
	mod->init_len=firstseql;
	rewind(f);
}

/**********************************************************/

align **Read_Seq_Fasta(option *io, model *mod){
  int i,pos;
  char *line;
  align **data;
  char c;
  data   = (align **)mCalloc(mod->n_otu,sizeof(align *));
  For(i,mod->n_otu){
     data[i]        = (align *)mCalloc(1,sizeof(align));
     data[i]->name  = (char *)mCalloc(T_MAX_NAME,sizeof(char));
     data[i]->state = (char *)mCalloc(mod->init_len*mod->state_len+1,sizeof(char));
     data[i]->is_ambigu = NULL;
     data[i]->len = 0;
  }

  int seq=0;
  int id=0;
  int otu=-1;
  int idpos=0;
  do{
	  c=(char)fgetc(mod->fp_in_align);
	  if(c == '>'){
		  id=1;
		  otu++;
		  idpos=0;
	  }else if(c == '\n' && id){
		  id=0;
	  }else if(c == EOF)break;

	  if(id){
		  if(c != '>'){
			  data[otu]->name[idpos]=c;
			  idpos++;
		  }
	  }else if(c!='\n' && otu >= 0){
		  Uppercase(&c);
		  data[otu]->state[data[otu]->len++]=c;
	  }
  }while(1);
  return data;
}

/**********************************************************/


align **Read_Seq_Interleaved(option *io, model *mod)
{
  int i,end,num_block;
  char *line;
  align **data;//!< array of align pointers.
  char *format;

  line   = (char *)mCalloc(T_MAX_LINE,sizeof(char));
  format = (char *)mCalloc(T_MAX_NAME, sizeof(char));
  data   = (align **)mCalloc(mod->n_otu,sizeof(align *));

  end = 0;
  For(i,mod->n_otu){
    data[i]        = (align *)mCalloc(1,sizeof(align));
    data[i]->name  = (char *)mCalloc(T_MAX_NAME,sizeof(char));
    data[i]->state = (char *)mCalloc(mod->init_len*mod->state_len+1,sizeof(char));

    data[i]->len       = 0;
    data[i]->is_ambigu = NULL;

    sprintf(format, "%%%ds", T_MAX_NAME);

    if(!fscanf(mod->fp_in_align,format,data[i]->name)) Exit("\n");

    if(!Read_One_Line_Seq(&data,i,mod->fp_in_align)){
      end = 1;
      if((i != mod->n_otu) && (i != mod->n_otu-1)){
	      PhyML_Printf("\n. Err: Problem with species %s's sequence.\n",data[i]->name);
	      PhyML_Printf("\n. Observed sequence length: %d, expected length: %d\n",data[i]->len, mod->init_len * mod->state_len);
	      Warn_And_Exit("");
	    }
      break;
    }
  }

  if(data[0]->len == mod->init_len * mod->state_len) end = 1;

  if(!end){
    end = 0;
    num_block = 1;
    do{
      num_block++;
      /* interblock */
      if(!fgets(line,T_MAX_LINE,mod->fp_in_align)) break;
      if(line[0] != 13 && line[0] != 10){
	      PhyML_Printf("\n. One or more missing sequences in block %d.\n",num_block-1);
	      Warn_And_Exit("");
	  }
      For(i,mod->n_otu) if(data[i]->len != mod->init_len * mod->state_len) break;
      if(i == mod->n_otu) break;
      For(i,mod->n_otu){
	      if(data[i]->len > mod->init_len * mod->state_len){
          PhyML_Printf("\n. Observed length=%d expected length=%d.\n",data[i]->len,mod->init_len * mod->state_len);
          PhyML_Printf("\n. Err: Problem with species %s's sequence.\n",data[i]->name);
          Warn_And_Exit("");
        }else if(!Read_One_Line_Seq(&data,i,mod->fp_in_align)){
          end = 1;
          if((i != mod->n_otu) && (i != mod->n_otu-1)){
            PhyML_Printf("\n. Err: Problem with species %s's sequence.\n",data[i]->name);
            PhyML_Printf("\n. Observed sequence length: %d, expected length: %d.\n",data[i]->len, mod->init_len * mod->state_len);
            Warn_And_Exit("");
          }
          break;
       }
	 }
    }while(!end);
  }

  For(i,mod->n_otu) data[i]->state[data[i]->len] = '\0';

  For(i,mod->n_otu){
    if(data[i]->len != mod->init_len * mod->state_len){
      PhyML_Printf("\n. Check sequence '%s' length...\n",data[i]->name);
      Warn_And_Exit("");
    }
  }

  free(format);
  free(line);

  return data;
}

/*********************************************************/

align **Read_Seq_Sequential(option *io, model *mod)
{
  int i;
  char *line;
  align **data;
  /*   char c; */
  char *format;


  format = (char *)mCalloc(T_MAX_NAME,sizeof(char));
  line   = (char *)mCalloc(T_MAX_LINE,sizeof(char));
  data   = (align **)mCalloc(mod->n_otu,sizeof(align *));

  For(i,mod->n_otu){
    data[i]        = (align *)mCalloc(1,sizeof(align));
    data[i]->name  = (char *)mCalloc(T_MAX_NAME,sizeof(char));
    data[i]->state = (char *)mCalloc(mod->init_len*mod->state_len+1,sizeof(char));

    data[i]->is_ambigu = NULL;
    data[i]->len = 0;

    sprintf(format, "%%%ds", T_MAX_NAME);
    if(!fscanf(mod->fp_in_align,format,data[i]->name)) Exit("\n");

    while(data[i]->len < mod->init_len * mod->state_len) Read_One_Line_Seq(&data,i,mod->fp_in_align);

    if(data[i]->len != mod->init_len * mod->state_len){
      PhyML_Printf("\n. Err: Problem with species %s's sequence (check the format).\n",data[i]->name);
      PhyML_Printf("\n. Observed sequence length: %d, expected length: %d\n",data[i]->len, mod->init_len * mod->state_len);
      Warn_And_Exit("");
    }
  }

  For(i,mod->n_otu) data[i]->state[data[i]->len] = '\0';

  free(format);
  free(line);

  return data;
}

/*********************************************************/

align **Get_Seq_Phylip(option *io, model* mod)
{
  Read_Ntax_Len_Phylip(mod->fp_in_align,&mod->n_otu,&mod->init_len);

  if(io->interleaved) mod->data = Read_Seq_Interleaved(io,mod);
  else                mod->data = Read_Seq_Sequential(io,mod);

  return mod->data;
}

/*********************************************************/

void Read_Ntax_Len_Phylip(FILE *fp ,int *n_otu, int *n_tax)
{
  char *line;
  int readok;

  line = (char *)mCalloc(T_MAX_LINE,sizeof(char));

  readok = 0;
  do{
    if(fscanf(fp,"%s",line) == EOF){
      free(line);
      PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
      Warn_And_Exit("");
    }else{
      if(strcmp(line,"\n") && strcmp(line,"\r") && strcmp(line,"\t")){
	      sscanf(line,"%d",n_otu);
	      if(*n_otu <= 0) Warn_And_Exit("\n. The number of taxa cannot be negative.\n");

	      if(!fscanf(fp,"%s",line)) Exit("\n");
	      sscanf(line,"%d",n_tax);
	      if(*n_tax <= 0) Warn_And_Exit("\n. The sequence length cannot be negative.\n");
	      else readok = 1;
	    }
    }
  }while(!readok);

  free(line);
}

/*********************************************************/

int Read_One_Line_Seq(align ***data, int num_otu, FILE *in)
{
  char c = ' ';
  int nchar = 0;

  while(1){
    if((c == 13) || (c == 10)){
      if(!nchar){
	      c=(char)fgetc(in);
	      continue;
	  }else{
	      break;
	  }
    }else if(c == EOF){
      break;
    }else if((c == ' ') || (c == '\t') || (c == 32)){
      c=(char)fgetc(in);
      continue;
    }

    nchar++;
    Uppercase(&c);

    if(c == '.'){
      c = (*data)[0]->state[(*data)[num_otu]->len];
      if(!num_otu)
        Warn_And_Exit("\n. Err: Symbol \".\" should not appear in the first sequence\n");
    }
    (*data)[num_otu]->state[(*data)[num_otu]->len]=c;
    (*data)[num_otu]->len++;
    c = (char)fgetc(in);
    if(c == ';') break;
  }

  if(c == EOF) return 0;
  else return 1;
}

/*********************************************************/

FILE* GetTreeFile(option *io) {
    FILE *treefile = NULL;
    char *filename;
    char *filenr;
    if(io->logtree == 0) {
        Warn_And_Exit("(in GetTreeFile), should not happen");
    } else if(io->logtree == 1) {
        treefile = openOutputFile(io->out_boot_tree_file, "_igphyml_logtree", ".txt", io);
    } else if(io->logtree == 2) {
        filename = mCalloc(30, sizeof(char));
        strcat(filename, "_igphyml_logtree");
        filenr = mCalloc(8, sizeof(char));
        sprintf(filenr, "_%d", io->treecounter);
        strcat(filename, filenr);
        treefile = openOutputFile(io->out_boot_tree_file, filename, ".txt", io);
        io->treecounter++;
        free(filename);
        free(filenr);
    }
    return(treefile);
}

/*********************************************************/

t_tree *Read_Tree(char *s_tree){
  char **subs;
  int i,n_ext,n_int,n_otu;
  t_tree *tree;
  int degree;
  t_node *root_node;

  n_int = n_ext = 0;

  n_otu=0;
  For(i,(int)strlen(s_tree)) if(s_tree[i] == ',') n_otu++;
  n_otu+=1;

  tree = (t_tree *)Make_Tree(n_otu);
  Init_Tree(tree,tree->n_otu);
  Make_All_Tree_Nodes(tree);
  Make_All_Tree_Edges(tree);
  Make_Tree_Path(tree);

  subs = Sub_Trees(s_tree,&degree);
  Clean_Multifurcation(subs,degree,3);
  if(degree == 2){
    root_node      = tree->noeud[2*n_otu-2];
    root_node->num = 2*n_otu-2;
    tree->n_root   = root_node;
    n_int         -= 1;
  }else{
    root_node      = tree->noeud[n_otu];
    root_node->num = n_otu;
    tree->n_root   = NULL;
  }

  root_node->tax = 0;

  tree->has_branch_lengths = 0;
  tree->num_curr_branch_available = 0;
  For(i,degree) R_rtree(s_tree,subs[i],root_node,tree,&n_int,&n_ext);

  if(tree->n_root){
    tree->e_root = tree->t_edges[tree->num_curr_branch_available];

    tree->n_root->v[0]->v[0] = tree->n_root->v[1];
    tree->n_root->v[1]->v[0] = tree->n_root->v[0];

    Connect_One_Edge_To_Two_Nodes(tree->n_root->v[0],
                                  tree->n_root->v[1],
                                  tree->e_root,
                                  tree);
    tree->e_root->l = tree->n_root->l[0] + tree->n_root->l[1];
    tree->n_root_pos = tree->n_root->l[0] / tree->e_root->l;
  }

  For(i,NODE_DEG_MAX) free(subs[i]);
  free(subs);
  return tree;
}

/*********************************************************/
/* 'a' in t_node a stands for ancestor. 'd' stands for descendant */
void R_rtree(char *s_tree_a, char *s_tree_d, t_node *a, t_tree *tree, int *n_int, int *n_ext)
{
  int i;
  t_node *d;
  int n_otu = tree->n_otu;

  if(strstr(s_tree_a," ")){
    PhyML_Printf("\n. [%s]",s_tree_a);
    Warn_And_Exit("\n. Err: the tree must not contain a ' ' character\n");
  }

  if(s_tree_d[0] == '('){
    char **subs;
    int degree;

    (*n_int)+=1;
    d      = tree->noeud[n_otu+*n_int];
    d->num = n_otu+*n_int;
    d->tax = 0;

    Read_Branch_Label(s_tree_d,s_tree_a,tree->t_edges[tree->num_curr_branch_available]);
    Read_Branch_Length(s_tree_d,s_tree_a,tree);

    For(i,3){
      if(!a->v[i]){
	      a->v[i]=d;
	      d->l[0]=tree->t_edges[tree->num_curr_branch_available]->l;
	      a->l[i]=tree->t_edges[tree->num_curr_branch_available]->l;
	      break;
	    }
    }
    d->v[0]=a;

    if(a != tree->n_root){
      Connect_One_Edge_To_Two_Nodes(a,d,tree->t_edges[tree->num_curr_branch_available],tree);
      tree->num_curr_branch_available++;
    }

    subs=Sub_Trees(s_tree_d,&degree);
    Clean_Multifurcation(subs,degree,2);
    R_rtree(s_tree_d,subs[0],d,tree,n_int,n_ext);
    R_rtree(s_tree_d,subs[1],d,tree,n_int,n_ext);
    For(i,NODE_DEG_MAX) free(subs[i]);
    free(subs);
  }else{
    int i;
    d      = tree->noeud[*n_ext];
    d->tax = 1;
    Read_Node_Name(d,s_tree_d,tree);
    Read_Branch_Label(s_tree_d,s_tree_a,tree->t_edges[tree->num_curr_branch_available]);
    Read_Branch_Length(s_tree_d,s_tree_a,tree);

    For(i,3){
      if(!a->v[i]){
        a->v[i]=d;
        d->l[0]=tree->t_edges[tree->num_curr_branch_available]->l;
        a->l[i]=tree->t_edges[tree->num_curr_branch_available]->l;
        break;
      }
    }
    d->v[0]=a;

    if(a != tree->n_root){
      Connect_One_Edge_To_Two_Nodes(a,d,tree->t_edges[tree->num_curr_branch_available],tree);
      tree->num_curr_branch_available++;
    }

    d->num=*n_ext;
    (*n_ext)+=1;
  }
}

/*********************************************************/

void Read_Branch_Label(char *s_d, char *s_a, t_edge *b)
{
  char *sub_tp;
  char *p;
  int i,pos;

  sub_tp = (char *)mCalloc(3+(int)strlen(s_d)+1,sizeof(char));

  sub_tp[0] = '(';
  sub_tp[1] = '\0';
  strcat(sub_tp,s_d);
  strcat(sub_tp,"#");
  p = strstr(s_a,sub_tp);
  if(!p){
    sub_tp[0] = ',';
    sub_tp[1] = '\0';
    strcat(sub_tp,s_d);
    strcat(sub_tp,"#");
    p = strstr(s_a,sub_tp);
  }

  i = 0;
  b->n_labels = 0;
  if(p){
    if(!(b->n_labels%BLOCK_LABELS)) Make_New_Edge_Label(b);
    b->n_labels++;

    pos = 0;
    do{
      b->labels[b->n_labels-1][pos] = p[i+strlen(s_d)+1];
      i++;
      pos++;
      if(p[i+strlen(s_d)+1] == '#'){
	      b->labels[b->n_labels-1][pos] = '\0';
	      b->n_labels++;
	      if(!(b->n_labels%BLOCK_LABELS)) Make_New_Edge_Label(b);
	      i++;
	      pos=0;
	    }
    }
    while((p[i+strlen(s_d)+1] != ':') &&
          (p[i+strlen(s_d)+1] != ',') &&
          (p[i+strlen(s_d)+1] != '('));

    b->labels[b->n_labels-1][pos] = '\0';
  }

  if(p){
    if(b->n_labels == 1)
      PhyML_Printf("\n. Found label '%s' on t_edge %3d.",b->labels[0],b->num);
    else{
      PhyML_Printf("\n. Found labels ");
      For(i,b->n_labels) PhyML_Printf("'%s' ",b->labels[i]);
      PhyML_Printf("on t_edge %3d.",b->num);
    }
  }
  free(sub_tp);
}

/*********************************************************/

void Read_Branch_Length(char *s_d, char *s_a, t_tree *tree)
{
  char *sub_tp;
  char *p;
  t_edge *b;
  int i;

  b = tree->t_edges[tree->num_curr_branch_available];

  sub_tp = (char *)mCalloc(4+(int)strlen(s_d)+1,sizeof(char));

  For(i,b->n_labels){
    strcat(s_d,"#");
    strcat(s_d,b->labels[i]);
  }

  sub_tp[0] = '(';
  sub_tp[1] = '\0';
  strcat(sub_tp,s_d);
  strcat(sub_tp,":");
  p = strstr(s_a,sub_tp);
  if(!p){
    sub_tp[0] = ',';
    sub_tp[1] = '\0';
    strcat(sub_tp,s_d);
    strcat(sub_tp,":");
    p = strstr(s_a,sub_tp);
  }

  if(p){
    b->l = atof((char *)p+(int)strlen(sub_tp));
    tree->has_branch_lengths = 1;
  }
  free(sub_tp);
}

/*********************************************************/

void Read_Node_Name(t_node *d, char *s_tree_d, t_tree *tree)
{
  int i;

  if(!tree->t_edges[tree->num_curr_branch_available]->n_labels){
    strcpy(d->name,s_tree_d);
  }else{
    i = 0;
    do{
      d->name[i] = s_tree_d[i];
      i++;
    }
    while(s_tree_d[i] != '#');
    d->name[i] = '\0';
  }
  strcpy(d->ori_name , d->name);//!< Changed by Marcelo.

}
/*********************************************************/

void Clean_Multifurcation(char **subtrees, int current_deg, int end_deg)
{

  if(current_deg <= end_deg) return;
  else{
    char *s_tmp;
    int i;

    s_tmp = (char *)mCalloc(6+
                            (int)strlen(subtrees[0])+1+
                            (int)strlen(subtrees[1])+1,
                            sizeof(char));

    strcat(s_tmp,"(\0");
    strcat(s_tmp,subtrees[0]);
    strcat(s_tmp,",\0");
    strcat(s_tmp,subtrees[1]);
    strcat(s_tmp,")\0");
    free(subtrees[0]);
    subtrees[0] = s_tmp;

    for(i=1;i<current_deg-1;i++) strcpy(subtrees[i],subtrees[i+1]);

    Clean_Multifurcation(subtrees,current_deg-1,end_deg);
  }
}

/*********************************************************/

char **Sub_Trees(char *tree, int *degree){
  char **subs;
  int posbeg,posend;
  int i;

  if(tree[0] != '(') {*degree = 1; return NULL;}

  subs=(char **)mCalloc(NODE_DEG_MAX,sizeof(char *));

  For(i,NODE_DEG_MAX) subs[i]=(char *)mCalloc((int)strlen(tree)+1,sizeof(char));

  posbeg=posend=1;
  (*degree)=0;
  do{
    posbeg = posend;
    if(tree[posend] != '('){
      while((tree[posend] != ',' ) &&
            (tree[posend] != ':' ) &&
            (tree[posend] != '#' ) &&
            (tree[posend] != ')' )){
	      posend++ ;
	    }
      posend -= 1;
    }
    else posend=Next_Par(tree,posend);

    while((tree[posend+1] != ',') &&
          (tree[posend+1] != ':') &&
          (tree[posend+1] != '#') &&
          (tree[posend+1] != ')')) {posend++;}


    strncpy(subs[(*degree)],tree+posbeg,posend-posbeg+1);
    subs[(*degree)][posend-posbeg+1]='\0'; /* Thanks to Jean-Baka Domelevo-Entfellner */

    posend += 1;
    while((tree[posend] != ',') &&
          (tree[posend] != ')')) {posend++;}
    posend+=1;


    (*degree)++;
    if((*degree) == NODE_DEG_MAX){
      For(i,(*degree))
	    PhyML_Printf("\n. Subtree %d : %s\n",i+1,subs[i]);

      PhyML_Printf("\n. The degree of a t_node cannot be greater than %d\n",NODE_DEG_MAX);
      Warn_And_Exit("\n");
    }
  }
  while(tree[posend-1] != ')');

  return subs;
}


/*********************************************************/
int Next_Par(char *s, int pos)
{
  int curr;

  curr=pos+1;

  while(*(s+curr) != ')'){
    if(*(s+curr) == '(') curr=Next_Par(s,curr);
    curr++;
  }

  return curr;
}

/*********************************************************/

void Print_Tree(FILE *fp, t_tree *tree)
{
  char *s_tree;
  int i;

  s_tree = (char *)Write_Tree(tree);

  if(OUTPUT_TREE_FORMAT == 0) PhyML_Fprintf(fp,"%s\n",s_tree);
  else if(OUTPUT_TREE_FORMAT == 1){
    PhyML_Fprintf(fp,"#NEXUS\n");
    PhyML_Fprintf(fp,"BEGIN TREES;\n");
    PhyML_Fprintf(fp,"\tTRANSLATE\n");
    For(i,tree->n_otu) PhyML_Fprintf(fp,"\t%3d\t%s,\n",i+1,tree->noeud[i]->name);
    PhyML_Fprintf(fp,"\tUTREE PAUP_1=\n");
    PhyML_Fprintf(fp,"%s\n",s_tree);
    PhyML_Fprintf(fp,"ENDBLOCK;");
  }
  free(s_tree);
}

/*********************************************************/

char *Write_Tree(t_tree *tree)
{
  char *s;
  int i,available;
  int pos;

  s=(char *)mCalloc((int)T_MAX_NAME,sizeof(char));
  available = (int)T_MAX_NAME-1;


  s[0]='(';
  pos = 1;

#ifdef PHYML
  tree->n_root = NULL;
  tree->e_root = NULL;
#endif

  if(tree->mod->whichrealmodel <= HLP17){ //added by Ken 16/2/2017 to output as a tree ith rooted branch length of zero
	  t_node* r = Make_Node_Light(-1);
	  tree->noeud[tree->mod->startnode]->v[1] = r;
	  r->pstate=tree->noeud[tree->mod->startnode]->pstate;
	  t_edge* blank = (t_edge *)mCalloc(1,sizeof(t_edge));
	  blank->l=0.00001;
	  blank->left=tree->noeud[tree->mod->startnode];
	  blank->rght=r;
	  tree->noeud[tree->mod->startnode]->b[1] = blank;
	  r->b[0] = blank;
	  r->tax=1;
	  tree->n_root=r;
	  tree->e_root=blank;
	  tree->n_root_pos=tree->mod->startnode;
	  //r->name = (char*)mCalloc(T_MAX_OPTION,sizeof(char*));
	  strcpy(r->name,tree->noeud[tree->mod->startnode]->name);
	  R_wtree(tree->noeud[tree->mod->startnode],tree->noeud[tree->mod->startnode]->v[0],&available,&s,&pos,tree);
	  R_wtree(tree->noeud[tree->mod->startnode],tree->noeud[tree->mod->startnode]->v[1],&available,&s,&pos,tree);
	  /*if(tree->io->mod->ASR){
		  int num = tree->noeud[tree->mod->startnode]->anc_edge->num;
		  tree->mod->mlCodon[num][tree->mod->init_len]='\0';
		  fprintf(tree->io->fp_out_seqs,">%s\n%s\n",
		  	    	tree->mod->rootname,tree->mod->mlCodon[num]);
	  }*/

	  /*free(r->name);
	  free(r);*/
    Free_Node(r);
	  free(blank);
	  tree->noeud[tree->mod->startnode]->b[1]=NULL;
	  tree->noeud[tree->mod->startnode]->v[1]=NULL;
	  tree->n_root=NULL;
	  tree->e_root=NULL;
  }else{
	  if(!tree->n_root){
	  i = 0;
    	while((!tree->noeud[tree->n_otu+i]->v[0]) ||
    		(!tree->noeud[tree->n_otu+i]->v[1]) ||
			(!tree->noeud[tree->n_otu+i]->v[2])) i++;

    	R_wtree(tree->noeud[tree->n_otu+i],tree->noeud[tree->n_otu+i]->v[0],&available,&s,&pos,tree);
    	R_wtree(tree->noeud[tree->n_otu+i],tree->noeud[tree->n_otu+i]->v[1],&available,&s,&pos,tree);
    	R_wtree(tree->noeud[tree->n_otu+i],tree->noeud[tree->n_otu+i]->v[2],&available,&s,&pos,tree);
  	  }else{
  		R_wtree(tree->n_root,tree->n_root->v[0],&available,&s,&pos,tree);
    	R_wtree(tree->n_root,tree->n_root->v[1],&available,&s,&pos,tree);
  	  }
  }

  s[pos-1]=')';
  s[pos]=';';
  s[pos+1]='\0';

  return s;
}

/*********************************************************/
//modified by ken
void R_wtree(t_node *pere, t_node *fils, int *available, char **s_tree, int *pos, t_tree *tree)
{
  int i,p,ori_len;

  p = -1;
  if(fils->tax){
    ori_len = *pos;

    if(OUTPUT_TREE_FORMAT == 0){
      if(tree->io->long_tax_names){
	      strcat(*s_tree,tree->io->long_tax_names[fils->num]);
	      (*pos) += (int)strlen(tree->io->long_tax_names[fils->num]);
	    }else{
	      strcat(*s_tree,fils->name);
	      (*pos) += (int)strlen(fils->name);
	    }
    }else{
      (*pos) += sprintf(*s_tree+*pos,"%d",fils->num+1);
    }
    if((fils->b) && (fils->b[0]) && (fils->b[0]->l > -1.)){
      if(tree->print_labels){
	      if(fils->b[0]->n_labels < 10)
          For(i,fils->b[0]->n_labels){
          //(*pos) += sprintf(*s_tree+*pos,"#%s",fils->b[0]->labels[i]);
	      //(*pos) += sprintf(*s_tree+*pos,"[%d_%s]",tree->mod->num,fils->b[0]->labels[i]);
	    	  int num = atoi(fils->b[0]->labels[i]);
	    	  tree->mod->mlCodon[num][tree->mod->init_len]='\0';
	    	  //(*pos) += sprintf(*s_tree+*pos,"[%s]",tree->mod->mlCodon[num]);
	    	  if(!fils->b[0]->des_node->tax){
	    		  (*pos) += sprintf(*s_tree+*pos,"%d_%s",tree->mod->num,fils->b[0]->labels[i]);
	    	  }
	      //printf("%d\t%lf\th\n",fils->b[0]->num,fils->b[0]->l);
	      }else{
	    	  (*pos) += sprintf(*s_tree+*pos,"#%d_labels",fils->b[0]->n_labels);
	      }
	      if(tree->io->precon && tree->io->precon%10==0){
	    	  (*pos) += sprintf(*s_tree+*pos,"[%s]",tree->chars[fils->pstate]);
	      }
	    }
      strcat(*s_tree,":");
      (*pos)++;

#ifndef MC
      if(!tree->n_root){
	      (*pos) += sprintf(*s_tree+*pos,"%.10f",fils->b[0]->l);
	    }else{
	      if(pere == tree->n_root){
	    	  phydbl root_pos = (fils == tree->n_root->v[0])?(tree->n_root_pos):(1.-tree->n_root_pos);
	    	  (*pos) += sprintf(*s_tree+*pos,"%.10f",tree->e_root->l * root_pos);
	      }else{
          (*pos) += sprintf(*s_tree+*pos,"%.10f",fils->b[0]->l);
	      }
	    }
#else
      if(!tree->n_root){
	      (*pos) += sprintf(*s_tree+*pos,"%.10f",fils->b[0]->l);
	    }else{
	      (*pos) += sprintf(*s_tree+*pos,"%.10f",tree->rates->cur_l[fils->num]);
	    }
#endif
    }
    strcat(*s_tree,",");
    (*pos)++;

    (*available) = (*available) - (*pos - ori_len);

    if(*available < 0){
      PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
      Warn_And_Exit("");
    }

    if(*available < (int)T_MAX_NAME/2){
      (*s_tree) = (char *)mRealloc(*s_tree,*pos+(int)T_MAX_NAME,sizeof(char));
      (*available) = (int)T_MAX_NAME;
    }
  }else{
    (*s_tree)[(*pos)]='(';
    (*s_tree)[(*pos)+1]='\0';
    (*pos)++;
    (*available)--;

    if(tree->n_root){
      For(i,3){
	      if((fils->v[i] != pere) && (fils->b[i] != tree->e_root))
          R_wtree(fils,fils->v[i],available,s_tree,pos,tree);
	      else p=i;
	    }
    }else{
      For(i,3){
	      if(fils->v[i] != pere)
          R_wtree(fils,fils->v[i],available,s_tree,pos,tree);
	      else p=i;
	    }
    }

    ori_len = *pos;

    (*s_tree)[(*pos)-1] = ')';
    (*s_tree)[(*pos)]   = '\0';

    if((fils->b) && (fils->b[0]->l > -1.)){
      if(tree->print_boot_val){
	      (*pos) += sprintf(*s_tree+*pos,"%d",fils->b[p]->bip_score);
	  }else if(tree->print_alrt_val){
	      (*pos) += sprintf(*s_tree+*pos,"%.10f",fils->b[p]->ratio_test);
	  }
     if(tree->print_labels){
	      if(fils->b[p]->n_labels < 10)

        For(i,fils->b[p]->n_labels){
	    	  //(*pos) += sprintf(*s_tree+*pos,"[%d_%s]",tree->mod->num,fils->b[p]->labels[i]);
	    	   int num = atoi(fils->b[p]->labels[i]);
	    	  tree->mod->mlCodon[num][tree->mod->init_len]='\0';
	    	  //(*pos) += sprintf(*s_tree+*pos,"[%s]",tree->mod->mlCodon[num]);
	    	  if(!fils->b[p]->des_node->tax){
	    		  (*pos) += sprintf(*s_tree+*pos,"%d_%s",tree->mod->num,fils->b[p]->labels[i]);
	    	  }
	      }else{
	    	  (*pos) += sprintf(*s_tree+*pos,"#%d_labels",fils->b[p]->n_labels);
	      }
	      if(tree->io->precon && tree->io->precon%10==0){
	    	  (*pos) += sprintf(*s_tree+*pos,"[%s]",tree->chars[fils->pstate]);
	      }
     }

      strcat(*s_tree,":");
      (*pos)++;

#ifndef MC
      if(!tree->n_root){
	      (*pos) += sprintf(*s_tree+*pos,"%.10f",fils->b[p]->l);
	  }else{
	      if(pere == tree->n_root){
	    	  phydbl root_pos = (fils == tree->n_root->v[0])?(tree->n_root_pos):(1.-tree->n_root_pos);
          	  (*pos) += sprintf(*s_tree+*pos,"%.10f",tree->e_root->l * root_pos);
	      }else{
	    	  (*pos) += sprintf(*s_tree+*pos,"%.10f",fils->b[p]->l);
	      }
	  }
#else
      if(!tree->n_root){
	      (*pos) += sprintf(*s_tree+*pos,"%.10f",fils->b[p]->l);
	  }else{
	      (*pos) += sprintf(*s_tree+*pos,"%.10f",tree->rates->cur_l[fils->num]);
	  }
#endif
    }
    strcat(*s_tree,",");
    (*pos)++;
    (*available) = (*available) - (*pos - ori_len);

    if(*available < 0){
      PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
      Warn_And_Exit("");
    }

    if(*available < (int)T_MAX_NAME/2){
      (*s_tree) = (char *)mRealloc(*s_tree,*pos+(int)T_MAX_NAME,sizeof(char));
      (*available) = (int)T_MAX_NAME;
    }
  }
}


/********************************************************/
//Print model parameters
void Print_Params_Tab(model* mod, FILE* f){
	int i;
	if(mod->kappaci){
		fprintf(f,"\t%.4f\t%.4f\t%.4f",mod->kappa,mod->kappalci,mod->kappauci);
	}else if(mod->io->mod->kappaci){
		fprintf(f,"\t%.4f\t%.4f\t%.4f",mod->io->mod->kappa,mod->io->mod->kappalci,mod->io->mod->kappauci);
	}else{
		fprintf(f,"\t%.4f",mod->kappa);
	}
	For(i,mod->nomega_part){
		if(mod->omega_part_ci[i]){
			fprintf(f,"\t%.4f\t%.4f\t%.4f",mod->omega_part[i],mod->omega_part_lci[i],mod->omega_part_uci[i]);
		}else if(mod->io->mod->omega_part_ci[i]){
			fprintf(f,"\t%.4f\t%.4f\t%.4f",mod->io->mod->omega_part[i],mod->io->mod->omega_part_lci[i],mod->io->mod->omega_part_uci[i]);
		}else{
			fprintf(f,"\t%.4f",mod->omega_part[i]);
		}
	}
	For(i,mod->nmotifs){
		if(mod->io->mod->hoptci[mod->io->mod->motif_hotness[i]] == 2){
			fprintf(f,"\t%.4f\t%.4f\t%.4f",mod->hotness[mod->motif_hotness[i]],mod->hoptlci[mod->motif_hotness[i]],mod->hoptuci[mod->motif_hotness[i]]);
		}else if(mod->io->mod->hoptci[mod->io->mod->motif_hotness[i]] == 1){
			fprintf(f,"\t%.4f\t%.4f\t%.4f",mod->io->mod->hotness[mod->motif_hotness[i]],
					mod->io->mod->hoptlci[mod->motif_hotness[i]],mod->io->mod->hoptuci[mod->motif_hotness[i]]);
		}else{
			fprintf(f,"\t%.4f",mod->hotness[mod->motif_hotness[i]]);
		}
	}
}

/********************************************************/
void Print_Param_Header(char* str,int ci, FILE* f){
	if(ci){
		fprintf(f,"\t%s_MLE\t%s_LCI\t%s_UCI",str,str,str);
	}else{
		fprintf(f,"\t%s_MLE",str);
	}
}

/********************************************************/
//Tabular output function
void Print_Tab_Out(option *io){
	FILE* f = io->fp_out_stats;
	int i,j;
	int nseq = 0;
	phydbl treel = 0.0;
	fprintf(f,"CLONE\tNSEQ\tNSITE\tTREE_LENGTH\tLHOOD");
	Print_Param_Header("KAPPA",io->mod->kappaci,f);
	if(io->mod->nomega_part==1)Print_Param_Header("OMEGA",io->mod->omega_part_ci[0],f);
	else if(io->mod->nomega_part==2){
		Print_Param_Header("OMEGA_FWR",io->mod->omega_part_ci[0],f);
		Print_Param_Header("OMEGA_CDR",io->mod->omega_part_ci[1],f);
	}
	For(i,io->mod->nmotifs){
		Print_Param_Header(io->mod->motifs[i],io->mod->hoptci[io->mod->motif_hotness[i]],f);
	}

	fprintf(f,"\tTREE");

	//print repertoire-wide information
	phydbl msite = 0.0;
	For(i,io->ntrees){
		nseq += io->mod_s[i]->n_otu - 1;
		treel += Get_Tree_Size(io->tree_s[i]);
		msite += io->tree_s[i]->n_pattern;
	}
	treel = treel/(io->ntrees*1.0);
	msite = msite/(io->ntrees*1.0);

	fprintf(f,"\nREPERTOIRE\t%d\t%.4f\t%.4f\t%.4f",nseq,msite,treel,io->replnL);
	Print_Params_Tab(io->mod,f);
	fprintf(f,"\t`%s`",io->command); //command is repertoire-wide "tree"

	//print information for all subtrees
	For(i,io->ntrees){
		model* mod = io->mod_s[i];
		char* root = strdup(mod->rootname);
		char* cloneid = mCalloc(T_MAX_LABEL,sizeof(char));
		char* fragment = strsep(&root, "_");
		strcat(cloneid,fragment);
		fragment = strsep(&root, "_");
		while(strcmp(fragment,"GERM") != 0){
			strcat(cloneid,"_");
			strcat(cloneid,fragment);
			fragment = strsep(&root, "_");
		}
		fprintf(f,"\n%s\t%d\t%d\t%.4f\t%.4f",cloneid,mod->n_otu-1,mod->tree->n_pattern,Get_Tree_Size(io->tree_s[i]),mod->tree->c_lnL);
		Print_Params_Tab(mod,f);
		char* s_tree = (char *)Write_Tree(mod->tree);
		fprintf(f,"\t%s",s_tree);
		free(cloneid);
		fflush(f);

	}
	fprintf(f,"\n");

	if(io->mod->ASR == 1){
		//printf("\n%s",io->out_seqs_file);
		FILE* fastaout = io->fp_out_seqs;
		For(i,io->ntrees){
			model* mod = io->mod_s[i];
			t_tree* tree = io->tree_s[i];
			For(j,mod->nedges+1){
				if(j<mod->nedges){
					if(tree->t_edges[j]->des_node->tax)continue;
					if(!tree->t_edges[j]->des_node->name) fprintf(fastaout,">%d_%d\n",mod->num,j);
					else if(strcmp(tree->t_edges[j]->des_node->name,"")!=0)fprintf(fastaout,">%s\n",tree->t_edges[j]->des_node->name);
					else fprintf(fastaout,">%d_%d\n",mod->num,j);
			}else{
				fprintf(fastaout,">%s\n",tree->mod->rootname);
			}
			mod->mlCodon[j][mod->init_len]='\0';
			fprintf(fastaout,"%s\n",mod->mlCodon[j]);
			}
		}
		fclose(fastaout);
	}
}


/********************************************************/
//replacement output function
void Print_IgPhyML_Out(option* io){
	FILE* f = io->fp_out_stats;
	int i,j,k;

	fprintf(f, "oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo\n");
	fprintf(f,"%s%s ---\n", "                      --- IgPhyML ", VERSION);
	fprintf(f, "               Kenneth B Hoehn, Gerton Lunter, Oliver G Pybus\n");
	fprintf(f,"\n\n");
	fprintf(f,"                         Based off of codonPhyML\n");
	fprintf(f,"Marcelo S Zanetti, Stefan Zoller, Manuel Gil, Louis du Plessis, Maria Anisimova\n");
	fprintf(f,"                 http://sourceforge.net/projects/codonphyml/\n");
	fprintf(f,"oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo\n");
	fprintf(f,"                                 Summary\n");
	fprintf(f,"     See doc/IgPhyML_Manual.pdf for further detail on interpreting results\n");
	fprintf(f,"oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo\n");
	fprintf(f,". Command: %s\n",io->command);
  	fprintf(f,". Data sets: %d\n",io->ntrees);
	fprintf(f,". Model name:\t%s\n",io->mod->modelname);
	fprintf(f,". Hotspots:\t%s\n",io->mod->motifstring);
	fprintf(f,". h optimization:\t%s\n",io->mod->hotnessstring);
	if(io->mod->s_opt->opt_topo) {
		switch( io->mod->s_opt->topo_search){
		case NNI_MOVE:
			fprintf(f,". Tree topology search:\tNNIs\n");
			break;
		case SPR_MOVE:
			fprintf(f,". Tree topology search:\tSPRs\n");
			break;
		default:
			Lazy_Exit("Topology option output option not supported",__FILE__,__LINE__);
			break;
		}
	}else{
		fprintf(f,". Tree topology search:\tfixed\n");
	}
	fprintf(f,". Combined log-likelihood: 	%.4lf\n",io->replnL);
	phydbl treel=0.0;
	For(i,io->ntrees)treel+=(Get_Tree_Size(io->tree_s[i])*io->mod_s[i]->init_len/3);
	fprintf(f,". Estimated substitutions in repertoire:\t%.4f\n",treel);

	treel = 0.0;
	For(i,io->ntrees)treel+=Get_Tree_Size(io->tree_s[i]);
	treel = treel/(io->ntrees*1.0);
	fprintf(f,". Mean tree length:\t%.4f\n",treel);

	if(io->mod->optKappa != 2)fprintf(f,". Transition/transversion ratio: 	%.5lf",io->mod->kappa); //!< Kappa
	else fprintf(f,". Transition/transversion ratio: 	%s","Submodel_optimized"); //!< Kappa
	if(io->mod->kappaci)fprintf(f,"\t(%.5f, %.5f)",io->mod->kappalci,io->mod->kappauci);
	fprintf(f,"\n");

	int omegai; //Added by Ken 23/8
		for(omegai=0;omegai<io->mod->nomega_part;omegai++){
		   char* buf = mCalloc(T_MAX_OPTION,sizeof(char));
		   int temp = asprintf(&buf, ". Omega %d", omegai);//, io->mod->partNames[omegai]);//changed by Ken 12/1/2017
		   if(io->mod->omega_part_opt[omegai] != 2)fprintf(f,"%s:\t%.4f",buf,io->mod->omega_part[omegai]);
		   else fprintf(f,"%s:\t%s",buf,"Submodel_optimized");
		   if(io->mod->omega_part_ci[omegai])fprintf(f,"\t(%.4f, %.4f)",io->mod->omega_part_lci[omegai],io->mod->omega_part_uci[omegai]);
		   fprintf(f,"\n");
		}

	if(io->modeltypeOpt<=HLP17){
		fprintf(f,". Hotspot model\t\th_index\toptimized?\th_value\n");
	      int mot;
	      for(mot=0;mot<io->mod->nmotifs;mot++){
	    	  	  if(io->mod->hoptindex[io->mod->motif_hotness[mot]]!=2){
	    	  		  fprintf(f,"\tMotif:\t%s\t\t%d\t\t%d\t\t%.8lf",io->mod->motifs[mot],
	    	  				  io->mod->motif_hotness[mot],io->mod->hoptindex[io->mod->motif_hotness[mot]],
							  io->mod->hotness[io->mod->motif_hotness[mot]]);
	    	  	  }else{
	    	  		  fprintf(f,"\tMotif:\t%s\t\t%d\t\t%d\t\t%s",io->mod->motifs[mot],
	    	  				  io->mod->motif_hotness[mot],io->mod->hoptindex[io->mod->motif_hotness[mot]],"Submodel_optimized");
	    	  	  }
	    	  	  if(io->mod->hoptci[io->mod->motif_hotness[mot]]){
	    	  		  fprintf(f,"\t(%.8lf, %.8lf)",io->mod->hoptlci[io->mod->motif_hotness[mot]],
	    	  				  io->mod->hoptuci[io->mod->motif_hotness[mot]]);
	    	  	  }
	    	  	  fprintf(f,"\n");
	      }
	}
	if(io->mod->freq_model < ROOT){
	fprintf(f,". Nucleotides frequencies: \n");
	fprintf(f,"\tPosition 1:\tf(%s)=%.8lf\tf(%s)=%.8lf\tf(%s)=%.8lf\tf(%s)=%.8lf\n","T1",io->mod->base_freq[0],
			"C1",io->mod->base_freq[1],"A1",io->mod->base_freq[2],"G1",io->mod->base_freq[3]);
	fprintf(f,"\tPosition 2:\tf(%s)=%.8lf\tf(%s)=%.8lf\tf(%s)=%.8lf\tf(%s)=%.8lf\n","T2",io->mod->base_freq[4],
			"C2",io->mod->base_freq[5],"A2",io->mod->base_freq[6],"G2",io->mod->base_freq[7]);
	fprintf(f,"\tPosition 3:\tf(%s)=%.8lf\tf(%s)=%.8lf\tf(%s)=%.8lf\tf(%s)=%.8lf\n","T3",io->mod->base_freq[8],
			"C3",io->mod->base_freq[9],"A3",io->mod->base_freq[10],"G3",io->mod->base_freq[11]);

	fprintf(f,". Codon frequencies\n");
	For(i,64){
		char *buf;
        char s1[4];
		Sprint_codon(s1, i);
		phydbl val1=0.0;
		if(!io->stopCodons[i]){
			val1=io->mod_s[0]->pi[io->indexSenseCodons[i]];
		}
		fprintf(f,"\tf(%s)=%lf",s1,val1);
		if((i-3)%4==0){
			fprintf(f,"\n");
		}
	}
	}
	 div_t hour = div((int)(io->t_current-io->t_beg),3600);
	 div_t min  = div((int)(io->t_current-io->t_beg),60  );
	 min.quot -= hour.quot*60;
	 fprintf(f,". Time used:\t%dh%dm%ds\n",hour.quot,min.quot,(int)(io->t_current-io->t_beg)%60);
	 fprintf(f,". Seconds:\t%d\n",(int)io->t_current-(int)io->t_beg);


	 fprintf(f,"\noooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo\n");
	 fprintf(f,"                                 Submodels\n");
	 fprintf(f,"oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo\n");

	 fprintf(f,"Ind\tSeq\tTreeL\tLikelihood\tKappa");
	 For(omegai,io->mod->nomega_part)fprintf(f,"\tOmega%d",omegai);
	 For(omegai,io->mod->nmotifs)fprintf(f,"\t%s",io->mod->motifs[omegai]);
	 fprintf(f,"\tRootID\tSeq_File\tTree_File\n");
	 For(i,io->ntrees){
		 model* mod=io->mod_s[i];
		 t_tree* tree=io->tree_s[i];
		 fprintf(f,"%d\t%d\t%.2lf\t%.2lf\t%.4lf",i,mod->n_otu,Get_Tree_Size(tree),tree->c_lnL,mod->kappa);
		 if(mod->kappaci==1)fprintf(f,"(%.4lf,%.4lf",mod->kappalci,mod->kappauci);
		 For(omegai,io->mod->nomega_part){
			 fprintf(f,"\t%.4lf",mod->omega_part[omegai]);
			 if(mod->omega_part_ci[omegai]==1)fprintf(f,"(%.4lf,%.4lf)",mod->omega_part_lci[omegai],mod->omega_part_uci[omegai]);
		 }
		 For(omegai,io->mod->nmotifs){
			 fprintf(f,"\t%.4lf",mod->hotness[mod->motif_hotness[omegai]]);
			 if(mod->hoptci[mod->motif_hotness[omegai]]==1)fprintf(f,"(%.4lf,%.4lf)",
					 mod->hoptlci[mod->motif_hotness[omegai]],mod->hoptuci[mod->motif_hotness[omegai]]);
		 }
		 fprintf(f,"\t%s\t%s\t%s\n",io->rootids[i],io->datafs[i],io->treefs[i]);
	 }
	 /*fprintf(f,"oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo\n");
	 fprintf(f,"                                 Subtrees\n");
	 fprintf(f,"oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo\n");
	 For(i,io->ntrees){
		 fprintf(f,"%d\t%s\n",i,Write_Tree(io->tree_s[i]));
	 }*/


	if(io->mod->ASR){
	fprintf(f,"oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo\n");
	fprintf(f,"                     Ancestral Sequence Reconstruction\n");
	fprintf(f,"Marginal ML codon predictions for internal nodes. Open specified tree files using\n");
	fprintf(f,"R phylotate::read_annotated() to view sequence placement on tree.\n");
	fprintf(f,"Codon assignments use ambiguous nucleotides to summarize the minimal set of codons\n");
	fprintf(f,"making up 95%% relative probability. Check out\n");
	fprintf(f,"bioinformatics.org/sms/iupac.html for a reference on nucleotide ambiguity codes.\n");
	fprintf(f,"             See Manual for more information on this section.\n");
	fprintf(f,"oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo");

	  For(i,io->ntrees){
		 t_tree* tree=io->tree_s[i];
		 model* mod = io->mod_s[i];
		 if(io->append_run_ID)fprintf(f,"\n# Lineage %d Reconstructions\n# Sequence file: %s%s%s%s\n",i,io->datafs[i],"_igphyml_asr_",io->run_id_string,".nex");
		 else fprintf(f,"\n# Lineage %d Reconstructions\n# Sequence file: %s%s%s\n",i,io->datafs[i],"_igphyml_asr",".nex");
		 if(io->append_run_ID)fprintf(f,"# Tree file: %s%s%s%s\n",io->datafs[i],"_igphyml_asr_",io->run_id_string,".nex");
		 else fprintf(f,"# Tree file: %s%s\n",io->datafs[i],"_igphyml_asr.nex");

		 //output individual ASR fasta files
		 char fastafile[T_MAX_FILE];
		 strcpy(fastafile,io->datafs[i]);
		 strcat(fastafile,"_igphyml_asr");
		 if(io->append_run_ID){
			 strcat(fastafile,"_");
			 strcat(fastafile,io->run_id_string);
			 strcat(fastafile,".fasta");
		 }else strcat(fastafile,".fasta");
		 char buf[T_MAX_FILE];
		 snprintf(buf, sizeof(buf), "%s", fastafile);
		 FILE* fastaout = fopen(buf,"w");
		 For(j,mod->nedges+1){
			if(j<mod->nedges){
				 if(tree->t_edges[j]->des_node->tax)continue;
			 	 if(!tree->t_edges[j]->des_node->name) fprintf(fastaout,">%d_%d\n",mod->num,j);
			 	 else if(strcmp(tree->t_edges[j]->des_node->name,"")!=0)fprintf(fastaout,">%s\n",tree->t_edges[j]->des_node->name);
			 	 else fprintf(fastaout,">%d_%d\n",mod->num,j);
			}else{
				fprintf(fastaout,">%s\n",tree->mod->rootname);
			}
			mod->mlCodon[j][mod->init_len]='\0';
			fprintf(fastaout,"%s\n",mod->mlCodon[j]);
		   }
		 fclose(fastaout);
		 }
	   }

	 fprintf(f,"\n\n#\tIf you use IgPhyML, please cite:\n");
	 fprintf(f,"#\tK.B. Hoehn, G Lunter, O.G. Pybus\n");
	 fprintf(f,"#\tA phylogenetic codon substitution model for antibody lineages\n");
	 fprintf(f,"#\tGenetics. 2017 May; 206(1): 417–427\n");
	 fprintf(f,"#\n");
	 fprintf(f,"#\tand\n");
	 fprintf(f,"#\n");
	 fprintf(f,"#\tM. Gil, M.S. Zanetti, S. Zoller and M. Anisimova 2013\n");
	 fprintf(f,"#\tCodonPhyML: Fast Maximum Likelihood Phylogeny Estimation under Codon Substitution Models\n");
	 fprintf(f,"#\tMolecular Biology and Evolution, pages 1270-1280, volume 30, number 6\n");
	 fprintf(f,"#\n");
	 fprintf(f,"#\tIf you use aBayes branch supports please cite:\n");
	 fprintf(f,"#\tM. Anisimova, M. Gil, J.F. Dufayard, C. Dessimoz and O. Gascuel 2011\n");
	 fprintf(f,"#\tSurvey of branch support methods demonstrates accuracy, power, and robustness of fast likelihood-based approximation schemes. Syst Biol 60:685-699\n");


	 FILE* orep;
	 if(io->outrepspec){
		 orep=Openfile(io->outrep, 1 );
		 fprintf(orep,"%d\n",io->ntrees);
	 }
	 For(i,io->ntrees){
		 t_tree* tree=io->tree_s[i];
		 char fout[T_MAX_FILE];
		 strcpy(fout,io->datafs[i]);
		 if(io->mod->ASR)strcat(fout,"_igphyml_asr");
		 else strcat(fout,"_igphyml_tree");
		 if(io->append_run_ID){
			 strcat(fout, "_");
			 strcat(fout, io->run_id_string);
		 }
		 if(io->mod->ASR)strcat(fout,".nex");
		 else strcat(fout,".txt");
		 //potentially output new repertoire file
		 if(io->outrepspec)fprintf(orep,"%s\t%s\t%s\t%s\n",io->datafs[i],fout,io->rootids[i],io->partfs[i]);
		 FILE* treeout = Openfile(fout, 1 );
		 if(io->mod->ASR){
			 fprintf(treeout,"#NEXUS\nBegin taxa;\n");
			 fprintf(treeout,"\tDimensions ntax=%d;\n",tree->mod->n_otu);
			 fprintf(treeout,"\tTaxlabels\n");
			 For(j,(tree->mod->n_otu-1)*2){
			 	 if(tree->noeud[j]->name){
				 	 if(strcmp(tree->noeud[j]->name,"")!=0){
					 	 fprintf(treeout,"\t%s\n",tree->noeud[j]->name);
				 	 }
			 	 }
		 	 }
		 	 fprintf(treeout,"\t;\nEnd;\n");
		 	 fprintf(treeout,"Begin trees;Tree TREE1 = [&R] ");
		 }
		 Print_Tree(treeout,io->tree_s[i]);
		 if(io->mod->ASR){
			 fprintf(treeout,"End;");
			 /*fprintf(treeout,"begin figtree;\nset nodeLabels.colorAttribute=\"User selection\";\nset "
					 "nodeLabels.displayAttribute=\"Num\";\nset nodeLabels.fontName=\"sansserif\";"
					 "\nset nodeLabels.fontSize=9;\nset nodeLabels.fontStyle=0;\nset nodeLabels."
					 "isShown=true;\nset nodeLabels.significantDigits=4;\nset trees.order=true;"
					 "\nset trees.orderType=\"increasing\";\nend;");*/
		 }
		 fclose(treeout);
		 if(io->precon){
			 //free(io->mod_s[0]->fp_in_tree);
			//io->mod_s[i]->fp_in_tree = Openfile(fout,0);
		 }
	 }
}
