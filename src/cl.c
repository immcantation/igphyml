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
#include "cl.h"
#include "io.h"

extern struct option longopts[];

/*********************************************************/
/**
 * Fill the Option fields, with the argc array
 */ 
void Read_Command_Line( option *io, int argc, char **argv ){
    int c,i,j,switchResult;

    io->command=mCalloc(T_MAX_LINE,sizeof(char));
    For(i,argc){
    	strcat(io->command,argv[i]);
    	strcat(io->command," ");
    }
    printf("COMMAND: %s\n",io->command);

    while((c = getopt_long_only(argc,argv,"qi:d:g:m:b:n:w:t:f:v:c:a:u:ho:s:p",longopts,NULL)) != -1){
        if(c == 139 || c == 143) {
        	Warn_And_Exit("\n\nIgPhyML Doesn't take config files right now.\n\n");
        }else{
            // command line interface
            switchResult = mainOptionSwitch( c, optarg, io );
            if( c == '?' ) {
                if( isprint( optopt ) ) {
                    PhyML_Printf("\n. Unknown option `-%c'.\n", optopt );
                } else {
                    PhyML_Printf("\n. Unknown option character `\\x%x'.\n", optopt );
                }
                Warn_And_Exit( "" );
            }
            
            if( switchResult != 0 ) {
                Usage();
            }
        }
    }

    finishOptions( io );
	return;
}

/*********************************************************/
