/*
 * io.h
 *
 *  Created on: Jan 9, 2018
 *      Author: kenneth
 */

#ifndef SRC_IO_H_
#define SRC_IO_H_

#include "utilities.h"

FILE *Openfile(char *filename,int mode);
int Filexists(char *filename);
void Print_Settings(option *input,model* mod);


//sequence reading methods
align **Get_Seq(option *input, model *mod);
align **Read_Seq_Sequential(option *io, model *mod);
align **Read_Seq_Interleaved(option *io, model *mod);
align **Get_Seq_Phylip(option *input,model* mod);
void Read_Ntax_Len_Phylip(FILE *fp ,int *n_otu, int *n_tax);
int Read_One_Line_Seq(align ***data,int num_otu,FILE *in);

align** Read_Seq_Fasta(option*,model*);
void scanFasta(model*,FILE*);

//tree reading methods
FILE* GetTreeFile(option *io);
t_tree *Read_Tree(char *s_tree);
void R_rtree(char *s_tree_a, char *s_tree_d, t_node *a, t_tree *tree, int *n_int, int *n_ext);
void Read_Branch_Label(char *sub_part, char *full_part, t_edge *b);
void Read_Branch_Length(char *s_d, char *s_a, t_tree *tree);
void Read_Node_Name(t_node *d, char *s_tree_d, t_tree *tree);
void Clean_Multifurcation(char **subtrees,int current_deg,int end_deg);
char **Sub_Trees(char *tree,int *degree);
int Next_Par(char *s,int pos);
char *Write_Tree(t_tree *tree);
void R_wtree(t_node *pere,t_node *fils,int *available,char **s_tree, int *pos,t_tree *tree);


void Print_IgPhyML_Out(option* io);

void Print_Fp_Out(FILE *fp_out, time_t t_beg, time_t t_end, t_tree *tree, option *input, int n_data_set, int num_rand_tree, model* mod);
void Print_Fp_Out_Lines(FILE *fp_out,time_t t_beg,time_t t_end,t_tree *tree,option *input,int n_data_set,model* mod);


#endif /* SRC_IO_H_ */
