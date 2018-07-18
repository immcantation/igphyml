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

#ifndef PARS_H
#define PARS_H

#include "utilities.h"
#include "lk.h"
#include "optimiz.h"
#include "models.h"
#include "free.h"

void Make_Tree_4_Pars(t_tree *tree, int n_site);
int  Pars(t_tree *tree);
void Post_Order_Pars(t_node *a, t_node *d, t_tree *tree);
void Pre_Order_Pars(t_node *a, t_node *d, t_tree *tree);
void Get_Partial_Pars(t_tree *tree, t_edge *b_fcus, t_node *a, t_node *d);
void Site_Pars(t_tree *tree);
void Init_Ui_Tips(t_tree *tree);
void Update_P_Pars(t_tree *tree, t_edge *b_fcus, t_node *n);
int Pars_At_Given_Edge(t_edge *b, t_tree *tree);
void Get_All_Partial_Pars(t_tree *tree, t_edge *b_fcus, t_node *a, t_node *d);
int Update_Pars_At_Given_Edge(t_edge *b_fcus, t_tree *tree);
void Init_P_Pars_Tips(t_tree *tree);
void Get_Step_Mat(t_tree *tree);
int Pars_Core(t_edge *b, t_tree *tree);
int One_Pars_Step(t_edge *b,t_tree *tree);
void Init_Class_Tips(t_tree* tree);
int Fill_Sankoff(t_node* d,t_tree* tree, int root);
void Get_First_Path(t_node* d,int index,t_tree* tree, int root);
void Get_Rand_Path(t_node* d,int index,t_tree* tree, int root);
void printTreeState(t_node *d, t_tree* tree, int root);
void Clean_Tree(t_tree* tree);
void Set_Pars_Counters(t_node *d, t_tree *tree,int root);
void Copy_Sankoff_Tree(t_tree* tree1,t_tree* tree2);
int Get_All_Paths(t_node *d, int index, t_tree *tree, t_tree** btrees, int root,int maxtrees,int treeindex,int repindex);
int NNI_Pars_Edge(t_tree *tree, t_edge *b_fcus, int do_swap);
int Resolve_Polytomies_Pars(t_tree* tree,phydbl minbl);
int NNI_ParsSwaps(t_node *a, t_node *b, t_node *c, t_node *d, t_tree *tree);
int NNI_Pars_Search(t_node *c, t_node *d,t_edge* c_fcus,t_edge* d_fcus, int pars0, phydbl thresh,t_tree* tree);
void Get_Pars_Stats(t_tree** trees, int ntrees, int index, FILE* out);
void Fill_Pars_Stats(t_node* d,t_tree* tree, phydbl* switches, phydbl* classl, int root);
#endif





