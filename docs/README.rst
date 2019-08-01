IgPhyML - B cell phylogenetic inference package
================================================================================

IgPhyML is a program designed to build phylogenetic trees and test evolutionary
hypotheses regarding B cell affinity maturation. 

The biology of B cell somatic
hypermutation (SHM) violates important assumptions in most standard phylogenetic 
substitution models; further, while most phylogenetics programs are designed
to analyze single lineages, B cell repertoires typically contain thousands of
lineages. IgPhyML addresses both of these issues by implementing substitution
models that correct for the context-sensitive nature of SHM, and combines 
information from multiple lineages to give more precisely estimated repertoire-wide
model parameter estimates.

IgPhyML is integrated into the Immcantation framework, meaning users can easily 
build and visualize B cell lineage trees and test hypotheses about mutation and 
selection in B cell repertoires.
