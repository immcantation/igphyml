Using Dowser
===============================================================================

Currently the easiest way to use IgPhyML is through the R package 
`Dowser <https://dowser.readthedocs.io>`__. Dowser supports multiple methods for 
building B cell lineage trees, including IgPhyML. 

To build trees with Dowser, follow `this tutorial <https://dowser.readthedocs.io/en/latest/vignettes/Building-Trees-Vignette/#build-lineage-trees>`__. 
For building IgPhyML trees specifically, follow `this link <https://dowser.readthedocs.io/en/latest/vignettes/Building-Trees-Vignette/#build-igphyml-b-cell-trees>`__
 once formatting steps are complete. 

 In either case, IgPhyML will need to either have already been compiled from
 source code or contained within the Docker container.

 While :ref:`BuildTrees <BuildTrees-processing>` currently offers more direct options, 
 and masks codons split by V-gene reference alignment, these features under in development
 within Dowser.