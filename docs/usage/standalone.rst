.. _igphyml-standalone:


.. note::
   We now recommend using IgPhyML through `Dowser <http://dowser.readthedocs.io>`__.
   Dowser is our newer R package for phylogenetic analysis with B cells and supports
   building/visualizing/analyzing trees using IgPhyML and multiple other methods.


Standalone operation
=============================================================================== 

This section details how IgPhyML may be run directly as a standalone program. As before,
this section requires IgPhyML to be installed, with the executable in your ``PATH`` variable.
This is already done in the :ref:`Docker image <docker-image>`.

:ref:`Model specification<parameter-specification>` and :ref:`confidence interval estimation<ci-estimation>` are the same as when running IgPhyML from 
:ref:`BuildTrees<BuildTrees-processing>`. The means of specifying input, however,
are different.

To view all options for IgPhyML directly, run the command::

 igphyml --help

.. _single-lineage:

Analyzing a single lineage
-------------------------------------------------------------------------------

In its most basic and flexible operation, IgPhyML operates on single lineage. Its only required
input is an in-frame multiple sequence alignment, without any stop codons, in FASTA
or PHYLIP format. The model must be specified. The most basic is GY94, which is fast
but doesn't account for SHM context sensitivity. The HLP19 model corrects for SHM
context sensitivity but is slower requires the root sequence be specified::

    cd examples

    #Estimate parameters and topology using GY94 model
    igphyml -i example.fasta -m GY --run_id gy

    #Estimate parameters and topology using HLP19 model
    igphyml -i example.fasta -m HLP --root V4-59 --run_id hlp

It is also possible to use a fixed tree topology::
    
    igphyml -i example.fasta -m HLP --root V4-59 -u example.fasta_igphyml_tree_gy.txt

Or estimate separate :math:`\omega` values for CDR and FWR partitions::

    igphyml -i example.fasta -m HLP --root V4-59 --partfile part.example.txt

Generally, all the following features of repertoire-wide phylogenetic analysis
can also be performed on single lineages by specifying input files in the above
manner. 


Analyzing repertoires
-------------------------------------------------------------------------------

These commands should work as a first pass on many reasonably sized
datasets, but if you really want to understand what’s going on or make
sure what you’re doing makes sense, please check out the rest of the
manual.
 
**Convert Change-O files into IgPhyML inputs**
 
Move to the ``examples`` subfolder and run, in order::

    BuildTrees.py -d example.tsv --outname ex --log ex.log --collapse

This will create the directory ``ex`` and the file
``ex_lineages.tsv``. Each ``ex/<clone ID>.fasta`` contains the IMGT
mutliple sequence alignemt for a particular clone, and each
``ex/<clone ID>.part.txt`` file contains information about V and J
germline assignments, as well as IMGT unique numbering for each site.
The file ``ex.log`` will contain information about whether or not each
sequence was included in the analysis. The file ``ex_lineages.tsv`` is
the direct input to IgPhyML. Each line represents a clone and shows
the multiple sequence alignment, starting tree topology (N if
ignored), germline sequence ID in alignment file, and partition file
(N if ignored). These repertoire files start with the number of
lineages in the repertoire, and lineages are arranged from most to
least number of sequences. Here, it is important to not use 
``--igphyml`` or ``--clean all`` to prevent IgPhyML from being run 
by BuildTrees, and to prevent these intermediate files from being deleted.
 
**Build lineage trees using the GY94 model**

This option is fast and makes good starting topologies, but doesn't correct
for SHM mutation biases. Use the ``--outrep`` option to make a modified ::

 igphyml --repfile ex_lineages.tsv -m GY --outrep ex_lineages_gy.tsv --run_id gy

**Build lineage trees using the HLP19 model with GY94 starting trees** 

This option is slower but corrects for mutation biases of SHM::

 #HLP topology, branch lengths, and parameters (slow)
 igphyml --repfile ex_lineages_gy.tsv -m HLP --threads 1

 #GY94 topology, HLP19 branch lengths and parameters (faster)
 igphyml --repfile ex_lineages_gy.tsv -m HLP --optimize lr --threads 1

Both of these can be parallelized by modifying the
``--threads <thread count>`` option. Trees files are listed as
``ex/<clone id>.fasta_igphyml_tree.txt``, and can be viewed with most
tree viewers (I recommend
`FigTree <http://tree.bio.ed.ac.uk/software/figtree/>`__). Parameter
estimates are in ``ex_lineages.tsv_igphyml_stats.txt``. 


Controlling output format
------------------------------------------------------------------------

Alternatively, run using ``--oformat tab`` to create input readable by 
the 
`readIgphyml <https://alakazam.readthedocs.io/en/stable/topics/readIgphyml>`__ 
function of 
`Alakazam <https://alakazam.readthedocs.io>`__.::

 #Output can be read using readIgphyml function
 igphyml --repfile ex_lineages_gy.tsv -m HLP --optimize lr --threads 1 --oformat tab

Open an ``R`` session, and run the following commands. Note the results are the same as in the :ref:`quickstart example<igphyml-quickstart>` ::

 library(alakazam)
 library(igraph)
 
 db = readIgphyml("ex_lineages_gy.tsv_igphyml_stats.tab")

 #plot largest lineage tree
 plot(db$trees[[1]],layout=layout_as_tree)

 #show HLP10 parameters
 print(t(db$param[1,]))
 CLONE         "REPERTOIRE"
 NSEQ          "4"         
 NSITE         "107"       
 TREE_LENGTH   "0.286"     
 LHOOD         "-290.7928" 
 KAPPA_MLE     "2.266"     
 OMEGA_FWR_MLE "0.5284"    
 OMEGA_CDR_MLE "2.3324"    
 WRC_2_MLE     "4.8019"    
 GYW_0_MLE     "3.4464"    
 WA_1_MLE      "5.972"     
 TW_0_MLE      "0.8131"    
 SYC_2_MLE     "-0.99"     
 GRS_0_MLE     "0.2583"

.. figure:: ../_static/t1.png
   :scale: 25 %
   :align: center
   :alt: map to buried treasure

   Lineage tree of example clone.

