Installation
================================================================================

.. _docker-image: 

Docker Image (recommended)
--------------------------------------------------------------------------------

The simplest way to use IgPhyML is via the 
`Immcantation Docker image <https://immcantation.readthedocs.io/en/stable/docker/intro.html>`__.

Briefly, all the example commands on rest of this site can be run by first installing Docker and
downloading the Immcantation Docker image. In a terminal, enter::

 # pull Immcantation development Docker image
 docker pull kleinstein/immcantation:devel

 # clone IgPhyML repository to get example files
 git clone https://bitbucket.org/kbhoehn/igphyml

Then, move to the examples directory and load it into the Docker image depending on your operating system::
 
 # move to examples directory
 cd igphyml/examples

 # load Docker, Linux/ Mac OS X:
 docker run -it --workdir /data -v $(pwd):/data:z kleinstein/immcantation:devel bash

 # or load Docker, Windows:
 docker run -it --workdir /data -v %cd%:/data:z kleinstein/immcantation:devel bash

Note for some operating systems it may be necessary to have 
`Docker Desktop <https://hub.docker.com/editions/community/docker-ce-desktop-windows>`__
running before entering these commands. Once inside the container, check everything is properly configured::

 # should give example.fasta  example.tab  part.example.txt
 ls

 # should be IgPhyML 1.0.7 201902.19
 igphyml --version

More generally, use this command to load the Docker image on the current directory of a Linux/Max OS X system::

 docker run -it --workdir /data -v $(pwd):/data:z kleinstein/immcantation:devel bash

or for a Windows Command Prompt::

 docker run -it --workdir /data -v %cd%:/data:z kleinstein/immcantation:devel bash

For further information on using the Immcantation Docker image, see 
`Immcantation Docker image <https://immcantation.readthedocs.io/en/stable/docker/intro.html>`__.

Building from source
--------------------------------------------------------------------------------
If using the Docker image is not possible or preferable, the 
source code of the current development version can be downloaded using git::

    > git clone https://bitbucket.org/kbhoehn/igphyml
    > cd igphyml

Requirements
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

+ GNU autoconf
+ GNU automake
+ (optional) `Change-O 0.4.6 <https://changeo.readthedocs.io>`__
+ (optional) `Alakazam 0.3.0 <https://alakazam.readthedocs.io>`__
+ (optional) OpenMP-enabled C compiler (e.g. gcc or LLVM)
+ (optional) BLAS and LAPACK optimization libraries

Linux
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

On Linux operating systems, you can usually just run::

    ./make_phyml_omp

If you have BLAS and LAPACK installed,
which provide libraries for faster and more accurate matrix exponentiation
operations. In Ubuntu Linux, these are provided in the packages
``libblas-dev`` and ``liblapack-dev``. Other distros probably have
similar package names. To compile IgPhyML with BLAS and LAPACK 
support, run::
 
    ./make_phyml_blas_omp
 
If OpenMP is also not available, IgPhyML may be compiled without it,
but this will substantially reduce performance on multicore machines::
 
    ./make_phyml

Once everything is compiled, add the igphyml/src directory to your
``PATH`` variable, and IgPhyML should work. Importantly, some directory
files are hardcoded in at compilation, so re-compile IgPhyML if you move
the installation directory. Alternatively, you can set the ``IGPHYML_PATH``
environment variable to the location of the igphyml/src/motifs folder for
the same effect.

Mac OS X
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Installation on Mac OS X is trickier, but possible. The primary issue
is gaining OpenMP support, and installing some GNU command line tools.
The best way is to just install the latest version of ``llvm``
available through ``homebrew``, as well as ``autoconf`` and
``automake``. To do these you’ll need to install
`homebrew <http://brew.sh/index.html>`_. If it’s already installed be
sure it’s at the latest version (``brew update``). You may need to install
Xcode as well. Next, install ``autoconf``, ``automake``, and ``llvm``::

    brew install autoconf
    brew install automake
    brew install llvm

Specify the ``llvm`` version of ``clang`` in ``Makefile.am`` and
``src/Makefile.am`` by adding the line ``CC=<path to llvm clang>``
to the beginning of both files. You will also need to add
``MACOMP=<path to omp.h>`` and ``MACLLVM=<path to llvm lib>`` to
``src/Makefile.am``. For instance, if you’ve install ``llvm 3.9.1``
via homebrew, you will likely need to add the line
``CC=/usr/local/Cellar/llvm/3.9.1/bin/clang``
to ``Makefile.am`` and the lines::

    CC=/usr/local/Cellar/llvm/3.9.1/bin/clang
    MACOMP=/usr/local/Cellar/llvm/3.9.1/lib/clang/3.9.1/include/omp.h
    MACLLVM=/usr/local/Cellar/llvm/3.9.1/lib

to ``src/Makefile.am``.
Your specific path may look different, but you can check locations
of these files and folders by looking around in
``/usr/local/Cellar/llvm/``. The directory structure should be
similar. Run ``./make_blas_phyml_omp``, or other versions, as desired, and add
the ``src`` folder to your ``PATH`` variable.

On some versions of OS X it may be necessary to install XCode command
line tools using::

    xcode-select --install
    cd /Library/Developer/CommandLineTools/Packages/
    open macOS_SDK_headers_for_macOS_<OS X version>.pkg