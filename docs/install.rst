
.. _install:

Installation
================================================================================

.. _docker-image: 

Docker Image (recommended for non-Linux systems)
--------------------------------------------------------------------------------

The simplest way to use IgPhyML is via the 
`Immcantation Docker image <https://immcantation.readthedocs.io/en/stable/docker/intro.html>`__. 
``4.5.0`` is used for these instructions, but you should use the latest version on the Immcantation site.

Briefly, all the example commands on rest of this site can be run by first installing Docker and
downloading the Immcantation Docker image. Note for some operating systems it may be necessary to have 
`Docker Desktop <https://hub.docker.com/editions/community/docker-ce-desktop-windows>`__
running before entering these commands.
In a terminal, enter::

 # pull Immcantation Docker image
 docker pull immcantation/suite:4.5.0

 # clone IgPhyML repository to get example files
 git clone https://github.com/immcantation/igphyml

You may also need to install ``git`` to complete the last command. Then, move to the examples directory and load it into the Docker image depending on your operating system::
 
 # move to examples directory
 cd igphyml/examples

 # load Docker, Linux/ Mac OS X:
 docker run -it --workdir /data -v $(pwd):/data:z immcantation/suite:4.5.0 bash

 # or load Docker, Windows:
 docker run -it --workdir /data -v %cd%:/data:z immcantation/suite:4.5.0 bash

Once inside the container, check everything is properly configured::

 # should give example.fasta  example.tab  part.example.txt
 ls

 # should be IgPhyML 2.0.0 092223
 igphyml --version

More generally, use this command to load the Docker image on the current directory of a Linux/Max OS X system::

 docker run -it --workdir /data -v $(pwd):/data:z immcantation/suite:4.5.0 bash

or for a Windows Command Prompt::

 docker run -it --workdir /data -v %cd%:/data:z immcantation/suite:4.5.0 bash

For further information on using the Immcantation Docker image, see 
`Immcantation Docker image <https://immcantation.readthedocs.io/en/stable/docker/intro.html>`__.

Compiling from source (recommended for Linux)
--------------------------------------------------------------------------------
If using Linux, or if the Docker image is not possible or preferable, the 
source code of the current development version can be downloaded using git and compiled::

    git clone https://github.com/immcantation/igphyml
    cd igphyml
    ./make_phyml_omp 

If compilation (the last step) was successful, you should see the executable file ``src/igphyml``.
If this step fails, there is still hope. See below for more detailed requirements and
options.

Requirements
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

+ GNU autoconf
+ GNU automake
+ (optional) `Dowser >1.0.0 <https://dowser.readthedocs.io>`__
+ (optional) OpenMP-enabled C compiler (e.g. gcc or LLVM)
+ (optional) BLAS and LAPACK optimization libraries

For Ubuntu systems, you can install automake tools and BLAS/LAPACK packages using::

    apt-get install automake autoconf libblas-dev liblapack-dev libatlas-base-dev

This make require `sudo` permission. Other Linux distros should have similar pacakges.

Linux
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

On Linux operating systems with the above requirements installed, you can usually just run::

    ./make_phyml_omp

If you have BLAS and LAPACK installed,
which provide libraries for faster and more accurate matrix exponentiation
operations. In Ubuntu Linux, these are provided in the packages
``libblas-dev``, ``liblapack-dev``, and ``libatlas-base-dev``. Other distros probably have
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

Specify the ``llvm`` version of ``clang`` in both ``Makefile.am`` and
``src/Makefile.am`` by adding the line ``CC=<path to llvm clang>``
to the beginning of both files. You will also need to add
``MACOMP=<path to omp.h>`` and ``MACLLVM=<path to llvm lib>`` to
``src/Makefile.am``. 

**Edit Makefile.am**:
For example, if you’ve installed ``llvm 19.1.4`` add this line to ``Makefile.am``::

    CC=/usr/local/Cellar/llvm/19.1.4/bin/clang

**Then edit src/Makefile.am**:
For example, if you’ve installed ``llvm 19.1.4`` add these the lines to ``src/Makefile.am``::

    CC=/usr/local/Cellar/llvm/19.1.4/bin/clang
    MACOMP=/usr/local/Cellar/llvm/19.1.4/lib/clang/19/include/omp.h
    MACLLVM=/usr/local/Cellar/llvm/19.1.4/lib

A version of these are already shown as comments in these files.
Your specific paths may look different, but you can check locations
of these files and folders by looking around in
``/usr/local/Cellar/llvm/``. The directory structure should be
similar. 

Once you've edited both files, run::

    ./make_phyml_omp

or other compilation options as desired and add
the ``src`` folder to your ``PATH`` variable.

You can check if compilation was successful by running::

    ./src/igphyml -version

On some versions of OS X it may be necessary to install XCode command
line tools using::

    xcode-select --install
    cd /Library/Developer/CommandLineTools/Packages/
    open macOS_SDK_headers_for_macOS_<OS X version>.pkg