# README for MG #

This README documents the steps necessary to get your MG installation up and running.

This repository only contains mg_torch, the supernova module of mg.

### What is MG? ###

* MG is a C based compressible finite volume CFD code, equipped with adaptive mesh refinement (AMR), and can be deployed on high performance computing (HPC) clusters using message passing interface (MPI).
* 2026 Version
* [Bitbucket link](https://bitbucket.org/bmaxw005/mg2026/src/main/)
* [YouTUbe:  MG download using git and install using Bash](https://youtu.be/UiuRLH5mtxE)
* [YouTUbe:  MG Navigation](https://youtu.be/KGIHW5S3C7k)


### MG solvers included ###

* <b>mg_torch</b><br>
MG with supernova equations of state and fusion reactions.  It has 8-isotope, 13-isotope, and the torch reaction model.  Compile with double prescision only.

### Dependencies ###

The following packages should be installed on your system prior to installation:

* tcsh
* make
* gcc
* g++
* gfortran
* openmpi-bin
* libopenmpi-dev
* libx11-dev
* libreadline-dev
* libgsl-dev
* libncurses5-dev
* libncursesw5-dev

### HPC modules ###

At present, to install on the Digital Alliance of Canada computers (Rorqual, Fir, Nibi, Narval) you should only need to load the following modules (i.e. with module load)

* gsl

For the trillium system you will need to load the following:

* gcc/12.3
* StdEnv/2023
* openmpi
* gsl


### Installation ###

1. To install bash version, first add the following lines to your $HOME/.bashrc file (assuming mg has been downloaded to ~/mg):<br>
` export mg=~/mg `<br>
` source $mg/_bashrc `
* Source the .bashrc by opening a new window, or typing<br>
` cd $HOME `<br>
` source .bashrc `
* Edit the `$mg/_bashrc` file if needed.  Note, there is also a _cshrc if you plan
to install using tcsh (or mac?)
* Compile with:<br>
` cd $mg `<br>
` Createlib `
* Finally, navigate to the MG program you wish to compile:<br>
` cd $mg/mg_g `<br>
` makemg mpi `


### Documentation ###

* Consult the [guide.html](./docs/guide.html) for basic documentation.
* A full installation manual is currently being developed in Overleaf.  A draft version of the [MG_Installation_and_Training_Manual.pdf](./docs/MG_Installation_and_Training_Manual.pdf) is packaged with this code.


### Tests ###

* <b>mg_torch:</b><br>


### Support ###

* This version of MG is supported by Brian Maxwell
* [bmaxwell@uottawa.ca](bmaxwell@uottawa.ca)

