# README for MG #

This README documents the steps necessary to get your MG installation up and running.


### What is MG? ###

* MG is a C based compressible finite volume CFD code, equipped with adaptive mesh refinement (AMR), and can be deployed on high performance computing (HPC) clusters using message passing interface (MPI).
* 2026 Version
* [Bitbucket link](https://bitbucket.org/bmaxw005/mg2026/src/main/)
* [YouTUbe:  MG download using git and install using Bash](https://youtu.be/UiuRLH5mtxE)
* [YouTUbe:  MG Navigation](https://youtu.be/KGIHW5S3C7k)


### MG solvers included ###

The different subfolders (mg programs) included in this version of MG are as follows:

* <b>mg_g</b><br>
This is the original mg_g code (calorically perfect gas equations) but allows internal boundaries based on fluid variable iqbd, and both a RANS and LES (LKM) model, and a few other features.  Creating new mg versions are usually best to be based off this version.  This  solver can be compiled in single or double prescision mode.

* <b>mg_lem</b><br>
This is a calorically perfect gas model but with one, and 3 step combustion.  This is the only version that has the CLEM-LES model (only for one-step combustion).  It can be compiled in single or double prescision.

* <b>mg_cti</b><br>
This is the detailed chemistry version.  It reads legacy cantera cti files, and allowes variable temperature dependent thermodynamic properties (Cv,Cp).  It uses integer-based flags for internal geometry.  Compile with double prescision only.

* <b>mg_torch</b><br>
MG with supernova equations of state and fusion reactions.  It has 8-isotope, 13-isotope, and the torch reaction model.  Compile with double prescision only.

* <b>mg_2p</b><br>
Two phase flow version using the NASG equation of state.  It also has a reduced 4-step combustion model (Boivin) for hydrogen chemistry, but uses constant Cv,Cp in the NASG framework. Compile with double prescision only.



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

The following is a list of example scripts that can be run manually in each respective folder to demonstrate various features of MG:

* <b>mg_g:</b><br>
` sod.dat `<br>
` sodxy.dat `<br>
` circle.dat `<br>
` obstacle.dat `<br>
* <b>mg_lem:</b><br>
` box.dat `<br>
` wallflame.dat `<br>
* <b>mg_cti:</b><br>
` sod.dat `<br>
` sodxy.dat `<br>
` blast.dat `<br>
` equil.dat `<br>
` reactor.dat `<br>
` reactorKonnov.dat `<br>
` reactor_h2no2.dat `<br>
` diffusion2d_mgimp.dat `<br>
` diffusion2d_adi.dat `<br>
` shearlayer.dat `<br>
` H2jet.dat `<br>
` 1Ddetonation.dat `<br>
` 2Ddetonation.dat `<br>
* <b>mg_2p:</b><br>
` sod.dat `<br>
` sod2.dat `<br>
` reactor.dat `<br>
* <b>mg_torch:</b><br>


### HPC examples ###

Examples of larger simulations, requiring HPC resources, are found in the `hpc_examples` subfolder.  The following examples are included:


* `ch4det_lem`<br>
This is an example of a methane-oxygen detonation using the linear eddy model in the `mg_lem` framework.

* `h2det_cti`<br>
This is an example of a hydrogen-oxygen-argon detonation using detailed chemsitry in the `mg_cti` framework.



### Support ###

* This version of MG is supported by Brian Maxwell
* [bmaxwell@uottawa.ca](bmaxwell@uottawa.ca)

