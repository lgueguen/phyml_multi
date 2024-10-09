# PHYML_MULTI


## How to compile PhyML_multi ?
**********************


### Linux or Solaris (with gcc)

. 'cd' to the sources directory,
. type 'make'


### Windows (with Microsoft Visual C++ version 6)

. click on the .\exe\phyml.dsp file
. click on Build->Rebuild all
. click on the file 'phyml.exe' that has been 
  created in the folder .\Release


### OSX (Jaguar or Panther)

. 'cd' to the sources directory
. remove the option '-static' from the CFLAGS
. type 'make'


## Usage

Posterior analysis of phyml_multi reconstruction needs either to run
PartitioningHMM.py for HMM or PartitioningMM.py for maximum
predictive partitioning.

These files need the installation of SARMENT libraries,
[there](https://github.com/lgueguen/SARMENT).


