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


### Man

> NAME
> 	phyml_multi, derived from PhyML,
> 	A simple, fast, and accurate algorithm to estimate
> 	large phylogenies by maximum likelihood.
> 
> 	Stephane Guindon and Olivier Gascuel,
> 	Systematic Biology 52(5):696-704, 2003.
> 	Please cite this paper if you use this software in your publications.
> 
> COMMAND-LINE USE
> 	phyml [ sequences data_type format data_sets bootstrap_sets model [kappa] invar nb_categ alpha tree opt_topology opt_lengths use_HMM nb of trees ]
> 
> 	You can use phyml with no arguments, in this case change the value of
> 	a parameter by typing its corresponding character as shown on screen.
> 
> 	You can alternatively use phyml with the following arguments :
> 
> 	sequence_file	DNA or Amino-Acids sequence filename (PHYLIP format)
> 
> 	data type	0 = DNA | 1 = Amino-Acids
> 
> 	format		i = interleaved sequence format | s = sequential
> 
> 	data_sets	number of data sets to analyse (ex:3)
> 
> 	bootstrap_sets	number of bootstrap data sets to generate (ex:2)
> 			only works with one data set to analyse
> 
> 	model		substitution model name
> 			JC69 | K2P | F81 | HKY | F84 | TN93 (DNA)
> 			JTT | MtREV | Dayhoff | WAG (Amino-Acids)
> 
> 	kappa		transition/transversion ratio, only for DNA sequences,
> 			a fixed value (ex:4.0) | e to get the maximum likelihood estimate
> 
> 	invar		proportion of invariable sites,
> 			a fixed value (ex:0.0) | e to get the maximum likelihood estimate
> 
> 	nb_categ	number of relative substitution rate categories (ex:4)
> 
> 	alpha		gamma distribution parameter,
> 			a fixed value (ex:1.0) | e to get the maximum likelihood estimate
> 
> 	tree		starting tree filename (Newick format),
> 			your tree filename | BIONJ for a distance-based tree
> 
> 	opt_topology	optimise tree topology ? y | n
> 
> 	opt_lengths	optimise branch lengths and rate parameters ? y | n
> 
> 	use_HMM		use HMM ? n | y
> 
> 	number of trees		How many trees (ex:4) ? 
> 
> Examples
> 
> .DNA sequences, no HMM :   ./phyml_multi seqs1 0 i 2 0 HKY 4.0 e 1 1.0 BIONJ y n n 2
> 
> .AA sequences, HMM :    ./phyml_multi seqs2 1 i 1 5 JTT 0.0 4 1.0 BIONJ n n y 3
> 

## Example

The file simSeq.phy contains an example alignment in which a recombination event occured. This is one of the simulated alignments as found in the 
"Evolutionary Bioinformatics" article. The breakpoint should be found at position 599.
To use phyml_multi, type: 

For use with the phylo-HMM: 

phyml_multi simSeqs.phy 0 i 1 0 HKY e e 4 e BIONJ y y y 2


and then to see the segmentation (with SARMENT library):

python3 PartitioningHMM.py simSeqs.phy_phyml_siteLks.txt 0.997973
or : 
python3 PartitioningMM.py simSeqs.phy_phyml_siteLks.txt 


For use with the mixture model: 
phyml_multi simSeqs.phy 0 i 1 0 HKY e e 4 e BIONJ y y n 2

and then to see the segmentation:
python3 PartitioningMM.py simSeqs.phy_phyml_siteLks.txt


This should give you the right answer, in a matter of a few minutes on a single desktop computer.

For any question, contact bastien.boussau@univ-lyon1.fr

## Reference

Boussau B, Guéguen L, Gouy M. A mixture model and a hidden markov
model to simultaneously detect recombination breakpoints and
reconstruct phylogenies. Evol Bioinform Online. 2009 Jun 25;5:67-79.
doi: 10.4137/ebo.s2242. PMID: 19812727; PMCID: PMC2747125.

