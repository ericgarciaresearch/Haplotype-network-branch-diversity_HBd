# Haplotype network branch diversity (*HBd*)

Estimate the complexity of haplotype networks

A common way of illustrating phylogeographic results is through the use of haplotype networks.
 While these networks help to visualize relationships between individuals, populations, and species, evolutionary studies often only quantitatively analyze genetic diversity among haplotypes and ignore other network properties. 

[HapNetComplexity.R](https://github.com/ericgarciaresearch/Haplotype-network-branch-diversity_HBd/blob/main/HapNetComplexity.R) estimates complexity of haplotype networks (*HBd*) by combining haplotype (*Hd*) and topological (*Bd*) diversity

---
Citation:
[Garcia E, Wright D, Gatins R, Roberts MB, Pinheiro HT, Salas E, et al. (2021) Haplotype network branch diversity, a new metric combining genetic and topological diversity to compare the complexity of haplotype networks. PLoS ONE 16(6): e0251878. https://doi.org/10.1371/journal.pone.0251878](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0251878#sec011)

---

Use this script to build haplotype networks and calculate:

*	Nucleotide diversity (Pi) = a measurement of the genetic distance between sequences
*	Haplotype diversity (Hd) = a measurement of the genetic diversity in a population
*	Branch diversity (Bd) = a new measurement of the topological diversity in a haplotype network or the diversity of interrelationships among the observed haplotypes in a population
*	Haplotype Network Branch diversity (HBd) = a new measurement of the complexity in a haplotype network

Additional calculations:

*	Number of individuals (n)
*	Number of haplotypes (nH)
*	Number of haplotype classes (nHc)
*	Frequency of haplotype classes (niHc)

See [paper](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0251878#sec011) for definitions

---

**Usage**

Place this script in the same directory as:

1. A fasta file with your sequences
2. A file with site info (optional).
	* The site file consists of a .csv file with the sample names in the first column and site names in the second.
	* Call these columns "sample" and "site", respectively. Save the file as "sites.csv"
	* All metrics can still be calculated without any site information (see script)

Open the script to set your working directory, load alignments and specify or create a sites file.

Follow further directions within the `HapNetComplexity.R` script


Example model datasets of increasing complexity:

[S1_File_testing_A.fasta]() 
[S6_File_testing_F.fasta]()
[S12_File_testing_L.fasta]()

Sites file for all datasets

[S13_File_sites_for_all_Fig1.csv]()

**All example datasets contain 21 individuals and range from 6 to 21 haplotypes. See our [paper](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0251878#sec011) for more model and empirical examples**