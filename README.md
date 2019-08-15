[![DOI](https://zenodo.org/badge/186807430.svg)](https://zenodo.org/badge/latestdoi/186807430)

## Sheehan et al. (2019) Evolution of L-DOPA 4,5-dioxygenase activity allows for recurrent specialisation to betalain pigmentation in Caryophyllales 

All code in this repository is provided as-is and comes with _absolutely no warranty_. Feel free to use (at your own risk) and share, but please cite the paper above when you do. It is now out in early release, with DOI 10.1111/nph.16089. 

Homolog scripts come courtesy of Ya Yang and Stephen Smith, citation is Lopez-Nieves et al 2018: Relaxation of tyrosine pathway regulation underlies the evolution of betalain pigmentation in Caryophyllales and are sourced from https://bitbucket.org/yangya/adh_2016/

misc_scripts contains various utilities for the analysis, including R scripts for ancestral reconstruction 

alignments contains all the alignments used in this study:
* Caryophyllales_3524_outaln.cn: alignment from pyphlawd used for species tree
* Caryophyllales_3524_outpart: partitioning file 
* DODA_full: all DODA sequences used to produce the tree in Figure 6 and Figure S5
* DODA_alpha: the reduced DODAa sequences used to produce the tree in Figure S16
* DODA_asr: the subsampled DODAa alignment used to produce the tree in Fig Figure 7 and Figure S12 

trees contains the newick files of the inferred DODA trees, inferred from the cleaned (-cln) version of each alignment. Labels are Rapid Boostrap Support and SH-aLRT support (RBS/SH). The Caryophyllales tree was inferred from the alignment as is, and branch labels are RBS only.

If you have any questions with these scripts or data, don't hesitate to open an issue or contact me directly at nw387@cam.ac.uk 
