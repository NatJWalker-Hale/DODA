[![DOI](https://zenodo.org/badge/186807430.svg)](https://zenodo.org/badge/latestdoi/186807430)

## Sheehan et al. (2019) Evolution of L-DOPA 4,5-dioxygenase activity allows for recurrent specialisation to betalain pigmentation in Caryophyllales 

All code in this repository is provided as-is and comes with _absolutely no warranty_. Feel free to use (at your own risk) and share, but please cite the paper above when you do. It is available [here](https://nph.onlinelibrary.wiley.com/doi/full/10.1111/nph.16089).

The full citation is:

> Sheehan, H., Feng, T., Walker‐Hale, N., Lopez‐Nieves, S., Pucker, B., Guo, R., Yim, W. C., Badgami, R., Timoneda, A., Zhao, L., Tiley, H., Copetti, D., Sanderson, M. J., Cushman, J. C., Moore, M. J., Smith, S. A., & Brockington, S. F. (2020). Evolution of l-DOPA 4,5-dioxygenase activity allows for recurrent specialisation to betalain pigmentation in Caryophyllales. New Phytologist, 227(3), 914–929. https://doi.org/10.1111/nph.16089

Homolog scripts come courtesy of Ya Yang and Stephen Smith, citation is [Lopez-Nieves et al 2018](https://nph.onlinelibrary.wiley.com/doi/epdf/10.1111/nph.14822) and are sourced from https://bitbucket.org/yangya/adh_2016/

misc_scripts contains various utilities for the analysis, including R scripts for ancestral reconstruction 

alignments contains all the alignments used in this study:
* Caryophyllales_3524_outaln.cn: alignment from pyphlawd used for species tree
* Caryophyllales_3524_outpart: partitioning file 
* DODA_full: all DODA sequences used to produce the tree in Figure 6 and Figure S5
* DODA_alpha: the reduced DODAa sequences used to produce the tree in Figure S16
* DODA_asr: the subsampled DODAa alignment used to produce the tree in Fig Figure 7 and Figure S12
  * The alignment file with the suffix .name contains sequences with taxon codes replaced with full taxon names according to the correspondence in /homolog_search/taxon_table

trees contains the newick files of the inferred DODA trees, inferred from the cleaned (-cln) version of each alignment. Labels are Rapid Boostrap Support and SH-aLRT support (RBS/SH). The Caryophyllales tree was inferred from the alignment as is, and branch labels are RBS only.

If you have any questions with these scripts or data, don't hesitate to open an issue or contact me directly at nw387@cam.ac.uk 
