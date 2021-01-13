## MetWAMer - Translation initiation site prediction in eukaryotes

MetWAMer is a bioinformatics application used to predict translation initiation sites (TISs) in eukaryotic protein-coding genes of non-viral origin. More generally, however, it is a software framework that can be used to assist with developing TIS prediction methods.

The software was described in Sparks, M.E. and Brendel, V. (2008) MetWAMer: eukaryotic translation initiation site prediction. *BMC Bioinformatics*. **9**:381 (https://doi.org/10.1186/1471-2105-9-381). Supporting data from that report, as well as prior versions of the code, are available at http://brendelgroup.org/SB08B/.

This repository provides a "rolling snapshot" of the package and, as of 13 January 2021, includes demonstrations of how an end user might train the system using publicly available data from *Bombyx mori*, and then use the results for predicting TISs in *Manduca sexta* datasets.

MetWAMer has been tested most recently on a 64-bit GNU/Linux system and it compiles cleanly under at least GCC versions 10.1.0 and 7.5.0. The package requires the libxml2 library and was observed to work properly under version 2.9.4. It is otherwise self-contained, although host systems are encouraged to provide the GNU implementation of awk, gawk, if possible. Please see the manual page for usage details. This software is provided "**AS IS**".
