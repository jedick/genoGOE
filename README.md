# genoGOE

This repository holds the code and data for "Evolutionary oxidation of proteins in Earth's history".

## Requirements

- R - download from [CRAN](https://cran.r-project.org)
- R packages - install with `install.packages(c("canprot", "CHNOSZ", "plotrix", "beanplot"))`
  - canprot - provides `Zc()`, `read_fasta()`, and `Cost()`
  - CHNOSZ - provides `protein.formula()`, `as.chemical.formula()`, and `expr.species()`
  - plotrix - provides `smoothColors()` and `gradient.rect()`
  - beanplot - provides `beanplot()`
  - ggplot2 - provides plotting functions used in `genoGOE_2()`
  - patchwork - provides `plot_layout()`

## Making the plots

Set `pdf = TRUE` to make PDF files or leave as FALSE (the default) to make plots on screen.

```R
source("genoGOE.R")
genoGOE_1(pdf = TRUE)
genoGOE_2(pdf = TRUE)
genoGOE_3(pdf = TRUE)
genoGOE_4(pdf = TRUE)
genoGOE_5(pdf = TRUE)
# Needed for thermodynamic calculations
library(CHNOSZ)
genoGOE_6(pdf = TRUE)
```

## Data sources

[UniProt/UP000000625_83333.csv.xz](UniProt/UP000000625_83333.csv.xz) - *E. coli* reference proteome from [UniProt](https://uniprot.org)

[GTDB](GTDB): Processed data from GTDB with additional features

- [methanogen_genomes.csv](GTDB/methanogen_genomes.csv): Genome IDs and taxonomy in GTDB r220 for 19 Class I and 19 Class II methanogen species selected from Fig. 1 of [Lyu and Lu (2018)](https://doi.org/10.1038/ismej.2017.173)
- [process_GTDB.R](GTDB/process_GTDB.R): Script to obtain DNA and protein sequences for 53 archaeal marker genes in GTDB and amino acid compositions for all proteins in 38 methanogen genomes
- [ar53_msa_marker_info_r220_XHZ+06.csv](GTDB/ar53_msa_marker_info_r220_XHZ+06.csv): List of archaeal marker genes from GTDB with protein abundance information for *Methanococcus maripaludis* added from [Xia et al. (2006)](https://doi.org/10.1074/mcp.M500369-MCP200)

[methanogen](methanogen): Sequences and amino acid compositions generated using `process_GTDB.R`

- [marker/fna](methanogen/marker/fna): Nucleotide sequences of marker genes
- [marker/faa](methanogen/marker/faa): Amino acid sequences of marker genes
- [aa](methanogen/aa): Amino acid compositions of all proteins

[MCK+23](MCK+23): Data for genomes with sulfur-cycling genes from [Mateos et al. (2023)](https://doi.org/10.1126/sciadv.ade4847)

- [get_genomes.R](MCK+23/get_genomes.R): Script to identify genomes from bootstrap files, download genomes from NCBI, and make the following CSV files
- [genomes.csv](MCK+23/genomes.csv): Table of genomes with binary labels denoting presence of specific S-cycling gene
- [genomes_aa.csv](MCK+23/genomes_aa.csv): Amino acid composition for each genome with availability in NCBI
- [sulfur_genomes.xlsx](MCK+23/sulfur_genomes.xlsx): Spreadsheet listing genomes that exclusively contain one S-cycling gene or gene cluster

[GMKK20](GMKK20): Data for extant and reconstructed ancestral nitrogenase sequences taken from [Garcia et al. (2020)](https://doi.org/10.1111/gbi.12381)

- [nitrogenase_aa.csv](GMKK20/nitrogenase_aa.csv): Amino acid compositions computed from the [Extant-MLAnc_Align.fasta](https://github.com/kacarlab/AncientNitrogenase/blob/master/Extant-MLAnc_Align.fasta) file in the [kacarlab](https://github.com/kacarlab) GitHub repo

[PIZ+11](PIZ+11): Data for reconstructed ancestral thioredoxin sequences derived from [Perez-Jimenez et al. (2011)](https://doi.org/10.1038/nsmb.2020)

- [thioredoxin.fasta](PIZ+11/thioredoxin.fasta): Protein sequences obtained from the RCSB PDB using the accessions listed in the next file.
- [DAAD19.csv](PIZ+11/DAAD19.csv): RCSB PDB IDs and ages for proteins listed by [Del Galdo et al. (2019)](https://doi.org/10.1007/s00239-019-09894-4).

[CDY+25/IPMDH.fasta](CDY+25/IPMDH.fasta): Data for reconstructed ancestral 3-isopropylmalate dehydrogenases sequences taken from Supporting Information of [Cui et al. (2025)](https://doi.org/10.1002/pro.70071)

[LMM16](LMM16): Scripts and processed data files for consensus gene ages from [Liebeskind et al. (2016)](https://doi.org/10.1093/gbe/evw113), modified from the files used by [Dick (2022)](https://doi.org/10.1007/s00239-022-10051-7)

- [mkaa.R](LMM16/mkaa.R): *script*: sum amino acid compositions of proteins in each gene age category
- [reference_proteomes.csv](LMM16/reference_proteomes.csv): *data*: IDs of UniProt reference proteomes for 31 organisms
- [modeAges_names.csv](LMM16/modeAges_names.csv): *output file*: Names of gene age categories for each organism
- [modeAges_aa.csv](LMM16/modeAges_aa.csv): *output file*: Summed amino acid composition for proteins in each gene age category

[KHAB17/rubisco.fasta](KHAB17/rubisco.fasta): Reconstructed ancestral Rubisco sequences taken from [Kaçar et al. (2017)](https://doi.org/10.1111/gbi.12243)

[DBCS23](DBCS23): Amino acid compsitions of methanogen and *Nitrososphaeria* genomes used by [Dick et al. (2023)](https://doi.org/10.1111/gbi.12532)

## BacDive analysis

See [bacdive/README.md](bacdive/README.md)
