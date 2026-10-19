# WGCNA of TCGA-PRAD by Gleason grade group

Co-expression analysis of TCGA-PRAD primary prostate tumours, relating gene modules to Gleason grade group (GG1–GG5). 

See `WGCNA_notes.md` for writeup and plots. 

## Requirements

- R packages: `TCGAbiolinks`, `SummarizedExperiment`, `edgeR`, `WGCNA`
- System libraries for building WGCNA's dependencies (Ubuntu): fonts and the usual `-dev` packages for `Hmisc`/`htmlTable`


## Running it

```bash
Rscript wgcna_mytest.R 
```
Note that I haven't inccluded the rds file (too large for github without LFS) but you can easily query it yourself. 


## Output

Will give a bunch of pdfs, related to filtering, normalisation and the threshold of the network. 
Also the dendrogram and the module boxplots split by Gleason Grade. 
