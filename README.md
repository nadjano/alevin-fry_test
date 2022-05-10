
# snRNA-mapping-rate: Pipeline to compare mapping rates accross different tools for single-nuclei RNA sequencing experiments.

## Requirements
Install Dockker image for alevin-fry: 
singularity pull docker://combinelab/usefulaf:latest

## How to run
```
nextflow run 'https://github.c/ebi-gene-expression-group/snRNA-mapping-rate.git' -r main \
--sdrf [srdf.txt] \
--resultsRoot [resultsDir] \
--referenceGenome [DNA.fa] \
--referencecDNA [cDNA.fa.gz] --referenceGtf [File.gtf] \
--config
```