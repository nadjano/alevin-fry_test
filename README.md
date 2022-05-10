
# snRNA-mapping-rate: Pipeline to compare mapping rates accross different tools for single-nuclei RNA sequencing experiments.

## Requirements
Install Dockker image for alevin-fry: 
singularity pull docker://combinelab/usefulaf:latest

## How to run
```
nextflow run 'https://github.c/ebi-gene-expression-group/snRNA-mapping-rate.git' -r main \
--sdrf snRNA-mapping-rate/E-CURD-90.sdrf.txt \
--resultsRoot results \
--referenceGenome E-CURD_90/reference/Drosophila_melanogaster.BDGP6.32.dna.toplevel.fa \
--referencecDNA E-CURD_90/reference/Drosophila_melanogaster.BDGP6.32.cdna.all.fa.gz --referenceGtf E-CURD_90/reference/Drosophila_melanogaster.BDGP6.32.106.gtf \
--config
```