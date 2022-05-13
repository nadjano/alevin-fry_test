
# snRNA-mapping-rate: Pipeline to compare mapping rates accross different tools for single-nuclei RNA sequencing experiments.
Runs Alevin, Alevin-fry, kb-tools and STARsolo on snRNA-Seq data with intron containing and intron free references.
## Requirements
Nextflow (version 21.10.6)

## How to run
```
nextflow run 'https://github.c/ebi-gene-expression-group/snRNA-mapping-rate.git' -r main \
--sdrf [srdf.txt] \
--resultsRoot [resultsDir] \
--referenceGenome [DNA.fa] \
--referencecDNA [cDNA.fa.gz] --referenceGtf [File.gtf] \
-config (change sequncing type to uppercase (10xv2 --> 10XV2))
```

## Output

|           | MR1 | MR2 | MR3 |      
| --------- | ----| ----| --- |       
| Alevin    |     |     |     |      
| Alevin-fry|     |     |     |
| kb-tools  |     |     |     |
| STARsolo  |     |     |     |

 MR1: Mapping rate with standard single-cell options
 MR2: Mapping rate with custom gtf file where 'trancript' is replaced with 'exon'
 MR3: Mapping rate with single-nucleus options 