#!/usr/bin/env nextflow
resultsRoot = params.resultsRoot
referenceFasta = params.referenceFasta
referenceGtf = params.referenceGtf
fastq_path = params.fastq_path



REFERENCE_GENOME = Channel.fromPath(referenceFasta, checkIfExists: true ).first()
REFERENCE_GTF = Channel.fromPath(referenceGtf, checkIfExists: true ).first()

REFERENCE_GTF.into {
    REFERENCE_GTF_FOR_CDNA
    REFERENCE_GTF_FOR_T2G
}
// generate splici index if it does not exits
process make_cDNA_from_Genome {

    cache 'lenient'
   
    conda "${baseDir}/envs/gff_read.yml"

    memory { 10.GB * task.attempt }

    errorStrategy { task.exitStatus !=2 && (task.exitStatus == 130 || task.exitStatus == 137 || task.attempt < 3)  ? 'retry' : 'ignore' }
    maxRetries 10

    input:
        path referenceGenome from REFERENCE_GENOME
        path referenceGtf from REFERENCE_GTF_FOR_CDNA

    output:
        path("transcriptome.fa") into CUSTOM_CDNA
       
    """
    gffread -F -w transcriptome.fa -g ${referenceGenome}  ${referenceGtf} 
    """
}

process make_splici {
    publishDir "t2g_alevin_fry/${species}", mode: 'copy', overwrite: true

    cache 'lenient'
   
    conda "${baseDir}/envs/pyroe.yml"

    memory { 10.GB * task.attempt }

    errorStrategy { task.exitStatus !=2 && (task.exitStatus == 130 || task.exitStatus == 137 || task.attempt < 3)  ? 'retry' : 'ignore' }
    maxRetries 10

    input:
        path referenceGenome from REFERENCE_GENOME
        path referenceGtf from REFERENCE_GTF

    output:
        path("splici_out/splici_fl*.fa") into SPLICI_FASTA_FOR_FRY
        path("splici_out/splici_fl*.tsv") into T2G_3_FOR_FRY

     
    """
    pyroe make-splici ${referenceGenome} ${referenceGtf} 90 splici_out
    """
}

process make_t2g {

    input:
    path referenceGtf from REFERENCE_GTF_FOR_T2G

    output:
    path "t2g_transcriptome.txt" into TRANSCRIPT_TO_GENE

    """
    cat ${referenceGtf} | grep -vE "^#" | awk '\$3=="transcript" {split(\$0, array, "transcript_id"); print array[2]}' | awk '{split(\$0, array, ";"); print array[1]}' | sed 's/"//g' | sed 's/^ *//g' | sed 's/transcript://g' > t.txt

    cat ${referenceGtf} | grep -vE "^#" | awk '\$3=="transcript" {split(\$0, array, "gene_id"); print array[2]}' | awk '{split(\$0, array, ";"); print array[1]}' | sed 's/"//g' | sed 's/^ *//g'| sed 's/gene://g'  > g.txt

    paste t.txt g.txt > t2g_transcriptome.txt
    """

}

process index_for_alevin {
    
    cache 'deep'
    memory { 10.GB * task.attempt }
    errorStrategy { task.exitStatus !=2 && (task.exitStatus == 130 || task.exitStatus == 137 || task.attempt < 3)  ? 'retry' : 'ignore' }
    maxRetries 10
    conda "${baseDir}/envs/alevin.yml"

    input:
        path reference from CUSTOM_CDNA
        
    output:
        path "alevin_index" into ALEVIN_INDEX
    

    """
    salmon index --transcript ${reference}  -i alevin_index
    """

 }

process index_for_alevin_fry {
    
    cache 'deep'
    memory { 10.GB * task.attempt }
    errorStrategy { task.exitStatus !=2 && (task.exitStatus == 130 || task.exitStatus == 137 || task.attempt < 3)  ? 'retry' : 'ignore' }
    maxRetries 10
    conda "${baseDir}/envs/alevin-fry_2.yml"

    input:
        path reference from SPLICI_FASTA_FOR_FRY
        
    output:
        path "alevin_index_splici" into ALEVIN_FRY_INDEX_SPLICI
    
  
    """
    salmon index --transcript ${reference}  -i alevin_index_splici
    """

 }

 process alevin_fry {
    
    cache 'lenient'
    cpus 8
    memory { 20.GB * task.attempt }
    errorStrategy { task.exitStatus !=2 && (task.exitStatus == 130 || task.exitStatus == 137 || task.attempt < 3)  ? 'retry' : 'ignore' }
    maxRetries 10
    conda "${baseDir}/envs/alevin-fry_2.yml"
    input:
        path "alevin_index_splici" from ALEVIN_FRY_INDEX_SPLICI
        path "t2g_cDNA.txt" from T2G_3_FOR_FRY
       
    output:
        set val(runId), path("test_ALEVIN_fry_quant"), path("test_ALEVIN_fry_quant/featureDump.txt") into ALEVIN_FRY_RESULTS_SPLICI
      

    """
    salmon alevin -l ISR --chromiumV3 --sketch -1 \$(ls ${params.fastq_path}/barcodes*.fastq.gz | tr '\\n' ' ') -2 \$(ls ${params.fastq_path}cdna*.fastq.gz | tr '\\n' ' ') \
        -i alevin_index_splici -p 12 -o test_ALEVIN_fry_map 
   
    alevin-fry generate-permit-list --input test_ALEVIN_fry_map -d fw --output-dir test_ALEVIN_fry_quant  --force-cells 50000
 
    alevin-fry collate -i test_ALEVIN_fry_quant -r test_ALEVIN_fry_map -t 16
    alevin-fry quant -i test_ALEVIN_fry_quant -m t2g_cDNA.txt -t 16 -r cr-like-em -o test_ALEVIN_fry_quant --use-mtx

        
    """
}

ALEVIN_FRY_RESULTS_SPLICI
    .into{
        ALEVIN_RESULTS_FOR_QC
        ALEVIN_RESULTS_FOR_PROCESSING
        ALEVIN_RESULTS_FOR_OUTPUT
    }

process mtx_alevin_fry_to_mtx {

    conda "${baseDir}/envs/parse_alevin_fry.yml"

    memory { 10.GB * task.attempt }
   

    input:
    set val("alevin-fry"), path("test_ALEVIN_fry_quant"), file(rawBarcodeFreq) from ALEVIN_RESULTS_FOR_PROCESSING

    output:
    set val("alevin-fry"), path("counts_mtx_test") into ALEVIN_FRY_MTX
    // file("counts_mtx_${protocol}") into PROTOCOL_COUNT_MATRICES


    """
    alevinFryMtxTo10x.py --cell_prefix test- test_ALEVIN_fry_quant counts_mtx_test
    """      
}
ALEVIN_FRY_MTX
    .into{
        ALEVIN_FRY_MTX_FOR_QC
        ALEVIN_FRY_MTX_FOR_EMPTYDROPS
        ALEVIN_FRY_MTX_FOR_OUTPUT
        ALEVIN_FRY_MTX_FOR_MERGE

    }



 process alevin {

    conda "${baseDir}/envs/alevin.yml"
    
    cache 'lenient'

    memory { 20.GB * task.attempt }
    cpus 12

    errorStrategy { task.exitStatus !=2 && (task.exitStatus == 130 || task.exitStatus == 137 || task.attempt < 3)  ? 'retry' : 'ignore' }
    maxRetries 10

    input:
        path(transcriptToGene) from TRANSCRIPT_TO_GENE
        path(index) from ALEVIN_INDEX

    output:
        set  file("test"),  file("test/alevin/raw_cb_frequency.txt") into ALEVIN_RESULTS

    """
    salmon alevin -l ISR --chromiumV3 --sketch -1 \$(ls ${params.fastq_path}/barcodes*.fastq.gz | tr '\\n' ' ') -2 \$(ls ${params.fastq_path}cdna*.fastq.gz | tr '\\n' ' ') \
        -i ${index} -p 12 -o test_tmp --tgMap ${transcriptToGene} --dumpFeatures --keepCBFraction 1 \
        --freqThreshold 10 --dumpMtx
 
    mv test_tmp test
    """
}

ALEVIN_RESULTS
    .into{
        ALEVIN_RESULTS_FOR_QC
        ALEVIN_RESULTS_FOR_PROCESSING
        ALEVIN_RESULTS_FOR_OUTPUT
    }

process alevin_to_mtx {

    conda "${baseDir}/envs/parse_alevin_fry.yml"

    // conda "/nfs/production/irene/ma/users/nnolte/conda/envs/parse_alevin_fry"
    
    memory { 10.GB * task.attempt }
    errorStrategy { task.exitStatus == 130 || task.exitStatus == 137 ? 'retry' : 'finish' }
    maxRetries 20

    input:
        set file(alevinResult), file(rawBarcodeFreq) from ALEVIN_RESULTS_FOR_PROCESSING
        val('alevin')

    output:
        set val('alevin'), file("test_counts_mtx") into ALEVIN_MTX

    """
    alevinMtxTo10x.py --cell_prefix test- $alevinResult test_counts_mtx
    """ 
}

ALEVIN_MTX
    .into{
        ALEVIN_MTX_FOR_QC
        ALEVIN_MTX_FOR_EMPTYDROPS
        ALEVIN_MTX_FOR_OUTPUT
        ALEVIN_MTX_FOR_MERGE
    
    }




// // Convert Alevin output to MTX. There will be one of these for every run, or
// // technical replicate group of runs




// Make a diagnostic plot



process droplet_qc_plot{
    
    conda "${baseDir}/envs/alevin.yml"
    
    memory { 10.GB * task.attempt }
    errorStrategy { task.exitStatus == 130 || task.exitStatus == 137 ? 'retry' : 'finish' }
    maxRetries 20

    input:
        set val(type), path(mtx) from ALEVIN_MTX_FOR_QC.join(ALEVIN_FRY_MTX_FOR_QC)

    output:
        set file("${type}.png") into ALEVIN_QC_PLOTS

    """
    dropletBarcodePlot.R --mtx-matrix $mtx --output ${type}.png
    """ 
}

// // Remove empty droplets from Alevin results


process remove_empty_drops {
    
    conda "${baseDir}/envs/dropletutils.yml"

    memory { 10.GB * task.attempt }
    errorStrategy { task.exitStatus == 130 || task.exitStatus == 137 ? 'retry' : 'ignore' }
    maxRetries 20
   
    input:
        set val(type), file(countsMtx) from ALEVIN_MTX_FOR_EMPTYDROPS

    output:
        set val(type), file("${type}_nonempty.rds") into NONEMPTY_RDS

    """
        dropletutils-read-10x-counts.R -s counts_mtx -c TRUE -o matrix.rds
        dropletutils-empty-drops.R -i matrix.rds --lower ${params.emptyDrops.lower} --niters ${params.emptyDrops.nIters} --filter-empty ${params.emptyDrops.filterEmpty} \
        --filter-fdr ${params.emptyDrops.filterFdr} --ignore ${params.minCbFreq} -o ${type}_nonempty.rds -t nonempty.txt
    """
}

// // Convert R matrix object with filtered cells back to .mtx

process rds_to_mtx{

    conda "${baseDir}/envs/dropletutils.yml"

    memory { 10.GB * task.attempt }
    errorStrategy { task.exitStatus == 130 || task.exitStatus == 137 ? 'retry' : 'finish' }
    maxRetries 20
   
    input:
        set val(type), file(rds) from NONEMPTY_RDS

    output:
        set val(type), file("${type}_counts_mtx_nonempty") into NONEMPTY_MTX

    """ 
        #!/usr/bin/env Rscript
        
        suppressPackageStartupMessages(require(DropletUtils))

        counts_sce <- readRDS('$rds')
        write10xCounts(assays(counts_sce)[[1]], path = '${type}_counts_mtx_nonempty', barcodes = colData(counts_sce)\$Barcode, gene.id = rownames(counts_sce))
    """
}





// // Compile raw results with raw and emptyDrops-filtered MTX

// ALEVIN_RESULTS_FOR_OUTPUT
//     .join(ALEVIN_MTX_FOR_OUTPUT)
//     .join(NONEMPTY_MTX)
//     .join(ALEVIN_QC_PLOTS)
//     .set{ COMPILED_RESULTS }

// process compile_results{

//     publishDir "$resultsRoot/alevin", mode: 'copy', overwrite: true
    
//     input:
//         set val(runId), file('raw_alevin'), file(rawBarcodeFreq), file(countsMtx), file(countsMtxNonempty), file(qcPlot) from COMPILED_RESULTS

//     output:
//         set val(runId), file("$runId") into RESULTS_FOR_COUNTING

//     """
//         mkdir -p raw_alevin/alevin/mtx
//         cp -P $countsMtx $countsMtxNonempty raw_alevin/alevin/mtx 
//         mkdir -p raw_alevin/alevin/qc
//         cp -P $qcPlot raw_alevin/alevin/qc
//         cp -P raw_alevin $runId
//     """
// }

// // Check the total number of runs we have 

