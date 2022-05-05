#!/usr/bin/env nextflow

sdrfFile = params.sdrf
resultsRoot = params.resultsRoot
referenceGenome = params.referenceGenome
referencecDNA = params.referencecDNA
referenceGtf = params.referenceGtf
protocol = params.protocol
outdir = "out_dir"
ref_type = ['splici', 'cDNA']


REFERENCE_CDNA = Channel.fromPath( referencecDNA)
REFERENCE_GTF = Channel.fromPath( referenceGtf)


manualDownloadFolder =''
if ( params.containsKey('manualDownloadFolder')){
    manualDownloadFolder = params.manualDownloadFolder
}

fastqProviderConfig = ''
if ( params.containsKey('fastqProviderConfig')){
    fastqProviderConfig = params.fastqProviderConfig
}
Channel
    .fromPath(sdrfFile, checkIfExists: true)
    .splitCsv(header:true, sep:"\t")
    .filter{ row -> (! row.containsKey(params.fields.quality)) || ( row["${params.fields.quality}"].toLowerCase() != 'not ok') }
    .into {
        SDRF_FOR_FASTQS
        SDRF_FOR_STRAND
        SDRF_FOR_TECHREP
        SDRF_FOR_COUNT
    }

//  Read URIs from SDRF, generate target file names, and barcode locations

SDRF_FOR_FASTQS
    .map{ row-> 
      controlled_access='no'
      if (  params.fields.containsKey('controlled_access')){
        controlled_access=row["${params.fields.controlled_access}"]
      }
      tuple(row["${params.fields.run}"], row["${params.fields.cdna_uri}"], row["${params.fields.cell_barcode_uri}"], file(row["${params.fields.cdna_uri}"]).getName(), file(row["${params.fields.cell_barcode_uri}"]).getName(), row["${params.fields.cell_barcode_size}"], row["${params.fields.umi_barcode_size}"], row["${params.fields.end}"], row["${params.fields.cell_count}"], controlled_access) 
    }    
    .set { FASTQ_RUNS }

process make_t2g_file {

    input:
        path reference from REFERENCE_CDNA

    output:
        file("t2g_cDNA.txt")

    """
    cat ${reference} | awk '{if(\$1~/>/)print \$1"\t"\$4"\t"}' \\
     > t2g_cDNA.txt; sed -i 's/>//g' t2g_cDNA.txt; sed -i 's/gene://g' t2g_cDNA.txt; \\
     sed -i 's/gene_symbol://g' t2g_cDNA.txt
    """
}

process download_fastqs {
    
    conda "${baseDir}/envs/atlas-fastq-provider.yml"
    
    maxForks params.maxConcurrentDownloads
    time { 10.hour * task.attempt }
    memory { 6.GB * task.attempt }

    errorStrategy { task.attempt<=10 & task.exitStatus != 4 ? 'retry' : 'finish' } 
    
    input:
        set runId, cdnaFastqURI, barcodesFastqURI, cdnaFastqFile, barcodesFastqFile, val(barcodeLength), val(umiLength), val(end), val(cellCount), val(controlledAccess) from FASTQ_RUNS

    output:
        set val(runId), file("${cdnaFastqFile}"), file("${barcodesFastqFile}"), val(barcodeLength), val(umiLength), val(end), val(cellCount) into DOWNLOADED_FASTQS

    """
        if [ -n "$manualDownloadFolder" ] && [ -e $manualDownloadFolder/${cdnaFastqFile} ] && [ -e $manualDownloadFolder/${barcodesFastqFile} ]; then
           ln -s $manualDownloadFolder/${cdnaFastqFile} ${cdnaFastqFile}
           ln -s $manualDownloadFolder/${barcodesFastqFile} ${barcodesFastqFile}
        elif [ -n "$manualDownloadFolder" ] && [ -e $manualDownloadFolder/${cdnaFastqFile} ] && [ ! -e $manualDownloadFolder/${barcodesFastqFile} ]; then
            echo 'cDNA file $cdnaFastqFile is available locally, but barcodes file $barcodesFastqFile is not 1>&2
            exit 2    
        elif [ -n "$manualDownloadFolder" ] && [ ! -e $manualDownloadFolder/${cdnaFastqFile} ] && [ -e $manualDownloadFolder/${barcodesFastqFile} ]; then
            echo 'cDNA file $cdnaFastqFile is not available locally, but barcodes file $barcodesFastqFile is 1>&2
            exit 3 
        elif [ "$controlledAccess" = 'yes' ]; then
            echo "One or both of ${cdnaFastqFile}, ${barcodesFastqFile} are not available at $manualDownloadFolder/ for this controlled access experiment" 1>&2
            exit 4   
        else
            confPart=''
            if [ -n "$fastqProviderConfig" ] && [ -e "$fastqProviderConfig" ]; then
                confPart=" -c $fastqProviderConfig"
            fi 
            # Stop fastq downloader from testing different methods -assume the control workflow has done that 
            export NOPROBE=1
        
            fetchFastq.sh -f ${cdnaFastqURI} -t ${cdnaFastqFile} -m ${params.downloadMethod} \$confPart
            
            # Allow for the first download also having produced the second output already
            if [ ! -e ${barcodesFastqFile} ]; then
                fetchFastq.sh -f ${barcodesFastqURI} -t ${barcodesFastqFile} -m ${params.downloadMethod} \$confPart
            fi
        fi
    """
}

// Group read files by run name, or by technical replicate group if specified

if ( params.fields.containsKey('techrep')){

    // If technical replicates are present, create a channel containing that info 

    SDRF_FOR_TECHREP
        .map{ row-> tuple(row["${params.fields.run}"], row["${params.fields.techrep}"]) }
        .groupTuple()
        .map{ row-> tuple( row[0], row[1][0]) }
        .set{ TECHREPS }

    // The target set of results will now be the technical replicate group number

    SDRF_FOR_COUNT
        .map{ row-> tuple(row["${params.fields.techrep}"]) }
        .unique()
        .count()
        .set { TARGET_RESULT_COUNT }
    
    // Now add the tech rep group to the run info, group by it, and create a
    // tuple of files keyed by techrep group

    TECHREPS.join( DOWNLOADED_FASTQS )
        .groupTuple(by: 1)
        .map{ row-> tuple( row[1], row[2].flatten(), row[3].flatten(), row[4][0], row[5][0], row[6][0], row[7][0]) }
        .set{
            FINAL_FASTQS
        }
}else{
    DOWNLOADED_FASTQS.set{ FINAL_FASTQS }
    
    SDRF_FOR_COUNT
      .map{ row-> tuple(row["${params.fields.run}"]) }
      .unique()
      .count()
      .set { TARGET_RESULT_COUNT }
}

FINAL_FASTQS.into{
    FINAL_FASTQS_FOR_CONFIG
    FINAL_FASTQS_FOR_ALEVIN
}


// make splici transcript 
process build_splici {
   
   conda "${baseDir}/envs/pyroe.yml"

   memory { 2.GB * task.attempt }
    cpus 4

    errorStrategy { task.exitStatus !=2 && (task.exitStatus == 130 || task.exitStatus == 137 || task.attempt < 3)  ? 'retry' : 'ignore' }
    maxRetries 10

    input:
        path(referenceGenome) 
        path(referenceGtf) 

    output:
        // publishDir "${outdir}"
        path("${outdir}/splici_fl45.fa") into splici_fasta
        path("${outdir}/splici_fl45*.tsv") into T2G_3
        
    """
    pyroe make-splici  ${referenceGenome}   ${referenceGtf}  50 ${outdir}

    """
}


process index_alevin {

    conda "${baseDir}/envs/alevin.yml"

    input:
        path reference from splici_fasta
        
    output:
        path("alevin_index_splici")

    """
    salmon index --transcript $reference   -i alevin_index_splici
    """

 }

 process index_alevin_cDNA {

    conda "${baseDir}/envs/alevin.yml"

    input:
        path reference from referencecDNA
        
    output:
        path("alevin_index_cDNA") 

    """
    salmon index --transcript $reference   -i alevin_index_cDNA
    """

 }

process t2g_splici{
    input:
        file("${outdir}/splici_fl45*.tsv") from T2G_3
      
    
    output:
        file("t2g_splici.txt") 

    """
    cat "${outdir}/splici_fl45*.tsv" | awk  '{print\$1"\t"\$1}'  > "t2g_splici.txt"
    """

}



process alevin_config {

    input:
        set val(runId), file("cdna*.fastq.gz"), file("barcodes*.fastq.gz"), val(barcodeLength), val(umiLength), val(end), val(cellCount) from FINAL_FASTQS_FOR_CONFIG

    output:
        set val(runId), stdout into ALEVIN_CONFIG
    
    script:

        def barcodeConfig = ''

        if ( params.containsKey(protocol) ){

            canonicalProtocol = params.get(protocol)
            alevinType = canonicalProtocol.alevinType

            // Non-standard barcode config is supplied as a custom method

            if ( alevinType == 'custom' || "${canonicalProtocol.barcodeLength}" != barcodeLength || "${canonicalProtocol.umiLength}" != umiLength || "${canonicalProtocol.end}" != end ){
                barcodeConfig = "--barcodeLength ${barcodeLength} --umiLength ${umiLength} --end ${end}" 

            }else{
                barcodeConfig = "--$alevinType"
            }
            barcodeConfig = "-l ${canonicalProtocol.libType} $barcodeConfig" 
        }

        """
        if [ -z "$barcodeConfig" ]; then
            echo Input of $protocol results is misconfigured 1>&2
            exit 1
        fi
        # Also check barcode read lengths and return non-0 if they're not what they should be
        targetLen=\$(($umiLength + $barcodeLength))
        barcodesGood=0
        set +e
        while read -r l; do
            checkBarcodeRead.sh -r \$(readlink -f \$l) -b $barcodeLength -u $umiLength -n 1000000 1>&2
            if [ \$? -ne 0 ]; then
                barcodesGood=1
            fi
        done <<< "\$(ls barcodes*.fastq.gz)"
        set -e
        
        echo -n "$barcodeConfig"
        exit \$barcodesGood
        """
}



ALEVIN_INDEX = Channel.fromPath('alevin_index_*')
T2G = Channel.fromPath('t2g_*')


process alevin {

    conda "${baseDir}/envs/alevin.yml"
    
    cache 'deep'

    memory { 2.GB * task.attempt }
    cpus 4

    errorStrategy { task.exitStatus !=2 && (task.exitStatus == 130 || task.exitStatus == 137 || task.attempt < 3)  ? 'retry' : 'ignore' }
    maxRetries 10

    input:
        set val(runId), file("cdna*.fastq.gz"), file("barcodes*.fastq.gz"), val(barcodeLength), val(umiLength), val(end), val(cellCount), val(barcodeConfig) from FINAL_FASTQS_FOR_ALEVIN.join(ALEVIN_CONFIG)
        path index_dir from ALEVIN_INDEX
        path t2g from T2G

    output:
        set val(runId), file("${runId}"),  file("${runId}/alevin/raw_cb_frequency.txt") into ALEVIN_RESULTS

    """
    salmon alevin ${barcodeConfig} -1 \$(ls barcodes*.fastq.gz | tr '\\n' ' ') -2 \$(ls cdna*.fastq.gz | tr '\\n' ' ') \
        -i ${index_dir} -p ${task.cpus} -o ${runId}_tmp --tgMap ${t2g.tsv} --dumpFeatures --keepCBFraction 1 \
        --freqThreshold ${params.minCbFreq} --dumpMtx
    min_mapping=\$(grep "percent_mapped" ${runId}_tmp/aux_info/meta_info.json | sed 's/,//g' | awk -F': ' '{print \$2}' | sort -n | head -n 1)   
    if [ "\${min_mapping%.*}" -lt "${params.minMappingRate}" ]; then
        echo "Minimum mapping rate (\$min_mapping) is less than the specified threshold of ${params.minMappingRate}" 1>&2
        exit 1 
    fi
 
    mv ${runId}_tmp ${runId}
    """
}

// process index_star {

//     conda "${baseDir}/envs/star.yml"

//     input:
//         val genomeDir
//         path(referenceGenome) 
//         path(referenceGtf) 
//     output:
//         file("${genomeDir}")
        

// }

// index for STARSolo
// environment STARsolo
// STAR --runMode genomeGenerate --genomeDir STAR_index --genomeFastaFiles Drosophila_melanogaster.BDGP6.32.dna.toplevel.fa --sjdbGTFfile Drosophila_melanogaster.BDGP6.32.106.gtf --genomeSAindexNbases 12

// // Call the download script to retrieve run fastqs

