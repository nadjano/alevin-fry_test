#!/usr/bin/env nextflow

sdrfFile = params.sdrf
resultsRoot = params.resultsRoot
referenceGenome = params.referenceGenome
referencecDNA = params.referencecDNA
referenceGtf = params.referenceGtf
protocol = params.protocol
outdir = "out_dir"

type = ('no_introns, introns')

REFERENCE_GENOME = Channel.fromPath( referenceGenome, checkIfExists: true ).first()
REFERENCE_GTF = Channel.fromPath( referenceGtf, checkIfExists: true ).first()
REFERENCE_CDNA = Channel.fromPath( referencecDNA, checkIfExists: true ).first()


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

// Read URIs from SDRF, generate target file names, and barcode locations

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
        path referencecDNA

    output:
        file("t2g.txt") into TRANSCRIPT_TO_GENE

    """
    cat Solanum_lycopersicum.SL3.0.cdna.all.fa.gz | awk '{if(\$1~/>/)print \$1"\t"\$4"\t"}' \\
     > t2g.txt; sed -i 's/>//g' t2g.txt; sed -i 's/gene://g' t2g.txt; \\
     sed -i 's/gene_symbol://g' t2g.txt
    """
}

// process download_fastqs {
    
//     conda "${baseDir}/envs/atlas-fastq-provider.yml"
    
//     maxForks params.maxConcurrentDownloads
//     time { 10.hour * task.attempt }
//     memory { 6.GB * task.attempt }

//     errorStrategy { task.attempt<=10 & task.exitStatus != 4 ? 'retry' : 'finish' } 
    
//     input:
//         set runId, cdnaFastqURI, barcodesFastqURI, cdnaFastqFile, barcodesFastqFile, val(barcodeLength), val(umiLength), val(end), val(cellCount), val(controlledAccess) from FASTQ_RUNS

//     output:
//         set val(runId), file("${cdnaFastqFile}"), file("${barcodesFastqFile}"), val(barcodeLength), val(umiLength), val(end), val(cellCount) into DOWNLOADED_FASTQS

//     """
//         if [ -n "$manualDownloadFolder" ] && [ -e $manualDownloadFolder/${cdnaFastqFile} ] && [ -e $manualDownloadFolder/${barcodesFastqFile} ]; then
//            ln -s $manualDownloadFolder/${cdnaFastqFile} ${cdnaFastqFile}
//            ln -s $manualDownloadFolder/${barcodesFastqFile} ${barcodesFastqFile}
//         elif [ -n "$manualDownloadFolder" ] && [ -e $manualDownloadFolder/${cdnaFastqFile} ] && [ ! -e $manualDownloadFolder/${barcodesFastqFile} ]; then
//             echo 'cDNA file $cdnaFastqFile is available locally, but barcodes file $barcodesFastqFile is not 1>&2
//             exit 2    
//         elif [ -n "$manualDownloadFolder" ] && [ ! -e $manualDownloadFolder/${cdnaFastqFile} ] && [ -e $manualDownloadFolder/${barcodesFastqFile} ]; then
//             echo 'cDNA file $cdnaFastqFile is not available locally, but barcodes file $barcodesFastqFile is 1>&2
//             exit 3 
//         elif [ "$controlledAccess" = 'yes' ]; then
//             echo "One or both of ${cdnaFastqFile}, ${barcodesFastqFile} are not available at $manualDownloadFolder/ for this controlled access experiment" 1>&2
//             exit 4   
//         else
//             confPart=''
//             if [ -n "$fastqProviderConfig" ] && [ -e "$fastqProviderConfig" ]; then
//                 confPart=" -c $fastqProviderConfig"
//             fi 
//             # Stop fastq downloader from testing different methods -assume the control workflow has done that 
//             export NOPROBE=1
        
//             fetchFastq.sh -f ${cdnaFastqURI} -t ${cdnaFastqFile} -m ${params.downloadMethod} \$confPart
            
//             # Allow for the first download also having produced the second output already
//             if [ ! -e ${barcodesFastqFile} ]; then
//                 fetchFastq.sh -f ${barcodesFastqURI} -t ${barcodesFastqFile} -m ${params.downloadMethod} \$confPart
//             fi
//         fi
//     """
// }

// // Group read files by run name, or by technical replicate group if specified

// if ( params.fields.containsKey('techrep')){

//     // If technical replicates are present, create a channel containing that info 

//     SDRF_FOR_TECHREP
//         .map{ row-> tuple(row["${params.fields.run}"], row["${params.fields.techrep}"]) }
//         .groupTuple()
//         .map{ row-> tuple( row[0], row[1][0]) }
//         .set{ TECHREPS }

//     // The target set of results will now be the technical replicate group number

//     SDRF_FOR_COUNT
//         .map{ row-> tuple(row["${params.fields.techrep}"]) }
//         .unique()
//         .count()
//         .set { TARGET_RESULT_COUNT }
    
//     // Now add the tech rep group to the run info, group by it, and create a
//     // tuple of files keyed by techrep group

//     TECHREPS.join( DOWNLOADED_FASTQS )
//         .groupTuple(by: 1)
//         .map{ row-> tuple( row[1], row[2].flatten(), row[3].flatten(), row[4][0], row[5][0], row[6][0], row[7][0]) }
//         .set{
//             FINAL_FASTQS
//         }
// }else{
//     DOWNLOADED_FASTQS.set{ FINAL_FASTQS }
    
//     SDRF_FOR_COUNT
//       .map{ row-> tuple(row["${params.fields.run}"]) }
//       .unique()
//       .count()
//       .set { TARGET_RESULT_COUNT }
// }

// FINAL_FASTQS.into{
//     FINAL_FASTQS_FOR_CONFIG
//     FINAL_FASTQS_FOR_ALEVIN
// }


// make splici transcript 
process build_splici {
   
   conda "${baseDir}/envs/pyroe.yml"

   memory { 2.GB * task.attempt }
    cpus 4

    errorStrategy { task.exitStatus !=2 && (task.exitStatus == 130 || task.exitStatus == 137 || task.attempt < 3)  ? 'retry' : 'ignore' }
    maxRetries 10

    input:
        file(referenceGenome) from REFERENCE_GENOME
        file(referenceGtf) from REFERENCE_GTF

    output:
        // publishDir "${outdir}"
        file("${outdir}/splici_fl45.fa") 
        file("${outdir}/splici_fl45*.tsv") into SPLICI_T2G_3
        

    """
    pyroe make-splici  ${referenceGenome}   ${referenceGtf}  50 ${outdir}

    """

}

process index_alevin{

    conda "${baseDir}/envs/alevin.yml"

    input:
        path reference from Channel.fromPath(["${outdir}/splici_fl45*.tsv", referenceGenome,])
        val x from Channel.from('splici', 'cDNA')
        // file("${outdir}/splici_fl45.fa") from SPLICI_REFERENCE
        // val(referenceType_splic)
    
    output:
        path("alevin_index_$x") into INDEX_ALEVIN

    """
    salmon index --transcript $reference   -i alevin_index_$x
    """

}

process t2g_splici{
    input:
        file("${outdir}/splici_fl45*.tsv") from SPLICI_T2G_3
      
    
    output:
        file("${outdir}/splici_t2g.tsv") into SPLICI_T2G

    """
    cat "${outdir}/splici_fl45*.tsv" | awk  '{print\$1"\t"\$1}'  > "${outdir}/splici_t2g.tsv"
    """

}

// process cDNA_index_alevin{

//     conda "${baseDir}/envs/alevin.yml"

//     input:
//         file(referencecDNA) from REFERENCE_CDNA
    
//     output:
//         path("alevin_index_cDNA") into INDEX_ALEVIN_CDNA

//     """
//     salmon index --transcript ${referencecDNA}   -i alevin_index_cDNA
//     """

// }

// ALEVIN_INDEX = Channel.from(['cDNA', "alevin_index_cDNA", TRANSCRIPT_TO_GENE],['splici', "alevin_index_splici", "${outdir}/splici_t2g.tsv"] )

// process alevin_config {

//     input:
//         set val(runId), file("cdna*.fastq.gz"), file("barcodes*.fastq.gz"), val(barcodeLength), val(umiLength), val(end), val(cellCount) from FINAL_FASTQS_FOR_CONFIG

//     output:
//         set val(runId), stdout into ALEVIN_CONFIG
    
//     script:

//         def barcodeConfig = ''

//         if ( params.containsKey(protocol) ){

//             canonicalProtocol = params.get(protocol)
//             alevinType = canonicalProtocol.alevinType

//             // Non-standard barcode config is supplied as a custom method

//             if ( alevinType == 'custom' || "${canonicalProtocol.barcodeLength}" != barcodeLength || "${canonicalProtocol.umiLength}" != umiLength || "${canonicalProtocol.end}" != end ){
//                 barcodeConfig = "--barcodeLength ${barcodeLength} --umiLength ${umiLength} --end ${end}" 

//             }else{
//                 barcodeConfig = "--$alevinType"
//             }
//             barcodeConfig = "-l ${canonicalProtocol.libType} $barcodeConfig" 
//         }

//         """
//         if [ -z "$barcodeConfig" ]; then
//             echo Input of $protocol results is misconfigured 1>&2
//             exit 1
//         fi
//         # Also check barcode read lengths and return non-0 if they're not what they should be
//         targetLen=\$(($umiLength + $barcodeLength))
//         barcodesGood=0
//         set +e
//         while read -r l; do
//             checkBarcodeRead.sh -r \$(readlink -f \$l) -b $barcodeLength -u $umiLength -n 1000000 1>&2
//             if [ \$? -ne 0 ]; then
//                 barcodesGood=1
//             fi
//         done <<< "\$(ls barcodes*.fastq.gz)"
//         set -e
        
//         echo -n "$barcodeConfig"
//         exit \$barcodesGood
//         """
// }


// ALEVIN_INDEX = Channel.of(['cDNA', "alevin_index_cDNA", TRANSCRIPT_TO_GENE],['splici', "alevin_index_splici","${outdir}/splici_t2g.tsv" ] )


// process alevin {

//     conda "${baseDir}/envs/alevin.yml"
    
//     cache 'deep'

//     memory { 2.GB * task.attempt }
//     cpus 4

//     errorStrategy { task.exitStatus !=2 && (task.exitStatus == 130 || task.exitStatus == 137 || task.attempt < 3)  ? 'retry' : 'ignore' }
//     maxRetries 10

//     input:
//         set val(runId), file("cdna*.fastq.gz"), file("barcodes*.fastq.gz"), val(barcodeLength), val(umiLength), val(end), val(cellCount), val(barcodeConfig) from FINAL_FASTQS_FOR_ALEVIN.join(ALEVIN_CONFIG)
//         tuple val(ref_type), path('index_dir'), path('t2g.tsv') from ALEVIN_INDEX

//     output:
//         set val(runId), file("${runId}"),  file("${runId}/alevin/raw_cb_frequency.txt") into ALEVIN_RESULTS

//     """
//     salmon alevin ${barcodeConfig} -1 \$(ls barcodes*.fastq.gz | tr '\\n' ' ') -2 \$(ls cdna*.fastq.gz | tr '\\n' ' ') \
//         -i ${index_dir} -p ${task.cpus} -o ${runId}_tmp --tgMap ${t2g.tsv} --dumpFeatures --keepCBFraction 1 \
//         --freqThreshold ${params.minCbFreq} --dumpMtx
//     min_mapping=\$(grep "percent_mapped" ${runId}_tmp/aux_info/meta_info.json | sed 's/,//g' | awk -F': ' '{print \$2}' | sort -n | head -n 1)   
//     if [ "\${min_mapping%.*}" -lt "${params.minMappingRate}" ]; then
//         echo "Minimum mapping rate (\$min_mapping) is less than the specified threshold of ${params.minMappingRate}" 1>&2
//         exit 1 
//     fi
 
//     mv ${runId}_tmp ${runId}
//     """
// }



// // Call the download script to retrieve run fastqs

