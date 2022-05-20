#!/usr/bin/env nextflow

sdrfFile = params.sdrf
resultsRoot = params.resultsRoot
referenceGenome = params.referenceGenome
referencecDNA = params.referencecDNA
referenceGtf = params.referenceGtf
protocol = params.protocol


REFERENCE_CDNA = Channel.fromPath(referencecDNA,checkIfExists: true).first()
REFERENCE_GTF = Channel.fromPath( referenceGtf,checkIfExists: true ).first()
REFERENCE_GENOME = Channel.fromPath( referenceGenome,checkIfExists: true ).first()

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
        path "t2g_cDNA.txt" into T2G_CDNA
        path "t2g_cDNA.txt" into T2G_FOR_FRY


    """
    zcat ${reference} | awk '{if(\$1~/>/)print \$1"\t"\$4}' \\
     > t2g_cDNA.txt; sed -i 's/>//g' t2g_cDNA.txt; sed -i 's/gene://g' t2g_cDNA.txt; \\
     sed -i 's/gene_symbol://g' t2g_cDNA.txt
    """
}
 
T2G_CDNA.into {
    T2G_CDNA_FOR_ALEVIN
    T2G_CDNA_FOR_ALEVIN_FRY
}

process download_fastqs {
    
    conda "${baseDir}/envs/atlas-fastq-provider.yml"
    
    maxForks params.maxConcurrentDownloads
    time { 10.hour * task.attempt }
    memory { 20.GB * task.attempt }
    cache 'lenient'

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
    FINAL_FASTQS_FOR_ALEVIN_SPLICI
    FINAL_FASTQS_FOR_ALEVIN_CDNA
    FINAL_FASTQS_FOR_STAR
    FINAL_FASTQS_FOR_KB_TOOLS
    FINAL_FASTQS_FOR_KB_TOOLS_SPLICI
    FINAL_FASTQS_FOR_KB_TOOLS_PRERNA
    FINAL_FASTQS_FOR_ALEVIN_FRY
    FINAL_FASTQS_FOR_ALEVIN_FRY_CDNA
    FINAL_FASTQS_FOR_ALEVIN_FRY_TRANSCRIPTOME
}


// make splici transcript (spliced transcript + introns) 
// see https://github.com/COMBINE-lab/alevin-fry
process build_splici {
    cache 'lenient'
   
    conda "${baseDir}/envs/pyroe.yml"

    memory { 20.GB * task.attempt }
    cpus 4

    errorStrategy { task.exitStatus !=2 && (task.exitStatus == 130 || task.exitStatus == 137 || task.attempt < 3)  ? 'retry' : 'ignore' }
    maxRetries 10

    input:
        path referenceGenome from REFERENCE_GENOME
        path referenceGtf from REFERENCE_GTF

    output:
        path("splici_out/splici_fl*.fa") into SPLICI_FASTA
        path("splici_out/splici_fl*.fa") into SPLICI_FASTA_FOR_FRY
        path("splici_out/splici_fl*.tsv") into T2G_3
        path("splici_out/splici_fl*.tsv") into T2G_3_FOR_FRY
        
    """
    pyroe make-splici ${referenceGenome} ${referenceGtf} 90 splici_out
    """
}



// build index for alevin with splici transcript
process index_alevin_MR2 {
    cache 'lenient'

    memory { 40.GB * task.attempt }
    cpus 4
    errorStrategy { task.exitStatus !=2 && (task.exitStatus == 130 || task.exitStatus == 137 || task.attempt < 3)  ? 'retry' : 'ignore' }
    maxRetries 10


    conda "${baseDir}/envs/alevin.yml"

    input:
        path reference from SPLICI_FASTA
        
    output:
        path "alevin_index_splici" into ALEVIN_INDEX_SPLICI

    """
    salmon index --transcript ${reference}   -i alevin_index_splici
    """

 }

// build index for alevin with cdDNA 
 process index_alevin_MR1 {
    cache 'lenient'
    memory { 40.GB * task.attempt }
    cpus 4
    errorStrategy { task.exitStatus !=2 && (task.exitStatus == 130 || task.exitStatus == 137 || task.attempt < 3)  ? 'retry' : 'ignore' }
    maxRetries 10


    conda "${baseDir}/envs/alevin.yml"

    input:
        path reference from REFERENCE_CDNA
        
    output:
        path "alevin_index_cDNA" into ALEVIN_INDEX_CDNA

    """
    salmon index --transcript ${reference}   -i alevin_index_cDNA
    """

 }

// extract first two collumns from t2g file produced by pyroe
process t2g_splici{
    cache 'lenient'
    memory { 20.GB * task.attempt }
    cpus 4
    input:
        file("splici_out/splici_fl*.tsv") from T2G_3
      
    
    output:
        path "t2g_splici.txt" into T2G_SPLICI

    """
    cat splici_out/splici_fl*.tsv | awk  '{print\$1"\t"\$2}'  > t2g_splici.txt
    """

}
// Derive Alevin barcodeconfig

process alevin_config {
    cache 'lenient'
    memory { 20.GB * task.attempt }
    cpus 4

    input:
        set val(runId), file("cdna*.fastq.gz"), file("barcodes*.fastq.gz"), val(barcodeLength), val(umiLength), val(end), val(cellCount) from FINAL_FASTQS_FOR_CONFIG

    output:
        set val(runId), stdout into CONFIG
    
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

CONFIG
    .into{
        ALEVIN_CONFIG
        ALEVIN_CONFIG_SPLICI
        ALEVIN_CONFIG_CDNA
        STAR_CONFIG
        KB_CONFIG
        KB_CONFIG_SPLICI
        KB_CONFIG_PRERNA
        ALEVIN_FRY_CONFIG
        ALEVIN_FRY_CONFIG_CDNA
        ALEVIN_FRY_CONFIG_TRANSCRIPTOME
    }

// run Alevin for splici 
process alevin_MR2 {
    cache 'lenient'
    memory { 100.GB * task.attempt }
    cpus 4
    errorStrategy { task.exitStatus !=2 && (task.exitStatus == 130 || task.exitStatus == 137 || task.attempt < 3)  ? 'retry' : 'ignore' }
    maxRetries 10


    conda "${baseDir}/envs/alevin.yml"
    
    // cache 'deep'

    // memory { 2.GB * task.attempt }
    cpus 4

    // errorStrategy { task.exitStatus !=2 && (task.exitStatus == 130 || task.exitStatus == 137 || task.attempt < 3)  ? 'retry' : 'ignore' }
    // maxRetries 10

    input:
        set val(runId), file("cdna*.fastq.gz"), file("barcodes*.fastq.gz"), val(barcodeLength), val(umiLength), val(end), val(cellCount), val(barcodeConfig) from FINAL_FASTQS_FOR_ALEVIN_SPLICI.join(ALEVIN_CONFIG_SPLICI)
        path "alevin_index_splici" from ALEVIN_INDEX_SPLICI
        path "t2g_splici.txt" from T2G_SPLICI

    output:
        // publishDir path "${runId}_ALEVIN"
        set stdout, val(runId), file("${runId}_splici_ALEVIN") into ALEVIN_RESULTS_SPLICI 
        set val(runId), stdout into ALEVIN_SPLICI_MAPPING
        path ".command.log" into MEM_ALEVIN_MR2
    

    """
    salmon alevin ${barcodeConfig} -1 \$(ls barcodes*.fastq.gz | tr '\\n' ' ') -2 \$(ls cdna*.fastq.gz | tr '\\n' ' ') \
        -i alevin_index_splici -p ${task.cpus} -o ${runId}_splici_ALEVIN_tmp --tgMap t2g_splici.txt --dumpFeatures --keepCBFraction 1 \
        --freqThreshold ${params.minCbFreq} --dumpMtx > /dev/null
    mapping_rate=\$(grep "mapping_rate" ${runId}_splici_ALEVIN_tmp/aux_info/alevin_meta_info.json | \
    sed 's/,//g' | awk -F': ' '{print \$2}' | sort -n | head -n 1 | cut -c 1-4)
    
    echo -n "\$mapping_rate"

    
    
  
    mv ${runId}_splici_ALEVIN_tmp ${runId}_splici_ALEVIN
    """
}

// run alevin for cDNA
process alevin_MR1 {
    cache 'lenient'
    memory { 100.GB * task.attempt }
    cpus 4
    errorStrategy { task.exitStatus !=2 && (task.exitStatus == 130 || task.exitStatus == 137 || task.attempt < 3)  ? 'retry' : 'ignore' }
    maxRetries 10


    conda "${baseDir}/envs/alevin.yml"
    
    // cache 'deep'

    // memory { 2.GB * task.attempt }
    cpus 4

    // errorStrategy { task.exitStatus !=2 && (task.exitStatus == 130 || task.exitStatus == 137 || task.attempt < 3)  ? 'retry' : 'ignore' }
    // maxRetries 10

    input:
        set val(runId), file("cdna*.fastq.gz"), file("barcodes*.fastq.gz"), val(barcodeLength), val(umiLength), val(end), val(cellCount), val(barcodeConfig) from FINAL_FASTQS_FOR_ALEVIN_CDNA.join(ALEVIN_CONFIG_CDNA)
        path "alevin_index_cDNA" from ALEVIN_INDEX_CDNA
        path "t2g_cDNA.txt" from T2G_CDNA_FOR_ALEVIN

    output:
        // publishDir path "${runId}_ALEVIN"
        set stdout, val(runId), file("${runId}_cdna_ALEVIN") into ALEVIN_RESULTS_CDNA
        set val(runId), stdout into ALEVIN_CDNA_MAPPING
        path ".command.log"  into MEM_ALEVIN_MR1
            


    """
    salmon alevin ${barcodeConfig} -1 \$(ls barcodes*.fastq.gz | tr '\\n' ' ') -2 \$(ls cdna*.fastq.gz | tr '\\n' ' ') \
        -i alevin_index_cDNA -p ${task.cpus} -o ${runId}_cdna_ALEVIN_tmp --tgMap t2g_cDNA.txt --dumpFeatures --keepCBFraction 1 \
        --freqThreshold ${params.minCbFreq} --dumpMtx > /dev/null
    mapping_rate=\$(grep "mapping_rate" ${runId}_cdna_ALEVIN_tmp/aux_info/alevin_meta_info.json |\
     sed 's/,//g' | awk -F': ' '{print \$2}' | sort -n | head -n 1 | cut -c 1-4)
    echo -n "\$mapping_rate" 
    
    
    mv ${runId}_cdna_ALEVIN_tmp ${runId}_cdna_ALEVIN
    """
}


// build index to runSTARSolo
process index_star {
    cache 'lenient'
    memory { 100.GB * task.attempt }
    cpus 4
    errorStrategy { task.exitStatus !=2 && (task.exitStatus == 130 || task.exitStatus == 137 || task.attempt < 3)  ? 'retry' : 'ignore' }
    maxRetries 10


    conda "${baseDir}/envs/star.yml"

    input:
        path(referenceGenome) from REFERENCE_GENOME
        path(referenceGtf) from REFERENCE_GTF
    output:
        path("STAR_index") into STAR_INDEX
    
    """
    STAR --runMode genomeGenerate --genomeDir STAR_index --genomeFastaFiles ${referenceGenome} --sjdbGTFfile ${referenceGtf} --genomeSAindexNbases 12 

    """
}

// run STARSolo 

process run_STARSolo {
   
    cache 'lenient'

    memory { 100.GB * task.attempt }
    cpus 10
    errorStrategy { task.exitStatus !=2 && (task.exitStatus == 130 || task.exitStatus == 137 || task.attempt < 3)  ? 'retry' : 'ignore' }
    maxRetries 10


    conda "${baseDir}/envs/star.yml"


    input:
    set val(runId), file("cdna*.fastq.gz"), file("barcodes*.fastq.gz"), val(barcodeLength), val(umiLength), val(end), val(cellCount), val(barcodeConfig) from FINAL_FASTQS_FOR_STAR.join(STAR_CONFIG)
    path("STAR_index") from STAR_INDEX

    output:
    set val(runId), path("${runId}_STAR_tmpSolo.out") into STAR_RESULTS
    path ".command.log"   into  MEM_STAR
    

    script:
    if( barcodeConfig == '10XV3' )
        """
        STAR --genomeDir STAR_index --readFilesIn \$(ls cdna*.fastq.gz | tr '\\n' ', '| sed  's/,*\$//g') \$(ls barcodes*.fastq.gz | tr '\\n' ', '| sed  's/,*\$//g') --soloType Droplet --soloCBwhitelist '${baseDir}/whitelist/3M-february-2018.txt.gz' --soloUMIlen ${umiLength} --soloCBlen ${barcodeLength} --soloUMIstart \$(($barcodeLength+1)) --soloCBstart 1 --runThreadN 12 --soloFeatures Gene GeneFull --outFileNamePrefix ${runId}_STAR_tmp --readFilesCommand zcat --soloBarcodeReadLength 0

        mapping_rate=\$(grep "Uniquely mapped reads %" ${runId}_STAR_tmpLog.final.out |\
         awk '{split(\$0, array, "|"); print array[2]}')
        echo  "\${mapping_rate}"

        
        
        """
    else if( barcodeConfig == '10XV2' )
  
        """
        STAR --genomeDir STAR_index --readFilesIn \$(ls cdna*.fastq.gz | tr '\\n' ', '| sed  's/,*\$//g') \$(ls barcodes*.fastq.gz | tr '\\n' ', '| sed  's/,*\$//g') --soloType Droplet --soloCBwhitelist '${baseDir}/whitelist/737K-august-2016.txt' --soloUMIlen ${umiLength} --soloCBlen ${barcodeLength} --soloUMIstart \$(($barcodeLength+1)) --soloCBstart 1 --runThreadN 12 --soloFeatures Gene GeneFull --outFileNamePrefix ${runId}_STAR_tmp --readFilesCommand zcat --soloBarcodeReadLength 0

        mapping_rate=\$(grep "Uniquely mapped reads %" ${runId}_STAR_tmpLog.final.out |\
         awk '{split(\$0, array, "|"); print array[2]}')
        echo  "\${mapping_rate}"

        
        
        """
    else
        """
        STAR --genomeDir STAR_index --readFilesIn \$(ls cdna*.fastq.gz | tr '\\n' ', '| sed  's/,*\$//g') \$(ls barcodes*.fastq.gz | tr '\\n' ', '| sed  's/,*\$//g') --soloType Droplet --readFilesCommand zcat --soloCBwhitelist None --soloUMIlen ${umiLength} --soloCBlen ${barcodeLength} --soloUMIstart \$(($barcodeLength+1)) --soloCBstart 1 --runThreadN 12 --soloFeatures Gene GeneFull --outFileNamePrefix ${runId}_STAR_tmp --soloBarcodeReadLength 0

        mapping_rate=\$(grep "Uniquely mapped reads %" ${runId}_STAR_tmpLog.final.out | \
        awk '{split(\$0, array, "|"); print array[2]}')
        echo  "\${mapping_rate}"


        
        
        """
}

// extract mapping rates from star and turn them into percentages
process get_STAR_mapping {
    cache 'lenient'
    input:
    set val(runId), path("${runId}_STAR_tmpSolo.out") from STAR_RESULTS
    
    output:
    set val(runId), env(GENE_PERCENT) into STAR_MAPPING_GENE
    set val(runId), env(GENEFULL_PERCENT) into STAR_MAPPING_GENEFULL

    """
    GENE=\$(grep "Reads Mapped to Gene: Unique Gene" ${runId}_STAR_tmpSolo.out/Gene/Summary.csv | \
    awk '{split(\$0, array, ","); print array[2]}' | cut -c 1-4)

    GENE_PERCENT=\$(echo "scale=2;((\$GENE * 100))"|bc)

    GENEFULL=\$(grep "Reads Mapped to GeneFull: Unique GeneFull" ${runId}_STAR_tmpSolo.out/GeneFull/Summary.csv |\
     awk '{split(\$0, array, ","); print array[2]}' | cut -c 1-4)

    GENEFULL_PERCENT=\$(echo "scale=2;((\$GENEFULL * 100))"|bc)
    """
}

// build cDNA index for kb-tools 
process index_kb_MR1 {
    cache 'lenient'
    memory { 100.GB * task.attempt }
    errorStrategy { task.exitStatus !=2 && (task.exitStatus == 130 || task.exitStatus == 137 || task.attempt < 3)  ? 'retry' : 'ignore' }
    maxRetries 10
    cpus 4

    conda "${baseDir}/envs/kb-tools.yml"
    
    input:
        path(referenceGenome) from REFERENCE_GENOME
        path(referenceGtf) from REFERENCE_GTF
       

    output:
        set file("kb_index_cDNA"), file("t2g_kb.txt") into KB_INDEX_CDNA
    
       
    """
    kb ref -i kb_index_cDNA -g t2g_kb.txt -f1 cDNA.fa ${referenceGenome} ${referenceGtf} 
    """
}  
// run kb tools count for cDNA reference
process kb_count_MR1 {
    cache 'lenient'
    memory { 100.GB * task.attempt }
    errorStrategy { task.exitStatus !=2 && (task.exitStatus == 130 || task.exitStatus == 137 || task.attempt < 3)  ? 'retry' : 'ignore' }
    maxRetries 10
    conda "${baseDir}/envs/kb-tools.yml"


    input:
        set file("kb_index_cDNA"), file("t2g_kb") from KB_INDEX_CDNA
        set val(runId), file("cdna*.fastq.gz"), file("barcodes*.fastq.gz"), val(barcodeLength), val(umiLength), val(end), val(cellCount), val(barcodeConfig) from FINAL_FASTQS_FOR_KB_TOOLS.join(KB_CONFIG)
        val protocol
    output:
        set val(runId), stdout into KB_CDNA_MAPPING
        path ".command.log"   into MEM_KB_MR1


    """
    kb count -i ${kb_index_cDNA} -t 2 -g ${t2g_kb} -x $protocol \
    -c1 cDNA_cDNA.fa \$(ls barcodes*.fastq.gz | tr '\\n' ' ') \$(ls cdna*.fastq.gz | tr '\\n' ' ') \
     -o "${runId}_out_kb_cDNA"

    mapping_rate=\$(grep "p_pseudoaligned" ${runId}_out_kb_cDNA/run_info.json |sed 's/,//g' | \
    awk '{split(\$0, array, ":"); print array[2]}' | sed 's/^ *//g' | cut -c 1-4) 
    echo -n "\$mapping_rate"
    
    

    """

}

// build index from splici for kb-tools
process index_kb_MR3 {
    cache 'lenient'
    memory { 100.GB * task.attempt }
    cpus 4
    errorStrategy { task.exitStatus !=2 && (task.exitStatus == 130 || task.exitStatus == 137 || task.attempt < 3)  ? 'retry' : 'ignore' }
    maxRetries 10

    conda "${baseDir}/envs/kb-tools.yml"
    
    input:
        path(referenceGenome) from REFERENCE_GENOME
        path(referenceGtf) from REFERENCE_GTF
       
    output:
        set file("kb_index_splici"), file("t2g_kb_splici.txt"), file("cDNA_kb.txt"), file("intron_kb.txt") into KB_INDEX_SPLICI
    
    """
    kb ref -i kb_index_splici -g t2g_kb_splici.txt -f1 cDNA.fa \
        -f2 intron.fa -c1 cDNA_kb.txt -c2  intron_kb.txt \
         ${referenceGenome} ${referenceGtf}  --workflow nucleus
        
    """
}  
// run kb tools count for splici reference
process kb_count_MR3 {
    cache 'lenient'
    memory { 100.GB * task.attempt }
    errorStrategy { task.exitStatus !=2 && (task.exitStatus == 130 || task.exitStatus == 137 || task.attempt < 3)  ? 'retry' : 'ignore' }
    maxRetries 10
    conda "${baseDir}/envs/kb-tools.yml"

    input:
        set file("kb_index_splici"), file("t2g_kb_splici"),file("cDNA_kb.txt"), file("intron_kb.txt") from KB_INDEX_SPLICI
        set val(runId), file("cdna*.fastq.gz"), file("barcodes*.fastq.gz"), val(barcodeLength), val(umiLength), val(end), val(cellCount), val(barcodeConfig) from FINAL_FASTQS_FOR_KB_TOOLS_SPLICI.join(KB_CONFIG_SPLICI)
        val protocol
    output:
        set val(runId), stdout into KB_SPLICI_MAPPING
        path ".command.log"  into MEM_KB_MR3

    """
    kb count -i ${kb_index_splici} -t 2 -g ${t2g_kb_splici} -x $protocol \
    -c1 cDNA.fa \$(ls barcodes*.fastq.gz | tr '\\n' ' ') \$(ls cdna*.fastq.gz | tr '\\n' ' ')  -o "${runId}_out_kb_splici" \
    --workflow nucleus -c1 cDNA_kb.txt -c2 intron_kb.txt

    mapping_rate=\$(grep "p_pseudoaligned" ${runId}_out_kb_splici/run_info.json |sed 's/,//g' | \
    awk '{split(\$0, array, ":"); print array[2]}'| sed 's/^ *//g' | cut -c 1-4)

    echo -n  "\$mapping_rate"
    
    

    """

}
// index preRNA (from adapted gtf file) for kb_tools
process index_kb_MR2 {
    cache 'lenient'
    memory { 100.GB * task.attempt }
    errorStrategy { task.exitStatus !=2 && (task.exitStatus == 130 || task.exitStatus == 137 || task.attempt < 3)  ? 'retry' : 'ignore' }
    maxRetries 10
    cpus 4
    conda "${baseDir}/envs/gff_read.yml"

    input:
    path(referenceGenome) from REFERENCE_GENOME
    path(referenceGtf) from REFERENCE_GTF

    output:
    set file("kb_index_preRNA"), file("t2g_kb_preRNA.txt"), file("cDNA_preRNA.fa") into KB_INDEX_PRERNA
    
    """
    awk 'BEGIN{FS="\t"; OFS="\t"} \$3 == "transcript"{ \$3="exon"; print}' ${referenceGtf}  > preRNA_referenceGtf.gtf

    kb ref -i kb_index_preRNA -g t2g_kb_preRNA.txt -f1 cDNA_preRNA.fa ${referenceGenome} preRNA_referenceGtf.gtf
    """     
}
// run kb count for preRNA 
process kb_count_MR2 {
    cache 'lenient'
    memory { 100.GB * task.attempt }
    errorStrategy { task.exitStatus !=2 && (task.exitStatus == 130 || task.exitStatus == 137 || task.attempt < 3)  ? 'retry' : 'ignore' }
    maxRetries 10
    cpus 4
    conda "${baseDir}/envs/kb-tools.yml"

    input:
        set file("kb_index_preRNA"), file("t2g_kb_preRNA"), file("cDNA_preRNA.fa") from KB_INDEX_PRERNA
        set val(runId), file("cdna*.fastq.gz"), file("barcodes*.fastq.gz"), val(barcodeLength), val(umiLength), val(end), val(cellCount), val(barcodeConfig) from FINAL_FASTQS_FOR_KB_TOOLS_PRERNA.join(KB_CONFIG_PRERNA)
        val protocol
    output:
        set val(runId), stdout into KB_PRERNA_MAPPING
        path ".command.log"  into  MEM_KB_MR2
        


    """
    kb count -i ${kb_index_preRNA} -t 2 -g ${t2g_kb_preRNA} -x $protocol \
    -c1 cDNA_preRNA.fa \$(ls barcodes*.fastq.gz | tr '\\n' ' ')  \$(ls cdna*.fastq.gz | tr '\\n' ' ') \
    -o "${runId}_out_kb_preRNA"

    mapping_rate=\$(grep "p_pseudoaligned" ${runId}_out_kb_preRNA/run_info.json |sed 's/,//g'| awk '{split(\$0, array, ":"); print array[2]}' | sed 's/^ *//g' | cut -c 1-4) 
    echo -n "\$mapping_rate"
    
    

    """

}
// build salmon index for alevin-fry with splici
process index_alevin_fry_MR3 {
    cache 'lenient'
    memory { 100.GB * task.attempt }
    cpus 4
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
// run alevin-fry for quantification with splici index
 process alevin_fry_MR3 {
    cache 'lenient'
    memory { 100.GB * task.attempt }
    errorStrategy { task.exitStatus !=2 && (task.exitStatus == 130 || task.exitStatus == 137 || task.attempt < 3)  ? 'retry' : 'ignore' }
    maxRetries 10
    conda "${baseDir}/envs/alevin-fry_2.yml"
    input:
        set val(runId), file("cdna*.fastq.gz"), file("barcodes*.fastq.gz"), val(barcodeLength), val(umiLength), val(end), val(cellCount), val(barcodeConfig) from FINAL_FASTQS_FOR_ALEVIN_FRY.join(ALEVIN_FRY_CONFIG)
        path "alevin_index_for_fry" from ALEVIN_FRY_INDEX_SPLICI
        path "t2g_cDNA.txt" from T2G_3_FOR_FRY
       
    output:
        // publishDir path "${runId}_ALEVIN"
        set val(runId), file("${runId}_ALEVIN_fry_quant") into ALEVIN_FRY_RESULTS_SPLICI
        set val(runId), env(FRY_MAPPING) into ALEVIN_FRY_MAPPING_SPLICI
        path ".command.log"  into MEM_ALEVIN_FRY_MR3

    """
    salmon alevin ${barcodeConfig} --sketch -1 \$(ls barcodes*.fastq.gz | tr '\\n' ' ') -2 \$(ls cdna*.fastq.gz | tr '\\n' ' ') \
        -i alevin_index_for_fry -p ${task.cpus} -o ${runId}_ALEVIN_fry_map t2g_cDNA.txt 

    if (${barcodeConfig} == "--chromium")
    then
        alevin-fry generate-permit-list --unfiltered-pl '${baseDir}/whitelist/737K-august-2016.txt' --input ${runId}_ALEVIN_fry_map -d fw --output-dir ${runId}_ALEVIN_fry_quant -k --min-reads 1
    elif (${barcodeConfig} == "--chromiumV3")
    then
        alevin-fry generate-permit-list --unfiltered-pl '${baseDir}/whitelist/3M-february-2018.txt.gz' t --input ${runId}_ALEVIN_fry_map -d fw --output-dir ${runId}_ALEVIN_fry_quant -k --min-reads 1
    else
        alevin-fry generate-permit-list --input ${runId}_ALEVIN_fry_map -d fw --output-dir ${runId}_ALEVIN_fry_quant -k --min-reads 1
    fi

    alevin-fry collate -i ${runId}_ALEVIN_fry_quant -r ${runId}_ALEVIN_fry_map -t 16
    alevin-fry quant -i ${runId}_ALEVIN_fry_quant -m t2g_cDNA.txt -t 16 -r cr-like -o ${runId}_ALEVIN_fry_quant --use-mtx

    TOTAL=\$(grep "num_processed" ${runId}_ALEVIN_fry_map/aux_info/meta_info.json |  awk '{split(\$0, array, ": "); print array[2]}'| sed 's/,//g')

    MAPPED=\$(grep "num_mapped" ${runId}_ALEVIN_fry_map/aux_info/meta_info.json |  awk '{split(\$0, array, ": "); print array[2]}'| sed 's/,//g')

    FRY_MAPPING=\$(echo "scale=2;((\$MAPPED * 100) / \$TOTAL)"|bc)

    
    
    
    """
}

process index_alevin_fry_MR2 {
    cache 'lenient'
    memory { 100.GB * task.attempt }
    errorStrategy { task.exitStatus !=2 && (task.exitStatus == 130 || task.exitStatus == 137 || task.attempt < 3)  ? 'retry' : 'ignore' }
    maxRetries 10
    cpus 4

    conda "${baseDir}/envs/alevin-fry_2.yml"

    input:
        path reference from REFERENCE_GENOME
        path referenceGtf from REFERENCE_GTF
        
    output:
        path "alevin_index_for_fry_transcriptome" into ALEVIN_INDEX_FOR_FRY_TRANSCRIPTOME
        path "t2g_transcriptome.txt" into T2G_FOR_FRY_TRANSCRIPTOME

    """
    awk 'BEGIN{FS="\\t"; OFS="\\t"} \$3 == "transcript"{ \$3="exon"; print}' ${referenceGtf} > preRNA_referenceGtf.gtf

    gffread -w transcriptome -g ${reference} preRNA_referenceGtf.gtf

    sed -i 's/transcript://g' transcriptome
    
    cat ${referenceGtf} | grep -vE "^#" | awk '\$3=="transcript" {split(\$0, array, "transcript_id"); print array[2]}' | awk '{split(\$0, array, ";"); print array[1]}' | sed 's/"//g' | sed 's/^ *//g' | sed 's/transcript://g' > t.txt

    cat ${referenceGtf} | grep -vE "^#" | awk '\$3=="transcript" {split(\$0, array, "gene_id"); print array[2]}' | awk '{split(\$0, array, ";"); print array[1]}' | sed 's/"//g' | sed 's/^ *//g'| sed 's/gene://g'  > g.txt

    paste t.txt g.txt > t2g_transcriptome.txt
    salmon index --transcript transcriptome  -i alevin_index_for_fry_transcriptome -k 19
    """

 }
// cat preRNA_referenceGtf.gtf | awk  '{print \$10"\\t"\$12}' | awk  '{print \$2"\\t"\$1}' > t2g_transcriptome.txt
  
 process alevin_fry_MR2 {
    cache 'lenient'
    memory { 20.GB * task.attempt }
    errorStrategy { task.exitStatus !=2 && (task.exitStatus == 130 || task.exitStatus == 137 || task.attempt < 3)  ? 'retry' : 'ignore' }
    maxRetries 10
    conda "${baseDir}/envs/alevin-fry_2.yml"
    input:
        set val(runId), file("cdna*.fastq.gz"), file("barcodes*.fastq.gz"), val(barcodeLength), val(umiLength), val(end), val(cellCount), val(barcodeConfig) from FINAL_FASTQS_FOR_ALEVIN_FRY_TRANSCRIPTOME.join(ALEVIN_FRY_CONFIG_TRANSCRIPTOME)
        path "alevin_index_for_fry" from ALEVIN_INDEX_FOR_FRY_TRANSCRIPTOME
        path "t2g_transcriptome.txt" from T2G_FOR_FRY_TRANSCRIPTOME
       
    output:
        // publishDir path "${runId}_ALEVIN"
        set val(runId), file("${runId}_ALEVIN_fry_quant") into ALEVIN_FRY_RESULTS_TRANSCRIPTOME
        set val(runId), env(FRY_MAPPING) into ALEVIN_FRY_MAPPING_TRANSCRIPTOME
        path ".command.log"  into MEM_ALEVIN_FRY_MR2

    """
    salmon alevin ${barcodeConfig} --sketch -1 \$(ls barcodes*.fastq.gz | tr '\\n' ' ') -2 \$(ls cdna*.fastq.gz | tr '\\n' ' ') \
        -i alevin_index_for_fry -p ${task.cpus} -o ${runId}_ALEVIN_fry_map 

    if (${barcodeConfig} == "--chromium")
    then
        alevin-fry generate-permit-list --unfiltered-pl '${baseDir}/whitelist/737K-august-2016.txt' --input ${runId}_ALEVIN_fry_map -d fw --output-dir ${runId}_ALEVIN_fry_quant -k --min-reads 1
    elif (${barcodeConfig} == "--chromiumV3")
    then
        alevin-fry generate-permit-list --unfiltered-pl '${baseDir}/whitelist/3M-february-2018.txt.gz' t --input ${runId}_ALEVIN_fry_map -d fw --output-dir ${runId}_ALEVIN_fry_quant -k --min-reads 1
    else
        alevin-fry generate-permit-list --input ${runId}_ALEVIN_fry_map -d fw --output-dir ${runId}_ALEVIN_fry_quant -k --min-reads 1
    fi

    alevin-fry collate -i ${runId}_ALEVIN_fry_quant -r ${runId}_ALEVIN_fry_map -t 16
    alevin-fry quant -i ${runId}_ALEVIN_fry_quant -m t2g_transcriptome.txt -t 16 -r cr-like -o ${runId}_ALEVIN_fry_quant --use-mtx

    TOTAL=\$(grep "num_processed" ${runId}_ALEVIN_fry_map/aux_info/meta_info.json |  awk '{split(\$0, array, ": "); print array[2]}'| sed 's/,//g')

    MAPPED=\$(grep "num_mapped" ${runId}_ALEVIN_fry_map/aux_info/meta_info.json |  awk '{split(\$0, array, ": "); print array[2]}'| sed 's/,//g')

    FRY_MAPPING=\$(echo "scale=2;((\$MAPPED * 100) / \$TOTAL)"|bc)

    
    
    """
}

// build salmon index for alevin-fry
 process index_alevin_fry_MR1 {
    cache 'lenient'
    memory { 100.GB * task.attempt }
    errorStrategy { task.exitStatus !=2 && (task.exitStatus == 130 || task.exitStatus == 137 || task.attempt < 3)  ? 'retry' : 'ignore' }
    maxRetries 10
    cpus 4

    conda "${baseDir}/envs/alevin-fry_2.yml"

    input:
        path reference from REFERENCE_CDNA
        
        
    output:
        path "alevin_index_for_fry" into ALEVIN_INDEX_FOR_FRY_CDNA

    """
    salmon index --transcript ${reference}  -i alevin_index_for_fry -k 19
    """
 }

// run alevin-fry for cdna
process alevin_fry_MR1 {
    cache 'lenient'
    memory { 20.GB * task.attempt }
    errorStrategy { task.exitStatus !=2 && (task.exitStatus == 130 || task.exitStatus == 137 || task.attempt < 3)  ? 'retry' : 'ignore' }
    maxRetries 10
    conda "${baseDir}/envs/alevin-fry_2.yml"
    input:
        set val(runId), file("cdna*.fastq.gz"), file("barcodes*.fastq.gz"), val(barcodeLength), val(umiLength), val(end), val(cellCount), val(barcodeConfig) from FINAL_FASTQS_FOR_ALEVIN_FRY_CDNA.join(ALEVIN_FRY_CONFIG_CDNA)
        path "alevin_index_for_fry" from ALEVIN_INDEX_FOR_FRY_CDNA
        path "t2g_cDNA.txt" from T2G_CDNA_FOR_ALEVIN_FRY
       
    output:
        // publishDir path "${runId}_ALEVIN"
        set val(runId), file("${runId}_ALEVIN_fry_quant") into ALEVIN_FRY_RESULTS_CDNA
        set val(runId), env(FRY_MAPPING) into ALEVIN_FRY_MAPPING_CDNA
        path ".command.log" into MEM_ALEVIN_FRY_MR1

    """
    salmon alevin ${barcodeConfig} --sketch -1 \$(ls barcodes*.fastq.gz | tr '\\n' ' ') -2 \$(ls cdna*.fastq.gz | tr '\\n' ' ') \
        -i alevin_index_for_fry -p ${task.cpus} -o ${runId}_ALEVIN_fry_map t2g_cDNA.txt 

    if (${barcodeConfig} == "--chromium")
    then
        alevin-fry generate-permit-list --unfiltered-pl '${baseDir}/whitelist/737K-august-2016.txt' --input ${runId}_ALEVIN_fry_map -d fw --output-dir ${runId}_ALEVIN_fry_quant -k --min-reads 1
    elif (${barcodeConfig} == "--chromiumV3")
    then
        alevin-fry generate-permit-list --unfiltered-pl '${baseDir}/whitelist/3M-february-2018.txt.gz' t --input ${runId}_ALEVIN_fry_map -d fw --output-dir ${runId}_ALEVIN_fry_quant -k --min-reads 1
    else
        alevin-fry generate-permit-list --input ${runId}_ALEVIN_fry_map -d fw --output-dir ${runId}_ALEVIN_fry_quant -k --min-reads 1
    fi

    alevin-fry collate -i ${runId}_ALEVIN_fry_quant -r ${runId}_ALEVIN_fry_map -t 16
    alevin-fry quant -i ${runId}_ALEVIN_fry_quant -m t2g_cDNA.txt -t 16 -r cr-like -o ${runId}_ALEVIN_fry_quant --use-mtx

    TOTAL=\$(grep "num_processed" ${runId}_ALEVIN_fry_map/aux_info/meta_info.json |  awk '{split(\$0, array, ": "); print array[2]}'| sed 's/,//g')

    MAPPED=\$(grep "num_mapped" ${runId}_ALEVIN_fry_map/aux_info/meta_info.json |  awk '{split(\$0, array, ": "); print array[2]}'| sed 's/,//g')

    FRY_MAPPING=\$(echo "scale=2;((\$MAPPED * 100) / \$TOTAL)"|bc)

    
    

    
    """
}
// make a table for the mapping rates for different tools
process write_table {
    publishDir "$resultsRoot", mode: 'copy', overwrite: true
   
    input:
    set val(runId), mr1, mr2, mr3, mr4, mr5, mr8, mr9, mr10 from ALEVIN_CDNA_MAPPING.join(ALEVIN_SPLICI_MAPPING).join(KB_CDNA_MAPPING).join(KB_PRERNA_MAPPING).join(KB_SPLICI_MAPPING).join(ALEVIN_FRY_MAPPING_CDNA).join(ALEVIN_FRY_MAPPING_TRANSCRIPTOME).join(ALEVIN_FRY_MAPPING_SPLICI)
    set val(key), mr6, mr7 from STAR_MAPPING_GENE.join(STAR_MAPPING_GENEFULL)
    
    output:
    file("*_${runId}.txt") into RESULTS_FOR_COUNTING
    
    """
    echo "tool\tMPR1\tMPR2\tMPR3\nAlevin (%)\t${mr1}\t${mr2}\tNA\nAlevin-fry (%)\t${mr8}\t${mr9}\t${mr10}\nkb-tools (%)\t${mr3}\t${mr4}\t${mr5}\nSTARSolo (%)\t${mr6}\tNA\t${mr7}\n" > ${params.name}_${runId}.txt
         
    """
}
MEM = [MEM_ALEVIN_MR1, MEM_ALEVIN_MR2, MEM_ALEVIN_FRY_MR1, MEM_ALEVIN_FRY_MR2, MEM_ALEVIN_FRY_MR3, MEM_KB_MR1, MEM_KB_MR2, MEM_KB_MR3, MEM_STAR]
// MEM = MEM_ALEVIN_MR1.join(MEM_ALEVIN_MR2).join(MEM_ALEVIN_FRY_MR1).join(MEM_ALEVIN_FRY_MR2).join(MEM_ALEVIN_FRY_MR3).join(MEM_KB_MR1).join(MEM_KB_MR2).join(MEM_KB_MR3).join(MEM_STAR)
// MEM.view()
process parse_command_log {

    input: 
    file("log_file_*") from MEM
    output:
    env AVG_MEM into AVG_MEMORIES
    env RUN_TIME into RUN_TIMES
    
    """

    AVG_MEM=\$(grep "Average Memory : " log_file_* | awk '{split(\$0, array, ":"); print array[2]}' | sed 's/^ *//g' |sed 's/ MB//g' )
    RUN_TIME=\$(grep "Run time : " log_file_* | awk '{split(\$0, array, ":"); print array[2]}' | sed 's/^ *//g' |sed 's/ sec.//g' )

    """

}

// MEM=MEM_ALEVIN_MR1.join(MEM_ALEVIN_MR2).join(MEM_ALEVIN_FRY_MR1).join(MEM_ALEVIN_FRY_MR2).join(MEM_ALEVIN_FRY_MR3).join(MEM_KB_MR1).join(MEM_KB_MR2).join(MEM_KB_MR3).join(MEM_STAR)
// TIME=TIME_ALEVIN_MR1.join(TIME_ALEVIN_MR2).join(TIME_ALEVIN_FRY_MR1).join(TIME_ALEVIN_FRY_MR2).join(TIME_ALEVIN_FRY_MR3).join(TIME_KB_MR1).join(TIME_KB_MR2).join(TIME_KB_MR3).join(TIME_STAR)

process write_table_benchmark {
    publishDir "$resultsRoot/memory", mode: 'copy', overwrite: true
   
    input:
    set  mr1, mr2, mr3, mr4, mr5, mr6, mr7, mr8, mr9, mr10 from AVG_MEMORIES
    output:
    file("*_memory.txt") into RESULTS_MEMORY
 
    """
    echo "memory\tMPR1\tMPR2\tMPR3\nAlevin\t${mr1}\t${mr2}\tNA\nAlevin-fry\t${mr3}\t${mr4}\t${mr5}\nkb-tools\t${mr6}\t${mr7}\t${mr8}\nSTARSolo\t${mr9}\tNA\t${mr9}\n" > ${params.name}_memory.txt    
    """
}