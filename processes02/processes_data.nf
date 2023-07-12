/* 
 * pipeline input parameters 
 */
nextflow.enable.dsl=2
// params.reads = "$baseDir/data/gut_{1,2}.fq"
params.reads = "$baseDir/../data/*_{1,2}.fq"
params.transcript = "$baseDir/../data/transcriptome.fa"

log.info """\
         R N A S E Q - N F   P I P E L I N E    
         ===================================
         transcriptome: ${params.transcript}
         reads        : ${params.reads}
         outdir       : ${params.outdir}
         """
         .stripIndent()

 
/* 
 * define the `index` process that create a binary index 
 * given the transcriptome file
 */
process index {
    
    input:
    path transcriptome
     
    output:
    path 'index' , emit: index

    script:       
    """
    salmon index --threads $task.cpus -t $transcriptome -i index
    """
}


/*
 * Run Salmon to perform the quantification of expression using
 * the index and the matched read files
 */
process quantification {
     
    input:
    path index 
    tuple val(pair_id), path(reads) 
 
    output:
    path(pair_id)
 
    script:
    """
    salmon quant --threads $task.cpus --libType=U -i $index -1 ${reads[0]} -2 ${reads[1]} -o $pair_id
    """
}

workflow {
    publishDir = "/Users/ebuehler/Documents/GitHub/nextflow-training/results/"
    // Index file 
    index(params.transcript)

    // add local Data to Channel 
    all_pairs = params.reads
    Channel.fromFilePairs( all_pairs, checkIfExists:true ).set { all_pairs } 

    // add online hosted data to Channel
    lung_pairs = [
        ['lung', 
        [params.test_data_msk['lung_reads']['reads']['lung_1'], params.test_data_msk['lung_reads']['reads']['lung_2']]
        ]
    ]

    Channel.fromList( lung_pairs).set { lung_pairs } 
    
    // combine
    all_pairs = all_pairs.concat(lung_pairs)

    // Quantification 
    quantification(index.out.index, all_pairs)



    emit: 
        quantification = quantification.out
}