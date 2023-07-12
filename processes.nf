/* 
 * pipeline input parameters 
 */
nextflow.enable.dsl=2
// params.reads = "$baseDir/data/gut_{1,2}.fq"
params.reads = "$baseDir/data/*_{1,2}.fq"
params.transcript = "$baseDir/data/transcriptome.fa"
params.outdir = "results/"

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

    // add reads to Channel 
    Channel
    .fromFilePairs( params.reads, checkIfExists:true )
    .set { read_pairs_ch } 
    // read_pairs_ch.view()

    // Example 1: Quantification 
    quantification(index.out.index, read_pairs_ch)

    // Example 2: Quantification
    // index_dir = "$baseDir/nextflow-training/data/index"
    // quantification(index_dir, read_pairs_ch)


    emit: 
        quantification = quantification.out
}