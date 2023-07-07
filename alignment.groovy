// import
// bwa2 for extra alignment option
include { BWA_MEM      } from './modules/bwa_mem.groovy'
include { PICARD_ADDORREPLACEREADGROUPS      } from './modules/picard_addorreplacereadgroups.groovy'

workflow ALIGNMENT {

    take:
    fastqs // channel: [ val(meta), [ bam ] ]
    reference
    bwa

    main:

    versions = Channel.empty()

    BWA_MEM ( fastqs, reference, true ).bam.map {
        meta, bam ->
            new_id = 'aligned_bam'
            [[id: new_id], bam ]
    }.set {aligned_bam}
    versions = versions.mix(BWA_MEM.out.versions.first())
    
    PICARD_ADDORREPLACEREADGROUPS(aligned_bam).bam.map {
        meta, bam ->
            new_id = 'grouped_aligned_bam'
            [[id: new_id], bam ]
    }.set {grouped_bam}
    versions = versions.mix(PICARD_ADDORREPLACEREADGROUPS.out.versions.first())

    // final output
    emit:

    bam      = PICARD_ADDORREPLACEREADGROUPS.out.bam           // channel: [ val(meta), [ bam ] ]

    versions = versions                     // channel: [ versions.yml ]
}

workflow test_alignment {

    // channels enable parallel: https://www.nextflow.io/docs/latest/faq.html?highlight=parallel
    // test data 
    fastqs = [
    [[id:'gene', single_end:false], [params.test_data_msk['uncollapsed_bam_generation']['merged_fastq']['merged_1'], params.test_data_msk['uncollapsed_bam_generation']['merged_fastq']['merged_2']]]
    ]
    reference = [
        [id:'reference'], 
        file('test_nucleo/reference/')
    ]
    fastqs = ch_fastq = Channel.fromList(fastqs)

    // workflow 
    ALIGNMENT ( fastqs, reference, 1)
}