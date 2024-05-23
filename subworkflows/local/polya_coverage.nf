include { SAMTOOLS_MERGE } from '../../modules/nf-core/samtools/merge/main.nf'
include {
    BEDTOOLS_GENOMECOV as BEDTOOLS_GENOMECOV_POS;
    BEDTOOLS_GENOMECOV as BEDTOOLS_GENOMECOV_NEG;
    BEDTOOLS_GENOMECOV as BEDTOOLS_GENOMECOV_POS_REVERSED;
    BEDTOOLS_GENOMECOV as BEDTOOLS_GENOMECOV_NEG_REVERSED
} from '../../modules/nf-core/bedtools/genomecov/main.nf'
include { RANDOM_PRIMING } from '../../modules/local/random_priming.nf'
include { ANNOTATE_POLYA_SITES } from '../../modules/local/annotate_polya_sites.nf'

workflow POLYA_COVERAGE {
    take:
    read_list
    gtf
    quantseq_rev //boolean

    main:

    SAMTOOLS_MERGE(
        read_list,
        [] // reference fasta - Not needed here, and no, I don't know why it needs to be a [] and not '' or false
    )

    if (quantseq_rev) {
        BEDTOOLS_GENOMECOV_POS_REVERSED(
            SAMTOOLS_MERGE.out.bam.map{ tuple -> [tuple[0], tuple[1], 1.0] }, // tuple val(meta), path(intervals), val(scale)
            [],                                                               // sizes
            'bedgraph'                                                        // extension
        )
        ch_pos_genomecov = BEDTOOLS_GENOMECOV_POS_REVERSED.out.genomecov
        BEDTOOLS_GENOMECOV_NEG_REVERSED(
            SAMTOOLS_MERGE.out.bam.map{ tuple -> [tuple[0], tuple[1], 1.0] },
            [],
            'bedgraph'
        )
        ch_neg_genomecov = BEDTOOLS_GENOMECOV_NEG_REVERSED.out.genomecov
    } else {
        BEDTOOLS_GENOMECOV_POS(
            SAMTOOLS_MERGE.out.bam.map{ tuple -> [tuple[0], tuple[1], 1.0] }, // tuple val(meta), path(intervals), val(scale)
            [],                                                               // sizes
            'bedgraph'                                                        // extension
        )
        ch_pos_genomecov = BEDTOOLS_GENOMECOV_POS.out.genomecov

        BEDTOOLS_GENOMECOV_NEG(
            SAMTOOLS_MERGE.out.bam.map{ tuple -> [tuple[0], tuple[1], 1.0] },
            [],
            'bedgraph'
        )
        ch_neg_genomecov = BEDTOOLS_GENOMECOV_NEG.out.genomecov
    }

    RANDOM_PRIMING(
        ch_pos_genomecov,
        ch_neg_genomecov
    )

    ANNOTATE_POLYA_SITES(
        RANDOM_PRIMING.out.filteredunique_bed,
        gtf
    )

    emit:
    polya_bed = ANNOTATE_POLYA_SITES.out.bed
}
