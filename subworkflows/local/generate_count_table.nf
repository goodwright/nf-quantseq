include { CUTADAPT } from '../../modules/nf-core/cutadapt/main.nf'
include { STAR_ALIGN } from '../../modules/nf-core/star/align/main.nf'
include { SAMTOOLS_INDEX as SAMTOOLS_INDEX} from '../../modules/nf-core/samtools/index/main.nf'
include { SAMTOOLS_INDEX as SAMTOOLS_INDEX_DEDUPE} from '../../modules/nf-core/samtools/index/main.nf'
include { UMICOLLAPSE } from '../../modules/nf-core/umicollapse/main.nf'
include { STAR_INPUTALIGNMENTSFROMBAM } from '../../modules/local/star_inputalignmentsfrombam.nf'
include {
    BEDTOOLS_SORT as BEDTOOLS_SORT_POS;
    BEDTOOLS_SORT as BEDTOOLS_SORT_NEG
} from '../../modules/nf-core/bedtools/sort/main.nf'
include {
    UCSC_BEDGRAPHTOBIGWIG as UCSC_BEDGRAPHTOBIGWIG_POS;
    UCSC_BEDGRAPHTOBIGWIG as UCSC_BEDGRAPHTOBIGWIG_NEG
} from '../../modules/nf-core/ucsc/bedgraphtobigwig/main.nf'
include { BEDTOOLS_BAMTOBED } from '../../modules/nf-core/bedtools/bamtobed/main.nf'
include { BEDTOOLS_WINDOW } from '../../modules/local/bedtools_window.nf'
include { BEDTOOLS_WINDOW as BEDTOOLS_WINDOW_REVERSED } from '../../modules/local/bedtools_window.nf'
include { AWK as BEDGRAPHCONVERT_AWK } from '../../modules/local/awk.nf'

workflow GENERATE_COUNT_TABLE {
    take:
    reads
    star_index
    gtf
    fai
    polya_bed
    quantseq_rev
    bc_pattern

    main:

    reads
        .map{ tuple ->
            def new_meta = tuple[0].clone()
            new_meta["id"] += ".polya_trimmed"
            [new_meta, tuple[1]]
        }
        .set{ ch_cutadapt_input }

    CUTADAPT(
        ch_cutadapt_input
    )

    STAR_ALIGN(
        CUTADAPT.out.reads,
        star_index,
        gtf,
        true,  // star_ignore_sjdbgtf - Required for the GTF to be used to detect splice junctions
        false, // seq_platform
        false  // seq_center
    )
    ch_bam = STAR_ALIGN.out.bam_sorted

    SAMTOOLS_INDEX(
        STAR_ALIGN.out.bam_sorted
    )
    ch_bai = SAMTOOLS_INDEX.out.bai

    if(bc_pattern) {
        ch_bam_bai = ch_bam
            .map { row -> [row[0].id, row ].flatten()}
            .join ( ch_bai.map { row -> [row[0].id, row ].flatten()} )
            .map { row -> [row[1], row[2], row[4]] }

        UMICOLLAPSE(
            ch_bam_bai,
            "bam"
        )
        ch_bam = UMICOLLAPSE.out.bam
        SAMTOOLS_INDEX_DEDUPE(
            UMICOLLAPSE.out.bam
        )
        ch_bai = SAMTOOLS_INDEX_DEDUPE.out.bai
    }

    STAR_INPUTALIGNMENTSFROMBAM(
        ch_bam
    )

    STAR_INPUTALIGNMENTSFROMBAM.out.unique_coverage_pos
        .map{ tuple ->
            def new_meta = tuple[0].clone()
            new_meta["id"] += "_pos"
            [new_meta, tuple[1]] 
        }
        .set{ bedgraph_pos }

    STAR_INPUTALIGNMENTSFROMBAM.out.unique_coverage_neg
        .map{ tuple ->
            def new_meta = tuple[0].clone()
            new_meta["id"] += "_neg"
            [new_meta, tuple[1]] 
        }
        .set{ bedgraph_neg }

    BEDTOOLS_SORT_POS(
        bedgraph_pos,
        "sorted.bedgraph"
    )

    BEDTOOLS_SORT_NEG(
        bedgraph_neg,
        "sorted.bedgraph"
    )

    UCSC_BEDGRAPHTOBIGWIG_POS(
        BEDTOOLS_SORT_POS.out.sorted,
        fai.collect()
    )

    UCSC_BEDGRAPHTOBIGWIG_NEG(
        BEDTOOLS_SORT_NEG.out.sorted,
        fai.collect()
    )

    // Create count table

    BEDTOOLS_BAMTOBED(
        ch_bam
    )

    if(quantseq_rev) {
        BEDTOOLS_WINDOW_REVERSED(
            polya_bed.collect(),
            BEDTOOLS_BAMTOBED.out.bed
        )
        ch_window_counts = BEDTOOLS_WINDOW_REVERSED.out.overlap
    } else {
        BEDTOOLS_WINDOW(
            polya_bed.collect(),
            BEDTOOLS_BAMTOBED.out.bed
        )
        ch_window_counts = BEDTOOLS_WINDOW.out.overlap
    }

    BEDGRAPHCONVERT_AWK(
        ch_window_counts,
        '{{OFS="\t"}}{{if($6 == "+") {{print $1, $2, $3, $10}} else {{print $1, $2, $3, -$10}}}}',
        '| sort -k1,1 -k2,2n'
    )

    emit:
    bam = ch_bam
    bai = ch_bai
}
