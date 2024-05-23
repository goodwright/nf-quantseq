include { SAMTOOLS_FAIDX } from '../modules/nf-core/samtools/faidx/main.nf'
include { SAMTOOLS_INDEX } from '../modules/nf-core/samtools/index/main.nf'
include { FASTQC  } from '../modules/nf-core/fastqc/main.nf'
include { CUTADAPT as CUTADAPT_ADAPTERS } from '../modules/nf-core/cutadapt/main.nf'
include { CUTADAPT as CUTADAPT_ADAPTERS_TSO } from '../modules/nf-core/cutadapt/main.nf'
include { UMITOOLS_EXTRACT } from '../modules/nf-core/umitools/extract/main.nf'
include { GET_POLYA_READS } from '../subworkflows/local/get_polya_reads.nf'
include { STAR_GENOMEGENERATE } from '../modules/nf-core/star/genomegenerate/main.nf'
include { STAR_ALIGN } from '../modules/nf-core/star/align/main.nf'
include { UMICOLLAPSE } from '../modules/nf-core/umicollapse/main.nf'
include { GUNZIP } from '../modules/nf-core/gunzip/main.nf'
include { POLYA_COVERAGE } from '../subworkflows/local/polya_coverage.nf'
include { GENERATE_COUNT_TABLE } from '../subworkflows/local/generate_count_table.nf'

workflow QUANTSEQ {
    // Input stuff, kind of messy

    Channel.fromPath( params.input )
        .map{ path -> [
            [
                'id': path.getSimpleName(),
                'single_end': true
            ],
            path
        ]}
        .set{ ch_fastq }

    params.fasta = WorkflowMain.getGenomeAttribute(params, 'fasta')

    // WORKFLOW

    SAMTOOLS_FAIDX(
        Channel.fromPath(params.fasta).map{ path -> [[], path] }
    )
    
    SAMTOOLS_FAIDX.out.fai.map{ tuple -> tuple[1] }
        .set{ fai }

    FASTQC(
        ch_fastq
    )

    if(params.umitools_bc_pattern) {
        UMITOOLS_EXTRACT(
            ch_fastq
        )
        ch_fastq = UMITOOLS_EXTRACT.out.reads
    }

    if(params.tso) {
        CUTADAPT_ADAPTERS_TSO(
            ch_fastq
        )
        ch_fastq = CUTADAPT_ADAPTERS_TSO.out.reads
    } else {
    CUTADAPT_ADAPTERS(
        ch_fastq
    )
    ch_fastq = CUTADAPT_ADAPTERS.out.reads
    }

    // Decompress the gtf / gff
    GUNZIP(
        Channel.fromPath(params.gtf).map{
            path -> [['id': 'gtf'], path]
        }
    )

    GUNZIP.out.gunzip
        .map{ tuple -> tuple[1] }
        .collect()
        .set{ gtf_gunzip }

    STAR_GENOMEGENERATE(
        params.fasta,
        gtf_gunzip
    )
    
    STAR_GENOMEGENERATE.out.index
        .collect()
        .set{ star_index }

    if(!params.polya_bed){
        GET_POLYA_READS(
            ch_fastq
        )

        STAR_ALIGN(
            GET_POLYA_READS.out.reads,
            star_index,
            gtf_gunzip,
            true,  // star_ignore_sjdbgtf - Required for the GTF to be used to detect splice junctions
            false, // seq_platform
            false  // seq_center
        )

        SAMTOOLS_INDEX ( 
            STAR_ALIGN.out.bam_sorted 
            )
        ch_bam = STAR_ALIGN.out.bam_sorted
        ch_bai = SAMTOOLS_INDEX.out.bai

        if(params.umitools_bc_pattern) {
            ch_bam_bai = ch_bam
                .map { row -> [row[0].id, row ].flatten()}
                .join ( ch_bai.map { row -> [row[0].id, row ].flatten()} )
                .map { row -> [row[1], row[2], row[4]] }

            UMICOLLAPSE(
                ch_bam_bai,
                "bam"
            )
            ch_bam = UMICOLLAPSE.out.bam
        }

        ch_bam
            .map{ tuple -> tuple[1] }
            .collect()
            .map{ paths -> [['id': 'merged_polya'], paths] }
            .set{ collected_reads }

        POLYA_COVERAGE(
            collected_reads,
            gtf_gunzip,
            params.quantseq_rev
        )
    ch_polya_bed = POLYA_COVERAGE.out.polya_bed
    } else {
        ch_polya_bed = Channel.fromPath(params.polya_bed).map{
            path -> [['id': 'polya_bed'], path]
        }
    
    }

    GENERATE_COUNT_TABLE(
        ch_fastq,
        star_index,
        gtf_gunzip,
        fai,
        ch_polya_bed,
        params.quantseq_rev,
        params.umitools_bc_pattern
    )

}
