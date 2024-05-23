include { CUTADAPT_UNTRIMMED } from '../../modules/local/cutadapt_untrimmed.nf'
include { CUTADAPT } from '../../modules/nf-core/cutadapt/main.nf'

workflow GET_POLYA_READS {
    take:
    reads

    main:
    CUTADAPT_UNTRIMMED(
        reads
    )

    CUTADAPT(
        CUTADAPT_UNTRIMMED.out.reads
    )

    emit:
    reads = CUTADAPT.out.reads
}
