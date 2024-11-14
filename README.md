# nf-quantseq

## Usage
You must have Nextflow and either Docker/Singularity installed.

An example usage, FASTQ must be single end and provided as a reference to a directory.
The organism provided as "org" must be one of: "hg38", "mm10", "mm39", "rn6", "rn7".

```
nextflow run ulelab/nf-quantseq -r feat-userInputSites \
    -profile crick \
    -resume \
    --input '/camp/lab/ulej/home/users/griffil2/RNA_seq/CDK11i_seq/fastq_r1/*.fastq.gz' \
    --outdir ./results \
    --fasta $REFDIR/GRCh38.primary_assembly.genome.fa \
    --gtf /camp/home/griffil2/home/RNA_seq/gencode.v45.annotation.gtf.gz \
    --org hg38 \
    --minreads 10
```


## Additional params
### Provide PolyA sites from a previous run
If you have a set of PolyA sites from a previous pipeline run that you would like to quantify against, provide the bed with the `--polya_bed` parameter. It will be in the polyaclusters folder and named `merged_polya.filteredunique.annotated.bed`. You can also provide a bed file that wasn't produced by this pipeline, ensure the polyA sites are single nucleotide positions.

### TSO sequencing
If you sequenced using TSO protocol you will have NNNNNGGG at the start of your reads. Use param `--tso` for removing the GGG and param `--umitools_bc_pattern 'NNNNN'` to extract the UMI to the read header and trigger PCR duplicate removal before PolyA site atlas generation and quantification.

### Quantseq Reverse
The pipeline is defaulted to dealing with libraries that produce forward stranded reads. Quantseq REV produces reverse stranded reads and so read counting must be on the opposite strand. To do this provide the option `-quantseq_rev` at command line or `quantseq_rev=true` in config file.

*!!Note!! currently you must provide a polya_bed with this option.*

