# nf-quantseq

**drimseq.R script**
Required input files:
- poly(A) atlas : merged_polya.filteredunique.annotated.bed from /polyaclusters main pipeline output folder
- bed count files : e.g CDK11i_1_r1.polya_trimmed_window_merged_polya.bed from /counts main pipeline output folder
- metadata.txt file, e.g:
      SampleName	Condition
      DMSO_1    	DMSO
      DMSO_2    	DMSO
      DMSO_3    	DMSO
      CDK11i_1	  CDK11i
      CDK11i_2  	CDK11i
      CDK11i_3	  CDK11i
