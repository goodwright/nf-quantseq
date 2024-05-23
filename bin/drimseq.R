library(data.table)
library(DRIMSeq)
library(BiocParallel)
library(tidyverse)
library(stageR)

setwd("/Users/llywelyngriffith/Documents/AZ_postdoc/CDK11i_3seq/Nobby_APA_analysis/")

#define PAS thresholds using PAS atlas
pas.dt = fread('merged_polya.filteredunique.annotated.bed', col.names = c("seqnames", "start", "end", "id", "score", "strand", "ensg", "hgnc", "region"))
pas.dt = pas.dt[region == "UTR3"]
pas.dt[, total_score := sum(score), by = ensg]
pas.dt[, percent_5 := 0.05 * total_score]
pas.dt[, percent_1 := 0.01 * total_score]
pas.dt[, multi := .N, by = ensg]

# Load and merge all count files
all.count.files = list.files("counts", full.names = TRUE, pattern = ".bed")
count.list = lapply(1:length(all.count.files), function(i) fread(all.count.files[[i]], select = c(4, 7, 8, 10), col.names = c("PAS", "ensg", "hgnc", gsub("_r1.polya_trimmed_window_merged_polya.bed", "", basename(all.count.files)[i]))))
count.dt = Reduce(function(x, y) merge(x = x, y = y, by = c("PAS", "ensg", "hgnc")), count.list)

# Remove intergenic PAS and those genes with only 1 PAS
count.dt = count.dt[ensg != "intergenic"]
count.dt = count.dt[PAS %in% pas.dt[multi > 1]$id]

#make columns called gene_id and feature_id for drimseq
counts.df = as.data.frame(count.dt[, `:=` (gene_id = paste0(hgnc, "_", ensg),
                                           feature_id = PAS)])

counts.df=counts.df[,-c(1:3)]

#make metadata dataframe
quantseq_metadata = read.table("quantseq_metadata.txt",
                               header = TRUE, as.is = TRUE)

quantseq_samples = data.frame(sample_id = quantseq_metadata$SampleName,
                              group = quantseq_metadata$Condition)


##run drim-seq
quantseq_samples$group = factor(quantseq_samples$group, levels = c("DMSO","CDK11i"))

#create drimseq object
d = dmDSdata(counts = counts.df, samples = quantseq_samples)
head(counts(d), 3)
plotData(d)

# Require gene expression in at least 75% of all samples, and PAS expression in 75% of samples for either DMSO or inhibition (whichever has fewest)
d = dmFilter(d, 
             min_samps_gene_expr = floor(nrow(quantseq_samples) * 0.75), 
             min_samps_feature_expr = floor(min(table(quantseq_samples$group)) * 0.75), 
             min_gene_expr = 10, 
             min_feature_expr = 5)

#drimseq design
design_full = model.matrix(~ group, data = DRIMSeq::samples(d))
set.seed(42)
d = dmPrecision(d, design = design_full, verbose = 1)
d = dmFit(d, design = design_full, verbose = 1)
plotPrecision(d)
head(coefficients(d))

#run drimseq test
d = dmTest(d, coef = "groupCDK11i", verbose = 1)
View(results(d))
View(results(d,level="feature"))

#visualise results
res_genes = results(d)
res_pA_sites = results(d,level="feature")
res_genes = res_genes[order(res_genes$adj_pvalue, decreasing = FALSE), ]
res_pA_sites = res_pA_sites[order(res_pA_sites$adj_pvalue, decreasing = FALSE), ]
top_gene_id = res_genes$gene_id[1]
plotProportions(d, gene_id = top_gene_id, group_variable = "group")

#change NA values to 1s to prevent problems with two stage test (ref = https://f1000research.com/articles/7-952/v3)
no.na = function(x) ifelse(is.na(x), 1, x)
res_genes$pvalue = no.na(res_genes$pvalue)
res_pA_sites$pvalue = no.na(res_pA_sites$pvalue)

##two stage test
#assign gene-level pvalues to the screening stage
pScreen = res_genes$pvalue
names(pScreen) = res_genes$gene_id

#assign transcript level pvalues to the confirmation stage
pConfirmation = matrix(res_pA_sites$pvalue, ncol = 1)
rownames(pConfirmation) = res_pA_sites$feature_id

#create the gene-transcript mapping
tx2gene = res_pA_sites[,c("feature_id","gene_id")]

#create the stageRTx object and perform the stage-wise analysis
stageR0bj = stageRTx(pScreen = pScreen, 
                     pConfirmation = pConfirmation, 
                     pScreenAdjusted = FALSE, 
                     tx2gene = tx2gene)

stageR0bj = stageWiseAdjustment(object = stageR0bj, method = "dtu", alpha = 0.05)

getSignificantGenes(stageR0bj)
getSignificantTx(stageR0bj)

two_stage_padj=getAdjustedPValues(stageR0bj, order=FALSE, onlySignificantGenes = FALSE)

colnames(two_stage_padj)[1]='gene_id'
colnames(two_stage_padj)[2]='feature_id'
colnames(two_stage_padj)[3]='twostep_gene_padj'
colnames(two_stage_padj)[4]='twostep_transcript_padj'

#add column of difference between mean usage and only keep those with > 10% change
proportions=proportions(d)
proportions$mean_diff = (rowMeans(proportions[6:8]) - rowMeans(proportions[3:5]))
res_pA_sites_with_mean_diff = left_join(res_pA_sites,proportions, by = "feature_id")

#combine two-test table with mean diff table
two_step_pA_sites_with_mean_diff=inner_join(two_stage_padj,res_pA_sites_with_mean_diff,by='feature_id')
sig_two_step_pA_sites_with_mean_diff=two_step_pA_sites_with_mean_diff %>% filter(twostep_transcript_padj <= 0.05 & abs(mean_diff)>=0.1)

#write to file
write_csv(two_step_pA_sites_with_mean_diff,"two_step_CDK11i_pA_sites_with_mean_diff.csv")
write_csv(sig_two_step_pA_sites_with_mean_diff,"sig_two_step_CDK11i_pA_sites_with_mean_diff.csv")
