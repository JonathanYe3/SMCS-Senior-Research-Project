#load packages
pacman::p_load("curatedMetagenomicData", "EnrichmentBrowser", "SummarizedExperiment", "dplyr", "stringr", "EnhancedVolcano")

#Format data for SE
hmp <- readr::read_tsv("data/HMP2_dataset_mpa3.tsv")
hmp.metadata <- as.data.frame(t(hmp[1:135,])) %>% 
   janitor::row_to_names(row_number = 1)

assay <- as.matrix(hmp[136:nrow(hmp),])
rownames(assay) <- assay[,1]
assay <- assay[,-1]
mode(assay) <- "numeric"

hmp.se <- SummarizedExperiment(assays = list(assay), colData = hmp.metadata)

#get raw read counts
mode(hmp.se$number_reads) <- "numeric"
mode(assay(hmp.se)) <- "numeric"
hmp.counts = sweep(assay(hmp.se), 2, hmp.se$number_reads / 100, "*")
hmp.counts = round(hmp.counts)
hmp.se <- SummarizedExperiment(assays = hmp.counts, colData = hmp.metadata)

#location subset
hmp.se <- subset(hmp.se, , hmp.se$location == "Boston")

grp1 = ifelse(hmp.se$study_condition == "IBD", 1, 0)
hmp.se$GROUP = grp1

# keep only those microbes that have more than 430 read counts across all samples
keep <- rowSums(assay(hmp.se)) > 430
hmp.se <- hmp.se[keep,]

#differential expression analysis
assay(hmp.se) <- assay(hmp.se) + 1
hmp.se <- deAna(hmp.se, de.method = "DESeq2", filter.by.expr = F)

#reformat hmp microbe names, no MetaPhlAn string format
getGenus <- function(MetaPhlAn){
   temp <- gsub("s__", "", MetaPhlAn)
   gsub("_noname", "", temp)
   gsub("\\_.*","",temp)
}

rownames(hmp.se) <- vapply(rownames(hmp.se), getGenus,
                               character(1), USE.NAMES = FALSE)

#set-based enrichment analyses
bergeys.gmt <- qusage::read.gmt("data/physiologies4.gmt")

#ORA
ora_results <- sbea("ora", hmp.se, bergeys.gmt, perm = 0)
ora_results <- gsRanking(ora_results, signif.only = F) %>% 
      as.data.frame()
colnames(ora_results) <- sub("GENE", "MICROBE", colnames(ora_results))

#GSEA
gsea_results <- sbea("gsea", hmp.se, bergeys.gmt, perm = 1000)
gsea_results <- gsRanking(gsea_results, signif.only = F) %>% 
      as.data.frame()
colnames(gsea_results) <- sub("GENE", "MICROBE", colnames(gsea_results))

gdata::keep(hmp.se, ora_results, gsea_results, sure = TRUE)

#Volcano plot for DESeq2 results
png(filename = "data/results/HMP_Volcano.png", res = 300, height = 1800, width = 1800)
EnhancedVolcano(rowData(hmp.se),
                lab = rownames(hmp.se),
                pCutoff = 0.05,
                title = 'Differential Abundance[DESeq2]',
                x = 'FC',
                y = 'ADJ.PVAL',
                xlim = c(-20, 20),
                pointSize = 1.5)
dev.off()

#Results output
write.csv(ora_results, file = "data/results/HMP_ora_res.csv")
write.csv(gsea_results, file = "data/results/HMP_gsea_res.csv")
