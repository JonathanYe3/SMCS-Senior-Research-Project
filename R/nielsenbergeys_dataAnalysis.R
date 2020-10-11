pacman::p_load("curatedMetagenomicData", "EnrichmentBrowser", "SummarizedExperiment", "dplyr", "stringr", "EnhancedVolcano")

#get IBD data and format
nielsen <- curatedMetagenomicData("NielsenHB_2014.metaphlan_bugs_list.stool", dryrun = F, counts = T)
nielsen.eset <- nielsen[[1]]
nielsen.se <- as(nielsen.eset, "SummarizedExperiment")

#assign 1 to IBD, 0 to control
grp1 = ifelse(nielsen.se$study_condition == "control", 0, 1)
nielsen.se$GROUP = grp1

# keep only those microbes that have more than 148 read counts across all samples
keep <- rowSums(assay(nielsen.se)) > 148
nielsen.se <- nielsen.se[keep,]

#differential expression analysis
nielsen.se <- deAna(nielsen.se, de.method = "DESeq2", filter.by.expr = F)

#reformat nielsen microbe names, no MetaPhlAn string format
reg_name <- function(MetaPhlAn){
      temp <- sub(".*g__", "", MetaPhlAn)
      gsub("\\|s__.*", "", temp)
      gsub("_noname", "", temp)
}

.getLast <- function(n)
{   
   spl <- unlist(strsplit(n, "\\|"))
   spl[length(spl)]
}

rownames(nielsen.se) <- vapply(rownames(nielsen.se), reg_name,
                               character(1), USE.NAMES = FALSE)
rownames(nielsen.se) <- vapply(rownames(nielsen.se), .getLast,
                               character(1), USE.NAMES = FALSE)

#set-based enrichment analyses
bergeys.gmt <- qusage::read.gmt("data/physiologies4.gmt")

#ORA
ora_results <- sbea("ora", nielsen.se, bergeys.gmt, perm = 0)
ora_res.df <- gsRanking(ora_results, signif.only = F) %>% 
      as.data.frame()
colnames(ora_res.df) <- sub("GENE", "MICROBE", colnames(ora_res.df))

#GSEA
gsea_results <- sbea("gsea", nielsen.se, bergeys.gmt, perm = 1000)
gsea_res.df <- gsRanking(gsea_results, signif.only = F) %>% 
      as.data.frame()
colnames(gsea_res.df) <- sub("GENE", "MICROBE", colnames(gsea_res.df))

#Volcano plot for DESeq2 results
png(filename = "data/results/Volcano2.png", res = 300, height = 1800, width = 1800)
#pdf("data/results/Volcano2.pdf")
EnhancedVolcano(rowData(nielsen.se),
                lab = rownames(nielsen.se),
                pCutoff = 0.05,
                title = 'Differential Abundance[DESeq2]',
                x = 'FC',
                y = 'ADJ.PVAL',
                xlim = c(-40, 40),
                pointSize = 1.5)
dev.off()

#Results output
write.csv(ora_res.df, file = "data/results/B_ora_res.csv")
write.csv(gsea_res.df, file = "data/results/B_gsea_res.csv")
