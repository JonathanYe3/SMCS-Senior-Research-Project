pacman::p_load("curatedMetagenomicData", "EnrichmentBrowser", "SummarizedExperiment", "dplyr", "stringr")

#get IBD data and format
nielsen <- curatedMetagenomicData("NielsenHB_2014.metaphlan_bugs_list.stool", dryrun = F, counts = T)
nielsen.eset <- nielsen[[1]]
nielsen.se2 <- as(nielsen.eset, "SummarizedExperiment")

#assign 1 to IBD, 0 to control
grp1 = ifelse(nielsen.se2$study_condition == "control", 0, 1)
nielsen.se2$GROUP = grp1

# keep only those microbes that have more than 148 read counts across all samples
keep <- rowSums(assay(nielsen.se2)) > 148
nielsen.se2 <- nielsen.se2[keep,]

#differential expression analysis
nielsen.se2 <- deAna(nielsen.se2, de.method = "DESeq2", filter.by.expr = F)

#reformat nielsen microbe names, no MetaPhlAn string format
n <- rownames(nielsen.se2)

reg_name <- function(MetaPhlAn){
      temp <- sub(".*s__", "", MetaPhlAn)
      tolower(temp)
      #gsub("\\|s__", "_", temp)
}

for(i in 1L:length(n)){
      reg_names <- n %>% 
            reg_name()
}

rownames(nielsen.se2) <- reg_names

#set-based enrichment analysis
basicChars.gmt <- qusage::read.gmt("data/basic_characteristics_SPECIES.gmt")
results <- sbea("ora", nielsen.se2, basicChars.gmt, perm = 0)
gsRanking(results, signif.only = F)

results2 <- sbea("gsea", nielsen.se2, basicChars.gmt, perm = 1000)
gsRanking(results2, signif.only = F)

