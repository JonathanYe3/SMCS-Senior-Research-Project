pacman::p_load("curatedMetagenomicData", "EnrichmentBrowser", "SummarizedExperiment", "dplyr", "stringr", "qusage")

#get IBD data and format
nielson <- curatedMetagenomicData("NielsenHB_2014.metapfhlan_bugs_list.stool", dryrun = F, counts = T)
nielson.eset <- nielson[[1]]
nielson.se <- as(nielson.eset, "SummarizedExperiment")

#assign 1 to IBD, 0 to control
grp1 = ifelse(nielson.se$study_condition == "control", 0, 1)
nielson.se$GROUP = grp1

#differential expression analysis
nielson.se <- deAna(nielson.se, de.method = "DESeq2")

#reformat nielson microbe names, no MetaPhlAn string format
n <- rownames(nielson.se)

reg_name <- function(MetaPhlAn){
      temp <- sub(".*g__", "", MetaPhlAn)
      gsub("\\|s__.*", "", temp)
      gsub("_noname", "", temp)
}

for(i in 1L:length(n)){
      reg_names <- n %>% 
            reg_name()
}

rownames(nielson.se) <- reg_names

#set-based enrichment analysis
bergeys.gmt <- read.gmt("data/physiologies3.gmt")
results <- sbea("ora", nielson.se, bergeys.gmt, perm = 0)
gsRanking(results, signif.only = F)
