library(curatedMetagenomicData)
library(bugSetEnrichment)
library(EnrichmentBrowser)
library(SummarizedExperiment)
library(Biobase)
library(coin)

nielson <- curatedMetagenomicData("NielsenHB_2014.metaphlan_bugs_list.stool", dryrun = F)

#group 1, study_condition: 0 is control, 1 is IBD
grp1=ifelse(nielson[[1]]$study_condition == "control", 0, 1)
nielson[[1]]$GROUP=grp1

#group 2, age_category: 1 is adult,  0 is senior and schoolage
grp2=ifelse(nielson[[2]]$age_category == "adult", 1, 0)
nielson[[2]]$GROUP=grp2

#differential expression analysis

z=deAna(nielson[[1]], de.method="coin", filter.by.expr = FALSE)
