pacman::p_load("dplyr", "bugphyzz", "tidyr", "EnrichmentBrowser")

phys <- physiologies() %>% 
      lapply(as_tibble)

get_signatures <- function(df)
{
      temp <- split(df, df$Attribute)
      temp <- lapply(temp, function(x) subset(x, select = Taxon_name))
      temp <- lapply(temp, dplyr::pull, name = Taxon_name)
}

bugphyzz_GMT <- lapply(phys, get_signatures) %>% 
      unlist(recursive = F)
EnrichmentBrowser::writeGMT(bugphyzz_GMT, gmt.file = "data/bug_physiologies.gmt")
