#setup
pacman::p_load("tidyverse", "data.table")
load("data/physiologies3.Rda")

#change empty cells to na
phys <- physiologies3 %>% mutate_all(na_if, "")
rm(physiologies3)

#get gram and oxygen utilization data
gram <- select(phys, genus, gram_type) %>% 
      drop_na()
resp <- select(phys, genus, respiration) %>% 
      drop_na()
bugsID <- tidyr::unite(bugsPhys, col = "bugsID", "4", "1", sep = "NCBI:txid", remove = F)

#function for organizing gram stains by keyword-----------------

#keyword must be a string in quotes
sortgram <- function(dat, gram_type, keyword){
      filter(.data = dat, gram_type == keyword) %>% 
            subset(select = -c(gram_type)) %>% 
            as.vector()
}

#function for organizing oxygen utilization by keyword-----------------

#keyword must be a string in quotes
sortresp <- function(dat, respiration, keyword, keyword2 = "none"){
      filter(.data = dat, respiration == keyword | respiration == keyword2) %>% 
            subset(select = -c(respiration)) %>% 
            as.vector()
}

#extracting/separating the data----------------------

gram_negative <- sortgram(gram, gram_type, "negative") %>% 
      rename(gram_negative = genus)
gram_positive <- sortgram(gram, gram_type, "positive") %>% 
      rename(gram_positive = genus)

facultatively_aerobic <- sortresp(resp, respiration, "facultatively aerobic") %>% 
      rename(facultatively_aerobic = genus)
obligately_aerobic <- sortresp(resp, respiration, "obligately aerobic") %>% 
      rename(obligately_aerobic = genus)
aerobic <- sortresp(resp, respiration, "aerobic", "aerobic and facultatively anaerobic") %>% 
   rename(aerobic = genus)

facultatively_anaerobic <- sortresp(resp, respiration, "facultatively anaerobic", "aerobic and facultatively anaerobic") %>% 
      rename(facultatively_anaerobic = genus)
obligately_anaerobic <- sortresp(resp, respiration, "obligately anaerobic") %>% 
      rename(obligately_anaerobic = genus)
anaerobic <- sortresp(resp, respiration, "anaerobic") %>% 
      rename(anaerobic = genus)

#Combine data frames and write to .gmt--------------
gmt_phys <- c(gram_negative, gram_positive, 
              facultatively_aerobic, obligately_aerobic, aerobic,
              facultatively_anaerobic, obligately_anaerobic, anaerobic)

EnrichmentBrowser::writeGMT(gmt_phys, gmt.file = "data/physiologies4.gmt")
