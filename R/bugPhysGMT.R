pacman::p_load("dplyr", "googlesheets4", "data.table")

#specify google sheet and read into dataframe
webpage <- "https://docs.google.com/spreadsheets/d/1Vp5uVi_WhX-f33sR-I7azWcSnC5bkka75SnRsMfZf3U"
gs4_auth(email = "jonathanye2003@gmail.com")
getBugsPhys <- gs4_get(webpage)
bugsPhys <- read_sheet(webpage, col_names = F) %>% as.data.frame()

colnames(bugsPhys) <- 1:length(bugsPhys)
bugsPhys[bugsPhys == c("NULL", "NA")] <- NA
bugsPhys <- bugsPhys[rowSums(is.na(bugsPhys)) != ncol(bugsPhys),]
bugsID <- tidyr::unite(bugsPhys, col = "bugsID", "4", "1", sep = " NCBI:txid", remove = F)
bugsID <- tidyr::unite(bugsID, col = "moreBugsID", "3", "2", sep = ":",remove = F)
bugsID[bugsID == "NA:NA"] <- ""

bugsID <- tidyr::unite(bugsID, col = "bugsID", "bugsID", "moreBugsID", sep = " ")

#functions for sorting
sortfunc <- function(data, keyword, strict = FALSE, title = keyword){
      data <- rccmisc::lownames(data) %>% tidyr::drop_na()
      data[[2]] <- tolower(data[[2]])

      if(strict == TRUE){
            temp <- data[data[[2]] == keyword,]
            temp <- temp[1]
      }
      else{
            temp <- data[data[[2]] %like% keyword,]
            temp <- temp[1]
      }
      colnames(temp) <- title
      as.vector(temp)
}

addElements <- function(dat, traits, prefix, strict = F){
   for (i in 1L:length(traits)) {
      currTitle <- paste0(prefix ,traits[[i]])
      if(strict == T){
         newElement <- sortfunc(dat, traits[[i]], title = currTitle, strict = T)
      }
      else{
         newElement <- sortfunc(dat, traits[[i]], title = currTitle, strict = F)
      }
      bugs <- append(bugs, newElement)
   }
   bugs
}

#gram stain
bugs <- vector("list")
gram <- bugsPhys[c("4", "6")]
gram.traits <- c("positive", "negative", "var")
bugs <- addElements(dat = gram, gram.traits, "gram.")

#oxygen utilization/respiration
resp <- bugsPhys[c("4", "11")]
resp.strictTraits <- c("aerobic", "anaerobic")
resp.traits <- c("obligately aerobic", "microaerobic", "facultatively aerobic",
                 "obligately anaerobic", "facultatively anaerobic", "microaerophilic")
bugs <- addElements(dat = resp, resp.strictTraits, "resp.", strict = T)
bugs <- addElements(dat = resp, resp.traits, "resp.")

#cell Shape
shape <- bugsPhys[c("4", "27")]
shape.traits <- c("coccoid", "coccobacilli", "ovoid", "rod", "pleomorphic", "spir", 
                  "fila", "helical")
bugs <- addElements(dat = shape, shape.traits, "shape.")

#Cell arrangement
arrangement <- bugsPhys[c("4", "33")]
arrangement.traits <- c("aggregates", "pairs", "netlike", "single", "v-formation",
                        "branched", "filamentous", "forked", "chains", "palisades",
                        "nonaggregated", "clusters", "multicellular", "groups", "tetrads",
                        "rosettes")
arrangement.strictTraits <- c("encapsulated", "unencapsulated")
bugs <- addElements(dat = shape, arrangement.traits, "shape.")
bugs <- addElements(dat = shape, arrangement.strictTraits, "shape.", strict = T)


EnrichmentBrowser::writeGMT(bugs, gmt.file = "data/bugsPhys.gmt")
bugsPhys.gmt <- qusage::read.gmt("data/bugsPhys.gmt")


