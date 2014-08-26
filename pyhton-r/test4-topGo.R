library(topGO)

setwd("/work/bio-informatics/pyhton-r")

file <- 'gene_association.goa_human_noHeader'
sep <- "\t"
IDsep <- ","

a <- read.delim(file = file,
                header = FALSE,
                quote = "", 
                sep = sep,
                colClasses = "character")

map <- a[, 5]
str(map)
names(map) <- a[, 3]



?read.delim

a[3][1]


map <- lapply(map, function(x) gsub(" ", "", strsplit(x, split = IDsep)[[1]]))


str(map)



file <- "gene_association.goa_human_noHeader"
con  <- file(file, open = "r")

dataList <- list()
ecdfList <- list()

map <- list()

while (length(oneLine <- readLines(con, n = 1, warn = FALSE)) > 0) {
  
  linVec <- (strsplit(oneLine, "\t"))[[1]]
  
  if(map$linVec[3])
  
  map <- append ()
  

} 
close(con)



#######################


subset <- scan(file="subset.txt", what=character())
geneList <- factor(as.integer (all_genes %in% all_genes))
geneList <- setNames(geneList, all_genes)

