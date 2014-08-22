library(topGO)
library(qvalue)

setwd("/work/bio-informatics/pyhton-r")

file <- 'gene_association.goa_human_noHeader'
sep <- "\t"
IDsep <- ","

a <- read.delim(file = file, header = FALSE, quote = "", 
                sep = sep, colClasses = "character")

map <- a[, 5]
names(map) <- a[, 2]
map <- lapply(map, function(x) gsub(" ", "", strsplit(x, split = IDsep)[[1]]))


my.data <- read.csv('gene_data_hs.csv')
geneList <- my.data$pvalue
names(geneList) <- my.data$DB.Object.ID

topDiffGenes = function(allScore) {
  return (allScore < 0.01)
}

topDiffGenes(geneList)

sampleGOdataBP <- new("topGOdata",
                    description = "Simple session",
                    ontology = "BP", 
                    geneSelectionFun = topDiffGenes,
                    allGenes = geneList, 
                    annot = annFUN.gene2GO,
                    gene2GO = map
)
resultFisherBP <- runTest(sampleGOdataBP, algorithm = "classic", statistic = "fisher")

"""
s <- score(resultFisherBP)
pval <- c()
for (i in 1:length(s)){
  pval[i] <- s[i][1]
}
p.adjust(pval, method='fdr')
"""





