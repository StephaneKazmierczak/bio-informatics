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
names(map) <- a[, 3]
map <- lapply(map, function(x) gsub(" ", "", strsplit(x, split = IDsep)[[1]]))


my.data <- read.csv('gene_data_hs.csv')
geneList <- my.data$pvalue
names(geneList) <- my.data$DB.Object.Symbol

topDiffGenes = function(allScore) {
  return (allScore < 0.01)
}

GOdataBP <- new("topGOdata",
                    description = "Simple session",
                    ontology = "CC", 
                    geneSelectionFun = topDiffGenes,
                    allGenes = geneList, 
                    nodeSize = 1,
                    annot = annFUN.gene2GO,
                    gene2GO = map
)
resultFisherBP <- runTest(GOdataBP,
                          algorithm = "classic",
                          statistic = "fisher")



allRes <- GenTable(GOdataBP,
                   classic = resultFisherBP,
                   orderBy = "classic",
                   ranksOf = "classicFisher",
                   topNodes = 10)


"""
library(qvalue)
s <- score(resultFisherBP)
pval <- c()
for (i in 1:length(s)){
pval[i] <- s[i][1]
}
p.adjust(pval, method='fdr')
"""


