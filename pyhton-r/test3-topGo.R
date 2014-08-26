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
                    ontology = "MF", 
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

########################################################
library(topGO)
setwd("/work/bio-informatics/pyhton-r")

genes2go <- readMappings("gene_anno_custom")
str(head(genes2go))

testset <- scan(file='subset.txt', what=character())
testset <- sample(testset, length(testset)/2)
refset <- names(genes2go)

genes_of_interest <- factor(as.integer(refset %in% testset))
genes_of_interest <- setNames(genes_of_interest, refset)


GOdata = new("topGOdata",
               ontology="MF",
               annot=annFUN.gene2GO,
               allGenes=genes_of_interest,
               gene2GO=genes2go
)

#######################################################

library(topGO)
setwd("/work/bio-informatics/pyhton-r")

geneID2GO = readMappings("gene_anno_custom")

geneNames = names(geneID2GO)
myInterestingGenes = sample(geneNames, length(geneNames)/4)
geneList = factor(as.integer(geneNames %in% myInterestingGenes))
names(geneList) <- geneNames

GOdata = new("topGOdata",
             ontology="MF",
             allGenes=geneList,
             annot=annFUN.gene2GO,
             gene2GO=geneID2GO)

result <- runTest(GOdata,
                  algorithm = "classic",
                  statistic = "fisher")

result_table <- GenTable(GOdata,
                   classic = result,
                   orderBy = "classic",
                   ranksOf = "classicFisher",
                   topNodes = 10)


###################################################

geneID2GO <- readMappings(file = system.file("examples/geneid2go.map", package = "topGO"))
str(head(geneID2GO))

geneNames <- names(geneID2GO)
myInterestingGenes <- sample(geneNames, length(geneNames) / 80)
geneList <- factor(as.integer(geneNames %in% myInterestingGenes))
names(geneList) <- geneNames

GOdata <- new("topGOdata", ontology = "MF", allGenes = geneList,
              annot = annFUN.gene2GO, gene2GO = geneID2GO)





######################################################

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


geneNames <- scan(file='subset.txt', what=character())
myInterestingGenes = sample(geneNames, length(geneNames)/2)
geneList <- factor(as.integer(geneNames %in% myInterestingGenes))
names(geneList) <- geneNames

GOdata <- new("topGOdata",
                ontology = "MF", 
                allGenes = geneList, 
                annot = annFUN.gene2GO,
                gene2GO = map
)


result <- runTest(GOdata)




#####################################
# verif resultat 
# clean python code
# ensembl map
#
############


