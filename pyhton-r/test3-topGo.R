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
             ontology="BP",
             allGenes=geneList,
             annot=annFUN.gene2GO,
             gene2GO=geneID2GO,
             nodeSize = 10
            )

resultFisher <- runTest(GOdata,
                  algorithm = "classic",
                  statistic = "fisher")

resultKS <- runTest(GOdata,
                    algorithm = "classic",
                    statistic = "ks")

resultKS.elim <- runTest(GOdata,
                         algorithm = "elim",
                         statistic = "ks")

result_table <- GenTable(GOdata, classicFisher = resultFisher, topNodes = 100)

result_table <- GenTable(GOdata, classicFisher = resultFisher,
                         classicKS = resultKS, elimKS = resultKS.elim,
                         orderBy = "elimKS", ranksOf = "classicFisher", topNodes = 100)


head(score(resultFisher))

result <- score(resultFisher)

##########################################################

r <- score(resultFisher)
r <- r[r<0.05]

index <- order(r, decreasing = FALSE)
r <- r[index]

termStat(GOdata, 'GO:0031625')

######################################

GOdata <- list()
resultFisher <- list()
resultKS <- list()

for( o in c("MF","BP","CC")) {
  
  g = new("topGOdata",
                ontology=o,
                allGenes=geneList,
                annot=annFUN.gene2GO,
                gene2GO=geneID2GO,
                nodeSize = 10
  )
  l <- list(g)
  GOdata <- append(GOdata, setNames(l, o))
  
  r <- runTest(g,
          algorithm = "classic",
          statistic = "fisher")
  
  l <- list(r)
  resultFisher <- append(resultFisher, setNames(l, o))
    
  r <- runTest(g,
                      algorithm = "classic",
                      statistic = "ks")
  
  
  r.elim <- runTest(g,
                           algorithm = "elim",
                           statistic = "ks")
  
  l <- list(r)
  resultKS <- append(resultKS, setNames(l, o))  
 
  
}

result_table <- GenTable(GOdata$MF, 
                         classicFisher = resultFisher$MF,
                         classicKS = resultKS, elimKS = resultKS.elim,
                         orderBy = "elimKS", ranksOf = "classicFisher", topNodes = 100)




##################################################################


resultFisher <- runTest(GOdataBP,
                        algorithm = "classic",
                        statistic = "fisher")

pvalClassic <- score(resultFisher)
pvalClassicBP <- score(resultFisherBP)

pvalClassic[pvalClassic < 0.01]
pvalClassicBP[pvalClassicBP < 0.001]


# result_table <- GenTable(GOdata,
#                    classic = result,
#                    orderBy = "classic",
#                    ranksOf = "classicFisher",
#                    topNodes = 10)

pvalClassic <- score(resultFisher)
pvalClassic[pvalClassic < quantile(pvalClassic,0.05)]
pvalClassic[pvalClassic < 0.01]

pval <- score(resultKS.elim)[names(pvalClassic)]
pval[pval < 0.0001]

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


#####################################################################
# TODO 
#
# verif resultat 
# clean python code
# ensembl map
# http://davetang.org/muse/
#####################################

ontology = ontology(GOdata)
whichTerms = "GO:0005694"

qTerms <- paste(paste("'", whichTerms, "'", sep = ""), collapse = ",")
retVal <- dbGetQuery(GO_dbconn(), paste("SELECT term, go_id FROM go_term WHERE ontology IN",
                                        "('", ontology, "') AND go_id IN (", qTerms, ");", sep = ""))

whichTerms = "GO:0000279"
qTerms <- paste(paste("'", whichTerms, "'", sep = ""), collapse = ",")
retVal <- dbGetQuery(GO_dbconn(), paste("SELECT ontology, go_id, term, definition FROM go_term WHERE go_id IN (", qTerms, ");", sep = ""))
retVal


r <- score(resultKS.elim)
r <- score(resultKS)
r['GO:0000279']
r['GO:0005694']

library(qvalue)
padjusted <- p.adjust(r, method="fdr") 
padjusted['GO:0000279']



retVal <- dbGetQuery(GO_dbconn(), paste("SELECT term, go_id FROM go_term;", sep = ""))

####xx <- as.list(GOTERM)
retVal$go_id
retVal$definition
retVal$term
retVal$term




############################################## convert ensembl to hgnc

library(biomaRt)
mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")


test <- 'ENSG00000118473'
getBM(attributes=c('ensembl_gene_id', "hgnc_symbol"), filters = "ensembl_gene_id", values=test, mart=mart)

test <- c('CDC46', 'CDC46', 'CDCA8', 'RAD51', 'RRM2', 'FIGNL1', 'BUB1', 'CCNB1', 'AURKB', 'KNTC1', 'CDC18L')
test <- 'CDC46'
getBM(attributes=c('ensembl_gene_id', "hgnc_symbol"), filters = "hgnc_symbol", values=test, mart=mart)




