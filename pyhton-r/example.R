##################################
# Functional Enrichment Analysis #
# using topGO                    #
##################################
## Install topGO and Affymetrix Human Genome U133 Plus 2.0 Array annotation data
source("http://bioconductor.org/biocLite.R")
biocLite()
biocLite("topGO")
biocLite("hgu133plus2.db")
## boot the gaggle
library(gaggle)
gaggleInit()
## load libraries
library(topGO)
library(hgu133plus2.db)

####################
# Data preparation #
####################
### Initializing the analysis ###
# hgu133plus2ACCNUM is an R object that contains mappings between the manufacturers identifiers and gene names
# of Affymetrix Human Genome U133 Plus 2.0 Array.(Other annotation packages can be found on the 
# Bioconductor website: http://www.bioconductor.org/packages/release/data/annotation/ )
# all.genes: all background genes ( gene universe )
all.genes <- ls(hgu133plus2ACCNUM) 

### Make gene lists ###
# We will make a list that includes two sets of genes of interest
# Initialize the list:
glioblastoma.genes = list()
# Here we use gaggle to broadcast the genes from the following websites to R: 
# Go to http://baliga.systemsbiology.net/events/sysbio/content/bicluster-307
# broadcast 'bicluster 307 genes' to R
# then get the genelist from R, assign to the first element of glioblastoma.genes,
# and at the same time name this element as "bc307":
glioblastoma.genes[["bc307"]] = tolower(getNameList())
# then go to http://baliga.systemsbiology.net/events/sysbio/content/bicluster-353
# broadcast 'bicluster 353 genes' to R
glioblastoma.genes[["bc353"]] = tolower(getNameList())

#######################
# Make topGO object   #
#######################
# limit genes to those in bc353
relevant.genes <- factor(as.integer(all.genes %in% glioblastoma.genes[["bc353"]]))
names(relevant.genes) <- all.genes
## Make Biological Process GOData object ##
# ontology: character string specifying the ontology of interest 
#           ('BP':Biological process,'MF':molecular function or 'CC':cellular components)
# allGenes: named vector of type numeric or factor. The names attribute contains the genes identifiers. The
#           genes listed in this object define the gene universe.
# annotationFun: function that maps gene identifiers to GO terms. annFUN.db extracts the gene-to-GO mappings from the affyLib object
# affyLib: character string containing the name of the Bioconductor annotaion package for a specific microarray chip.
GOdata.BP <- new("topGOdata", ontology='BP', allGenes = relevant.genes, annotationFun = annFUN.db, affyLib = 'hgu133plus2.db')
############################
# Run enrichment analysis  #
############################
results <- runTest(GOdata.BP, algorithm = 'classic', statistic = 'fisher')

#######################
# Analysis of results #
#######################
## Biological Processes
# Basic information on input data can be accessed using the geneData function. The number of annotated
# genes, the number of significant genes (if it is the case), the minimal size of a GO category as well as the
# number of GO categories which have at least one signficant gene annotated are listed:

# generate a summary of the enrichment analysis
results.table <- GenTable(GOdata.BP, results, topNodes = length(results@score))
# How many GO terms were tested?
dim(results.table)[1]
# reduce results to GO terms passing Benjamini-Hochberg multiple hypothesis corrected pval <= 0.05, FDR <= 5%
results.table.bh <- results.table[which(p.adjust(results.table[,"result1"],method="BH")<=0.05),]
# How many terms are enriched?
dim(results.table.bh)[1]
# What are the top ten terms?
results.table.bh[1:10,]
# Get genes in most sig terms
# Get all the genes annotated to the top one GO term:
GOid.of.interest = results.table.bh[1,"GO.ID"]
all.term.genes = genesInTerm(GOdata.BP,GOid.of.interest)[[1]]
# Which of these genes is in the bicluster?
genes.of.interest <- intersect(glioblastoma.genes[["bc353"]],all.term.genes)
# print table with probe ID and gene symbol
toTable(hgu133plus2SYMBOL[genes.of.interest])

#######################
# Extension           #
#######################
#We can also look at multiple GO terms at the same time:
GOids.of.interest = results.table.bh[c(1:10),"GO.ID"]
all.term.genes = genesInTerm(GOdata.BP, GOids.of.interest)
# Which of these genes is in the bicluster?
genes.of.interest <- sapply(names(all.term.genes),function(x){intersect(all.term.genes[[x]],glioblastoma.genes[["bc353"]])})
# print table with probe ID and gene symbol
geneSynmol.of.interest <- lapply(names(genes.of.interest),function(x){toTable(hgu133plus2SYMBOL[genes.of.interest[[x]]])})
names(geneSynmol.of.interest)<- GOids.of.interest

#######################
# Automation          #
#######################
results <- list()
for( bc in names(glioblastoma.genes) ) {
  cat(paste("Computing functional enrichment for...",bc,"\n"))
  relevant.genes <- factor(as.integer(all.genes %in% glioblastoma.genes[[bc]]))
  names(relevant.genes) <- all.genes
  ## Make Biological Process GOData object ##
  # ontology: character string specifying the ontology of interest 
  #           ('BP':Biological process,'MF':molecular function or 'CC':cellular components)
  # allGenes: named vector of type numeric or factor. The names attribute contains the genes identifiers. The
  #           genes listed in this object define the gene universe.
  # annotationFun: function that maps gene identifiers to GO terms. annFUN.db extracts the gene-to-GO mappings from the affyLib object
  # affyLib: character string containing the name of the Bioconductor annotaion package for a specific microarray chip.
  GOdata.BP <- new("topGOdata", ontology='BP', allGenes = relevant.genes, annotationFun = annFUN.db, affyLib = 'hgu133plus2.db')
  # Run enrichment analysis 
  result.BP <- runTest(GOdata.BP, algorithm = 'classic', statistic = 'fisher')
  results[[bc]] <- GenTable(GOdata.BP,result.BP,topNodes = length(result.BP@score))
}

#####################################
# EXTRAS!!!
##################################### 

# Make GOdata object
relevant.genes <- factor(as.integer(all.genes %in% glioblastoma.genes[[1]]))
names(relevant.genes) <- all.genes
GOdata.BP <- new("topGOdata", ontology='BP', allGenes = relevant.genes, annotationFun = annFUN.db, affyLib = 'hgu133plus2.db')

#####################################
# working with the topGOdata object #
#####################################
## Getting familiar with the object ##
# Brief description of the GOdata object:
GOdata.BP
# obtaining all genes that can be used in the analysis:
a = genes(GOdata.BP)
# Check the structure of a:
str(a)
# number of genes that can be used in the analysis:
numGenes(GOdata.BP)
# The list of significant genes(our genes of interest) can also be accessed by:
sg = sigGenes(GOdata.BP)
str(sg)
# The GO terms are available for analysis:
ug = usedGO(GOdata.BP)
str(ug)
# Try selecting 10 random GO terms, count the number of annotated genes and obtain their annotation:
sel.terms <- sample(usedGO(GOdata.BP), 10)
# countGenesInTerm: to see the number of annotated genes in specified GO terms
num.ann.genes <- countGenesInTerm(GOdata.BP, sel.terms)
# Check results
num.ann.genes
# genesInTerm: list the annotated genes in specified GO terms
ann.genes <- genesInTerm(GOdata.BP, sel.terms)
# Check the above results
str(ann.genes)


##########################
# other useful functions #
##########################
# how the signicant GO terms are distributed over the GO graph
showSigOfNodes(GOdata.BP, score(r1.BP), firstSigNodes = 5, useInfo = 'all')
# output the above graph to a pdf file
printGraph(GOdata.BP, r1.BP, firstSigNodes = 5, useInfo = 'all', fn.prefix = paste("GOgraph_",geneSet,"_",names(glioblastoma.genes)[geneSet],sep=""),pdfSW=TRUE)

##############################
# other algorithms supported #
##############################
r1.BP = runTest(GOdata.BP, algorithm = 'classic', statistic = 'fisher')
# other algorithms that consider the GO structure-'elim' and 'weight':
r1.BP.elim = runTest(GOdata.BP, algorithm = 'elim', statistic = 'fisher')
r1.BP.weight = runTest(GOdata.BP, algorithm = 'weight', statistic = 'fisher')
# generate a summary of the results of the enrichment analysis, including the results from all 3 algorithms:
GenTable(GOdata.BP, classicFisher = r1.BP,elimFisher = r1.BP.elim,weightFisher = r1.BP.weight,weight01Fisher = r1.BP.weight01,orderBy = "elimFisher", ranksOf = "classicFisher", topNodes = 10)
# get the enrichment p-values from all three algorithms and compare
pValue.classic <- score(r1.BP)
pValue.elim <- score(r1.BP.elim)[names(pValue.classic)]
pValue.weight <- score(r1.BP.weight)[names(pValue.classic)]
# for plots of p-values: set plotting parameters:
gstat <- termStat(GOdata.BP, names(pValue.classic))
gSize <- gstat$Annotated / max(gstat$Annotated) * 4
# function for building color map for plotting
colMap <- function(x) {
  .col <- rep(rev(heat.colors(length(unique(x)))), time = table(x))
  return(.col[match(1:length(x), order(x))])
}
gCol <- colMap(gstat$Significant)
plot(pValue.classic, pValue.elim, xlab = "p-value classic", ylab = "p-value elim",pch = 19, cex = gSize, col = gCol)
plot(pValue.classic, pValue.elim, xlab = "p-value classic", ylab = "p-value elim",pch = 19, cex = gSize, col = gCol,log='xy')
plot(pValue.classic, pValue.weight, xlab = "p-value classic", ylab = "p-value weight",pch = 19, cex = gSize, col = gCol)