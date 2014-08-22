


library(topGO)



geneList
topDiffGenes


sampleGOdata <- new("topGOdata",
                    description = "Simple session",
                    ontology = "BP", 
                    allGenes = geneList, geneSel = topDiffGenes, 
                    nodeSize = 10,
                    annot = annFUN.db,
                    affyLib = affyLib
                    )

sampleGOdata

resultFisher <- runTest(sampleGOdata, algorithm = "classic", statistic = "fisher")
resultFisher


allRes <- GenTable(sampleGOdata, classicFisher = resultFisher, 
                   ranksOf = "classicFisher", topNodes = 10)

?GenTable
?runTest

