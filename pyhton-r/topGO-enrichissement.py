#!/usr/local/bin/python2.7
# -*- coding: utf-8 -*-
#
#       go-uniprot-to-topGO.py
#==============================================================================
import argparse
import sys
import os
from rpy2.robjects.packages import importr
from rpy2.robjects import r as R
#==============================================================================
#Command line options==========================================================
#==============================================================================
parser = argparse.ArgumentParser()
parser.add_argument("gene_association", type=str,
                    help="A Uniprot gene association file")
parser.add_argument("test_ids", type=str,
                    help="A file containing the IDs of interest (the subset)")
parser.add_argument("-a", "--algorithm", type=str, default="weight01",
                    help="The algorithm to use for the overrepresentation test")
parser.add_argument("-m", "--min", type=int, default=1,
                    help="Minimum number of genes per GO term for the GO term to be tested. default: 1")
parser.add_argument("-o", "--overrepresentationaction", action='store_true',
                    default=True, help="Test for overrepresentation")
parser.add_argument("-u", "--underrepresentation", action='store_true',
                    default=False, help="Test for underrepresentation")
if len(sys.argv) == 1:
    parser.print_help()
    sys.exit(1)
args = parser.parse_args()
#==============================================================================

def init_topGO():
    try:
        topGO = importr("topGO")
    except:
        print ("It looks like topGO is not installed. Trying to install topGO via"
               "Bioconductor...")
        try:
            R.source("http://bioconductor.org/biocLite.R")
            R.biocLite("topGO")
            topGO = importr("topGO")
        except:
            print "Problem installing topGO from Bioconductor!"
            print ("Please install manually from: "
                   "http://www.bioconductor.org/packages/2.13/bioc/html/topGO.html")
    return topGO


def adjust_pvalues(pvalueHash):
    importr("qvalue")
    terms, pvalues = [], []
    for k, v in pvalueHash.iteritems():
        terms.append(k)
        pvalues.append(v)
    padjusted = R["p.adjust"](pvalues, method="fdr")
    return terms, pvalues, padjusted


def go_enrichment(gene_association, subset, algo, minnode):
    significant = []
    geneID2GO = R.readMappings(file=gene_association)
    refset = R.names(geneID2GO)
    testset = R.scan(file=subset, what=R.character())
    # Double %% escapes % in strings1
    genes_of_interest = R("factor(as.integer(%s %%in%% %s))" % (refset.r_repr(), testset.r_repr()))
    # Use the setNames function instead of this R code,
    # names(genes_of_interest) <- refset
    # because Python cannot assign to values

    genes_of_interest = R.setNames(genes_of_interest, refset)
    for o in ["MF", "BP", "CC"]:
        GOdata = R.new("topGOdata",
                       ontology=o,
                       allGenes=genes_of_interest,
                       annot=R["annFUN.gene2GO"],
                       gene2GO=geneID2GO,
                       nodeSize=minnode)
        pvalueHash = R.score(R.runTest(GOdata, algorithm=algo, statistic="fisher"))
        terms, pvalues, padjusted = adjust_pvalues(pvalueHash)
        for i, t in enumerate(terms):
            if pvalues[i] < 0.05:
                significant.append([t, str(pvalues[i]), str(padjusted[i])])
    return significant


def main():
    topGO = init_topGO()
    significant = go_enrichment(args.gene_association,
                                args.test_ids,
                                args.algorithm,
                                args.min)
    name = os.path.splitext(args.test_ids)[0]+".ora"
    with open(name, "wb") as handle:
        handle.write("#GOTERM\tP-Value\tFDR-adjusted\n")
        for s in significant:
            line = "\t".join(s)+"\n"
            handle.write(line)
if __name__ == "__main__":
    main()