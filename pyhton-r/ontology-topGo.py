__author__ = 'sk'
import csv
import collections
from rpy2.robjects.packages import importr
from rpy2.robjects import r as R, DataFrame


__topGo = None
__biomaRt = None
__mart = None
__mart_db = "hsapiens_gene_ensembl"


def init_topGO():
    global __topGo

    if __topGo is None:
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
                       "http://www.bioconductor.org/packages/release/bioc/html/topGO.html")
        __topGo = topGO

def init_biomaRt():

    global __biomaRt
    global __mart
    if __biomaRt is None:

        try:
            biomaRt = importr("biomaRt")
        except:
            print ("It looks like biomaRt is not installed. Trying to install biomaRt via"
                   "Bioconductor...")
            try:
                R.source("http://bioconductor.org/biocLite.R")
                R.biocLite("biomaRt")
                biomaRt = importr("biomaRt")
            except:
                print "Problem installing biomaRt from Bioconductor!"
                print ("Please install manually from: "
                       "http://www.bioconductor.org/packages/release/bioc/html/biomaRt.html")

        __biomaRt = biomaRt
        __mart = R.useMart(biomart = "ensembl", dataset = __mart_db)

def main():

    init_topGO()

    nodeSize = 10
    algo = "classic"  # choice from classic, elim, weight
    input_csv = "gene_data_hs.csv"
    input_go_map = "gene_association.goa_human_noHeader"

    with open(input_csv) as input_csv_file:
        genes_list = parse_input_csv(input_csv_file)

    with open(input_go_map) as input_go_map_file:
        genes2go = parse_go_map(input_go_map_file, genes_list)

    go_enrichment(genes2go, genes_list, algo, nodeSize)



def go_enrichment(genes2go, genes_list, algo, nodeSize):
    print "start go enrichment func"


    # genes_of_interest = collections.defaultdict(list)
    # for key in genes2go.keys():
    #     if key in genes_list:
    #         genes_of_interest[key] = 1
    #     else:
    #         genes_of_interest[key] = 0
    #print genes_of_interest['BIRC5']

    genes2go = R.readMappings(file='gene_anno_custom')
    testset = R.scan(file='subset.txt', what=R.character())
    testset = R.sample(testset, R("length(%s)/10" % testset.r_repr()))

    refset = R.names(genes2go)



    genes_of_interest = R("factor(as.integer(%s %%in%% %s))" % (refset.r_repr(), testset.r_repr()))
    genes_of_interest = R.setNames(genes_of_interest, refset)

    significant = collections.defaultdict(float)
    for o in ["MF", "BP", "CC"]:
    #o = "BP"
        GOdata = R.new("topGOdata",
                   ontology=o,
                   annot=R["annFUN.gene2GO"],
                   allGenes=genes_of_interest,
                   gene2GO=genes2go,
                   nodeSize=nodeSize
                   )

        pvalueHash = R.score(R.runTest(GOdata, algorithm=algo, statistic="fisher"))

        for i in range(len(pvalueHash)):
            if pvalueHash[i] < 0.05:
                significant[pvalueHash.names[i]] = pvalueHash[i]

        #print significant

    GO2Pval = collections.OrderedDict(sorted(significant.items(), key=lambda t: t[1]))

    #print GO2Pval

    i=0
    for k,v in GO2Pval.items():
        if i < 100 :
            print str(k)+":"+str(v)
            i += 1
        else:
            break


    # print "go enrichment object created"
    #
    #results = R.runTest(GOdata, algorithm=algo, statistic="fisher")
    #print results

    #
    # scores = R.score(results)
    #results_table = R.GenTable(GOdata,
    #                           classic = results,
    #                           orderBy = "classic",
    #                           ranksOf = "classicFisher",
    #                           topNodes = 10)
    #print results
    #
    #print results_table




def __parse_go_map(input_go_map_file):

    """
    parse a go map file from http://geneontology.org/
    :param input_go_map_file: file
    :return: genes2go dictionary
    """

    genes2go = collections.defaultdict(list)

    for line in input_go_map_file:
        if line[0] != "!":
            sline = line.split("\t")
            gene_name = sline[2].strip()
            go_id = sline[4].strip()
            genes2go[gene_name].append(go_id)

    return dict(genes2go)


def mk_genes2go_file(gene_annotation_file):

    with open(gene_annotation_file) as f:
        genes2go = __parse_go_map(f)
    f.close()

    with open('gene_annotation_genes2go', 'w') as f:
        firstLine = True
        for key in genes2go:
            # ensembl_id = convert_hgnc2ensembl(key)
            # if ensembl_id is not None:
            fristItem = True
            for item in genes2go.get(key):
                if fristItem:
                    if firstLine:
                        f.write(key+"\t")
                        firstLine = False
                    else:
                        f.write("\n"+key+"\t")

                    f.write(item)
                    fristItem = False
                else:
                    f.write(","+item)
    f.close()


def convert_ensembl2hgnc(ensembl_id):
    init_biomaRt()

    v = R.c(ensembl_id)
    res = R.getBM(attributes=R.c("hgnc_symbol"), filters="ensembl_gene_id", values=v, mart=__mart)

    try:
        return R.get("hgnc_symbol", res)[0]
    except:
        print 'Error convert_ensembl2hgnc: '+str(ensembl_id)+' not found in database'
        return None

def convert_hgnc2ensembl(hgnc_id):
    init_biomaRt()

    v = R.c(hgnc_id)
    res = R.getBM(attributes=R.c("ensembl_gene_id"), filters="hgnc_symbol", values=v, mart=__mart)

    try:
        return R.get("ensembl_gene_id", res)[0]
    except:
        print 'Error convert_hgnc2ensembl: '+str(hgnc_id)+' not found in database'
        return None



def test():
    #print convert_hgnc2ensembl("BIRC5")
    #print convert_hgnc2ensembl("CDC46")

    #print convert_ensembl2hgnc("ENSG00000118473")
    #print convert_ensembl2hgnc("dfdfdf")

    mk_genes2go_file("gene_association.goa_human")

if __name__ == "__main__":
    #main()
    test()
