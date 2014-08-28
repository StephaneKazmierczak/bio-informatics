__author__ = 'sk'
import collections
import ntpath
import argparse
import csv
from rpy2.robjects.packages import importr
from rpy2.robjects import r as R
from datetime import datetime

__qvalue = None
__topGo = None
__biomaRt = None
__mart = None
__mart_dataset = "hsapiens_gene_ensembl"


parser = argparse.ArgumentParser()
parser.add_argument("-mkGenes2go", type=str, default=None,
                    help="Create a genes2go mapping file, needs a gene association file from geneontology.org")

args = parser.parse_args()


def init_topGO():
    global __topGo

    if __topGo is None:
        try:
            print "Importing topGO ..."
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
            print "Importing biomaRt ..."
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
        __mart = R.useMart(biomart = "ensembl", dataset = __mart_dataset)


def init_qvalue():
    global __qvalue

    if __qvalue is None:
        try:
            print "Importing qvalue ..."
            qvalue = importr("qvalue")
        except:
            print ("It looks like qvalue is not installed. Trying to install qvalue via"
                   "Bioconductor...")
            try:
                R.source("http://bioconductor.org/biocLite.R")
                R.biocLite("qvalue")
                qvalue = importr("qvalue")
            except:
                print "Problem installing qvalue from Bioconductor!"
                print ("Please install manually from: "
                       "http://www.bioconductor.org/packages/release/bioc/html/qvalue.html")
        __qvalue = qvalue


def __add_adjusted_pvalues(dictionary):


    pvalues = []

    for key, item in dictionary.items():
        pvalues.append(dictionary[key]['pval'])

    fdr = R["p.adjust"](pvalues, method="fdr")
    by = R["p.adjust"](pvalues, method="BY")

    i = 0
    for key, item in dictionary.items():
        item['FDR'] = fdr[i]
        item['BY'] = by[i]
        i += 1
    #
    #
    #
    #
    #     dictionary[key]['BH'] = bh
    #     dictionary[key]['FDR'] = fdr
    #
    #     print "pvalue "+str(dictionary[key]['pval'])
    #     print "fdr "+str(R["p.adjust"](dictionary[key]['pval'], method="fdr"))





    return dictionary


def __add_GO_info(dictionary):

    whichTerms = R.c(dictionary.keys())
    qTerms = R.paste(R.paste("'", whichTerms, "'", sep=""), collapse=",")
    retVal = R.dbGetQuery(R.GO_dbconn(), R.paste("SELECT ontology, go_id, term, definition FROM go_term WHERE go_id IN (", qTerms, ");", sep=""))


    for iter in retVal.iter_row():
        go_id = iter.rx2('go_id')[0]
        ontology = iter.rx2('ontology')[0]
        term = iter.rx2('term')[0]
        definition = iter.rx2('definition')[0]

        dictionary[go_id]['ontology'] = ontology
        dictionary[go_id]['term'] = term
        dictionary[go_id]['definition'] = definition

    return dictionary


def __save_results(dictionary):

    with open("ontology_results_"+datetime.now().strftime('%Y-%m-%d_%H-%M-%S')+".csv", 'w') as f:

        line="GO:ID\tTerm\tPvalue\tFDR\tBenjamini-Yekutieli"
        print line
        f.write(line+'\n')

        for k, v in dictionary.items():
            line = str(k)+"\t"+str(v['pval'])+"\t"+str(v['FDR'])+"\t"+str(v['BY'])+"\t"+str(v['term'])
            print line
            f.write(line+'\n')

    f.close()


def go_enrichment(genes2go, genes_list, algo, nodeSize):

    init_topGO()
    #init_qvalue()

    genes2go_map = R.readMappings(file=genes2go)
    subset = genes_list
    refset = R.names(genes2go_map)

    genes_of_interest = R("factor(as.integer(%s %%in%% %s))" % (refset.r_repr(), subset.r_repr()))
    genes_of_interest = R.setNames(genes_of_interest, refset)

    score = collections.defaultdict(dict)

    for o in ["MF", "BP", "CC"]:
    #for o in ["MF"]:
        GOdata = R.new("topGOdata",
                   ontology=o,
                   annot=R["annFUN.gene2GO"],
                   allGenes=genes_of_interest,
                   gene2GO=genes2go_map,
                   nodeSize=nodeSize
                   )

        scoreR = R.score(R.runTest(GOdata, algorithm=algo, statistic="fisher"))

        for i in range(len(scoreR)):
            if scoreR[i] < 0.05:
                score[scoreR.names[i]] = {"pval": scoreR[i]}

    score = collections.OrderedDict(sorted(score.items(), key=lambda t: t[1]))
    score = __add_adjusted_pvalues(score)
    score = __add_GO_info(score)

    __save_results(score)
    #results = R.runTest(GOdata, algorithm=algo, statistic="fisher")
    #print results

    # scores = R.score(results)
    #results_table = R.GenTable(GOdata,
    #                           classic = results,
    #                           orderBy = "classic",
    #                           ranksOf = "classicFisher",
    #                           topNodes = 10)
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
    print "Starting genes2go mapping"

    with open(gene_annotation_file) as f:
        genes2go = __parse_go_map(f)
    f.close()

    file_name = ntpath.basename(gene_annotation_file)

    with open(file_name+'_genes2go', 'w') as f:
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
    print "Genes2go mapping file "+file_name+"_genes2go completed"


def convert_list_ensembl2hgnc(ensembl_id_list):

    init_biomaRt()

    v = R.c(ensembl_id_list)
    res = R.getBM(attributes=R.c("hgnc_symbol"), filters="ensembl_gene_id", values=v, mart=__mart)

    try:
        return R.get("hgnc_symbol", res)
    except:
        print 'Error convert_ensembl2hgnc: '+str(ensembl_id)+' not found in database'
        return None


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

    #print convert_list_ensembl2hgnc(['ENSG00000118473','ENSG00000089685','ENSG00000100297','ENSG00000134690'])
    #print convert_ensembl2hgnc("dfdfdf")

    #mk_genes2go_file("gene_association.goa_human")
    pass


def main():

    init_topGO()

    nodeSize = 10
    algo = "classic"  # choice from classic, elim, weight
    input_gene_subset = "gene_subset.txt"
    input_genes2go_map = "gene_association.goa_human_genes2go"

    gene_list = R.scan(file=input_gene_subset, what=R.character())
    gene_list = convert_list_ensembl2hgnc(gene_list)


    #gene_list = R.sample(gene_list, R("length(%s)/10" % gene_list.r_repr()))

    go_enrichment(input_genes2go_map, gene_list, algo, nodeSize)



if __name__ == "__main__":

    if args.mkGenes2go is not None:
        mk_genes2go_file(args.mkGenes2go)
    else:
        main()

        # d = {'a': robjects.IntVector((1,2,3)), 'b': robjects.IntVector((4,5,6))}
        # dataf = robjects.DataFrame(d)
        # print(dataf)
        #
        # for iter in dataf.iter_row():
        #     print iter
        #     print iter.rx2('a')[0]

    #test()
