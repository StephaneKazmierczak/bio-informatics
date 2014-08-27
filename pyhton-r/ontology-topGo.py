__author__ = 'sk'
import csv
import collections
from rpy2.robjects.packages import importr
from rpy2.robjects import r as R, DataFrame


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




# def pyDict2List(d):
#
#     for key in d:
#         for item in d.get(key):

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




def parse_go_map(input_go_map_file, genes_list):

    genes2go = collections.defaultdict(list)

    for line in input_go_map_file:
        split_line = line.split("\t")

        gene_name = split_line[2].strip()
        go_id = split_line[4].strip()

        #print "line ="+str(line)
        #print "go_id ="+str(go_id)
        #if gene_name in genes_list:
        genes2go[gene_name].append(go_id)

    # with open('gene_anno_custom', 'w') as f:
    #
    #     firstLine = True
    #     for key in genes2go:
    #         fristItem = True
    #         for item in genes2go.get(key):
    #             if item[:2] == "GO":
    #                 if fristItem:
    #                     if firstLine:
    #                         f.write(key+"\t")
    #                         firstLine = False
    #                     else:
    #                         f.write("\n"+key+"\t")
    #
    #                     f.write(item)
    #                     fristItem = False
    #                 else:
    #                     f.write(","+item)
    #
    # f.close()

    return genes2go


def create_genes2go_file(gene_annotation_file, ):
    pass



def parse_input_csv(input_csv_file):
    # # extract the list of gene from a dummy file colum 2
    reader = csv.reader(input_csv_file)
    reader.next()  # remove header
    all_genes = list()
    for row in reader:
        all_genes.append(row[1])
    return all_genes

if __name__ == "__main__":
    main()