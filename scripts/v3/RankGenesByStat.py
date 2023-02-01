import sys
from collections import defaultdict as d
from optparse import OptionParser, OptionGroup
import math
import gzip

# Author: Martin Kapun v. 1.04 - 29/03/2017

#########################################################   HELP   #########################################################################
usage = """
"""
parser = OptionParser(usage=usage)
helptext = """    

H E L P :
_________
"""
group = OptionGroup(parser, helptext)
#########################################################   parameters   #########################################################################

parser.add_option("--input", dest="i", help="an input file")
parser.add_option("--gtf", dest="g", help="an input file")
parser.add_option("--perc", dest="p", help="top percent considered")

parser.add_option_group(group)
(options, args) = parser.parse_args()


################################### functions ######################################
def load_data(x):
    ''' import data either from a gzipped or or uncrompessed file or from STDIN'''
    import gzip
    if x == "-":
        y = sys.stdin
    elif x.endswith(".gz"):
        y = gzip.open(x, "r")
    else:
        y = open(x, "r")
    return y


def mean(x):
    return sum(x)/float(len(x))


def median(x):
    sortedLst = sorted(x)
    lstLen = len(x)
    index = (lstLen - 1) // 2

    if (lstLen % 2):
        return sortedLst[index]
    else:
        return (sortedLst[index] + sortedLst[index + 1])/2.0


def parse_gtf(x):
    '''parse a gtf 9th column'''
    items = x.replace("\"", "").split("; ")
    FBgn, ID = "NA", "NA"
    for pair in items:
        ident, symbol = pair.lstrip().split(" ")
        if ident == "gene_id":
            FBgn = symbol
        elif ident == "gene_symbol":
            ID = symbol[:-1]
    return ID, FBgn

# 3L      10004861        0.0277778       downstream_gene_variant FBgn0267537     CR45877 3994


Data = load_data(options.i)
GTF = load_data(options.g)
GeneH_raw = d(list)

GTH = d(str)
for l in GTF:
    a = l.rstrip().split("\t")
    if a[2] != "gene":
        continue
    C = a[0]
    S = str(min([int(a[3]), int(a[4])]))
    E = str(max([int(a[3]), int(a[4])]))
    L = str(abs(int(a[3])-int(a[4])))
    O = a[6]
    G, F = parse_gtf(a[8])
    GTH[F] = "\t".join([G, F, C, S, E, L, O])


for l in Data:
    Chr, Pos, Stat, Effect, Fbgn, Gene, Dist = l.rstrip("\n").split("\t")
    if Gene != "" and Stat != "nan":
        GeneH_raw[Fbgn].append(float(Stat))

GeneH = d(str)
CutH = d(str)
for Gene, FSTs in list(GeneH_raw.items()):
    cutoff = int(round(len(FSTs)*float(options.p), 0))
    if cutoff == 0:
        cutoff = 1
    F = sorted(FSTs, reverse=True)[:cutoff]
    GeneH[Gene] = median(F)
    CutH[Gene] = cutoff

for k, v in sorted(list(GeneH.items()), key=lambda x: (x[1]), reverse=True):
    print(GTH[k]+"\t"+"\t".join(map(str, [v, CutH[k]])))
