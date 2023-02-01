import sys
from rpy2.robjects import r
import rpy2.robjects as robjects
from optparse import OptionParser, OptionGroup
import collections
import random
import gzip
# Author: Martin Kapun

#########################################################   HELP   #########################################################################
usage = "\npython %prog --input output_file.consensus --ind 2,3,6,10 --subsample 500 --chromosome 2L --output output_2L"
parser = OptionParser(usage=usage)
group = OptionGroup(parser,
                    """
""")
#########################################################   CODE   #########################################################################

parser.add_option("--input", dest="input", help="Input consensus file")
parser.add_option("--output", dest="o", help="output file")
parser.add_option("--Inv", dest="Inv",
                  help="List of Ind with Inversion")
parser.add_option("--Std", dest="Std",
                  help="List of Ind with Standard arrangement")
parser.add_option("--ind", dest="i", help="indivduals used for the analysis")
parser.add_option("--chromosome", dest="c",
                  help="chromosome used for the analysis")
parser.add_option("--region", dest="e",
                  help="chromosome used for the analysis", default="NA")
parser.add_option("--N-cutoff", dest="u",
                  help="Not more N's than x percent, e.g. 0.1")
parser.add_option("--min-allele", dest="min",
                  help="Not more N's than x percent, e.g. 0.1")
parser.add_option("--color", dest="col",
                  help="Not more N's than x percent, e.g. 0.1", default="yellow,red,blue")

parser.add_option_group(group)
(options, args) = parser.parse_args()


def isolate(x, ind, inv, std):
    ''' get nucleotides for certain individuals'''
    nuc = ""
    for i in inv:
        nuc += x[ind.index(i)]
    for i in std:
        nuc += x[ind.index(i)]
    return nuc


def load_data(x):
    ''' import data either from a gzipped or or uncrompessed file or from STDIN'''
    import gzip
    if x == "-":
        y = sys.stdin
    elif x.endswith(".gz"):
        y = gzip.open(x, "rt", encoding="latin-1")
    else:
        y = open(x, "r", encoding="latin-1")
    return y


def Rsquared(pa, pb, pab):
    ''' Hill Robertson 1968'''
    D = pab-pa*pb
    if (pa*(1-pa)*pb*(1-pb)) == 0:
        return "NA"
    r2 = (D**2)/(pa*(1-pa)*pb*(1-pb))
    return r2


def AfHf(x, y):
    ''' calculate allele and haplotype frequencies for two biallelic loci x and y '''
    a = x.replace("N", "")[0]
    b = y.replace("N", "")[0]
    # remove N's and only retain individuals without ambiguous allels at both loci
    X, Y = ["".join(v) for v in zip(*[w for w in zip(x, y) if not "N" in w])]
    ab = ["".join(w) for w in zip(X, Y)]
    pa = X.count(a)/float(len(X))
    pb = Y.count(b)/float(len(Y))
    pab = ab.count(a+b)/float(len(ab))
    return pa, pb, pab


def ReadList(x):
    INDList = []
    for l in load_data(x):
        INDList.append(l.rstrip())
    return INDList


r.assign("col", robjects.vectors.StrVector(options.col.split(",")))
r('colormix=colorRampPalette(col)')
r('COL=colormix(100)')

individual = list(map(int, options.i.split(",")))
chromosome = options.c
cutoff = float(options.u)

CHR = {"X": 23542271, "2L": 23513712, "2R": 25286936,
       "3L": 28110227, "3R": 32079331, "4": 1348131}
ylim = CHR[chromosome]

r('pdf("'+options.o+'.pdf",width=20,height=4)')
r('plot(NULL, xlim=c(0,'+str(ylim/1000000)+'), ylim=c(0,1), ylab="", xlab="'+chromosome+'")')

INV = ReadList(options.Inv)
STD = ReadList(options.Std)

InvGeno = "1"*len(INV)+"0"*len(options.STD)

count = 1

for l in load_data(options.input):
    if count % 100000 == 0:
        print(count, "positions processed")
    count += 1
    C, P = l.split()[:2]
    D = l.rstrip().split()[-1]

    # exclude SNPs from other chromosomes
    if C != chromosome:
        continue

    if options.e != "NA":
        if int(P) < region[0] or int(P) > region[1]:
            continue

    DI = isolate(D, individual, INV, STD)
    # print l,DI

    # test if less N's than expected
    if DI.count("N")/float(len(DI)) > cutoff:
        continue

    # only use positions with more than 1 allele
    if len(set(DI.replace("N", ""))) == 1:
        continue

    a1 = DI.replace("N", "")[0]
    freq = DI.replace("N", "").count(a1)/float(len(DI.replace("N", "")))
    if freq < float(options.min) or freq > 1-float(options.min):
        continue

    # calculate Allele and haplotype frequencies
    pa, pb, pab = AfHf(DI, InvGeno)

    # calculate r-squared
    rsq = Rsquared(pa, pb, pab)
    if rsq == "NA":
        continue
    r('abline(v='+str(float(P)/1000000) +
      ',col=COL[as.integer('+str(rsq)+'*100)],lwd=0.5)')
r('rect(16.432209,0,24.744010,1,lwd=3,lty=2)')
r('dev.off()')
r('library(corrplot)')
r('pdf("'+options.o+'_legend.pdf",width=5,height=5)')
r('plot.new()')
r('colorlegend(COL, labels=seq(0,1,0.1))')
r('dev.off()')
