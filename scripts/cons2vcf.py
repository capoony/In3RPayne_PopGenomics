import sys
from optparse import OptionParser, OptionGroup
from collections import defaultdict as d
import random
import gzip
# Author: Martin Kapun

#########################################################   HELP   #########################################################################
usage = "\npython %prog --input output_file.consensus --ind 2,3,6,10 --subsample 500 --chromosome 2L --output output_2L"
parser = OptionParser(usage=usage)
group = OptionGroup(parser,
                    """
H E L P:
____________

The purpose of this script is to extract a subsample (e.g. --subsample 0.5) or a defined number of polymorphic SNPs (e.g. --subsample 500) in the Individuals defined with --pop from a consensus file produced with extract_consensus.py (--input). For this SNP-dataset the script will calculate all pairwise r^2 values (according to Hill and Robertson 1968). There will be two output files: The first is a tab delimited file containing a tabular representation of the distance matrix, i.e. the columns consist of Chromosome, Position of SNP1, position of SNP2, r^2, alleles of SNP1 and alleles of SNP2.
The second file is a visual representation of the distance matrix using the LDheatmap R package, which needs to be installed before. This figure will show the physical genomic positions of SNPs used r^2 values highlighted in colors and for chromosomes with cosmopolitan inversion the breakpoints and the name of the inversion. Due to memory constraints, this script can only process one chromosmome at a time, which needs to be defined with (--chromosome).
""")
#########################################################   CODE   #########################################################################

parser.add_option("--input", dest="input", help="Input consensus file")
parser.add_option("--output", dest="o", help="output file")
parser.add_option("--subsample", dest="r",
                  help="subsample certain percentage, e.g 0.1, or a defined number of SNPs, e.g. 500", default=1.0)
parser.add_option("--ind", dest="i", help="indivduals used for the analysis")
parser.add_option("--chromosome", dest="c",
                  help="chromosome used for the analysis", default="NA")
parser.add_option("--region", dest="e",
                  help="chromosome used for the analysis", default="NA")
parser.add_option("--N-cutoff", dest="u",
                  help="Not more N's than x percent, e.g. 0.1")
parser.add_option("--names", dest="name",
                  help="Not more N's than x percent, e.g. 0.1", default="yellow,red,blue")


parser.add_option_group(group)
(options, args) = parser.parse_args()


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


def isolate(x, ind):
    ''' get nucleotides for certain individuals'''
    nuc = ""
    for i in ind:
        nuc += x[i]
    return nuc


individual = map(int, options.i.split(","))
chromosome = options.c
cutoff = float(options.u)
# idhash=dict(zip(individual,positions))

if options.e != "NA":
    region = map(int, options.e.split(":"))

fullsnplist = []
count = 1

names = options.name.split(",")

indnames = []


for i in range(len(individual)):
    indnames.append(names[i])

for l in load_data(options.input):
    if count % 1000000 == 0:
        print count, "positions processed"
    count += 1
    C, P, D = l.rstrip().split()

    # exclude SNPs from other chromosomes
    if chromosome != "NA" and C != chromosome:
        print "wrong chromosome"
        continue

    if options.e != "NA":
        if int(P) < region[0] or int(P) > region[1]:
            continue

    DI = isolate(D, individual)
    # print l,DI

    # test if less N's than expected
    if DI.count("N")/float(len(DI)) > cutoff:
        # print "N's"
        continue

    # only use positions with more than 1 allele
    if len(set(DI.replace("N", ""))) == 1 or len(set(DI.replace("N", ""))) > 2:
        # print "no allele"
        continue

    a, b = set(DI.replace("N", ""))

    DI = DI.replace("N", ".")
    R = a
    A = b
    DI = DI.replace(a, "0")
    DI = DI.replace(b, "1")

    # use proportional subset for calculations. set to 1 to use all
    if float(options.r) < 1:
        if random.random() <= float(options.r):
            fullsnplist.append(C+"\t"+P+"\t"+DI+"\t"+R+"\t"+A)
    # use a defined number of SNPs.
    else:

        fullsnplist.append(C+"\t"+P+"\t"+DI+"\t"+R+"\t"+A)
newlist = fullsnplist
if float(options.r) > 1:
    newdict = random.sample(range(len(fullsnplist)), int(options.r))
    newfullsnps = []
    for i in sorted(newdict):
        newfullsnps.append(fullsnplist[i])
    newlist = newfullsnps

VCF = gzip.open(options.o+".vcf.gz", "w")
VCF.write('''##fileformat=VCFv4.0
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
''')
VCF.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" +
          "\t".join(indnames)+"\n")
pos = []
for l in newlist:
    GT = []
    C, P, D, R, A = l.split()
    sl = C+"\t"+P+"\t.\t"+R+"\t"+A+"\t.\tPASS\t.\tGT\t"
    for i in range(len(indnames)):
        if D[i] == "1":
            GT.append(D[i]+"/"+D[i])
        elif D[i] == "0":
            GT.append(D[i]+"/"+D[i])
        else:
            GT.append("./.")
    VCF.write(sl+"\t".join(GT)+"\n")
