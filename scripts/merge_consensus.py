import sys
from collections import defaultdict as d
from optparse import OptionParser, OptionGroup

# Author: Martin Kapun

#########################################################   HELP   #########################################################################
usage = "python %prog --input file --output file "
parser = OptionParser(usage=usage)
group = OptionGroup(parser, '< put description here >')

#########################################################   CODE   #########################################################################

parser.add_option("--SNPs", dest="SL", help="SNPs to be considered")
parser.add_option("--consensus", dest="CL",
                  help="List of Consensus files to merge")

(options, args) = parser.parse_args()
parser.add_option_group(group)


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


SL = d(lambda: d(str))
for l in load_data(options.SL):
    a = l.rstrip().split()
    SL[a[0]][int(a[1])]

CL = d(lambda: d(list))
for File in options.CL.split(","):
    for l in load_data(File):
        a = l.rstrip().split()
        if int(a[1]) in SL[a[0]]:
            CL[a[0]][int(a[1])].append(a[-1])

LEN = len(options.CL.split(","))

for C, v in sorted(CL.items()):
    for P, v1 in sorted(v.items()):
        if len(v1) != LEN:
            continue
        print(C, str(P), "".join(v1), sep="\t")
