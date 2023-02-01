import sys
from collections import defaultdict as d
from optparse import OptionParser, OptionGroup

# Author: Martin Kapun

#########################################################   HELP   #########################################################################
usage = "python %prog --input file --output file "
parser = OptionParser(usage=usage)
group = OptionGroup(parser, '< put description here >')

#########################################################   CODE   #########################################################################

parser.add_option("--input", dest="IN", help="Input file")
parser.add_option("--AFl", dest="AF", help="AFlist")
parser.add_option("--pops", dest="PO", help="PopList")
parser.add_option("--logical", dest="log",
                  help="logical parameter", action="store_true")
parser.add_option("--param", dest="param",
                  help="numerical parameter", default=1)

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


AFl = d(lambda: d(str))
POPS = []
for l in load_data(options.PO):
    POPS.append(l.rstrip)

for l in load_data(options.AF):
    a = l.rstrip().split()
    AFl[a[0]][a[1]] = a[2]

PopPos = []
for l in load_data(options.IN):
    if l.startswith("##"):
        continue
    elif l.startswith("#"):
        header = l.rstrip().split()[9:]
        # print(header)
        for i in range(len(header)):
            if header[i] in POPS:
                PopPos.append(i)
        continue

    a = l.rstrip().split()
    if a[1] not in AFl[a[0]]:
        continue
    pops = a[9:]
    Alleles = [a[3]]
    Inva = []
    Alleles.extend(a[4].split(","))
    for i in PopPos:
        if pops[i] == "./.":
            continue
        Inva.append(Alleles[int(pops[i][0])])
    AC = Inva.count(AFl[a[0]][a[1]])
    if len(Inva) == 0:
        continue
    print("\t".join([a[0], a[1], AFl[a[0]][a[1]], str(AC/len(Inva))]))
