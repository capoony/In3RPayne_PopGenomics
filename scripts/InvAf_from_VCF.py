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
parser.add_option("--inv", dest="inv")
parser.add_option("--std", dest="std")

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


def most_frequent(List):
    return max(set(List), key=List.count)


InvID, StdID = [], []

for l in open(options.inv, "rt"):
    InvID.append(l.rstrip())

for l in open(options.std, "rt"):
    StdID.append(l.rstrip())

Invh = d(list)
for l in load_data(options.IN):
    if l.startswith("##"):
        continue
    elif l.startswith("#"):
        header = l.rstrip().split()[9:]
        # print(header)
        for i in range(len(header)):
            if header[i] in InvID:
                Invh["I"].append(i)
            elif header[i] in StdID:
                Invh["S"].append(i)
        continue
    a = l.rstrip().split()
    pops = a[9:]
    Alleles = [a[3]]
    Alleles.extend(a[4].split(","))
    Inva = d(list)
    for k, v in Invh.items():
        for i in v:
            if pops[i] == "./.":
                continue
            Inva[k].append(Alleles[int(pops[i][0])])
    if len(Inva["I"]) < 5 or len(Inva["S"]) < 5:
        continue
    MA = ""
    for A in Alleles:
        if Inva["I"].count(A)/len(Inva["I"]) > Inva["S"].count(A)/len(Inva["S"]):
            MA = A
            break
    if MA == "":
        MA = most_frequent(Inva["I"])
    print("\t".join([a[0], a[1], MA, "".join(Inva["I"]), "".join(Inva["S"])]))
