import sys
from collections import defaultdict as d
from optparse import OptionParser, OptionGroup

# Author: Martin Kapun

#########################################################   HELP   #########################################################################
usage = "python %prog --input file --output file "
parser = OptionParser(usage=usage)
group = OptionGroup(parser, '< put description here >')

#########################################################   CODE   #########################################################################

parser.add_option("--karyo", dest="KA", help="Karyotype file")
parser.add_option("--input", dest="IN", help="Input file")
parser.add_option("--names", dest="NA", help="List of sample names")

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


def sync2string(x):
    ''' convert sync format to string of nucleotides  where x is a string in sync format '''
    string = ""
    alleles = ["A", "T", "C", "G"]
    ah = dict(list(zip(alleles, list(map(int, x.split(":")[:4])))))
    for k, v in list(ah.items()):
        string += v*k

    return string


names = options.NA.split(",")
snph = d(lambda: d(int))
kary = d(tuple)
for l in load_data(options.KA):
    C, P, I, S = l.split()
    kary[C+"_"+P] = (S, I)

for l in load_data(options.IN):
    C, P, R = l.split()[:3]
    Pops = l.split()[3:]
    if C+"_"+P not in kary:
        continue
    for p in range(len(names)):
        snph[names[p]]["S"] += sync2string(Pops[p]).count(kary[C+"_"+P][0])
        snph[names[p]]["I"] += sync2string(Pops[p]).count(kary[C+"_"+P][1])
print("Line\tStandard\tInverted\tkaryotype")
for k, v in sorted(snph.items()):
    if "Y" in k:
        K = "Y"
    elif v["S"]/float(sum(v.values())) > 0.9:
        K = "S"
    elif v["I"]/float(sum(v.values())) > 0.9:
        K = "I"
    else:
        K = "N"
    print(k+"\t"+str(v["S"]/float(sum(v.values()))) +
          "\t"+str(v["I"]/float(sum(v.values())))+"\t"+K)
