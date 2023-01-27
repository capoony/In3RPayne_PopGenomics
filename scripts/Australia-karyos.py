import sys
from collections import defaultdict as d


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


def sync2string(x):
    ''' convert sync format to string of nucleotides  where x is a string in sync format '''
    string = ""
    alleles = ["A", "T", "C", "G"]
    ah = dict(zip(alleles, map(int, x.split(":")[:4])))
    for k, v in ah.items():
        string += v*k

    return string


names = sys.argv[3].split(",")
snph = d(lambda: d(int))
kary = d(tuple)
for l in load_data(sys.argv[1]):
    C, P, I, S = l.split()
    kary[C+"_"+P] = (S, I)

for l in load_data(sys.argv[2]):
    C, P, R = l.split()[:3]
    Pops = l.split()[3:]
    if C+"_"+P not in kary:
        continue
    for p in range(len(names)):
        snph[names[p]]["S"] += sync2string(Pops[p]).count(kary[C+"_"+P][0])
        snph[names[p]]["I"] += sync2string(Pops[p]).count(kary[C+"_"+P][1])
print "Line\tStandard\tInverted\tkaryotype"
for k, v in sorted(snph.items()):
    if "Y" in k:
        K = "Y"
    elif v["S"]/float(sum(v.values())) > 0.9:
        K = "S"
    elif v["I"]/float(sum(v.values())) > 0.9:
        K = "I"
    else:
        K = "N"
    print k+"\t"+str(v["S"]/float(sum(v.values())))+"\t"+str(v["I"]/float(sum(v.values())))+"\t"+K
