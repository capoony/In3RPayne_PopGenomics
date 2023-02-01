import sys
import os
from collections import defaultdict as d


def median(x):
    mid = int(len(x)/2)
    sort = sorted(x)
    if len(x) == 0:
        return None
    if len(x) % 2 == 0:
        lower = sort[mid-1]
        upper = sort[mid]
        return (float(lower)+float(upper))/2.0
    else:
        return sort[mid]


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


G = d(str)
GTF = load_data(sys.argv[3])
for l in GTF:
    a = l.rstrip().split("\t")
    if a[2] != "gene":
        continue
    I, F = parse_gtf(a[8])
    G[I] = F

GeneC = d(lambda: d(int))
GeneL = d(lambda: d(list))

Dir = os.walk(sys.argv[1])

for f in Dir:
    Root, B, Files = f
    if len(Files) < 3:
        continue
    Name = Root.split("/")[-1]
    for l in open(Root+"/abundance.tsv", "r"):
        if l.startswith("target_id"):
            continue
        target_id, length, eff_length, est_counts, tpm = l.rstrip().split()
        gene = "-".join(target_id.split("-")[:-1])
        GeneL[gene][target_id] = float(eff_length)
        GeneC[Name][gene] += float(est_counts)

out1 = open(sys.argv[2]+"_counts.txt", "w")
out2 = open(sys.argv[2]+"_factors.txt", "w")

out1.write("Gene\tLength\t" +
           "\t".join([":".join(x.split("_")) for x in sorted(GeneC.keys())])+"\n")

for k, L in sorted(GeneL.items()):
    counts = []
    for n, v in sorted(GeneC.items()):
        counts.append(int(round(v[k], 0)))
    # out1.write(k+"\t"+str(median(L.values()))+"\t"+"\t".join(map(str,counts))+"\n")
    out1.write(G[k]+"\t"+str(median(list(L.values()))) +
               "\t"+"\t".join(map(str, counts))+"\n")


out2.write(
    "Treat=factor(c("+",".join(["'"+x.split("_")[0]+"'" for x in sorted(GeneC.keys())])+"))\n")
out2.write(
    "Temp=factor(c("+",".join(["'"+x.split("_")[1]+"'" for x in sorted(GeneC.keys())])+"))\n")
out2.write(
    "Rep=factor(c("+",".join(["'"+x.split("_")[2]+"'" for x in sorted(GeneC.keys())])+"))\n")
