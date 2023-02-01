import sys
from collections import defaultdict as d

def load_data(x):
	''' import data either from a gzipped or or uncrompessed file or from STDIN'''
	import gzip			
	if x=="-":
		y=sys.stdin
	elif x.endswith(".gz"):
		y=gzip.open(x,"r")
	else: 
		y=open(x,"r")
	return y

def parse_gtf(x):
    '''parse a gtf 9th column'''
    items=x.replace("\"","").split("; ")
    FBgn,ID="NA","NA"
    for pair in items:
        ident,symbol=pair.lstrip().split(" ")
        if ident=="gene_id":
            FBgn=symbol
        elif ident=="gene_symbol":
            ID=symbol[:-1]
    return ID,FBgn

Inp=load_data(sys.argv[1])
GTF=load_data(sys.argv[2])

head="GeneID\tFBgn\tChrom\tStart\tEnd\tLength\tOrientation"
unknown="\tNA\tNA\tNA\tNA\tNA\tNA\t"
Gene=d(str)
for l in GTF:
    a=l.rstrip().split("\t")
    if a[2]!="gene":
        continue
    C=a[0]
    S=str(min([int(a[3]),int(a[4])]))
    E=str(max([int(a[3]),int(a[4])]))
    L=str(abs(int(a[3])-int(a[4])))
    O=a[6]
    G,F=parse_gtf(a[8])
    Gene[F]="\t".join([G,F,C,S,E,L,O])


for l in Inp:
    if l.startswith("Gene"):
        print head+"\t"+"\t".join(l.rstrip().split()[1:])
        continue
    a=l.rstrip().split()
    if a[0] in Gene:
        print Gene[a[0]]+"\t"+"\t".join(a[1:])
    else:
        print a[0]+unknown+"\t".join(a[1:])
        



