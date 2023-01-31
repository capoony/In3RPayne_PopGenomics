import sys
import gzip
from collections import defaultdict as d
from optparse import OptionParser, OptionGroup

# Author: Martin Kapun

#########################################################   HELP   #########################################################################
usage = "python %prog --input file --output file "
parser = OptionParser(usage=usage)
group = OptionGroup(parser, '< put description here >')

#########################################################   CODE   #########################################################################

parser.add_option("--input", dest="IN", help="Input file")
parser.add_option("--names", dest="NA", help="Output file")
parser.add_option("--type", dest="TY", help="pi or TD")

(options, args) = parser.parse_args()
parser.add_option_group(group)

names = options.NA.split(",")
D = options.IN.split(",")
posh = d(lambda: d(lambda: d(list)))


for i in range(len(D)):
    data = gzip.open(D[i], "r")
    for l in data:
        if l.startswith("CHROM"):
            continue
        a = l.rstrip().split()
        posh[i][a[0]][float(a[1])] = a[-1]


def ValList(x, y):
    '''find listitem closest to value'''
    Val = min(x, key=lambda X: abs(X-y))
    return Val


print("Name\tOrigin\tInvStatus\tlocation\tC\tP\t"+options.TY)

for S, v1 in sorted(posh.items()):
    for C, v2 in sorted(v1.items()):
        for P, V in sorted(v2.items()):
            if V == "NA":
                continue
            if "Inv" in names[S]:
                Status = "Inv"
            else:
                Status = "Std"
            if C == "3R" and P > 14432209 and P < 26744010:
                Type = "inside"
            else:
                Type = "outside"
            print(names[S]+"\t"+names[S].split("-")[0]+"\t" +
                  Status+"\t"+Type+"\t"+C+"\t"+str(P)+"\t"+V)
