import sys
from collections import defaultdict as d
import math
from optparse import OptionParser, OptionGroup

# Author: Martin Kapun

#########################################################   HELP   #########################################################################
usage = "python %prog --input input.sync"
parser = OptionParser(usage=usage)
group = OptionGroup(parser,
                    """

""")
#########################################################   CODE   #########################################################################

parser.add_option("--input", dest="input", help="sync file with all SNPs")

(options, args) = parser.parse_args()
parser.add_option_group(group)


def FET(a, b, c, d):
    ''' Fisher's exact test'''
    import rpy2.robjects as robjects
    from rpy2.robjects import r
    r('counts<-matrix(c('+",".join(map(str, [a, b, c, d]))+'),ncol=2,byrow=F)')
    return str(list(r('fisher.test(counts)$p.value'))[0])


if "-F-" in options.input:
    stat = -1
else:
    stat = -2

CI, CO, NI, NO, CIBP, CIBO, NIBP, NIBO = 0, 0, 0, 0, 0, 0, 0, 0
name = options.input.split("/")[-1].split(".")[0]

for l in open(options.input, "r"):
    if l.startswith("GeneID"):
        continue
    a = l.split("\t")
    C = a[2]
    P = min(map(int, a[3:5]))
    pval = float(a[stat])
    if pval < 0.05:
        if C == "3R" and int(P) > 16232209 and int(P) < 24944010:
            CI += 1
            if (int(P) > 16232209 and int(P) < 18232209) or (int(P) > 22944010 and int(P) < 24944010):
                CIBP += 1
            else:
                CIBO += 1
        else:
            CO += 1
    elif C == "3R" and int(P) > 16232209 and int(P) < 24944010:
        NI += 1
        if (int(P) > 16232209 and int(P) < 18232209) or (int(P) > 22944010 and int(P) < 24944010):
            NIBP += 1
        else:
            NIBO += 1
    else:
        NO += 1


print(name+"\t"+"\t".join(map(str, [CI, CO, NI, NO]))+"\t"+FET(CI, CO, NI, NO) +
      "\t"+"\t".join(map(str, [CIBP, CIBO, NIBP, NIBO]))+"\t"+FET(CIBP, CIBO, NIBP, NIBO))
