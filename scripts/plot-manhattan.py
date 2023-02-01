import sys
import math
from optparse import OptionParser, OptionGroup
from collections import defaultdict as d
from rpy2.robjects import r
import rpy2.robjects as robjects

# Author: Martin Kapun

#########################################################   HELP   #########################################################################
usage = ""
parser = OptionParser(usage=usage)
group = OptionGroup(parser,
                    """

""")
#########################################################   CODE   #########################################################################

parser.add_option("--full", dest="full", help="")
parser.add_option("--color", dest="col", help="", default="blue")
parser.add_option("--ylim", dest="ylim", help="")
parser.add_option("--log", dest="log", help="", action="store_true")
parser.add_option("--CF", dest="cf", help="")
parser.add_option("--out", dest="out", help="")

(options, args) = parser.parse_args()
parser.add_option_group(group)


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


name = options.full.split("/")[-1].split(".")[0]
TH = int(options.cf)

# read data
posf, valf, posc, valc, mif, maf = d(list), d(list), d(list), d(list), 0.0, 0.0
for l in load_data(options.full):
    if l.startswith("GeneID"):
        continue
    a = l.split()
    if a[TH] in ["NA", "na", "nan"]:
        continue
    Chr = a[2]
    Pos = sum(map(int, a[3:5]))/2
    pval = float(a[-2])
    posf[Chr].append(float(Pos))

    if options.log:
        if float(a[TH]) == 0.0:
            val = 200
        else:
            val = -math.log10(float(a[TH]))
    else:
        val = float(a[TH])

    if val > maf:
        maf = val
    if val < mif:
        mif = val
    valf[Chr].append(val)
    if pval < 0.05:
        posc[Chr].append(float(Pos))
        valc[Chr].append(val)

code = ["X", "2L", "2R", "3L", "3R", "4"]
CHR = {"X": 23542271, "2L": 23513712, "2R": 25286936,
       "3L": 28110227, "3R": 32079331, "4": 1348131}

suml = 0
newCHR = d(list)
for x in code:
    if x not in posf:
        code.remove(x)
    suml += CHR[x]
S = 0.0
for x in code:
    E = S+(CHR[x]/float(suml))
    newCHR[x] = [S, E]
    S = E

r('pdf("'+options.out+"/"+name+'.pdf",width=30,height=8)')
r('par(mar = c(6, 0, 0, 0), oma = c(0, 6, 1, 0.5))')

for c in code:
    S, E = newCHR[c]
    r('par(fig=c('+str(S)+","+str(E)+',0,1), new=TRUE)')
    # print c,pih[c].values()
    r.assign('Position', robjects.vectors.FloatVector(posf[c]))
    r.assign('Pval', robjects.vectors.FloatVector(valf[c]))
    r('plot(Position/1000000,Pval,pch=16,ylim=c('+str(mif)+','+str(maf) +
      '),cex=1.5,col=rgb(0,0,0,0.1),cex.lab=2,axes = FALSE,xlab=expression(bold("'+c+'")))')
    if len(posc[c]) != 0:
        r.assign('P2', robjects.vectors.FloatVector(posc[c]))
        r.assign('Pval2', robjects.vectors.FloatVector(valc[c]))
        r('points(P2/1000000,Pval2,cex=1.5,pch=16,col="'+options.col+'")')
    r('box()')
    r('axis(1,at=seq(0,35,5),cex.axis=2)')
    if c == "X":
        r('axis(2,cex.axis=2)')
r('mtext(expression(bold("-log10(P)")), side = 2, outer = TRUE, cex = 2, line = 3.5)')
r('dev.off()')
