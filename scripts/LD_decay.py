import sys
from collections import defaultdict as d
from rpy2.robjects import r
import rpy2.robjects as robjects
import gzip
from optparse import OptionParser, OptionGroup
import random

# Author: Martin Kapun v. 1.04 - 29/03/2017

#########################################################   HELP   #########################################################################
usage = """python %prog \
      """
parser = OptionParser(usage=usage)
helptext = """    

H E L P :
_________

Calculate and plot pairwise LD among n combinations (defined with --combsample) based on m SNPs (defined with --SNPs) in a maximum distance of l (defined with --maxdist) within a genomic region (defined with --region) ignoring positions with a fraction k of N's (defined with --N-threshold). The analysis is separately carried out for two karyotypic groups (defined with --samples and --names)
"""
group = OptionGroup(parser, helptext)
#########################################################   parameters   #########################################################################

parser.add_option("--vcf", dest="v", help="")
parser.add_option("--samples", dest="i", help="")
parser.add_option("--names", dest="n", help="")
parser.add_option("--maxdist", dest="m", help="")
parser.add_option("--out", dest="o", help="")
parser.add_option("--region", dest="r", help="")
parser.add_option("--threshold", dest="t", help="")
parser.add_option("--N-threshold", dest="nt", help="")
parser.add_option("--SNPs", dest="sn", help="")
parser.add_option("--samplesize", dest="sa", help="", default="NA")
parser.add_option("--combsample", dest="c", help="", default="NA")


parser.add_option_group(group)
(options, args) = parser.parse_args()


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


def AfHf(x, y):
    ''' calculate allele and haplotype frequencies for two biallelic loci x and y '''
    if "".join(set(x)) == "N" or "".join(set(y)) == "N":
        return "NA"
    a = x.replace("N", "")[0]
    b = y.replace("N", "")[0]
    if len(zip(*[w for w in zip(x, y) if not "N" in w])) < 2:
        return "NA"
    # remove N's and only retain individuals without ambiguous allels at both loci
    X, Y = ["".join(v) for v in zip(*[w for w in zip(x, y) if not "N" in w])]
    ab = ["".join(w) for w in zip(X, Y)]
    pa = X.count(a)/float(len(X))
    pb = Y.count(b)/float(len(Y))
    pab = ab.count(a+b)/float(len(ab))
    return pa, pb, pab


def Rsquared(x):
    ''' Hill Robertson 1968'''
    pa, pb, pab = x
    D = pab-pa*pb
    if (pa*(1-pa)*pb*(1-pb)) == 0:
        return "NA"
    r2 = (D**2)/(pa*(1-pa)*pb*(1-pb))
    return r2


def parseVCF(AL, x):
    if x == ".":
        return "N"
    return AL[int(x)]


def testVCF(Pops, GT, th, nth, AL):
    Alleles = ""
    for i in Pops:
        Alleles += parseVCF(AL, GT[i][0])
    if "".join(set(Alleles)) == "N":
        return "NA"

    if Alleles.count("N")/float(len(Alleles)) > nth:
        return "NA"

    a = Alleles.replace("N", "")
    if len(set(a)) == 1:
        return "NA"
    A = "".join(set(Alleles))[0]
    Af = a.count(A)/float(len(a))
    if Af < th or Af > 1-th:
        return "NA"
    return Alleles


def decay(P, R, N, C):
    print P, R, N, C
    r('distance<-'+P)
    r('LD.data<-'+R)
    r('n<-'+N)
    r('HW.st<-c(C=0.00001)')

    r('HW.nonlinear<-nls(LD.data~((10+C*distance)/((2+C*distance)*(11+C*distance)))*(1+((3+C*distance)*(12+12*C*distance+(C*distance)^2))/(n*(2+C*distance)*(11+C*distance))),start=HW.st,control=nls.control(maxiter=1000,warnOnly=TRUE))')
    r('tt<-summary(HW.nonlinear)')
    r('new.C<-tt$parameters[1]')
    r('fpoints<-((10+new.C*distance)/((2+new.C*distance)*(11+new.C*distance)))*(1+((3+new.C*distance)*(12+12*new.C*distance+(new.C*distance)^2))/(n*(2+new.C*distance)*(11+new.C*distance)))')
    r('ld.df<-data.frame(distance,fpoints)')
    r('ld.df<-ld.df[order(ld.df$distance),]')
    r('lines(ld.df$distance,ld.df$fpoints,lty=2,lwd=6,col='+C+')')
    return str(list(r('tt$parameters[4]'))[0])


# get regional positions
Chr, Positions = options.r.split(":")
START, END = map(int, Positions.split("-"))

SsH = d(list)
SAMPLES = options.i.split(",")
NAMES = options.n.split(",")
for i in range(len(NAMES)):
    SsH[NAMES[i]] = [l.rstrip() for l in open(SAMPLES[i], "r")]

SH = d(list)
if options.sa != "NA":
    for k, v in SsH.items():
        SH[k] = random.sample(v, int(options.sa))
else:
    MinLen = min([len(x) for x in SsH.values()])
    for k, v in SsH.items():
        SH[k] = random.sample(v, MinLen)

PH = d(lambda: d(str))

c = 1
for l in load_data(options.v):
    if c % 100000 == 0:
        print c, "positions read"
    if l.startswith("##"):
        continue
    elif l.startswith("#"):
        header = l.rstrip().split()[9:]
        continue
    a = l.rstrip().split()
    if a[0] != Chr:
        continue
    if int(a[1]) < START or int(a[1]) > END:
        continue
    if len(a[4].split(",")) > 1:
        continue
    GT = dict(zip(header, a[9:]))
    AL = a[3]+a[4]
    c += 1
    for i in NAMES:
        AC = testVCF(SH[i], GT, float(options.t), float(options.nt), AL)
        if AC != "NA":
            PH[i][int(a[1])] = AC

print "VCF reading done"

RH = d(list)
Dist = d(list)

C1 = ["rgb(0,0.5,0.5)", "rgb(0.2,0.2,0.2)", "darkblue"]
C2 = ["rgb(0,1,1,0.01)", "rgb(0.2,0.2,0.2,0.01)", "rgb(0,0,1,0.01)"]

for i in NAMES:
    c = 1
    PP = random.sample(PH[i].keys(), int(options.sn))
    for P1 in PP:
        for P2 in PP:
            if P2 <= P1 or P2-P1 > int(options.m):
                continue
            if c % 10000 == 0:
                print c, "combinations calculated for ", i
            c += 1

            AFHF = AfHf(PH[i][P1], PH[i][P2])
            if AFHF == "NA":
                continue
            RS = Rsquared(AFHF)
            if RS == "NA":
                continue
            Dist[i].append(P2-P1)
            RH[i].append(RS)

Dist1 = d(list)
RH1 = d(list)
for i in NAMES:
    A = random.sample(range(len(RH[i])), int(options.c))
    Dist1[i] = [Dist[i][x] for x in A]
    RH1[i] = [RH[i][x] for x in A]

print "Rsquared calculations done"
r('pdf("'+options.o+'.pdf",width=12,height=6)')
for i in range(len(NAMES)):
    r.assign(NAMES[i]+"P", robjects.vectors.IntVector(Dist1[NAMES[i]]))
    r.assign(NAMES[i]+"R", robjects.vectors.FloatVector(RH1[NAMES[i]]))

    if i == 0:
        r('plot('+NAMES[i]+"P,"+NAMES[i]+"R,"+',ylim=c(0,1),xlim=c(0,' +
          str(options.m)+'),pch=16,xlab="Distance(bp)",ylab="r^2",col='+C2[i]+')')
    else:
        r('points('+NAMES[i]+"P,"+NAMES[i]+"R,"+',col='+C2[i]+',pch=16)')
for i in range(len(NAMES)):
    I = decay(NAMES[i]+"P", NAMES[i]+"R", str(len(SH[NAMES[i]])), C1[i])
out = open(options.o+".txt", "w")
out.write("Group\tGroupID\tdistance\tLD.data\tn\n")

for i in range(len(NAMES)):
    for j in range(len(Dist1[NAMES[i]])):
        out.write(str(i)+"\t"+NAMES[i]+"\t"+str(Dist1[NAMES[i]][j]) +
                  "\t"+str(RH1[NAMES[i]][j])+"\t"+str(len(SH[NAMES[i]]))+"\n")
# r('legend("topright", col=c("red","blue"),legend=c(paste("Inside: ","'+I+'",sep=""),paste("Outside: ","'+S+'",sep="")),lwd=2)')
r('dev.off()')
