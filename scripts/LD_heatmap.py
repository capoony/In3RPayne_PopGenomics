import sys
from rpy2.robjects import r
import rpy2.robjects as robjects
from optparse import OptionParser, OptionGroup
import collections
import random
import gzip
# Author: Martin Kapun

#########################################################   HELP   #########################################################################
usage = "\npython %prog --input output_file.consensus --ind 2,3,6,10 --subsample 500 --chromosome 2L --output output_2L"
parser = OptionParser(usage=usage)
group = OptionGroup(parser,
                    """
H E L P:
____________

The purpose of this script is to extract a subsample (e.g. --subsample 0.5) or a defined number of polymorphic SNPs (e.g. --subsample 500) in the Individuals defined with --pop from a consensus file produced with extract_consensus.py (--input). For this SNP-dataset the script will calculate all pairwise r^2 values (according to Hill and Robertson 1968). There will be two output files: The first is a tab delimited file containing a tabular representation of the distance matrix, i.e. the columns consist of Chromosome, Position of SNP1, position of SNP2, r^2, alleles of SNP1 and alleles of SNP2.
The second file is a visual representation of the distance matrix using the LDheatmap R package, which needs to be installed before. This figure will show the physical genomic positions of SNPs used r^2 values highlighted in colors and for chromosomes with cosmopolitan inversion the breakpoints and the name of the inversion. Due to memory constraints, this script can only process one chromosmome at a time, which needs to be defined with (--chromosome).
""")
#########################################################   CODE   #########################################################################

parser.add_option("--input", dest="input", help="Input consensus file")
parser.add_option("--output", dest="o", help="output file")
parser.add_option("--subsample", dest="r",
                  help="subsample certain percentage, e.g 0.1, or a defined number of SNPs, e.g. 500", default=1.0)
parser.add_option("--ind", dest="i", help="indivduals used for the analysis")
parser.add_option("--chromosome", dest="c",
                  help="chromosome used for the analysis")
parser.add_option("--region", dest="e",
                  help="chromosome used for the analysis", default="NA")
parser.add_option("--N-cutoff", dest="u",
                  help="Not more N's than x percent, e.g. 0.1")
parser.add_option("--min-allele", dest="min",
                  help="Not more N's than x percent, e.g. 0.1")
parser.add_option("--color", dest="col",
                  help="Not more N's than x percent, e.g. 0.1", default="yellow,red,blue")

parser.add_option_group(group)
(options, args) = parser.parse_args()


def load_data(x):
    ''' import data either from a gzipped or or uncrompessed file or from STDIN'''
    import gzip
    if x == "-":
        y = sys.stdin.decode('ASCII')
    elif x.endswith(".gz"):
        y = gzip.open(x, "rt", encoding="latin-1")
    else:
        y = open(x, "r", encoding="latin-1")
    return y


def isolate(x, ind):
    ''' get nucleotides for certain individuals'''
    nuc = ""
    for i in ind:
        nuc += x[i]
    return nuc


def AfHf(x, y):
    ''' calculate allele and haplotype frequencies for two biallelic loci x and y '''
    a = x.replace("N", "")[0]
    b = y.replace("N", "")[0]
    # remove N's and only retain individuals without ambiguous allels at both loci
    X, Y = ["".join(v) for v in zip(*[w for w in zip(x, y) if not "N" in w])]
    ab = ["".join(w) for w in zip(X, Y)]
    pa = X.count(a) / float(len(X))
    pb = Y.count(b) / float(len(Y))
    pab = ab.count(a + b) / float(len(ab))
    return pa, pb, pab


def Rsquared(pa, pb, pab):
    ''' Hill Robertson 1968'''
    D = pab - pa * pb
    if (pa * (1 - pa) * pb * (1 - pb)) == 0:
        return "NA"
    r2 = (D**2) / (pa * (1 - pa) * pb * (1 - pb))
    return r2


def inversionrect(chromo, steps, LO):
    ''' draw lines for inversionbreakpoints in the Heatmap!'''
    incoo = []
    if chromo == "2L":
        incoo.append(
            ((2225744 - LO) * steps, (13154180 - LO) * steps, "In(2L)t"))
    elif chromo == "3R":
        incoo.append(((20096867 - LO) * steps,
                      (32079331 - LO) * steps, "In(3R)C"))
        incoo.append(((21406917 - LO) * steps,
                      (29031297 - LO) * steps, "In(3R)Mo"))
        incoo.append(((16432209 - LO) * steps,
                      (24744010 - LO) * steps, "In(3R)Payne"))
    elif chromo == "2R":
        incoo.append(((15391154 - LO) * steps,
                      (20276334 - LO) * steps, "In(2R)Ns"))
    elif chromo == "3L":
        incoo.append(
            ((3173046 - LO) * steps, (16308841 - LO) * steps, "In(3L)P"))
    else:
        return "NA"
    return incoo


def execute(snplist, chr):
    out = gzip.open(options.o + ".dist.gz", "wt")
    count1, count2, count = 0, 0, 1
    rhash, alist, blist = [], [], []
    # loop through first SNP
    for i in snplist:
        if count % (float(options.r) / 10) == 0:
            print(count, "positions processed, now at position:", i.rstrip())
        count += 1
        pos = int(i.split()[1])
        # loop through second SNP
        for j in snplist[count1:]:
            pos1 = int(j.split()[1])
            # get allele frequencies
            if pos1 <= pos:
                continue
            nuc = i.split()[2]
            nuc1 = j.split()[2]
            pa, pb, pab = AfHf(nuc, nuc1)
            # print pa,pb,pab
            rsq = Rsquared(pa, pb, pab)
            if rsq == "NA":
                continue
            rhash.append(rsq)
            # store positions of SNP1
            alist.append(pos)
            # store positions of SNP2
            blist.append(pos1)
            # print inversion,rsq,pos,pos1
            out.write(
                "\t".join([str(x) for x in [chr, pos, pos1, rsq, nuc, nuc1]]) + "\n")
        # print count1,i.rstrip()
        count1 += 1

    # make x-axis labels based on the assumption that the labels go from 0-1 and now you need to scale the whole genome relative to these borders
    binlist, labellist = [], []
    # last SNP pos
    upper = max(alist + blist)
    # first SNP pos
    lower = min(alist + blist)
    # print upper,lower
    # make sure that there is at least ONE SNP
    if upper - lower != 0:
        # here caluclate the relative step of one basepair
        step = 1 / float(upper - lower)
        invcoo = inversionrect(chr, step, lower)
        # set step size to 2mb
        co = 2000000
        # this is the stepsize between the ticks
        stepsize = co * step
        # this is the start position
        start = (co - lower) * step
        # bin the steps in a list
        binlist.append(start)
        labellist.append(str(co / 1000000) + "mb")
        co += 2000000
        start += stepsize
        # append ticks until the step is larger than the position of the last SNP
        while (co < upper):
            labellist.append(str(co / 1000000) + "mb")
            binlist.append(start)
            co += 2000000
            start += stepsize

        # convert python to R
        cp = robjects.vectors.FloatVector(rhash)
        al = robjects.vectors.IntVector(alist)
        bl = robjects.vectors.IntVector(blist)
        bins = robjects.vectors.FloatVector(binlist)
        labels = robjects.vectors.StrVector(labellist)
        r.assign('values', cp)
        r.assign('al', al)
        r.assign('bl', bl)
        r.assign('bins', bins)
        r.assign('labels', labels)
        # open graphics device and load libraries
        r('library("LDheatmap")')
        r('png("' + options.o + "_" + chr + '.png",width=5000,height=5000)')
        # convert distance list to distance matrix
        r('x.names <- sort(unique(c(al, bl)))')
        r('x.dist <- matrix(0, length(x.names), length(x.names))')
        r('dimnames(x.dist) <- list(x.names, x.names)')
        r('x.ind <- rbind(cbind(match(al, x.names), match(bl, x.names)),cbind(match(bl, x.names), match(al, x.names)))')
        r('x.dist[x.ind] <- rep(values, 2)')
        # print r('t(arev(x.dist))')
        # make LDHeatmap grid object based on the r2 values. Use the topo color palette and put the Chromosome and Inversion in the title. Also print the number of SNPs used.Rotate the whole heatmap by 270 degrees.
        r('M<-LDheatmap(x.dist,sort(unique(c(al, bl)),decreasing=F),color=colormix(100),flip=T,geneMapLabelX=10000,title="")')
        # add an X-Axis above heatmap and use the labels generated above
        r('la<-LDheatmap.addGrob(M, grid.xaxis(label=labels,at=bins,main=F,gp=gpar(cex=10),name="axis"),height=0)')
        # add inversion breakpoints
        if invcoo != "NA":
            invcount = 0
            alphabet = ["a", "b", "c", "d", "e", "f", "g", "h"]
            for coord in invcoo:
                # print coord
                # add red line for the inversion boundaries
                r('l' + alphabet[invcount + 1] + '<-LDheatmap.addGrob(l' + alphabet[invcount] + ', grid.lines(x=c(' + str(coord[0]) + ',' + str(
                    coord[1]) + '),y=' + str(1.1 + (invcount / float(5))) + ',gp=gpar( lwd=8,col="red")),height=' + str(0.1 + (invcount / float(500))) + ')')
                # add label for the inversion
                r('l' + alphabet[invcount + 2] + '<-LDheatmap.addGrob(l' + alphabet[invcount + 1] + ', grid.text("' + str(coord[2]) + '",x=' + str(
                    coord[0]) + ',y=' + str(1.3 + (invcount / float(5))) + ',gp = gpar(cex = 5)),height=' + str(0.1 + (invcount / float(500))) + ')')
                invcount += 2
        # make everything white.
        r('grid.edit("axis", gp = gpar(col = "white"))')
        # and then just make the ticks and the labels black
        r('grid.edit(gPath("axis", "labels"), gp = gpar(col = "black"))')
        r('grid.edit(gPath("axis", "ticks"), gp = gpar(col = "black",lwd=4))')
        # resize the linewidth of the segments
        r('grid.edit(gPath("geneMap", "segments"), gp = gpar(lwd = 0.2))')
        # increae the size of the color key labels
        r('grid.edit("Key",  gp = gpar(cex = 8))')
        # increase the size of the title
        # r('grid.edit(gPath("heatMap", "title"), gp = gpar(cex=0))')

        r('dev.off()')


r.assign("col", robjects.vectors.StrVector(options.col.split(",")))
r('colormix=colorRampPalette(col)')

r('''
range01 <- function(x)(x-min(x))/diff(range(x))
cRamp <- function(x){
  cols <- colorRamp(color(100))(range01(x))
  apply(cols, 1, function(xt)rgb(xt[1], xt[2], xt[3], maxColorValue=255))
}
''')


individual = [int(x) for x in options.i.split(",")]
chromosome = options.c
cutoff = float(options.u)
# idhash=dict(zip(individual,positions))

if options.e != "NA":
    region = [int(x) for x in options.e.split(":")]

fullsnplist = []
count = 1

for l in load_data(options.input):
    if count % 10000000 == 0:
        print(count, "positions processed")
    count += 1
    a = l.rstrip().split()
    if len(a) == 3:
        C, P, D = l.split()
    else:
        C, P, Al, D = l.split()

    # exclude SNPs from other chromosomes
    if C != chromosome:
        continue

    if options.e != "NA":
        if int(P) < region[0] or int(P) > region[1]:
            continue

    DI = isolate(D, individual)
    # print l,DI

    # test if less N's than expected
    if DI.count("N") / float(len(DI)) > cutoff:
        continue

    # only use positions with more than 1 allele
    if len(set(DI.replace("N", ""))) == 1:
        continue

    a1 = DI.replace("N", "")[0]
    freq = DI.replace("N", "").count(a1) / float(len(DI.replace("N", "")))
    if freq < float(options.min) or freq > 1 - float(options.min):
        continue

    # use proportional subset for calculations. set to 1 to use all
    if float(options.r) < 1:
        if random.random() <= float(options.r):
            fullsnplist.append(C + "\t" + P + "\t" + DI)
    # use a defined number of SNPs.
    else:

        fullsnplist.append(C + "\t" + P + "\t" + DI)
if float(options.r) >= 1:
    newdict = random.sample(range(len(fullsnplist)), int(options.r))
    newfullsnps = []
    for i in sorted(newdict):
        newfullsnps.append(fullsnplist[i])

print("done")
# for each Chromosome loop through all possible combinations of SNPs
if float(options.r) < 1:
    execute(fullsnplist, chromosome)
else:
    execute(newfullsnps, chromosome)
