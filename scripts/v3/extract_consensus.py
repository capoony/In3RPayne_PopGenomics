import math
import collections
from optparse import OptionParser, OptionGroup


# Author: Martin Kapun


#########################################################   HELP   #########################################################################
usage = "\npython %prog --input individuals.sync --min-coverage 20 --max-coverage 0.05 --min-count 20 --CI-mode individual --output output_file\n\nor:\npython %prog --input individuals.sync --min-coverage 20 --max-coverage 100,100,100,100,100,90+200,200,200,199,100 --min-count 20 --CI-mode individual --output output_file"
parser = OptionParser(usage=usage)
group = OptionGroup(parser,
                    """
H E L P:
____________

The purpose of this script is to extract the sire allele from F1 hybrids sequenced as indivduals. The input has to be a sync file with the sequences of the dam (mel36) in the first column, followed by all indivduals. Note that only the Autosomes X,2L,2R,3L,3R and 4 will be considered. In contrast to the individuals, Mel36 has been sequenced as a pool of 10 females, therefore the mean allele frequency is not necessarily 50%. Similarily to the scripts in PopPoolation this script uses cummulative counts of alleles across all synced indivduals/populations to test for minimum count (--min-count). Additionally you have to set a minimum coverage (--min-coverage) threshold, which will be applied to all populations/individuals. 
You have two options for the maximum coverage thresholds (--max-coverage). You can either set the threshold in percent, e.g. 0.05. Then the script will calculate the cutoff based on the top 5 percentile. In contrast to PoPoolation this script calculates these thresholds for each chromosome separately. This is necessary, as male indivduals should only have half the coverage on the X. Then, the script will create an additional output with the cutoff threshold as defined by the input, the IDs ad the actual coverage cutoffs. 
The latter information can be alternatively used as a direct input for the maximum coverage threshold to avoid re-calculation of the percentiles. (In the format: Ind1_2L,
At every position, the script 1) determines, whether the reference (mel36) is within the coverage thresholds and 2) the reference is polymoprhic. I.e. the minor allele is above the min-count threshold and the minor allele frequency (AF) is larger than the lower boundary of the 90% binomial confidence interval (BCI) for a 5% allele frequency at the given coverage.
If the refernce does not fullfill the coverage criteria, the allele from the reference genome is used. If the reference is polyorphic, all indivduals will be set to "N", as the sire allele cannot inferred in this case. 
If the reference is not polymorphic, the individuals are tested for alleles different from the reference: Here you have two choices: 1) to calculate 90% BCI for an allele frequency of 50% for each indivdual separately and test whether the allele frequency is within the range (--CI-mode individual) or 2) to pool all individuals carrying the alternative allele and calculate the 95% BCI for an allele frequency of 50% for the pool and use the alternative allele for every individual if the cummulative AF is within the BCI boundaries (--CI-mode all). Tests have shown that method produces false positives. I would therefore strongly suggest method 1 (which is the default).
The output looks a bit like a pileup file. I have introduced a quality string, which allows to reconstruct, why at certain indivduals an "N" was used rather than an allele. See a description of the strings below:

1: sire allele similar to reference
2: sire allele different from reference
c: coverage too low or too high
r: reference is polymorphic
/: alternative AF outside boundaries of BCI for the particular individual
a: more than two alleles

See a template output below:

2L      6921    G       G       NGGNNNGTGNG     c11ccc121c1     T/G
2L      6922    A       A       AAANNNAAANA     111ccc111c1     A
2L      16933   T       T       TTTNTNTTTNT     111/1c111c1     A/T
2L      19795   T       A/T     NNNNNNNNNNN     rrrrrrrrrrr     A/T
2L      20008   A       A       AGGNANGAANA     122c1c211c1     A/G


Col1: Chromosome
Col2: position
Col3: allele from the reference
Col4: allele(s) from sequenced mel36
Col5: Alleles for each individual; in the same order as in the sync file
Col6: Quality strings for each individual; in the same order as in the sync file
Col7: all alleles above the min-count threshold

In addition, two more outputs will be produced:
1)The file with the extension *.af contains the allele frequencies of individuals with non-reference alleles full-filling all criteria. After columns for chromsome and position, all individuals are listed in the same order as in the sync file. Missing values are indicated as "NA".
2) The file with the extension *_ref.af lists the polymorphic sites in the reference. It contains columns for Chromosome, position, reference genome allele, alleles detected and the frequency of the minor allele.

""")
#########################################################   CODE   #########################################################################

parser.add_option("--input", dest="input", help="A synchronized input file")
parser.add_option("--min-coverage", dest="c",
                  help="minimum coverage threshold")
parser.add_option("--min-count", dest="m", help="minimum count threshold")
parser.add_option("--output", dest="o", help="output file(s)")
parser.add_option("--max-coverage", dest="a", help="maximum coverage cutoff")
parser.add_option("--CI-mode", dest="ci",
                  help="choose between 'all' and 'individual', see Help for details, default = 'individual' ", default="individual")
parser.add_option_group(group)
(options, args) = parser.parse_args()


def CI(x, limit):
    ''' calulate exact binomial confidence interval'''
    if 0.0 in x:
        return 0.0, 0.0
    else:
        cilim = str(float(1-limit)/2)
        n = sum(x)
        X = x[0]
        p = X/float(n)
        df1l = 2*(n-X+1)
        df2l = 2*X

        Fl = float(
            r('qf('+cilim+','+str(df1l)+","+str(df2l)+',lower.tail=F)')[0])
        df1u = df2l+2
        df2u = df1l-2

        Fu = float(
            r('qf('+cilim+','+str(df1u)+","+str(df2u)+',lower.tail=F)')[0])
        # lower limit
        Ll = X/(X+(n-X+1)*Fl)
        # upper limit
        Lu = (X+1)*Fu/(n-X+(X+1)*Fu)
        return Ll, Lu


def run():
    for pop, alleles in sorted(popalleles.items()):
        # test if coverage of individual is above threshold, else append "N"
        if sum(alleles.values()) < int(options.c):
            hapallele.append("N")
            qual.append("c")
        else:
            hapallele.append(refal)
            if refal == oldref:
                qual.append("1")
            else:
                qual.append("2")


def CI_ind():
    for pop, alleles in sorted(popalleles.items()):
        # print pop,alleles
        # print pop
        # test if coverage of individual is above threshold, else append "N"
        if sum(alleles.values()) < int(options.c):
            # print pop,"c"
            hapallele.append("N")
            qual.append("c")
            freqal.append("NA")
        elif sum(alleles.values()) > covdict[pop][chrom]:
            # print pop,"c"
            hapallele.append("N")
            qual.append("C")
            freqal.append("NA")
        else:
            # test if alleles found in individual is fullfilling the minimum count criteria else delete the allele
            for k, val in list(alleles.items()):
                if val <= 2:
                    del alleles[k]
            if len(alleles) == 2:
                if len(refal) > 1:
                    hapallele.append("N")
                    qual.append("r")
                    freqal.append("NA")
                    continue
                if refal not in list(alleles.keys()):
                    # print pop,"a"
                    hapallele.append("N")
                    qual.append("a")
                    freqal.append("NA")
                    continue
                total = sum(alleles.values())
                # test if CI fits the data:
                if CI([total/2, total/2], 0.9)[0] <= min(alleles.values())/float(sum(alleles.values())):
                    for k in list(alleles.keys()):
                        if k != refal[0]:
                            # test if alternative allele occurs more than options.m in cummulative sample
                            if fullalleles[k] >= int(options.m):
                                # print pop,"2"
                                hapallele.append(k)
                                if k == oldref:
                                    qual.append("1")
                                else:
                                    qual.append("2")
                                freqal.append(
                                    alleles[k]/float(sum(alleles.values())))
                            else:
                                # append "m" in qual string if count too low
                                # print pop,"a"
                                hapallele.append("N")
                                qual.append("m")
                                freqal.append("NA")
                else:
                    # print pop,"2"
                    hapallele.append("N")
                    qual.append("/")
                    freqal.append("NA")
            elif len(alleles) == 1:
                # print pop,"1"
                hapallele.append(list(alleles.keys())[0])
                if list(alleles.keys())[0] == oldref:
                    qual.append("1")
                else:
                    qual.append("2")
                freqal.append("NA")
            else:
                # print pop,"a"
                hapallele.append("N")
                qual.append("a")
                freqal.append("NA")


############################### max coverage #################################################

# individual=["mel36",52,53,80,100,117,132,136,150,168,85,106]
# iupac={"R":["A","G"],"Y":["C","T"],"S":["G","C"],"W":["A","T"],"K":["G","T"],"M":["C","T"]}
synccode = ["A", "T", "C", "G"]
chromo = ["2L", "2R", "3L", "3R", "4", "X"]
popstotest = list(range(3, len(open(options.input).readline().split()[3:])+3))

if "/" not in options.a:
    out_cov = open(options.o+".cov", "w")
    cutoff = float(options.a)
    covdict_full = collections.defaultdict(
        lambda: collections.defaultdict(lambda: []))
    covdict = collections.defaultdict(
        lambda: collections.defaultdict(lambda: 0))
    count = 0
    for line in open(options.input, "r"):
        a = line.split()
        popstotest = list(range(3, len(a[3:])+3))
        chrom = a[0]
        for pop in popstotest:
            if a[pop] != "-":
                # print line[pop], line, pop
                covdict_full[pop][chrom].append(
                    sum(map(int, a[pop].split(":"))))

        count += 1
        if count % 1000000 == 0:
            print(count, "positions processed")
    string = []
    chromlist = []
    for pop, hash1 in sorted(covdict_full.items()):
        string_chrom = []
        chrom_ID = []
        for chrom, values in sorted(hash1.items()):
            # make a dictionary of all unique coverages to later detect the next smaller coverage
            covclass = collections.defaultdict(lambda: 0)
            count = 1
            for item in sorted(set(values)):
                covclass[item] = count
                count += 1
            # make a reversed dictionary with the indices as keys and the coverages as values
            revcovclass = dict(list(zip(list(covclass.values()), list(covclass.keys()))))
            # determine the class with the top 2% coverages by sorting all coverages from largest to smallest and then slice the list at the position cutoff*length_list (round down to the lower integer)
            covpos = sorted(values, reverse=True)[
                int(math.floor(len(values)*cutoff))]

            # append the next smaller coverage to a list
            if covclass[covpos] < 2:
                covdict[pop][chrom] = 0
            else:
                covdict[pop][chrom] = revcovclass[covclass[covpos]-1]
            string_chrom.append(covdict[pop][chrom])
            chrom_ID.append(str(pop)+"_"+chrom)
        string.append(string_chrom)
        chromlist.append(chrom_ID)
    # write the ID's and coverages to an output file. the last line can be used as the input for another run if you do not want to calucalte everything again
    out_cov.write(" coverage cutoff of top "+str(cutoff*100)+"%:\n")
    print(" coverage cutoff of top "+str(cutoff*100)+"%:")
    out_cov.write("+".join([",".join(x) for x in chromlist])+"\n")
    print("+".join([",".join(x) for x in chromlist]))
    out_cov.write("+".join([",".join(map(str, x)) for x in string])+"\n")
    print("+".join([",".join(map(str, x)) for x in string]))

    covdict_full = 0

else:
    chromo = ["2L", "2R", "3L", "3R", "4", "X"]
    covdict = collections.defaultdict(
        lambda: collections.defaultdict(lambda: 0))
    cov = open(options.a, "r")
    cov.readline()
    header = cov.readline().rstrip().split("+")
    coverage = cov.readline().rstrip()

    for i in range(len(coverage.split("+"))):
        head = zip(*[x.split("_") for x in header[i].split(",")])[1]
        covdict[i+3] = dict(list(zip(head,
                            list(map(int, coverage.split("+")[i].split(","))))))
# mincov
# print covdict
############################### determine haplotype #########################################

# out=open(options.o+".consensus","a")
out2 = open(options.o+".af", "a")
out3 = open(options.o+"_ref.af", "a")
count = 1
for line in open(options.input, "r"):
    a = line.split()
    if len(a) < 2:
        continue
    popstotest = list(range(3, len(a[3:])+3))
    chrom = a[0]
    oldref = a[2]
    if count % 1000000 == 0:
        print(count, "lines processed")
    count += 1
    fullalleles = collections.defaultdict(lambda: 0)
    popalleles = collections.defaultdict(
        lambda: collections.defaultdict(lambda: 0))

    # extract the alleles with a count > 0 cummulatively for all populations (fullalleles) and for each population seperatly (popalleles)
    for pop in popstotest:
        popalleles[pop]

        # go through all nucleotides per populations and test whether larger than one, if yes, keep!
        for i in range(len(synccode)):
            if int(a[pop].split(":")[i]) > 0:

                # store counts cummulatively for all populations
                fullalleles[synccode[i]] += int(a[pop].split(":")[i])

                # store counts for each population separately
                popalleles[pop][synccode[i]] += int(a[pop].split(":")[i])

    # if overall coverage is too low just print N's for all pops and continue
    if len(fullalleles) == 0:
        print("\t".join(a[:3])+"\t"+a[2]+"\t"+"N"*(len(popstotest)-1)+"\t"+"c"*(len(popstotest)-1)+"\t-")
        # out.flush()
        continue

    ############################# extract reference allele ###################################

    reference = popalleles[3]

    for k, val in list(reference.items()):
        if val <= 2:
            del reference[k]

    # test if reference is above minimum and below maximum coverage
    if sum(reference.values()) < int(options.c) or sum(reference.values()) > covdict[3][chrom]:
        refal = a[2]
        del popalleles[3]

    # print "N's"	and exit, if reference is ambiguous
    elif len(reference) > 1:
        total = sum(reference.values())
        out3.write("\t".join(a[:3])+"\t"+"/".join(list(reference.keys())) +
                   "\t"+str(min(reference.values())/float(total))+"\n")
        out3.flush()

        # test if minor allele is within the range of 90% interval of a minimum AF of 0.05 (because reference is a pool of 10 females). If yes, print "N's" and exit, if reference is ambiguous
        if CI([total*0.05, total*0.95], 0.9)[0] <= min(reference.values())/float(total):
            refal = list(reference.keys())
            del popalleles[3]

        # elses use major allele as reference
        else:
            tar = dict((value, key) for key, value in list(reference.items()))
            del popalleles[3]
            refal = tar[max(tar.keys())]
    # if (for whatever reason) there is no allele, take the old reference
    elif len(reference) == 0:
        refal = a[2]
        del popalleles[3]
    # else use the only allele as the new reference allele
    else:
        refal = list(reference.keys())[0]
        del popalleles[3]

    ################# determine the sire allele in the indivduals ###########################

    freqal = []  # list of allelesfreqs from the haplotypes
    hapallele = []  # list of alleles from the haplotypes
    qual = []  # list of description symbols for each individual
    if len(fullalleles) > 1:
        if options.ci == "individual":
            CI_ind()
        else:
            CI_all()
    else:
        run()

    print("\t".join(a[:3])+"\t"+"/".join(refal)+"\t"+"".join(hapallele)+"\t"+"".join(qual)+"\t"+"/".join(list(fullalleles.keys())))
    # out.flush()
    if len(freqal) != 0 and list(set(freqal)) != ["NA"]:
        out2.write("\t".join(a[:2])+"\t"+"\t".join(map(str, freqal))+"\n")
        out2.flush()
