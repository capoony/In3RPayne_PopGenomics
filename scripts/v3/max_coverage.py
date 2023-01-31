import math
import collections
from optparse import OptionParser, OptionGroup

# Author: Martin Kapun


#########################################################   HELP   #########################################################################
usage = "\npython %prog --input individuals.sync --max-coverage 0.05 --output output_file"
parser = OptionParser(usage=usage)
group = OptionGroup(parser,
                    """
H E L P:
____________

calculates the maximum coverage threshold for a given threshold for every chromosome (2L,2R,3L,3R,4,X) and every individual. Can be used as the mx-coverage input for bcf2consensus.py
""")
#########################################################   CODE   #########################################################################

parser.add_option("--input", dest="input", help="A synchronized input file")
parser.add_option("--max-coverage", dest="a", help="maximum coverage cutoff")
parser.add_option("--output", dest="o", help="output file(s)")
parser.add_option_group(group)
(options, args) = parser.parse_args()


############################### max coverage #################################################

# individual=["mel36",52,53,80,100,117,132,136,150,168,85,106]
# iupac={"R":["A","G"],"Y":["C","T"],"S":["G","C"],"W":["A","T"],"K":["G","T"],"M":["C","T"]}
synccode = ["A", "T", "C", "G"]
chromo = ["2L", "2R", "3L", "3R", "4", "X"]
popstotest = list(range(3, len(open(options.input).readline().split()[3:])+3))

out_cov = open(options.o, "w")
cutoff = float(options.a)
covdict_full = collections.defaultdict(
    lambda: collections.defaultdict(lambda: []))
covdict = collections.defaultdict(lambda: collections.defaultdict(lambda: 0))
count = 0
for line in open(options.input, "r"):
    a = line.split()
    chrom = a[0]
    for pop in popstotest:
        if a[pop] != "-":
            # print line[pop], line, pop
            covdict_full[pop][chrom].append(sum(map(int, a[pop].split(":"))))

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
out_cov.write("+".join([",".join(x) for x in chromlist])+"\n")
out_cov.write("+".join([",".join(map(str, x)) for x in string])+"\n")
