import sys
import collections
from optparse import OptionParser, OptionGroup

#Author: Martin Kapun 

#########################################################   HELP   #########################################################################
usage="python %prog --input input.cons --names 0,1,2,3,4,5,6,7,8 --exclude 3,4 > out.nex"
parser = OptionParser(usage=usage)
group=OptionGroup(parser,
"""
H E L P:
____________

A very simple script, which converts a consensus file to the NEXUS file format. Here, all names of the individuals in the consensus file need to be provided. Additionally, the option --exclude allows to specify individuals, which should be excluded from the output.
""") 
#########################################################   CODE   #########################################################################

parser.add_option("--input", dest="input", help="A consensus file")
parser.add_option("--names", dest="i", help="The names need to be separated by a comma and in the same order as in the consensus string, (e.g. 0,1,2,3,4,5,6,7,8,9)")
parser.add_option("--population",dest="e",help="individuals, which should be excluded from the NEXUS output")

parser.add_option_group(group)
(options, args) = parser.parse_args()

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


data=load_data(options.input)
names=options.i.split(",")
population=map(int,options.e.split(","))
datahash=collections.defaultdict(str)


for l in data:
    CHR,POS,cons=l.split()
    # if len(set(cons.replace("N","")))!=2:
    #     continue
    alleles=list(set(cons))
    CHR,POS,cons=l.split()
    nuc=[]
    #print cons,len(cons)
    for i in population:
        #print i,cons[i],names[i]
        datahash[names[i]]+=cons[i]

for k,v in datahash.items():
    if len(set(v))==1:
        del(datahash[k])

print"#NEXUS"

print "Begin data;"
print "Dimensions ntax="+str(len(datahash))+" nchar="+str(len(datahash.values()[0]))+";"
print "Format datatype=dna symbols=\"ACTG\" missing=N gap=-;\nMatrix"

for k,v in sorted(datahash.items()):
    print k+"\t"+v

print ";\nEnd;"
