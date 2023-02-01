import sys
from collections import defaultdict as d
from optparse import OptionParser, OptionGroup
import math
import gzip

#Author: Martin Kapun v. 1.04 - 29/03/2017

#########################################################   HELP   #########################################################################
usage="""
"""
parser = OptionParser(usage=usage)
helptext="""    

H E L P :
_________
"""
group=OptionGroup(parser,helptext)
#########################################################   parameters   #########################################################################

parser.add_option("--inputs", dest="i", help="an input file")
parser.add_option("--genes", dest="g", help="an input file")
parser.add_option("--stats", dest="s", help="an input file")
parser.add_option("--type", dest="t", help="type of stat")

parser.add_option_group(group)
(options, args) = parser.parse_args()


################################### functions ######################################

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


RawH=d(lambda:d(str))
GeneH=d(str)
Inp=options.i.split(",")
Genes=list(map(int,options.g.split(",")))
Stats=list(map(int,options.s.split(",")))
Typ=options.t.split(",")
for i in range(len(Inp)):
    for l in load_data(Inp[i]):
        if l.startswith("GeneID"):
            continue
        a=l.rstrip().split()
        #print a
        G=a[Genes[i]]
        if Typ[i]=="P":
            S=-math.log10(float(a[Stats[i]]))
        else:
            S=float(a[Stats[i]])
        RawH[i][G]=S

RawH2=d(lambda:d(lambda:d(str)))

for i,v1 in list(RawH.items()):
    c=1
    for k,v in sorted(list(v1.items()),key=lambda x:x[1], reverse=True):
        RawH2[i][k]=(v,c)
        c+=1

print("genes1\tgenes2\trank1\trank2\tP1\tP2")
for k,v in sorted(RawH2[0].items()):
    if k not in RawH2[1]:
        continue
    print(k+"\t"+k+"\t"+str(v[1])+"\t"+str(RawH2[1][k][1])+"\t"+str(v[0])+"\t"+str(RawH2[1][k][0]))
    