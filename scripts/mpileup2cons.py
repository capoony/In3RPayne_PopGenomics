import sys
import collections
import re
from optparse import OptionParser, OptionGroup


#Author: Martin Kapun

#########################################################   HELP   #########################################################################
usage="python %prog --mpileup simulans_data.mpileup --library-names Kib32,Tam10 --fasta simulans_genome.fa --coverage-threshold 10,100 --base-quality-threshold 25 --output outputfolder/updated_genomes"
parser = OptionParser(usage=usage)
helptext="""	

H E L P :
_________

Description:
This script identifies non-reference SNP alleles in a mpileup file (--mpileup) and replaces the nucleotides at the corresponding position in the output FASTA file.
Only the most frequent alternative allele will be used to replace the reference nucleotide based on a reference genome (--fasta).
The script is able to handle multiple libraries in the mpileup file (--library-names) and will produce an outputfile for each library (--output).
You need to define base-quality and coverage thresholds (--base-quality-threshold, --coverage-threshold). Positions not fullfilling these criteria will be excluded from the analyses. 
"""
group=OptionGroup(parser,helptext)
#########################################################   parameters   #########################################################################

parser.add_option("-m","--mpileup", dest="m", help="A mpileup file")
parser.add_option("-c","--coverage-threshold", dest="c", help="a list containing the minimum and the maximum coverage thresholds separated by a comma: e.g. 10,80")
parser.add_option("-b","--base-quality-threshold", dest="b", help="The Base-quality threshold for Qualities encoded in Sanger format (Illumina 1.8 format)")
parser.add_option("-u","--cutoff", dest="u", help="output-file: Do not add any extension (such as .fasta) as the program will do this automatically")

parser.add_option_group(group)
(options, args) = parser.parse_args()

###################################### functions #############################################


def keywithmaxvalue(d):
    ''' This function resturns the key for the maximum value in a dictionary'''
    newhash=collections.defaultdict(list)
    for k,v in d.items():
        newhash[v].append(k)
    return newhash[max(newhash.keys())]


def splitter(l, n):
    ''' This generator function returns equally sized cunks of an list'''
    #credit: Meric Lieberman, 2012
    i = 0
    chunk = l[:n]
    while chunk:
        yield chunk
        i += n
        chunk = l[i:i+n]

def extract_indel(l,sign):
    ''' This function returns an Indel from a sequence string in a pileup'''
    position = l.index(sign)
    numb =""
    i = 0
    while True:
        if l[position+1+i].isdigit():
            numb+= l[position+1+i]
            i+=1
        else:
            break
        
    seqlength = int(numb)
    sequence = l[position:position+i+1+seqlength]
    indel=sequence.replace(numb,"")

    return sequence,indel

def load_data(x):
	import gzip			
	if x=="-":
		y=sys.stdin
	elif x.endswith(".gz"):
		y=gzip.open(x,"r")
	else: 
		y=open(x,"r")
	return y


##################################### data #################################################

data=load_data(options.m)
threshold=int(options.b)
cutoff=float(options.u)

if len(map(int,options.c.split(",")))>2:
    covmin=map(int,options.c.split(","))[0]
    covmax=map(int,options.c.split(","))[1:]
else:
    covmin,covmax=map(int,options.c.split(","))

##################################### main code #################################################

# parse mpileup and store alternative alleles:
allelehash=collections.defaultdict(lambda:collections.defaultdict(lambda:collections.defaultdict(str)))
for line in data:
    #print line
    k = line[:-1].split('\t')
    chrom,position,refbase = k[:3]
    div = list(splitter(k,3))
    libraries=div[1:]
    
    cons=[]
    # loop through libraries    
    for j in range(len(libraries)):
        
        nuc = libraries[j][1]
        qualities = libraries[j][2]
        alleles=collections.defaultdict(int)
        
        # test if seq-string is empty
        if nuc=="*":
            cons.append("N")
            continue
            
        # test whether position within coverage thresholds    
        cov = libraries[j][0]
        if isinstance(covmax, int):
            if int(cov)<covmin or int(cov)>covmax:
                cons.append("N")
                continue
        else:
            if int(cov)<covmin or int(cov)>covmax[j]:
                cons.append("N")
                continue
            
        # find and remove read indices and mapping quality string
        nuc = re.sub(r'\^.',r'',nuc)
        nuc = nuc.replace('$','')
        
        # find and remove InDels
        while "+" in nuc or "-" in nuc:
            if "+" in nuc:
                insertion,ins=extract_indel(nuc,"+")
                #alleles[ins.upper()]+=nuc.count(insertion)
                nuc=nuc.replace(insertion,"")
            else:
                deletion,dele=extract_indel(nuc,"-")
                #alleles[dele.upper()]+=nuc.count(deletion)
                nuc=nuc.replace(deletion,"")
        
        # read all alleles
        for i in range(len(nuc)):
            
            # test for base quality threshold (if below: ignore nucleotide)
            if ord(qualities[i])-33<threshold:
                continue
            # ignore single nucleotide deletions
            if nuc[i]=="*":
                continue
            # count nucleotides similar to reference base
            if nuc[i] =="," or nuc[i] == ".":
                alleles[refbase]+=1
                continue
            # count altenrative nucleotides
            alleles[nuc[i].upper()]+=1
        
        # ignore position if coverage after filtering for 1) InDels and base-quality below threshold
        if sum(alleles.values())<covmin:
            cons.append("N")
            continue
        
        af={}
        for k,v in alleles.items():
            af[v/float(sum(alleles.values()))]=k
        
        if max(af.keys())>cutoff:
            cons.append(af[max(af.keys())])
        else:
            cons.append("N")
    print chrom+"\t"+position+"\t"+"".join(cons)