import sys
from collections import defaultdict as d
import re
import gzip

header = d(lambda: d(lambda: d(str)))
# split by delimiter outside of quote symbol
p_equal = re.compile(r'''((?:[^="']|"[^"]*"|'[^']*')+)''')
p_semi = re.compile(r'''((?:[^;"']|"[^"]*"|'[^']*')+)''')
p_straight = re.compile(r'''((?:[^|"']|"[^"]*"|'[^']*')+)''')
p_comma = re.compile(r'''((?:[^,"']|"[^"]*"|'[^']*')+)''')


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


def isfloat(value):
    ''' try if value can be converted to a decimal '''
    try:
        float(value)
        return float(value)
    except ValueError:
        return value


def meta(meta, vcfh, subhead):
    ''' add info of first seven columns to the linehash '''
    for i in range(7):
        vcfh[subhead[i]] = isfloat(meta[i])


def info(INFO, vcfh):
    ''' add INFO field to the linehash '''
    vcfh["INFO"] = d(str)
    # First split the info field by a ";" and than by a "="
    # total=[p_equal.split(x)[1::2] for x in p_semi.split(INFO)[1::2]]
    total = [x.split("=") for x in INFO.split(";")]
    # loop through splitted Info field"
    for k, v in total:
        # test if infofield contains the annotation Annotation
        if k == "ANN":
            ann(v, vcfh)
            continue
        # test if info value is a number
        vcfh["INFO"][k] = isfloat(v)


def ann(anno, vcfh):
    ''' parse the annotation info in the INFO field and add to linehash '''
    # split annotation field by ","
    alleles = anno.split(",")
    # obtain the descriptions of the annotation field from the header, retain the string within single quotes, split this string by "|" and remove the white space
    anot_descr = [x.replace(" ", "") for x in p_straight.split(
        header["INFO"]["ANN"]["Description"].split("'")[1])[1::2]]
    # split
    annlist = []
    for allele in alleles:
        allele_anot = dict(zip(anot_descr, allele.split("|")))
        # print allele_anot

        if allele_anot["Annotation"] == "intergenic_region":
            annlist.append((allele_anot["Annotation"], "", "", ""))
        else:
            annlist.append((allele_anot["Annotation"], allele_anot["Gene_ID"],
                           allele_anot["Gene_Name"], allele_anot["Distance"]))
    vcfh["INFO"]["ANN"] = list(set(annlist))


def samples(sample, name, vcfh, FORMAT):
    ''' parse the annotation info in the INFO field and add to linehash '''
    # isolate all alleles:

    # loop through all samples
    for i in range(len(sample)):
        all_sample = sample[i].split(":")
        # return NA if no counts
        if all_sample[0] == "./.":
            vcfh[name[i]] = "NA"
            continue
        # determine the alleles in sample and remove the reference
        vcfh[name[i]] = d(lambda: d(str))

        alleleh = dict(zip(FORMAT.split(":"), all_sample))
        if "RD" in alleleh:
            vcfh[name[i]]["coverage"] = isfloat(
                alleleh["RD"])+isfloat(alleleh["AD"])
            vcfh[name[i]][vcfh["REF"]]["count"] = isfloat(alleleh["RD"])
            vcfh[name[i]][vcfh["REF"]]["freq"] = isfloat(
                alleleh["RD"])/(isfloat(alleleh["RD"])+isfloat(alleleh["AD"]))
            vcfh[name[i]][vcfh["ALT"]]["count"] = isfloat(alleleh["AD"])
            vcfh[name[i]][vcfh["ALT"]]["freq"] = isfloat(
                alleleh["AD"])/(isfloat(alleleh["RD"])+isfloat(alleleh["AD"]))
        elif "DP" in alleleh:
            vcfh[name[i]]["coverage"] = isfloat(alleleh["DP"])
            vcfh[name[i]][vcfh["REF"]]["count"] = isfloat(
                alleleh["DP"])-isfloat(alleleh["AD"])
            vcfh[name[i]][vcfh["REF"]]["freq"] = vcfh[name[i]
                                                      ][vcfh["REF"]]["count"]/isfloat(alleleh["DP"])
            vcfh[name[i]][vcfh["ALT"]]["count"] = isfloat(alleleh["AD"])
            vcfh[name[i]][vcfh["ALT"]]["freq"] = isfloat(
                alleleh["AD"])/isfloat(alleleh["DP"])


def vcfline(line, head, subhead):
    vcfh = d(str)
    linesp = line.split()
    meta(linesp[:7], vcfh, subhead)
    info(linesp[7], vcfh)
    FORMAT = linesp[8]
    samplenames = subhead[9:]
    samples(linesp[9:], samplenames, vcfh, FORMAT)
    return vcfh


code = ["X", "2L", "2R", "3L", "3R", "X"]
datah = {}
for l in load_data(sys.argv[2]):
    a = l.split()
    if a[0] not in code:
        continue
    ID = a[0]+"_"+a[1]
    datah[ID] = l.rstrip("\n")

for l in load_data(sys.argv[1]):
    # parse header
    if l.startswith("##"):
        new = l[2:-1]
        if new.split("=")[0] not in ["INFO", "FORMAT", "FILTER"]:
            continue
        ID, val = new.split("=<")
        desc = dict([p_equal.split(x)[1::2] for x in p_comma.split(val)[1::2]])
        for k, v in desc.items():
            header[ID][desc["ID"]][k] = v
        continue

    # parse Subheader
    if l.startswith("#"):
        subheader = l.rstrip()[1:].split()
        continue

    # test if SNP in dataset
    a = l.split()
    ID = a[0]+"_"+a[1]
    if ID not in datah:
        continue
    # generate lineobject
    result = vcfline(l, header, subheader)
    for annot in result["INFO"]["ANN"]:
        print datah[ID]+"\t"+"\t".join(map(str, annot))
    del datah[ID]

for k, v in sorted(datah.items()):
    print v+"\t\t\t\t"
