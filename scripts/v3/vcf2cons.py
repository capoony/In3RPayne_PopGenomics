import sys


def load_data(x):
    ''' import data either from a gzipped or or uncrompessed file or from STDIN'''
    import gzip
    if x == "-":
        y = sys.stdin
    elif x.endswith(".gz"):
        y = gzip.open(x, "rt", encoding="latin-1")
    else:
        y = open(x, "r", encoding="latin-1")
    return y

# linelist=sys.argv[2].split(",")


for l in load_data(sys.argv[1]):
    if l.startswith("#"):
        continue
    a = l.split()
    C, P = a[:2]
    pops = a[9:]
    R, A = a[3:5]
    if len(R) != 1 or len(A) != 1:
        continue
    pl = ""
    for x in range(len(pops)):
        if pops[x].split(":")[0] == "0/0" or pops[x].split(":")[0] == "0|0":
            pl += R
        elif pops[x].split(":")[0] == "1/1" or pops[x].split(":")[0] == "1|1":
            pl += A
        else:
            pl += "N"
    print(C+"\t"+P+"\t"+R+"\t"+pl)
