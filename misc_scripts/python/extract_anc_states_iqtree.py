import sys
import numpy
import pandas as pd

states = pd.read_csv(sys.argv[1],sep="\t",header=(8))

nodes = list(set([n for n in states["Node"]]))

seqdict = {} # key is node, value is list of codons
for n in nodes: 
    seqdict[n] = "".join([c for c in states.loc[states["Node"] == n]["State"]])

for k in seqdict.keys():
    print ">"+k
    print seqdict[k]



