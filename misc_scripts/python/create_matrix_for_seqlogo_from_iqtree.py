import sys
from Bio import Seq
import pandas as pd

# this is a script designed to produce CSVs to be read into matrices for the r package ggseqlogo, directly from iqtree .state files, for a given node
# runs in two modes: 1) by default, produce matrix for logo for every site in the alignment (for nuc, every site, for codon, every codon, for aa, every site). 2) alternatively, given a range of 
# (zero-indexed) sites in the alignment, return a matrix for just these sites.

def df_from_states(statefile):
    df = pd.read_csv(statefile,sep="\t",header=(8))
    return df

def rename_state_df(df):
    new = {}
    for n in df.columns[3:]:
        new[n] = n.lstrip("p_")
    renamed = df.rename(columns = new)
    return renamed

def get_node_df(df):
    nodedict = {} # key is node, value is pd df of node 
    nodes = list(set([n for n in df["Node"]]))
    for n in nodes:
        nodedict[n] = df.loc[df["Node"] == n]
    return nodedict

def process_and_transpose(df):
    processed = df.iloc[:,3:]
    transpose = processed.T
    return transpose

def translate_process_and_transpose(df):
    new = {}
    for n in df.columns[3:]:
        new[n] = Seq.translate(n)
    renamed = df.rename(columns = new)
    nonredundant = renamed.iloc[:,3:].groupby(renamed.iloc[:,3:].columns,1).sum() # because codons are redundant, sum across resultant translation (I think this is legit)
    transpose = nonredundant.T
    return transpose

if __name__ == "__main__":
    if len(sys.argv) < 2: 
        print "usage: python "+sys.argv[0]+" states node translate[y/n, do n for nuc or aa] site0 site1 site2 ...[optional]"
        sys.exit()

    df = df_from_states(sys.argv[1])
    renamed = rename_state_df(df)
    nodedf = get_node_df(renamed)
    if len(sys.argv) == 4: # first running mode      
        if sys.argv[3] == "y":
            print translate_process_and_transpose(nodedf[sys.argv[2]]).to_csv(header=False)
        else:
            print process_and_transpose(nodedf[sys.argv[2]]).to_csv(header=False)

    else: # second running mode
        sites = [int(s) for s in sys.argv[4:]]
        if sys.argv[3] == "y":
            print translate_process_and_transpose(nodedf[sys.argv[2]]).iloc[:,sites].to_csv(header=False)
        else:
            print process_and_transpose(nodedf[sys.argv[2]]).iloc[:,sites].to_csv(header=False)
        


