import sys,os,re

## this script takes a constraint tree backbone (e.g. a family-level tree, but it can be any) and a list of taxa with the corresponding higher level of taxonomy they come from (e.g. family for a 
## family-level tree. The taxon list must be each on a new line and take the form [higher_level_taxonomy]_[taxon], although the underscore can be replaced by any separator by changing the definition.
## this will produce a multifurcating constraint tree where the taxa in the file are place in the tree in a polytomy according to the position of the higher level taxon. 

separator = "_" # change this if you would like a different separator

def get_backbone_taxa(tr_str):
    tr_str = re.sub(",","\n",tr_str)
    taxa = re.sub("[\(\);]","",tr_str).split("\n")
    return taxa

def make_taxa_dict(higher_taxa_list,lower_taxa_list):
    taxa_dict = {}
    for f in higher_taxa_list:
        taxa_dict[f] = [s for s in lower_taxa_list if s.split(separator)[0] == f]
    return taxa_dict

def create_constraint(taxa_dict,tr_str):
    for f in list(taxa_dict.keys()):
        tr_str = re.sub(f,",".join(taxa_dict[f]),tr_str)
    return tr_str

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print "Usage: python "+sys.argv[0]+" constraint_tree taxon_file"
        sys.exit()
    tr_str = open(sys.argv[1],"r").read()
    lower_taxa_list = open(sys.argv[2],"r").read().splitlines()
    higher_taxa_list = get_backbone_taxa(tr_str)
    #print higher_taxa_list
    taxa_dict = make_taxa_dict(higher_taxa_list,lower_taxa_list)
    #print taxa_dict
    print create_constraint(taxa_dict,tr_str)

