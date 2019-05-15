import sys,os
import phylo3,newick3
import tree_utils

"""
to change sequence names with taxa names and make trees more readable,
and output a new file named infile.names

Create a tabular file that each line contains
code	taxon_name
separated by tab, and put it in the same dir as this script

if first character of the first line is >, treat as fasta file
otherwise treat as a tree file
"""

def get_long_id(label,DICT):
	if "@" in label: #for homolog
		spls = label.split("@")
		taxonID = spls[0]
		seqID = spls[1]
		try: return DICT[taxonID]+"@"+seqID
		except: return label
	else: #for ortholog
		if label in DICT:
			return DICT[label]
		else: return label

def taxon_name_subst(original,table=sys.path[0]+"/taxon_table"):
	DICT = {} # key is seq acronym, value is full taxon name, separated by tab
	with open(table, "r") as infile:
		for line in infile:
			spls = line.strip().split("\t")
			if len(spls) > 1:
				DICT[spls[0].replace("|","_")] = spls[1]
	
	with open(original,"r") as infile:
		line = infile.readline()
		is_fasta = True if line[0] == ">" else False
	
	if is_fasta: # for fasta files
		infile = open(original,"r")
		outfile = open(original+".name","w")
		for line in infile:
			if line[0] == ">":
				outfile.write('>'+get_long_id(line.strip()[1:],DICT)+"\n")
			else: outfile.write(line)
		infile.close()
		outfile.close()
	else: # tree file
		with open(original,"r") as infile:
			intree = newick3.parse(infile.readline())
		for i in intree.leaves():
			print i.label,
			i.label = get_long_id(i.label,DICT)
			print i.label
		with open(original+".name","w") as outfile:
			outfile.write(newick3.tostring(intree)+";\n")

if __name__ == "__main__":
	if len(sys.argv) != 2:
		print "python taxon_name_subst.py tree_or_fasta_file"
		sys.exit(0)
	
	original = sys.argv[1]
	taxon_name_subst(original)
