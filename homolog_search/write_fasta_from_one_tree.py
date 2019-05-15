import sys,os,newick3
import tree_utils
from Bio import SeqIO

if __name__ =="__main__":
	if len(sys.argv) != 3:
		print "usage: write_fasta_from_one_tree.py fasta tree"
		sys.exit()
	
	fasta,tree = sys.argv[1:]
	with open(sys.argv[2],"r") as infile:
		intree = newick3.parse(infile.readline())
	seqid_set = set(tree_utils.get_front_labels(intree))
	print len(seqid_set),"tips read"
	
	handle = open(sys.argv[1])
	outfile = open(tree+".fa","w")
	for s in SeqIO.parse(handle,"fasta"):
		if s.id in seqid_set:
			outfile.write(">"+str(s.id)+"\n"+str(s.seq)+"\n")
			seqid_set.remove(s.id)
			print s.id,"written to outfile"
	handle.close()
	
	if len(seqid_set) > 0:
		print seqid_set,"not found in fasta file"
