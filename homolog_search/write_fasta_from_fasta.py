import sys,os
from Bio import SeqIO

if __name__ =="__main__":
	if len(sys.argv) != 3:
		print "usage: write_fasta_from_one_tree.py all_fasta in_fasta"
		sys.exit()
	
	all_fasta,in_fasta = sys.argv[1:]
	seqids = []
	with open(in_fasta,"r") as infile:
		for line in infile:
			if line[0] == ">":
				seqids.append(line.strip()[1:])
	seqid_set = set(seqids)
	print len(seqid_set),"seqids read"
	
	handle = open(sys.argv[1])
	outfile = open(in_fasta+".fa","w")
	for s in SeqIO.parse(handle,"fasta"):
		if s.id in seqid_set:
			outfile.write(">"+str(s.id)+"\n"+str(s.seq)+"\n")
			seqid_set.remove(s.id)
			print s.id,"written to outfile"
	handle.close()
	
	if len(seqid_set) > 0:
		print seqid_set,"not found in fasta file"
