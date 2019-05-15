"""
back translation
input: peptide alignment and corresponding cds with exactly the same names
output: back translated cds alignments

pepDIR: DIR with peptide alignments
cdsFILE: dna seq matched to peptide tip names. cat them to one fasta file
all cds have to be in the plus direction otherwise pal2nal wouldn't recognize the patterns
outDIR: output slimed down pep alignment and back-translated dna alignment
"""

import sys,os
from Bio import SeqIO

PEP_ALN_FILE_ENDING = ".aln"

if __name__ =="__main__":
	if len(sys.argv) != 4:
		print "usage: python pal2nal_wrapper.py pepDIR cdsFILE outDIR"
		sys.exit()
	
	pepDIR = sys.argv[1]+"/"
	cdsFILE = sys.argv[2]
	outDIR = sys.argv[3]+"/"
	
	print "Reading in cds"
	cdsDICT = {} #key is taxon id, value is seqid
	handle = open(cdsFILE,"rU")
	for record in SeqIO.parse(handle,"fasta"):
		seqid,seq = str(record.id),str(record.seq)
		cdsDICT[seqid] = seq
	handle.close()
	
	#going through pep alignments and back translate
	len_file_ending = len(PEP_ALN_FILE_ENDING)
	for i in os.listdir(pepDIR):
		if i[-len_file_ending:] != PEP_ALN_FILE_ENDING: continue
		pepaln = pepDIR+i
		cds = outDIR+i.replace(PEP_ALN_FILE_ENDING,".cds")
		print "Reading in pep alignment from",pepaln
		outfile = open(cds,"w")
		handle = open(pepaln,"rU")
		for record in SeqIO.parse(handle,"fasta"):
			seqid = (str(record.id).replace("-","_")).replace(".","_")
			try:
				outfile.write(">"+seqid+"\n"+cdsDICT[seqid]+"\n")
			except:
				print seqid, "not found in fasta file"
		handle.close()
		outfile.close()
		
		#back translate using pal2nal
		outaln = outDIR+i.replace(PEP_ALN_FILE_ENDING,".cds.aln")
		cmd = "pal2nal.pl "+pepaln+" "+cds+" -output fasta >"+outaln
		#Use the paml format if doing codeml analysis next
		#cmd = "pal2nal.pl "+peppht+" "+dnafas+" -output paml >"+outaln
		print cmd
		os.system(cmd)

