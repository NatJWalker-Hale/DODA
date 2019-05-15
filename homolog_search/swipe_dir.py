"""
dependencies: put swipe and makeblastdb in the PATH

Also by default ignores hits with bitscore less than 5% of self hits.
"""

import sys,os
import seq
import ntpath

def get_filename_from_path(path):
	head, tail = ntpath.split(path)
	return tail or ntpath.basename(head)

def swipe(query_fasta,pepdir,outdir,num_cores,max_num_hits=10,min_bitscore=30.0):
	"""
	get the initial fasta files
	swipe query_fasta against each data set in DIR with peptide fasta files that
	either end with .pep.fa or.cdhit, swipe on each one
	write output in outdir
	take the top num_to_bait hits ranked by bitscore
	evalue set to 10
	ignore all hits with bitscore less than 0.1 of the highest hit from the query
	"""
	if pepdir[-1] != "/": pepdir += "/"
	if outdir[-1] != "/": outdir += "/"
	num_cores = str(num_cores)
	max_num_hits = int(max_num_hits)
	gene_name = get_filename_from_path(query_fasta).split(".")[0]
	
	# swipe with each taxon
	pepfiles = [i for i in os.listdir(pepdir) if (i.endswith(".pep.fa") or i.endswith(".cdhit"))]
	datasets = [i.split(".")[0] for i in pepfiles]
	print len(pepfiles),"input peptide files read"
	assert len(set(datasets)) == len(datasets),\
		"dataset name repeats. remove duplicated sets"
	temp_files = [] # will remove these later
	for i in os.listdir(pepdir):
		if i.endswith(".pep.fa") or i.endswith(".cdhit"):
			if not os.path.exists(pepdir+i+".psd"):
				os.system("makeblastdb -in "+pepdir+i+" -parse_seqids -dbtype prot -out "+pepdir+i)
			swipe_outname = pepdir+i+"."+get_filename_from_path(query_fasta).split(".")[0]+".swipe"
			temp_files.append(swipe_outname)
			temp_files.append(swipe_outname+".hits")
			if not os.path.exists(swipe_outname):
				cmd = "swipe -d "+pepdir+i+" -i "+query_fasta+" -a "+num_cores
				cmd += " -p blastp -o "+swipe_outname+" -m 8 -e 10"
				print cmd
				os.system(cmd)
			assert os.path.exists(swipe_outname), \
				"swipe did not finish correctly"
			"""
			swipe output colums are:
			Query id, Subject id, % identity, alignment length, mismatches, 
			gap openings, q. start, q. end, s. start, s. end, e-value, bit score
			"""
			# summarize the hit seq ids
			if not os.path.exists(swipe_outname+".hits"):
				hit_tuples = [] # a list of tuples (hit, bitscore)
				with open(swipe_outname,"r") as infile:
					for line in infile:
						if len(line) < 3: continue # skip empty lines
						spls = line.strip().split("\t")
						query,hit,bitscore = spls[0],spls[1].replace("lcl|",""),float(spls[-1])
						if query != hit and bitscore >= min_bitscore:
							hit_tuples.append((hit,bitscore))
				
				out = [] # unique hit ids
				highest = 0.0 # store the highest bitscore
				for hit, bitscore in sorted(hit_tuples,key=lambda x:x[1],reverse=True):
					if highest == 0.0: highest = bitscore # record the first hit
					if bitscore < 0.05*highest: break # stop recording
					if hit not in out:
						out.append(hit)
					if len(out) == max_num_hits: break
				if len(out) == 0: print "Warning: No hits found"
				with open(swipe_outname+".hits","w") as outfile:
					for hit in out:
						print hit
						outfile.write(hit+"\n")
	
	# write output summary files .rawswipe, .filetered_hits, .swipe.fa
	print "Writing output fasta"
	outfile1 = open(outdir+gene_name+".rawswipe","w")
	outfile2 = open(outdir+gene_name+".filetered_hits","w")
	outfile3 = open(outdir+gene_name+".swipe.fa","w")
	query_seqids = [] # avoid seq id repeats
	with open(query_fasta,"r") as infile:
		for line in infile:
			outfile3.write(line) # copy over query seqs
			if line[0] == ">":
				query_seqids.append(line.strip()[1:])
	for i in os.listdir(pepdir):
		if i.endswith(".pep.fa") or i.endswith(".cdhit"):
			# write the gene_name.rawswipe file to outdir
			with open(pepdir+i+"."+get_filename_from_path(query_fasta).split(".")[0]+".swipe","r") as infile:
				for line in infile:
					outfile1.write(line)
			seqDICT = {} # a dict for each taxon, key is seq name, value is seq
			for s in seq.read_fasta_file(pepdir+i):
				seqDICT[s.name] = s.seq
			with open(pepdir+i+"."+get_filename_from_path(query_fasta).split(".")[0]+".swipe.hits","r") as infile:
				for line in infile:
					line = line.strip()
					if len(line) == 0: continue
					outfile2.write(line+"\n")
					if line not in query_seqids:
						outfile3.write(">"+line+"\n"+seqDICT[line]+"\n")
	outfile1.close()
	outfile2.close()
	outfile3.close()
	
	# remove intermediate files
	for i in temp_files: os.remove(i)


if __name__ == "__main__":
	if len(sys.argv) != 6:
		print "python swipe_dir.py query_fasta pepdir outdir max_num_hits num_cores"
		sys.exit(0)
	
	query_fasta,pepdir,outdir,max_num_hits,num_cores = sys.argv[1:]
	swipe(query_fasta,pepdir,outdir,num_cores,max_num_hits)


