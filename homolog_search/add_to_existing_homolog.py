"""
Input:

- prepare a input fasta file that contains sequences used as bait:
	genename.pep.fa

- prepare a directory with candidate sequence pep fasta files named
	taxonid.pep.fa, OR
	taxonid.pep.fa.cdhit
  for each sequence, the sequence ids look like taxonid@seqid

swipe, taking the top output (ignoring taxon already in the existing homolog)

construct the homolog
and carry out two rounds of refinement
replace the seqids to long id for visualization in figtree
"""

import phylo3,newick3,os,sys
from swipe_dir import swipe
from seq import read_fasta_file
import fasta_to_tree
import ntpath
import refine_homolog
import taxon_name_subst



# do not mask tips for these genomes, 43 in total
GENOMES = ["Achn","Acoerulea","Athaliana","Atrichopoda","Beta",\
		   "Brapa","Bstricta","Caan","Cclementina","Cila",\
		   "Coca","Cpapaya","Crubella","Csativus","Csinensis",\
		   "Dica","Egrandis","Elgu","Esalsugineum","Fvesca",\
		   "Gmax","Graimondii","Lusitatissimum","Mesculenta","Mguttatus",\
		   "Mtruncatula","Musa","Nenu","Osativa","Pheq",\
		   "Phoe","Ppersica","Ptrichocarpa","Rara","Rcommunis",\
		   "Slycopersicum","Spolyrhiza","Spurpurea","Tcacao","Utri.plus",\
		   "Vvinifera","Zmays","Spol"]

min_bitscore=50.0

def get_filename_from_path(path):
	head, tail = ntpath.split(path)
	return tail or ntpath.basename(head)

if __name__ == "__main__":
	if len(sys.argv) != 7:
		print "python add_to_existing_homolog.py query_pep_fa pepDIR num_to_bait num_cores outdir existing_homolog_fasta"
		sys.exit(0)
	
	query_fasta,pepdir,num_to_bait,num_cores,out,existing_homolog = sys.argv[1:]
	# process input
	if pepdir[-1] != "/": pepdir += "/"
	if out[-1] != "/": out += "/"
	num_cores = str(num_cores)
	max_num_hits = int(num_to_bait)
	
	# the query fasta name is gene_name.pep.fa
	gene_name = get_filename_from_path(query_fasta).split(".")[0]
						
	outdir = out+"/"+gene_name+"_add/"
	try: os.stat(outdir)
	except: os.mkdir(outdir)
	if not os.path.exists(outdir+gene_name+".swipe.fa"):
		# get taxon ids appeared in the initial cluster and skip adding these
		old_taxa = []
		with open(existing_homolog,"r") as infile:
			for line in infile:
				if line[0] == ">": # fasta format
					taxonid = line[1:].strip().split("@")[0]
					if taxonid not in old_taxa:
						old_taxa.append(taxonid)
		print len(old_taxa),"taxa in the existing homolog:"
		print old_taxa
		"""
		get the initial fasta files
		wipe query_fasta against each data set in DIR with peptide fasta files that
		either end with .pep.fa or.cdhit, swipe on each one
		write output in outdir
		take the top num_to_bait hits ranked by bitscore
		evalue set to 10
		ignore all hits with bitscore less than 0.1 of the highest hit from the query
		"""
		# swipe with each taxon
		pepfiles = []
		for i in os.listdir(pepdir):
			if (i.endswith(".pep.fa") or i.endswith(".cdhit")):
				taxonid = i.split(".")[0]
				print taxonid
				if taxonid not in old_taxa:
					pepfiles.append(i)
		datasets = [i.split(".")[0] for i in pepfiles]
		print len(pepfiles),"input peptide files read"
		print pepfiles
		assert len(set(datasets)) == len(datasets),\
			"dataset name repeats. remove duplicated sets"
		temp_files = [] # will remove these later
		for i in pepfiles:
			if not os.path.exists(pepdir+i+".psd"): # make sure there's a blast database
				os.system("makeblastdb -in "+pepdir+i+" -parse_seqids -dbtype prot -out "+pepdir+i)
			swipe_outname = pepdir+i+"."+gene_name+".swipe"
			temp_files.append(swipe_outname)
			temp_files.append(swipe_outname+".hits")
			if not os.path.exists(swipe_outname):
				cmd = "swipe -d "+pepdir+i+" -i "+query_fasta+" -a "+num_cores
				cmd += " -p blastp -o "+swipe_outname+" -m 8 -e 10"
				print cmd
				os.system(cmd)
			assert os.path.exists(swipe_outname), \
				"swipe did not finish correctly"
			"""swipe output:
			Query id, Subject id, % identity, alignment length, mismatches, 
			gap openings, q. start, q. end, s. start, s. end, e-value, bit score"""
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
					if bitscore < 0.1*highest: break # stop recording
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
		seqids_written = [] # avoid seq id repeats
		for s in read_fasta_file(query_fasta):
			seqid,seq = str(s.name),str(s.seq)
			if seqid not in seqids_written:
				outfile3.write(">"+seqid+"\n"+seq+"\n")
				seqids_written.append(seqid)
		for s in read_fasta_file(existing_homolog):
			seqid,seq = str(s.name),str(s.seq)
			if seqid not in seqids_written:
				outfile3.write(">"+seqid+"\n"+seq+"\n")
				seqids_written.append(seqid)
		for i in pepfiles:
			# write the gene_name.rawswipe file to outdir
			with open(pepdir+i+"."+gene_name+".swipe","r") as infile:
				for line in infile:
					outfile1.write(line)
			seqDICT = {} # a dict for each taxon, key is seq name, value is seq
			for s in read_fasta_file(pepdir+i):
				seqDICT[s.name] = s.seq
			with open(pepdir+i+"."+gene_name+".swipe.hits","r") as infile:
				for line in infile:
					line = line.strip()
					if len(line) == 0: continue
					outfile2.write(line+"\n")
					if line not in seqids_written:
						outfile3.write(">"+line+"\n"+seqDICT[line]+"\n")
		outfile1.close()
		outfile2.close()
		outfile3.close()
		
		# remove intermediate files
		for i in temp_files: os.remove(i)

	# first round of refine
	refined_fastas = refine_homolog.refine(query_fasta=query_fasta,\
		start_fasta=outdir+gene_name+".swipe.fa",\
		deep_paralog_cutoff=1.0,
		num_cores=num_cores)
	print refined_fastas
	
	# second round of refine
	new_fastas = []
	if len(refined_fastas) > 0:
		for i in refined_fastas:
			new_fastas += refine_homolog.refine(query_fasta=query_fasta,\
				start_fasta=i,\
				deep_paralog_cutoff=0.5,
				num_cores=num_cores)
	print new_fastas
	refined_fastas = new_fastas
				
	# if last round was using fastree, refine one more time
	if len(refined_fastas) > 0 and os.path.exists(outdir+gene_name+"_1.fasttree.tre"):
		for i in refined_fastas:
			refine_homolog.refine(query_fasta=query_fasta,\
				start_fasta=i,\
				deep_paralog_cutoff=0.5,
				num_cores=num_cores)
	
	# change the names for all the .mm and .subtree files for visualization
	for i in os.listdir(outdir):
		if i.endswith(".mm") or i.endswith(".subtree"):
			taxon_name_subst.taxon_name_subst(original=outdir+i)
			
