"""
Input:

- prepare a input fasta file that contains sequences used as bait:
	genename.pep.fa

- prepare a directory with candidate sequence pep fasta files named
	taxonid.pep.fa, or
	taxonid.pep.fa.cdhit
  for each sequence, the sequence ids look like taxonid@seqid

swipe, taking the top 20 output, construct the homolog
and carry out two rounds of refinement
replace the seqids to long id for visualization in figtree
"""

import phylo3,newick3,os,sys
from swipe_dir import swipe
import fasta_to_tree
import ntpath
import refine_homolog
import taxon_name_subst

# do not mask tips for these high-quality genomes, 43 in total.
# still mask Dianthus, since it's low quality
GENOMES = ["Acoerulea","Acomosus","Alyrata","Atrichopoda","Bdistachyon","Bstricta","Bvulgaris","Cclementina","Cgigantea","Cgrandiflora","CorSFB","Cquinoa","Crubella","Csativus","Csinensis","Dcarota","Dica","Egrandis","Esalsugineum","Fvesca","Gmax","Graimondii","KecaSFB_gen","LimaeSFB_gen","Lusitatissimum","MacauSFB_gen","Macuminata","Mesculenta","Mguttatus","Osativa","PhaexSFB_gen","Ppersica","Ptrichocarpa","Pvirgatum","Rcommunis","Sbicolor","CVDF_genome","Soleracea","SpergArv","Spurpurea","Sviridis","Tcacao","Tpratense","Phallii","Ahypochondriacus","Acoerulea","Acomosus","Ahalleri","Ahypochondriacus","Ahypochondriacus","Alyrata","Aoccidentale","Aofficinalis","Athaliana","Atrichopoda","BdistachyonBd21","Bdistachyon","Bhybridum","Boleraceacapitata","BrapaFPsc","Bstacei","Bstricta","Bsylvaticum","Bvulgaris","Carietinum","Cclementina","Cgrandiflora","Cpapaya","Cquinoa","Crubella","Csativus","Csinensis","Dcarota","Egrandis","Esalsugineum","Fvesca","Ghirsutum","Gmax","Graimondii","Hannuus","Hvulgare","Kfedtschenkoi","Klaxiflora","Lsativa","Lusitatissimum","Macuminata","Mdomestica","Mesculenta","Mguttatus","Mpolymorpha","Msinensis","Mtruncatula","Oeuropaea","Osativa","Othomaeum","PdeltoidesWV94","PhalliiHAL","Phallii","Phallii","Ppatens","Ppersica","Ptrichocarpa","Ptrichocarpa","Pvirgatum","Pvirgatum","Rcommunis","SbicolorRio","Sbicolor","Sbicolor","Sfallax","Sitalica","Slycopersicum","Smoellendorffii","Soleracea","Spolyrhiza","Spurpurea","Stuberosum","Sviridis","Sviridis","Taestivum","Tcacao","Tintermedium","Tpratense","Vunguiculata","Vvinifera","Zmarina","Zmays","ZmaysPH207"] #"DGBStha"]


def get_filename_from_path(path):
	head, tail = ntpath.split(path)
	return tail or ntpath.basename(head)

if __name__ == "__main__":
	if len(sys.argv) != 7:
		print "python bait_homologs.py query_pep_fa pepDIR num_to_bait num_cores outdir deep_paralog_cutoff"
		sys.exit(0)
	
	query_fasta,pepdir,num_to_bait,num_cores,out,deep_paralog_cutoff = sys.argv[1:] 
	# the query fasta name is gene_name.pep.fa
	gene_name = get_filename_from_path(query_fasta).split(".")[0]
						
	"""
	get the initial fasta files
	swipe query_fasta against each data set in DIR
	write output in outdir
	take the top num_to_bait hits ranked by bitscore
	evalue set to 10
	"""
	outdir = out+"/" # +"/"+gene_name+"/"
	try: os.stat(outdir)
	except: os.mkdir(outdir)
	if not os.path.exists(outdir+gene_name+".swipe.fa"):
		swipe(query_fasta=query_fasta,\
			pepdir=pepdir,\
			outdir=outdir,\
			num_cores=num_cores,\
			max_num_hits=num_to_bait)
	
	# first round of refine
	refined_fastas = refine_homolog.refine(query_fasta=query_fasta,\
		start_fasta=outdir+gene_name+".swipe.fa",\
		deep_paralog_cutoff=float(deep_paralog_cutoff),
		num_cores=num_cores)
	print refined_fastas
	
	# second round of refine
	new_fastas = []
	if len(refined_fastas) > 0:
		for i in refined_fastas:
			new_fastas += refine_homolog.refine(query_fasta=query_fasta,\
				start_fasta=i,\
				deep_paralog_cutoff=float(deep_paralog_cutoff),
				num_cores=num_cores)
	print new_fastas
	refined_fastas = new_fastas
				
	# if last round was using fastree, refine one more time
	if len(refined_fastas) > 0 and os.path.exists(outdir+gene_name+"_1.fasttree.tre"):
		for i in refined_fastas:
			refine_homolog.refine(query_fasta=query_fasta,\
				start_fasta=i,\
				deep_paralog_cutoff=float(deep_paralog_cutoff),
				num_cores=num_cores)
	
	# change the names for all the .mm and .subtree files for visualization
	for i in os.listdir(outdir):
		if i.endswith(".mm") or i.endswith(".subtree"):
			taxon_name_subst.taxon_name_subst(original=outdir+i)
			
