"""
Input: a starting fasta, and the itinitial baits
Output: refined tree(s) and correcponding fasta file(s)
"""

import phylo3,newick3,os,sys
import tree_utils
from swipe_dir import swipe
import fasta_to_tree
import ntpath
import trim_tips
import mask_tips_by_taxonID_transcripts
import cut_long_internal_branches
import seq


# do not mask tips for these high-quality genomes, 43 in total.
# still mask Dianthus, since it's low quality


GENOMES = ["Acoerulea","Acomosus","Alyrata","Atrichopoda","Bdistachyon","Bstricta","Bvulgaris","Cclementina","Cgigantea","Cgrandiflora","CorSFB","Cquinoa","Crubella","Csativus","Csinensis","Dcarota","Dica","Egrandis","Esalsugineum","Fvesca","Gmax","Graimondii","KecaSFB_gen","LimaeSFB_gen","Lusitatissimum","MacauSFB_gen","Macuminata","Mesculenta","Mguttatus","Osativa","PhaexSFB_gen","Ppersica","Ptrichocarpa","Pvirgatum","Rcommunis","Sbicolor","CVDF_genome","Soleracea","SpergArv","Spurpurea","Sviridis","Tcacao","Tpratense","Phallii","Ahypochondriacus","Acoerulea","Acomosus","Ahalleri","Ahypochondriacus","Ahypochondriacus","Alyrata","Aoccidentale","Aofficinalis","Athaliana","Atrichopoda","BdistachyonBd21","Bdistachyon","Bhybridum","Boleraceacapitata","BrapaFPsc","Bstacei","Bstricta","Bsylvaticum","Bvulgaris","Carietinum","Cclementina","Cgrandiflora","Cpapaya","Cquinoa","Crubella","Csativus","Csinensis","Dcarota","Egrandis","Esalsugineum","Fvesca","Ghirsutum","Gmax","Graimondii","Hannuus","Hvulgare","Kfedtschenkoi","Klaxiflora","Lsativa","Lusitatissimum","Macuminata","Mdomestica","Mesculenta","Mguttatus","Mpolymorpha","Msinensis","Mtruncatula","Oeuropaea","Osativa","Othomaeum","PdeltoidesWV94","PhalliiHAL","Phallii","Phallii","Ppatens","Ppersica","Ptrichocarpa","Ptrichocarpa","Pvirgatum","Pvirgatum","Rcommunis","SbicolorRio","Sbicolor","Sbicolor","Sfallax","Sitalica","Slycopersicum","Smoellendorffii","Soleracea","Spolyrhiza","Spurpurea","Stuberosum","Sviridis","Sviridis","Taestivum","Tcacao","Tintermedium","Tpratense","Vunguiculata","Vvinifera","Zmarina","Zmays","ZmaysPH207"] #"DGBStha"]

def get_filename_from_path(path):
	if path[-1] == "/": path = path[:-1]
	head, tail = ntpath.split(path)
	if head == "": head = "."
	return head+"/",tail

def refine(query_fasta,start_fasta,deep_paralog_cutoff,num_cores):
	gene_name = get_filename_from_path(query_fasta)[1].split(".")[0]
	outdir,fasta = get_filename_from_path(start_fasta)
	#print outdir,fasta
	deep_paralog_cutoff = float(deep_paralog_cutoff)
	query_ids = [s.name for s in seq.read_fasta_file(query_fasta)]
	new_fasta = [] # list of output refined fasta files
	print outdir,fasta
	
	# make a tree from the start_fasta
	tree = fasta_to_tree.fasta_to_tree(outdir,fasta,num_cores,"aa")
	if tree == None: return []
	with open(tree,"r") as infile:
		intree = newick3.parse(infile.readline())
	root = trim_tips.trim(intree,relative_cutoff=deep_paralog_cutoff,absolute_cutoff=deep_paralog_cutoff*2)
	if os.path.exists(outdir+fasta+".pasta.aln-cln"):
		clnfile = outdir+fasta+".pasta.aln-cln"
	else: clnfile = outdir+fasta+".mafft.aln-cln"
	root = mask_tips_by_taxonID_transcripts.mask(root,\
		clnfile=clnfile,\
		para="y",
		ignore=GENOMES)
	if root != None:
		with open(tree+".tt.mm","w") as outfile:
			outfile.write(newick3.tostring(root)+"\n")
		subtrees = cut_long_internal_branches.cut_long_internal_branches(root,cutoff=deep_paralog_cutoff)
		count = 0
		base_name = fasta.split(".")[0]
		seqDICT = {} # key is seqid, value is seq
		for s in seq.read_fasta_file(start_fasta):
			seqDICT[s.name] = s.seq
		for tree in subtrees:
			if tree == None: continue
			label_set = set(tree_utils.get_front_labels(tree))
			if len(label_set) > 4 and len(label_set & set(query_ids)) > 0:
				count += 1
				with open(outdir+base_name+"_"+str(count)+".subtree","w") as outfile:
					outfile.write(newick3.tostring(tree)+";\n")
				with open(outdir+base_name+"_"+str(count)+".fa","w") as outfile:
					for seqid in tree_utils.get_front_labels(tree):
						try:
							outfile.write(">"+seqid+"\n"+seqDICT[seqid]+"\n")
						except:
							print seqid,"not found in fasta file"
				new_fasta.append(outdir+base_name+"_"+str(count)+".fa")
	
	return new_fasta

if __name__ == "__main__":
	if len(sys.argv) != 5:
		print "python refine_homolog.py query_fasta start_fasta deep_paralog_cutoff num_cores"
		print "Example: python src/refine_homolog.py DODA.pep.fa DODA/DODA.swipe.fa 1.0 4"
		sys.exit(0)
	
	query_fasta,start_fasta,deep_paralog_cutoff,num_cores = sys.argv[1:]
	refine(query_fasta,start_fasta,deep_paralog_cutoff,num_cores)

