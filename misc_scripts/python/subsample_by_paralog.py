from __future__ import division
import sys,os
from ete3 import Tree,PhyloTree

# script to return a list of subsampled sequences based on genera in paralogs, such that each genus is only represented once. All names must be in format family_genus_species@seqid and tip names must
# match alignment names exactly. Alignment in fasta format by default, this can easily be changed.

CANDIDATES = ["Amaranthaceae_Ptilotus_exaltatus_hyb@KF747354.1_cds_AIS23505.1_1","Cactaceae_Carnegia_gigantea@Cgig1_18485_CgDODAa1","Caryophyllaceae_Telephium_imperatii@DODAa","Chenopodiaceae_Beta_vulgaris@XM_010675993.2_cds_XP_010674295.1_1","Kewaceae_Kewa_caespitosa@DODAa2","Montiaceae_Parakeelya_mirabilis@KF747352.1_cds_AIS23503.1_1","Nyctaginaceae_Mirabilis_jalapa@c44693_g1_i1_97_870_minus","Portulacaceae_Portulaca_grandiflora@AJ580598.1_cds_CAE45178.1_1"]

def mask_monophyly(rootnode): # this will modify the given tree in place
    dtip = []
    species_list = list(set([ leaf.species for leaf in t.iter_leaves() ]))
    for s in species_list:
        masks = [ n for n in t.get_monophyletic(values=[s],target_attr="species") if len(n) > 1 ]
        nlist = []
        for n in masks:
            for l in n.iter_leaves():
                nlist.append((l.name,len(l.sequence.replace("-",""))))
            dtip += [ x[0] for x in sorted(nlist,reverse=True,key = lambda x:x[1])[1:] ]
            nlist = []
    tips = set( [ l.name for l in t.iter_leaves() ] )
    ktips = list(tips.symmetric_difference(set(dtip)))
    t.prune(ktips)

def sample_by_paralog(root,des_seq):
    dups_to_process = t.split_by_dups()
    keep_seqs = []
    for n in dups_to_process:
        orig_prop = len(n)/len(root)
        #print "original proportion is: "+str(orig_prop)
        des_num = int(round(orig_prop*int(des_seq)))
        #print "desired number is: "+str(des_num)
        fam_dict = {} # key is family, value is list of tuples, (seq,seq_len)
        tips = [ (l.name,len(l.sequence.replace("-",""))) for l in n.iter_leaves() ]
        for i in tips:
            if not i[0].split("_")[0] in fam_dict.keys():
                fam_dict[i[0].split("_")[0]] = [i]
            else:
                fam_dict[i[0].split("_")[0]].append(i)
        for j in fam_dict.keys():
            if len(fam_dict[j]) == 1:
                keep_seqs.append(fam_dict[j][0][0])
            else:
                fam_prop = len(fam_dict[j])/len(n)
                #print "proportion of family "+j+" in paralog is: "+str(fam_prop)
                des_fam = int(round(fam_prop*des_num))
                #print "desired number of sequences in family "+j+" in paralog is: "+str(des_fam)
                if des_fam <= 1:
                    keep_seqs.append(max(fam_dict[j],key=lambda x:x[1])[0])
                else:
                    keep_seqs += [ x[0] for x in sorted(fam_dict[j],reverse=True,key=lambda x:x[1])[:des_fam] ] 
    return keep_seqs     
    

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print "usage: python "+sys.argv[0]+" tree alignment final_number_of_sequences"
        sys.exit()

    t = PhyloTree(sys.argv[1],alignment=sys.argv[2],alg_format="fasta") # change format for different alignment formats
    t.set_species_naming_function(lambda node: node.name.split("@")[0].rstrip("_SRA").rstrip("_PAC") ) # change this id for different naming conventions, best to have family_genus_species@seqid.
    # Extensions  _SRA _PAC etc. are to allow for different sequencing efforts of the same taxa.
    mask_monophyly(t)
    print "\n".join(sample_by_paralog(t,sys.argv[3])+CANDIDATES)

