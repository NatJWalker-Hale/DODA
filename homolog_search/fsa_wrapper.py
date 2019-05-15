import os,sys

if __name__ == "__main__":
	if len(sys.argv) != 2:
		print "usage: fsa_wrapper faDIR"
		sys.exit()

	faDIR = sys.argv[1]+"/"
	for i in os.listdir(faDIR):
		if i.endswith(".fa"):
			cmd = "time fsa --fast "+faDIR+i+" >"+faDIR+i+".fsa.aln"
			print cmd
			os.system(cmd)
			
