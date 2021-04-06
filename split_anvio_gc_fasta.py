import sys, os.path
from Bio import SeqIO
import glob
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import argparse


parser = argparse.ArgumentParser(description="Split up anvio gene cluster fasta into individual fastas for each GC")
parser.add_argument('-f', help="The Gene-cluster fasta exported from anvio'anvi-get-sequences-for-gene-clusters'")

args = parser.parse_args()
arg_dict = vars(args)
fasta = arg_dict['f']

GC_list = []

for seq_record in SeqIO.parse(fasta, "fasta"):
			contig = str(seq_record.id)
			cuts = contig.rsplit('|',3)[1]
			GC = cuts.rsplit(":", 2)[1]
			if GC not in GC_list:
				GC_list.append(GC)
				with open("%s.faa" % (GC), "a") as outfile:
					ID = str(seq_record.id)		
					outfile.write(">%s\n" % (ID) + str(seq_record.seq)+"\n") 
			else:
				with open("%s.faa" % (GC), "a") as outfile:
					ID = str(seq_record.id)		
					outfile.write(">%s\n" % (ID) + str(seq_record.seq)+"\n") 
outfile.close()
