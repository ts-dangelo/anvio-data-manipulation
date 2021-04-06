import sys, os.path
from Bio import SeqIO
import glob
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import argparse
from collections import defaultdict
import pandas as pd

parser = argparse.ArgumentParser(description="Create table of protein family count per genome")
parser.add_argument('-genomes', help="Text file with genome identifiers that are the first field before the first underscore in the fasta headers")
parser.add_argument('-fam_path', help="Path to faa files of protein families")
parser.add_argument('-ext', help="fasta file extension")
parser.add_argument('-output', help="Output table file name")
args = parser.parse_args()
arg_dict = vars(args)

fpath = arg_dict['fam_path']
fpaths = glob.glob(os.path.join(fpath, "*%s" % arg_dict['ext']))
gen = arg_dict['genomes']

gen_fam = {} #list of genome names
file = open(gen,'r')
for line in file :
		line = line.rstrip()
		fam_count = defaultdict(int)
		for faa in fpaths:
				(dir, file) = os.path.split(faa)
				OG = file.rsplit('.',1)[0]
				for seq_record in SeqIO.parse(faa, "fasta"):
						header = str(seq_record.id)
						cuts = header.rsplit('|',3)[2]
						genid = cuts.rsplit(":",2)[1]
						if genid == line :
							fam_count[OG] += 1
						else:
							continue
		gen_fam[line] = fam_count
out = arg_dict['output']
with open("%s.txt" % (out), "w") as outfile:
        df = pd.DataFrame.from_dict(gen_fam, orient='index')
        df = df.fillna(0)               
        df.to_csv(outfile, index = True, header=True, sep='\t')
