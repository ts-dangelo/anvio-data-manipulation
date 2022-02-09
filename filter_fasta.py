import sys, os.path
from Bio import SeqIO
import glob
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import argparse
import multiprocessing as mp
from multiprocessing import Pool
from functools import partial
import os,shutil


parser = argparse.ArgumentParser(description="Remove Sequences from a pre-existing MSA or fasta file")
parser.add_argument('-input_dir', help="Alignment file, or fasta file")
parser.add_argument('-list', help="List of genomes you want kept in the files")
parser.add_argument('-out', help="name of output directory")

args = parser.parse_args()
arg_dict = vars(args)

genomes = arg_dict['list']
oname = arg_dict['out'] 

abs_path = os.path.abspath(os.getcwd()) 
wo_dir = os.path.join(abs_path, oname)
if not os.path.exists(wo_dir):
        os.makedirs(wo_dir)

gc_path = arg_dict['input_dir']
gc_paths = glob.glob(os.path.join(gc_path, "*"))

cut_list = []
file = open(genomes,'r')
for line in file :
        line = line.rstrip()
        cut_list.append(line)


def rename(wo_dir, gc_path):

        try:
                
                (dir, file) = os.path.split(gc_path)
                gc_name = file.rsplit('.')[0]
                gc_out = os.path.join(wo_dir, file)
                with open(gc_out, "a") as outfile:
                        for seq_record in SeqIO.parse(gc_path, "fasta"):
                                gen_id = str(seq_record.id)
                                id1 = gen_id.rsplit('|',3)[2]
                                genid = id1.rsplit(":",2)[1]
                                grid = genid.replace('_', '-')
                                id2 = gen_id.rsplit('|',3)[3]
                                orfid = id2.rsplit(":",2)[1]                            
                                if genid in cut_list:
                                        outfile.write(">%s_%s\n" % (grid, orfid) + str(seq_record.seq)+"\n") 
                                else:
                                        continue
        except:
                print("Error with " + gc_path)
        

pool = Pool(processes=96)
pool.map(partial(rename, wo_dir), gc_paths)
pool.close()
