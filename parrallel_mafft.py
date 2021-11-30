'''
Do a lot of alignments at once
'''


import multiprocessing as mp
from multiprocessing import Pool
from functools import partial
import subprocess 
import sys
import glob
import argparse
import os,shutil

parser = argparse.ArgumentParser(description="Parrallel MAFFT using ProcessPoolExecutor")
parser.add_argument('-i', help="Path to gene-cluster amino acid fastas")
parser.add_argument('-out_dir', help="Directory for output files")
args = parser.parse_args()
arg_dict = vars(args)

o_dir = arg_dict['out_dir'] 
abs_path = os.path.abspath(os.getcwd()) 
wo_dir = os.path.join(abs_path, o_dir)
if not os.path.exists(wo_dir):
        os.makedirs(wo_dir)
                
ipath = arg_dict['i']
gc_path = glob.glob(os.path.join(ipath, "*"))
        

def run_tree_workflow(treedir, gc_path):
    
        try:
                (dir, file) = os.path.split(gc_path)
                gc_name = file.rsplit(".")[0]
                aln_name = gc_name + ".aln"
                aln_out = os.path.join(wo_dir, aln_name)
                with open(aln_out, "w") as outfile:
                        subprocess.call(["mafft", "--retree","2", "--quiet", gc_path], stdout=outfile)
                

        except:
                print("Error with %s" %(gc_path))
                

pool = Pool(processes=96)
pool.map(partial(run_tree_workflow, wo_dir), gc_path)

pool.close()
