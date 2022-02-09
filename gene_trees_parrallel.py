'''
Using pool processing to make lots of trees in parrallel. Change subprocess parameters if needed
'''


import multiprocessing as mp
from multiprocessing import Pool
from functools import partial
import subprocess
import sys
import glob
import argparse
import os,shutil

parser = argparse.ArgumentParser(description="Make trees in Parrallel using ProcessPoolExecutor")
parser.add_argument('-i', help="Path to directory containing only the gene-cluster amino acid fastas")
parser.add_argument('-out_dir', help="Directory for output files")
args = parser.parse_args()
arg_dict = vars(args)

o_dir = arg_dict['out_dir']
abs_path = os.path.abspath(os.getcwd())
wo_dir = os.path.join(abs_path, o_dir)
if not os.path.exists(wo_dir):
  os.makedirs(wo_dir)

ipath = arg_dict['i']
aln_path = glob.glob(os.path.join(ipath, "*.trimmed"))


def run_tree_workflow(treedir, aln_path):

  try:
    subprocess.call(["iqtree", "-nt","AUTO", "-m", "TEST", "-st", "AA", "-bb", "1000", "-s", aln_path])
    gc_name_dir = aln_path.rsplit(".")[0]
    tree_files = glob.glob(gc_name_dir + ".aln.trimmed.*")
    for f in tree_files:
      tabs_path = os.path.abspath(f)
      tf_move = os.path.join(wo_dir, tabs_path)
      shutil.move(tf_move, wo_dir)
  except:
     print("Error with %s" %(aln_path))



pool = Pool(processes=96)
pool.map(partial(run_tree_workflow, wo_dir), aln_path)
pool.close()
