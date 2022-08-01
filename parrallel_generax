'''
Run GeneRax in Parrallel
'''

import subprocess 
import sys
import glob
import argparse
import os,shutil
import multiprocessing as mp
from multiprocessing import Pool
from functools import partial

parser = argparse.ArgumentParser(description="Run GeneRax in Parrallel using ProcessPoolExecutor")
parser.add_argument('-family', help="Path to directory containing family files")
parser.add_argument('-tree', help="Rooted species tree")
parser.add_argument('-results', help="Directory for GeneRax results subdirectories")
args = parser.parse_args()
arg_dict = vars(args)

abs_path = os.path.abspath(os.getcwd()) 

o_dir = arg_dict['results'] 
res_dir = os.path.join(abs_path, o_dir)
if not os.path.exists(res_dir):
	os.makedirs(res_dir)
		
fpath = arg_dict['family']
fam_paths = glob.glob(os.path.join(fpath, "*"))

species_tree = arg_dict['tree'] 
		
def run_generax(species_tree, res_dir, fam_path):
    
    try:
    	(dir, file) = os.path.split(fam_path)
        gc_id = file.rsplit('.')[0]
    	fam_res = res_dir + "/" + gc_id
    	subprocess.call(["generax", "-f ", fam_path , "-s", species_tree, "p", fam_res])
        
    except:
        print("Error with %s" %(fam_path))
        

pool = Pool(processes=96)
pool.map(partial(run_generax, species_tree, res_dir), fam_paths)
pool.close()
