'''

This finds the consensus eggNOG annotation for each gene cluster and combine with other GC metadata
exported from Anvio. This counts up the annotation info for all the sequences 
in a given gene cluster and picks the most frequent one as the consensus annotation for
the entire gene-cluster. It then concatenates that information with the gene cluster
statistics produce by anvio (anvi-summarize "--quick-summary" output : /misc_data_items/default.txt )

'''


import sys, os.path
import pandas as pd
from collections import defaultdict
import glob
import argparse


parser = argparse.ArgumentParser(description="Finds consensus eggNOG annotation for each gene cluster and combine with other GC metadata exported from anvio")
parser.add_argument('-i', help="Path to directory containing eggNOG '.annotations' files for each gene-cluster")
parser.add_argument('-meta', help="Anvio gene-cluster statistics file: This is file is produced via anvi-summarize --quick-summary and will be locoated here: /misc_data_items/default.txt")
parser.add_argument('-o', help="Name of output file")
args = parser.parse_args()
arg_dict = vars(args)

o_name = arg_dict['o'] 
abs_path = os.path.abspath(os.getcwd()) 
o_name_abs = os.path.join(abs_path, o_name)

md = arg_dict['meta']
gc_dat = pd.read_csv(md, sep = "\t")
gc_dat.rename(columns = {'items':'gene_cluster'}, inplace = True)

en = arg_dict['i']
en_files = glob.glob(os.path.join(en, "*"))

gc_cons_annot = pd.DataFrame(columns = ['gene_cluster', 'Preferred_name', 'KEGG_ko', 'COG', 'eggNOG free text desc.'])

for f in en_files:
	(dir, file) = os.path.split(f)
	gc = file.rsplit(".", 3)[0]
	KO = []
	pref_name = []
	eggnog = []
	COG = []   
	rfile = open(f, "r")
	for line in rfile.readlines():
		line = line.rstrip()
		if line.startswith("#"):
			continue
		else:
			tabs = line.rsplit("\t")
			KO.append(tabs[8])
			pref_name.append(tabs[5])
			eggnog.append(tabs[-1])
			COG.append(tabs[-2])
	if len(KO) > 0:
		ko_cons = max(set(KO), key = KO.count)
	if len(pref_name) > 0:
		name_cons = max(set(pref_name), key = pref_name.count)
	if len(eggnog) > 0:
		eggnog_cons = max(set(eggnog), key = eggnog.count)
	if len(COG) > 0:
		COG_cons = max(set(COG), key = COG.count)
	
	gc_cons_annot = gc_cons_annot.append({'gene_cluster' : gc, 'Preferred_name' : name_cons, 'KEGG_ko' : ko_cons, 'COG' : COG_cons,'eggNOG free text desc.' : eggnog_cons}, ignore_index=True)

gc_cons_w_meta = pd.merge(gc_cons_annot, gc_dat, on="gene_cluster")
gc_cons_w_meta.to_csv(o_name_abs, sep = "\t", header = True, index = False)
