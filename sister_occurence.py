''' 

This program calculate verticality and 'sister diversity' metrics for a set of rooted
phylogenetic trees.  The trees can be made via an ortholog identifying pipeline and a phylogeny
constructing tool.  Rooting by the Minimal Ancestor Deviation (M.A.D) method can be used. 

Tree tips must be formatted like this: GENOMENAME_GENEINFORMATION. The genome name must be
the first string of the tip ID, with the gene info being the last string in the tip ID seperated by an
underscore. 

Example: GCA_011056015_403 - GCA_011056015 is the genome name and 403 is the gene info.  

The genome name GCA_011056015 needs to be in the metadata file under column "ID" and the Phylum 
of the genome under a column named "Phylum"

'''


import sys, os.path
import pandas as pd
from collections import defaultdict
import glob
from scipy.stats import iqr
import numpy as np
from ete3 import Tree
import argparse
from skbio.diversity import alpha_diversity

parser = argparse.ArgumentParser(description="Calculating tree metrics")
parser.add_argument('-trees', help="Path of directory containing the trees, and nothing else - the trees must be rooted already.")
parser.add_argument('-metadata', help="Tab seperated text metadata file that needs to have at least the two columns 'ID' and 'Phylum'. See description")
parser.add_argument('-phyla_of_interest', help="The phyla you want sisters tallied for.  Same spelling as in metadata file")
parser.add_argument('-outname', help="Prefix for output files")
args = parser.parse_args()
arg_dict = vars(args)

tree_path = arg_dict['trees']
tree_paths = glob.glob(os.path.join(tree_path, "*"))

mdat = arg_dict['metadata']
mdata = pd.read_csv(mdat, sep = "\t")

poi = arg_dict['phyla_of_interest']

prefix = arg_dict['outname']
abs_path = os.path.abspath(os.getcwd()) 
out_pref = os.path.join(abs_path, prefix)


	
def sister_analysis(mdat, tree_paths, poi):

	phyla_l = list(set(mdata.Phylum.to_list()))
	phyla_dicts = defaultdict(list)
	for f in phyla_l:
		phy_list = []
		for row in mdata.itertuples():
			gen_id = row.ID
			phy = row.Phylum
			if phy is f:
				phy_list.append(gen_id)
			else:
				continue
		phyla_dicts[f] = phy_list
	
	phyla_sister_counts = {}
	
	for tree in tree_paths:
		(dir, file) = os.path.split(tree)
		t = Tree(tree)
		gc_id = file.rsplit('.',4)[0] #double check this
		
		print("calculating metrics for " + gc_id)
		
		genome_list=[]
		for leaf in t.iter_leaves():
			header = leaf.name
			splits = header.rsplit("_")[0:-1]
			nhead = '_'.join(splits)
			genome_list.append(nhead)
			for key, value in phyla_dicts.items():
				if nhead in value:
					leaf.add_feature("Phyla", key)
					
		phyla_list = []
		for genome in genome_list:
			for key, value in phyla_dicts.items():
				if genome in value:
					phyla_list.append(key)
					
		uniq_phyla = list(set(phyla_list))
		
		taxa_count = defaultdict(int)
		for node in t.get_monophyletic(values=[poi], target_attr="Phyla"):
			sisters = node.get_sisters()	
			for sister_node in sisters:
				for leaf in sister_node.iter_leaves():
					taxonomy = leaf.Phyla
					if taxonomy == poi: #skipping counting the same phyla as a sister of itself
						pass
					else:
						taxa_count[taxonomy] += 1
				
		
		phyla_sister_counts[gc_id] = taxa_count
					
	sister_df = pd.DataFrame.from_dict(phyla_sister_counts, orient='index')
	sister_df[sister_df > 0] = 1
	sister_df.loc['total_occurences_as_sister'] = sister_df.sum()
	
	return(sister_df)

def main():

	
	sister_df =  sister_analysis(mdat, tree_paths, poi)
	
	with open("%s_sister-analysis-occurence.txt" % (out_pref), "w") as outfile:
		sister_df.to_csv(outfile, index = True, header=True, sep='\t', na_rep='NA')
		print("Sister analysis results written to " + "%s_sister-analysis-occurence.txt" % (out_pref))
	
if __name__ == "__main__":
	main()	
	
