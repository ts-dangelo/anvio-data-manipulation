'''

This runs:

		1.	anvi-script-reformat-fasta
			
		2.	anvi-gen-contigs-databas
			
		3.	anvi-run-kegg-kofams (optional with --kegg TRUE)
			
	on an input directory of nucleotide fasta files.  It then merges them into a storage database using:

			anvi-gen-genomes-storage
			
	After that it runs the ANVIO pan-genomics pipeline via:
		
			anvi-pan-genome
			
	This will perform an all-to-all blastp via Diamond and use the similarity scores to create gene-clusters using MCL clustering (parameters flags -minbit -mcl-inflation can be set)
	
			
'''


def clean_up_fasta(fasta_paths, wo_dir):
	
	import os
	import subprocess
	import glob
	
	fod = "clean_fastas"
	fod_dir = os.path.join(wo_dir, fod) 
	if not os.path.exists(fod_dir):
		os.makedirs(fod_dir)
	
	
	for fasta in fasta_paths:
		(dir, file) = os.path.split(fasta)
		fname = file
		clean_out = os.path.join(fod_dir, fname)
		
		print ("Running anvi-script-reformat-fasta on "  + fname)		
		subprocess.call(["anvi-script-reformat-fasta", fasta, "-o", clean_out, "--simplify-names", "--seq-type", "NT"])
			
	cleaned = glob.glob(os.path.join(fod_dir, "*"))	
	return(cleaned)
	

def	gen_contig_dbs(cleaned_fastas, wo_dir):
	
	import os
	import subprocess
	import glob	
	
	external = "external_genomes.txt"
	extfile = os.path.join(wo_dir, external) 
	external_path = open(extfile, "w") 
	external_path.write("name" + "\t" "contigs_db_path" + "\n")


	for fasta in cleaned_fastas: 
		print(fasta)
		(dir, file) = os.path.split(fasta)
		genome = file.rsplit('.',1)[0]
		file_loc = fasta.rsplit('.',1)[0]
		
		db_nm = "dbs"
		db_dir = os.path.join(wo_dir, db_nm) 
		if not os.path.exists(db_dir):
			os.makedirs(db_dir)
			
		db_loc = os.path.join(db_dir, "%s.db" % (genome))
		external_path.write(genome + "\t" + db_loc+ "\n")
		
		print ("Generationg Contig DB for: "  + genome)	
		subprocess.call(["anvi-gen-contigs-database", "-f", fasta,  "-n", "%s_DB" % (genome), "-o", db_loc]) 
	
	external_path.close()
	
	return(db_dir, extfile)
	

def	annotate_kegg(db_dir, wo_dir):

	import os
	import subprocess
	import glob	
	
	dbs = glob.glob(os.path.join(db_dir, "*"))	
	for db in dbs:
		print("Running KOFAMSCAN on " + db) 
		subprocess.call(["anvi-run-kegg-kofams", "-c", db, "-T", "8"]) 
		
	kao = "kegg_annotations"
	kao_dir = os.path.join(wo_dir, kao)
	if not os.path.exists(kao_dir):
		os.makedirs(kao_dir)
	
	for db in dbs:
		(dir, file) = os.path.split(db)
		genome = file.rsplit('.',1)[0]
		db_ao = "%s_kegg.txt" % (genome)
		db_aout = os.path.join(kao_dir, db_ao)
		print("Annotation information for " + genome + " located here: " + db_aout)
		subprocess.call(["anvi-export-functions", "-c", db, "-o", db_aout])
	
	return
	
def merge_and_cluster(extfile, pan_name, minbit, infl, wo_dir):

	import os
	import subprocess
	import glob
	
	merge_fn = "merge_stdout.txt"
	merge_fnd = os.path.join(wo_dir, merge_fn)
	merge = open(merge_fnd, "w")
	
	gen_db = "%s-GENOMES.db" % (pan_name)
	gen_out = os.path.join(wo_dir, gen_db)
	
	print ("Merging Contig Profiles to create  "  + gen_db)
		
	subprocess.call(["anvi-gen-genomes-storage", "-e", extfile, "-o", gen_out], stdout=merge)
	
	merge.close()

	pan_fn = "pan_stdout.txt"
	pan_fnd = os.path.join(wo_dir, pan_fn)
	pan = open(pan_fnd, "w")
	
	cluster = "cluster_out"
	cluster_out = os.path.join(wo_dir, cluster) 
	if not os.path.exists(cluster_out):
		os.makedirs(cluster_out)
	
	pan_db = pan_name + "-PAN.db"
	pan_loc = os.path.join(cluster_out, pan_db)	

	print ("Running all-vs-all blast and MCL clustering on ") 
	subprocess.call(["anvi-pan-genome", "-g", gen_out, "--project-name", pan_name, "-o", cluster_out, "--num-threads", "8", "--minbit", minbit, "--mcl-inflation", infl], stdout=pan)
	
	pan.close()
		
	return(pan_loc, gen_out)
	
	
def	get_useful_info(pan_loc, gen_out, wo_dir, fasta_names):
	
	import os
	import subprocess
		
	gcs = "all_gene_clusters.faa"
	gc_faa = os.path.join(wo_dir, gcs)
		
	print ("Exporting gene-cluster sequences from "  + pan_loc)	
	subprocess.call(["anvi-get-sequences-for-gene-clusters", "-g", gen_out,  "-p", pan_loc,  "-o", gc_faa, "--min-num-genomes-gene-cluster-occurs", "2"])

	sum_files = "pangenome_summary_files"
	sum_loc = os.path.join(wo_dir, sum_files)	

	print ("Summarizing data and exporting gene-cluster annotations")
	subprocess.call(["anvi-script-add-default-collection", "-c", gen_out, "-p", pan_loc, "-C", "all_gcs", "-b",  "GC_collection"])
	subprocess.call(["anvi-summarize", "-p", pan_loc, "-g", gen_out, "-C", "all_gcs", "-o", sum_loc, "--quick-summary"])
                 
		
def main():

	import glob
	import argparse
	import pandas as pd
	import subprocess
	import os,shutil
	import os.path
	from collections import defaultdict
	
	parser = argparse.ArgumentParser(description="This wrapper script runs ANVIO scripts to automate pangenome analysis. It will: clean up fastas, translate to CDS, annotate CDS with KOFAMSCAN, merge genomes into genome db, cluster CDS to functional gene-clusters, export gene cluster sequences, summarize data.")
	parser.add_argument('-i', help="Path file folder containing genomes in nucelotide fasta format. Heed anvio's file name formatting requirements")
	parser.add_argument('-out_dir', help="Directory for output files. Doesn't need to exist already")
	parser.add_argument('-pandb_name', help="Name for 'Pangenome Database' in ANVIO")
	parser.add_argument('--kegg', help="TRUE to do KEGG annotation of CDS in each genome, default will not annote CDS")
	parser.add_argument('-minbit', help="Minbit value to pass to ANVIO pangenome blast results")
	parser.add_argument('-mcl_inflation', help="MCL inflation parameter to pass to Anvio")
	args = parser.parse_args()
	arg_dict = vars(args)
	
	pan_name = arg_dict['pandb_name']
		
	fpath = arg_dict['i']
	fasta_paths = glob.glob(os.path.join(fpath, "*"))
	
	fasta_names = []

	for fasta in fasta_paths:
		(dir, file) = os.path.split(fasta)
		fname = file
		fasta_names.append(fname)

	o_dir = arg_dict['out_dir'] #make the output/working dir.
	abs_path = os.path.abspath(os.getcwd()) 
	wo_dir = os.path.join(abs_path, o_dir) 
	if not os.path.exists(wo_dir):
		os.makedirs(wo_dir)
		
	minbit = arg_dict['minbit']
	infl = arg_dict['mcl_inflation']
	 
	cleaned_fastas = clean_up_fasta(fasta_paths, wo_dir)			
	db_dir, extfile = gen_contig_dbs(cleaned_fastas, wo_dir)
	
	ko_opt = arg_dict['kegg']
	if ko_opt == 'TRUE':
	
		annotate_kegg(db_dir, wo_dir)
		pan_loc, gen_out = merge_and_cluster(extfile, pan_name, minbit, infl, wo_dir)
		get_useful_info(pan_loc, gen_out, wo_dir, fasta_names)

	else:
		
		pan_loc, gen_out = merge_and_cluster(extfile, pan_name, minbit, infl, wo_dir)
		get_useful_info(pan_loc, gen_out, wo_dir, fasta_names)

	
if __name__ == "__main__":
	main()	
