## Most code within definitions excised and modified from: https://github.com/merenlab/anvio/blob/master/anvio/panops.py

def gen_mcl_input(blastall_results, min_percent_id, minbit_param, out_pref):

	import os
	import glob
       
	all_ids = set([])

        
	mapping = [str, str, float, int, int, int, int, int, int, int, float, float] # mapping for the fields in the blast output
	self_bit_scores = {}
	line_no = 1
    
	for line in open(blastall_results):
		fields = line.strip().split('\t')
		    	
		try:
			query_id, subject_id, perc_id, aln_length, mismatches, gaps, q_start, q_end, s_start, s_end, e_val, bit_score = [mapping[i](fields[i]) for i in range(0, len(mapping))]
			
		except Exception as e:
        
			raise ConfigError("Something went wrong while processing the blastall output file in line %d. "
                                   "Here is the error from the uppoer management: '''%s'''" % (line_no, e))
		line_no += 1
		all_ids.add(query_id)
		all_ids.add(subject_id)
		
		if query_id == subject_id:
			self_bit_scores[query_id] = bit_score
			
	ids_without_self_search = all_ids - set(self_bit_scores.keys())
	
	### Heuristic code removed here. Instead sequences without self bitscore will be removed prior to MCL and IDs printed to summary file.  
    
	abs_path = os.path.abspath(os.getcwd())
	mcl_input_file_path = os.path.join(abs_path, '%s-mcl-input.txt' % (out_pref))
	mcl_input = open(mcl_input_file_path, 'w')
	
	line_no = 1
	num_edges_stored = 0
	num_minbit_filtered = 0
	num_blast_filtered = 0
	
	for line in open(blastall_results):
		fields = line.strip().split('\t')
		
		query_id, subject_id, perc_id, aln_length, mismatches, gaps, q_start, q_end, s_start, s_end, e_val, bit_score = [mapping[i](fields[i]) for i in range(0, len(mapping))]
		
		if query_id in ids_without_self_search: # removing seqs that don't have self bit score instead of doing hueristic
			continue
		if subject_id in ids_without_self_search:
			continue
	
           
		if min_percent_id is not None:
			id_thresh = float(min_percent_id)
			if perc_id < id_thresh:
				num_blast_filtered += 1
				continue

    
		minbit = bit_score / min(self_bit_scores[query_id], self_bit_scores[subject_id])
		thresh = float(minbit_param)
		if minbit < thresh:
			num_minbit_filtered += 1
			continue
		
		mcl_input.write('%s\t%s\t%f\n' % (query_id, subject_id, perc_id / 100.0))
		num_edges_stored += 1

       
		
	mcl_input.close()
	
	return (mcl_input_file_path, num_edges_stored, num_minbit_filtered, num_blast_filtered, ids_without_self_search)
	

def cluster(mcl_input_file_path, inf, out_pref):
	
	import os
	import subprocess
	
	abs_path = os.path.abspath(os.getcwd())
	clusters_file_path = os.path.join(abs_path, '%s-mcl-cluster-output.txt' % (out_pref))
	mcl_stdout = os.path.join(abs_path, '%s-mcl-stdout.txt' % (out_pref))
	
	with open(mcl_stdout, 'w') as f:
		print("Clustering filtered blast results using MCL, inflation parameter set to " + inf)	
		subprocess.call(['mcl', mcl_input_file_path, '--abc', '-I', inf, '-o', clusters_file_path, '-te', '8'], stderr=f)
	
	return clusters_file_path
        
def get_clusters_dict(clusters_file_path):
        

	clusters = 0
	singletons = 0
	  
	for line in open(clusters_file_path).readlines():
	
		fields = line.rsplit("\t")
		
		if len(fields) > 1:
		
			clusters += 1
			
		if len(fields) == 1:
		
			singletons += 1
				
	num_clusters = clusters

	return (num_clusters, singletons)

        
def main():

	import glob
	import argparse
	import os
	import os.path
	from collections import defaultdict
	
	parser = argparse.ArgumentParser(description="Python definitions excised from https://github.com/merenlab/anvio/blob/master/anvio/panops.py to take Diamond or Blast results and do downstream Minbit filtering and MCL clustering. This was taken out of the anvio pangenomics workflow so that the influence of minbit and inflation parameters on cluster numbers could be investigated outside of the full workflow (to avoid having to do large blast-alls over and over again). You can't actually access gene sequences or clusters from the data produced here. It's only useful for investigating the graph clustering output. Will produce the MCL output files and a summary of Minbit and MCL filtering statistics")
	parser.add_argument('--blast', help="Output from all-vs-all Blast or Diamond output from anvio pangeomics workflow (or other source), Optional if you're passing previous minbit filtered results")
	parser.add_argument('--min_percent_id', help="Minimum percent identity between two sequences for incluesion in MCL clustering, Optional if you are passing previous minbit filtered results. As the anvio developers suggest, just use Minbit, don't use this")
	parser.add_argument('--minbit', help="Minbit parameter for filtering blast results, Optional if you are passing previous minbit filtered results")
	parser.add_argument('-inflation', help="Inflation parameter to pass to MCL")
	parser.add_argument('--minbit_out', help="Optional, minbit filtering output if you already have it")
	parser.add_argument('-out_pref', help="prefix for MCL output and log files")
	
	
	args = parser.parse_args()
	arg_dict = vars(args)
	
	minbit_param = arg_dict['minbit']
	min_percent_id = arg_dict['min_percent_id']
	blastall_results = arg_dict['blast']
	if blastall_results is not None:
		num_edges_raw = sum(1 for line in open(blastall_results))
	inf = arg_dict['inflation']
	minbit_out = arg_dict['minbit_out']
	out_pref = arg_dict['out_pref']
	
	abs_path = os.path.abspath(os.getcwd())
	res_file = os.path.join(abs_path, '%s-results-summary.txt' % (out_pref))
	
	with open(res_file, 'w') as rf:
	
		if minbit_out is not None:
		
			clusters_file_path = cluster(minbit_out, inf, out_pref)
			num_clusters, singletons = get_clusters_dict(clusters_file_path)
			num_edges_stored = sum(1 for line in open(minbit_out))
		
			rf.write("Connections passed to MCL after Minbit filtering \t " + str(num_edges_stored) + "\n")
			rf.write("Number of Gene (non-singleton) clusters identified by MCL \t " + str(num_clusters) + "\n")
			rf.write("Number of singleton clusters identified by MCL \t " + str(singletons) + "\n")
			
				
		else:
		
			mcl_input_file_path, num_edges_stored, num_minbit_filtered, num_blast_filtered, ids_without_self_search = gen_mcl_input(blastall_results, min_percent_id, minbit_param, out_pref)
			clusters_file_path = cluster(mcl_input_file_path, inf, out_pref)
			num_clusters, singletons = get_clusters_dict(clusters_file_path) 
		
			rf.write("Number of raw blast connections \t " + str(num_edges_raw) + "\n")
			rf.write("Weak connections filtered by Minbit \t " + str(num_minbit_filtered) + "\n")
			rf.write("Weak connections filtered by Percent ID \t " + str(num_blast_filtered) + "\n")
			rf.write("Connections passed to MCL after Minbit filtering \t " + str(num_edges_stored) + "\n")
			rf.write("Number of Gene (non-singleton) clusters identified by MCL \t " + str(num_clusters) + "\n")
			rf.write("Number of singleton clusters identified by MCL \t " + str(singletons) + "\n")
			rf.write("Sequences removed before MCL because they did not produce a self bitscore = " + str(ids_without_self_search))
		
			
if __name__ == "__main__":
	main()	      
 
