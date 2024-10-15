import sys, os.path
import pandas as pd
from collections import defaultdict
import glob
import numpy as np
import re
import glob
import time
import math
import argparse

parser = argparse.ArgumentParser(description="Create data files from GORG classifier results for downstream analysis.")
parser.add_argument('-i', help="Path to a directory containing the results files from running GORG Classifier")
args = parser.parse_args()
arg_dict = vars(args)
input_directory = arg_dict['i']

read_files = glob.glob(os.path.join(input_directory, "*")) # GORG classifier results

header= ["status", "sequence_id", "taxonomy_id", "length", "taxonomy_ids_lca", "sequence_ids_lca", "protein_sequence", "taxonomic_lineage", "prokka_EC_number", "prokka_product", "swissprot_gene", "swissprot_EC_number", "swissprot_eggNOG", "swissprot_KO", "swissprot_Pfam", "swissprot_CAZy", "swissprot_TIGRFAMs"]

abs_path = os.path.abspath(os.getcwd()) 

tax_dir = os.path.join(abs_path, "taxonomic_annotations") 
if not os.path.exists(tax_dir):
    os.makedirs(tax_dir)

print("summarizing taxonomy")

start_time = time.time()
for file in read_files:
    print(file)
    
    # extract the filename and modify it to use as the output file
    base=os.path.basename(file)
    base=base[:-7]
    out_n=(base+'_gorg_taxonomy_results.csv')
    out_name = os.path.join(tax_dir, out_n)
    
    chunksize = 10000
    counter=1
    check=1
    # read the file in chunks and group the reads by each gene/tax lineage for each chunk
    for chunk in pd.read_csv(file, sep='\t', chunksize=chunksize, skiprows=0, names=header):
        # tracker to keep track of the progress
        if counter % 10 == 0:
            print('the number of processed reads is',counter*chunksize)
            
        # read in the first chunk process it and store it as a summary df
        if check==1:
            # extract classified reads only
            DF_summary = chunk[(chunk['status'] == 'C')]
            if len(DF_summary)==0:
		
                
                check=1
                continue
                
            else:
                check=2


  # copy the df into a new df that I can use to generate noncompetitive counts of all the SAGs that a read was assigned to
		
                DF_explode=DF_summary.copy()
                DF_explode = DF_explode.fillna('')
                DF_explode=DF_explode.set_index(["status", "sequence_id", "taxonomy_id", "length", "taxonomy_ids_lca", "sequence_ids_lca", "protein_sequence", "taxonomic_lineage", "prokka_EC_number", "prokka_product", "swissprot_gene", "swissprot_EC_number", "swissprot_eggNOG", "swissprot_KO", "swissprot_Pfam", "swissprot_CAZy", "swissprot_TIGRFAMs"]).apply(lambda x: x.str.split(',').explode()).reset_index()
                DF_explode_summary=DF_explode.groupby(['taxonomic_lineage'], as_index=False, dropna=False)['status'].count()
                DF_explode_summary.rename(columns={'status':'number_of_reads'}, inplace=True)
                #DF_explode_summary[['Kingdom', 'Phylum', 'Class', 'Order', "Family", "Genus", "Species"]] = DF_explode_summary['taxonomic_lineage'].str.split('; ', n=6, expand=True)

 # read all other chunks into a temporary df and 
        else:
            # Same as above but store each successive chunk as a tempory_df
            tdf_summary = chunk[(chunk['status'] == 'C')]
            if len(tdf_summary)==0:
                continue
            
            
            else:
                # copy the df into a new df that I can use to generate noncompetitive counts of all the SAGs that a read was assigned to
                tDF_explode=tdf_summary.copy()

                tDF_explode=tDF_explode.set_index(["status", "sequence_id", "taxonomy_id", "length", "taxonomy_ids_lca", "sequence_ids_lca", "protein_sequence", "taxonomic_lineage", "prokka_EC_number", "prokka_product", "swissprot_gene", "swissprot_EC_number", "swissprot_eggNOG", "swissprot_KO", "swissprot_Pfam", "swissprot_CAZy", "swissprot_TIGRFAMs"]).apply(lambda x: x.str.split(',').explode()).reset_index()
                tDF_explode_summary=tDF_explode.groupby(['taxonomic_lineage'], as_index=False, dropna=False)['status'].count()
                tDF_explode_summary.rename(columns={'status':'number_of_reads'}, inplace=True)
                #tDF_explode_summary[['Kingdom', 'Phylum', 'Class', 'Order', "Family", "Genus", "Species"]] = tDF_explode_summary['taxonomic_lineage'].str.split('; ', n=6, expand=True)

                
                ### I need to merge the tmp explode df and sum the values and I think I might be done
                DF_explode_summary=DF_explode_summary.merge(tDF_explode_summary, on=['taxonomic_lineage'], how='outer')
                DF_explode_summary["number_of_reads_x"].fillna(0, inplace=True)
                DF_explode_summary["number_of_reads_y"].fillna(0, inplace=True)
        
                # add the read numbers, and drop the read number columns that were duplicated upon the merge
                DF_explode_summary['number_of_reads']=DF_explode_summary['number_of_reads_x']+DF_explode_summary['number_of_reads_y']
                DF_explode_summary.drop(columns=['number_of_reads_x', 'number_of_reads_y'], inplace=True)
                
        
        counter+=1
        
    DF_explode_summary[['Kingdom', 'Phylum', 'Class', 'Order', "Family", "Genus", "Species"]] = DF_explode_summary['taxonomic_lineage'].str.split('; ', n=6, expand=True)
    DF_explode_summary.drop(columns=['taxonomic_lineage'], inplace=True)
    DF_explode_summary.to_csv(out_name)
    

paths = glob.glob(os.path.join(tax_dir, "*")) # directory with files produced above 

tax_tab_dir = os.path.join(abs_path, "taxonomic_tables") 
if not os.path.exists(tax_tab_dir):
    os.makedirs(tax_tab_dir)
    
print("making taxonomy tables")

tax_list = ['Phylum', 'Class', 'Order', "Family", "Genus"]
for level in tax_list:

    master_dict = {}

    for txt in paths:

        (dir, file) = os.path.split(txt)
        sname = file.replace("_annotated_gorg_taxonomy_results.csv", "")

        df = pd.read_csv(txt, index_col= 0)
        df1 = df.groupby([level])['number_of_reads'].sum()
        dict1 = df1.to_dict()
        master_dict[sname] = dict1

    out_n=(level+'_count_table.txt')
    outname = os.path.join(tax_tab_dir, out_n)
    tax_table = pd.DataFrame.from_dict(master_dict, orient='index')
    tax_table = tax_table.fillna(0)
    tax_table.to_csv(outname, index = True, header=True, sep='\t')
 

## Summarizing functional annotations of the reads using the KO, eggNOG and Pfam annotations

func_dir = os.path.join(abs_path, "functional_annotations") 
if not os.path.exists(func_dir):
    os.makedirs(func_dir)
    
print("summarizing functional annotations")

## Kegg KO loop

ko_dir = os.path.join(func_dir, "KO") 
if not os.path.exists(ko_dir):
    os.makedirs(ko_dir)
    
print("KO annotations")

start_time = time.time()
for file in read_files:
    
    base=os.path.basename(file)
    base=base[:-7]
    print(base)
    out_n=(base+'_KO_results.csv')
    out_name = os.path.join(ko_dir, out_n)
    
    chunksize = 10000
    counter=1
    check=1

    for chunk in pd.read_csv(file, sep='\t', chunksize=chunksize, skiprows=0, names=header):
        
        if counter % 10 == 0:
            print('the number of processed reads is',counter*chunksize)
            
        if check==1:
    
            DF_summary = chunk[(chunk['status'] == 'C')]
            if len(DF_summary)==0:
                check=1
                continue
                
            else:
                check=2

                DF_explode=DF_summary.copy()
                DF_explode=DF_explode.set_index(["status", "sequence_id", "taxonomy_id", "length", "taxonomy_ids_lca", "sequence_ids_lca", "protein_sequence", "taxonomic_lineage", "prokka_EC_number", "prokka_product", "swissprot_gene", "swissprot_EC_number", "swissprot_eggNOG", "swissprot_KO", "swissprot_Pfam", "swissprot_CAZy", "swissprot_TIGRFAMs"]).apply(lambda x: x.str.split(',').explode()).reset_index()
                DF_explode_gene_annotations=DF_explode[['taxonomic_lineage']].copy()
                DF_explode['KO_ID'] =DF_explode['swissprot_KO'].str.split(",").str[0]
                DF_explode_summary=DF_explode.groupby(['KO_ID', "taxonomic_lineage"], as_index=False, dropna=False)['status'].count()
                DF_explode_summary.rename(columns={'status':'number_of_reads'}, inplace=True)

        else:
            
            tdf_summary = chunk[(chunk['status'] == 'C')]
            if len(tdf_summary)==0:
                continue
            
            else:
                
                tDF_explode=tdf_summary.copy()
                tDF_explode=tDF_explode.set_index(["status", "sequence_id", "taxonomy_id", "length", "taxonomy_ids_lca", "sequence_ids_lca", "protein_sequence", "taxonomic_lineage", "prokka_EC_number", "prokka_product", "swissprot_gene", "swissprot_EC_number", "swissprot_eggNOG", "swissprot_KO", "swissprot_Pfam", "swissprot_CAZy", "swissprot_TIGRFAMs"]).apply(lambda x: x.str.split(',').explode()).reset_index()
                tDF_explode_gene_annotations=tDF_explode[['taxonomic_lineage']].copy()
                
                DF_explode_gene_annotations=pd.concat([DF_explode_gene_annotations, tDF_explode_gene_annotations], ignore_index=True, axis=0)
                DF_explode_gene_annotations.drop_duplicates(inplace=True)
                
                tDF_explode['KO_ID'] =tDF_explode['swissprot_KO'].str.split(",").str[0]
                tDF_explode_summary=tDF_explode.groupby(['KO_ID', "taxonomic_lineage"], as_index=False, dropna=False)['status'].count()
                tDF_explode_summary.rename(columns={'status':'number_of_reads'}, inplace=True)

                DF_explode_summary=DF_explode_summary.merge(tDF_explode_summary, on=['KO_ID', "taxonomic_lineage"], how='outer')
                DF_explode_summary["number_of_reads_x"].fillna(0, inplace=True)
                DF_explode_summary["number_of_reads_y"].fillna(0, inplace=True)
        
                DF_explode_summary['number_of_reads']=DF_explode_summary['number_of_reads_x']+DF_explode_summary['number_of_reads_y']
                DF_explode_summary.drop(columns=['number_of_reads_x', 'number_of_reads_y'], inplace=True)
        
        counter+=1
        

    DF_explode_summary.to_csv(out_name)
    
## Pfam loop

pf_dir = os.path.join(func_dir, "Pfam") 
if not os.path.exists(pf_dir):
    os.makedirs(pf_dir)
    
print("Pfam Annotations")

for file in read_files:
    
    base=os.path.basename(file)
    base=base[:-7]
    print(base)
    out_n=(base+'_Pfam_results.csv')
    out_name = os.path.join(pf_dir, out_n)
    
    chunksize = 10000
    counter=1
    check=1
    
    for chunk in pd.read_csv(file, sep='\t', chunksize=chunksize, skiprows=0, names=header):
        
        if counter % 10 == 0:
            print('the number of processed reads is',counter*chunksize)
            
        if check==1:

            DF_summary = chunk[(chunk['status'] == 'C')]
            if len(DF_summary)==0:
                check=1
                continue
                
            else:
                check=2

                DF_explode=DF_summary.copy()
                DF_explode=DF_explode.set_index(["status", "sequence_id", "taxonomy_id", "length", "taxonomy_ids_lca", "sequence_ids_lca", "protein_sequence", "taxonomic_lineage", "prokka_EC_number", "prokka_product", "swissprot_gene", "swissprot_EC_number", "swissprot_eggNOG", "swissprot_KO", "swissprot_Pfam", "swissprot_CAZy", "swissprot_TIGRFAMs"]).apply(lambda x: x.str.split(',').explode()).reset_index()
                DF_explode_gene_annotations=DF_explode[['taxonomic_lineage']].copy()                
                DF_explode['Pfam'] =DF_explode['swissprot_Pfam'].str.split(",").str[0]
                DF_explode_summary=DF_explode.groupby(['Pfam', "taxonomic_lineage"], as_index=False, dropna=False)['status'].count()
                DF_explode_summary.rename(columns={'status':'number_of_reads'}, inplace=True)

        else:
            
            tdf_summary = chunk[(chunk['status'] == 'C')]
            if len(tdf_summary)==0:
                continue
            
            
            else:
                
                tDF_explode=tdf_summary.copy()
                tDF_explode=tDF_explode.set_index(["status", "sequence_id", "taxonomy_id", "length", "taxonomy_ids_lca", "sequence_ids_lca", "protein_sequence", "taxonomic_lineage", "prokka_EC_number", "prokka_product", "swissprot_gene", "swissprot_EC_number", "swissprot_eggNOG", "swissprot_KO", "swissprot_Pfam", "swissprot_CAZy", "swissprot_TIGRFAMs"]).apply(lambda x: x.str.split(',').explode()).reset_index()   
                tDF_explode_gene_annotations=tDF_explode[['taxonomic_lineage']].copy()
                
                DF_explode_gene_annotations=pd.concat([DF_explode_gene_annotations, tDF_explode_gene_annotations], ignore_index=True, axis=0)
                DF_explode_gene_annotations.drop_duplicates(inplace=True)
                
                tDF_explode['Pfam'] =tDF_explode['swissprot_Pfam'].str.split(",").str[0]
                tDF_explode_summary=tDF_explode.groupby(['Pfam', "taxonomic_lineage"], as_index=False, dropna=False)['status'].count()
                tDF_explode_summary.rename(columns={'status':'number_of_reads'}, inplace=True)

                DF_explode_summary=DF_explode_summary.merge(tDF_explode_summary, on=['Pfam', "taxonomic_lineage"], how='outer')
                DF_explode_summary["number_of_reads_x"].fillna(0, inplace=True)
                DF_explode_summary["number_of_reads_y"].fillna(0, inplace=True)
        
                DF_explode_summary['number_of_reads']=DF_explode_summary['number_of_reads_x']+DF_explode_summary['number_of_reads_y']
                DF_explode_summary.drop(columns=['number_of_reads_x', 'number_of_reads_y'], inplace=True)
        
        counter+=1
        


    DF_explode_summary.to_csv(out_name)
    
## eggNOG loop

og_dir = os.path.join(func_dir, "eggNOG") 
if not os.path.exists(og_dir):
    os.makedirs(og_dir)
    
print("eggNOG annotations")

for file in read_files:
    
    base=os.path.basename(file)
    base=base[:-7]
    print(base)
    out_n=(base+'_OG_results.csv')
    out_name = os.path.join(og_dir, out_n)
    
    chunksize = 10000
    counter=1
    check=1
    
    for chunk in pd.read_csv(file, sep='\t', chunksize=chunksize, skiprows=0, names=header):
        
        if counter % 10 == 0:
            print('the number of processed reads is',counter*chunksize)
            
        if check==1:
            
            DF_summary = chunk[(chunk['status'] == 'C')]
            if len(DF_summary)==0:
                check=1
                continue
                
            else:
                check=2

                DF_explode=DF_summary.copy()
                DF_explode=DF_explode.set_index(["status", "sequence_id", "taxonomy_id", "length", "taxonomy_ids_lca", "sequence_ids_lca", "protein_sequence", "taxonomic_lineage", "prokka_EC_number", "prokka_product", "swissprot_gene", "swissprot_EC_number", "swissprot_eggNOG", "swissprot_KO", "swissprot_Pfam", "swissprot_CAZy", "swissprot_TIGRFAMs"]).apply(lambda x: x.str.split(',').explode()).reset_index()
                DF_explode_gene_annotations=DF_explode[['taxonomic_lineage']].copy()
                DF_explode['eggNOG'] =DF_explode['swissprot_eggNOG'].str.split(",").str[0]
                DF_explode_summary=DF_explode.groupby(['eggNOG', "taxonomic_lineage"], as_index=False, dropna=False)['status'].count()
                DF_explode_summary.rename(columns={'status':'number_of_reads'}, inplace=True)

        else:
            
            tdf_summary = chunk[(chunk['status'] == 'C')]
            if len(tdf_summary)==0:
                continue
            
            
            else:
                
                tDF_explode=tdf_summary.copy()
                tDF_explode=tDF_explode.set_index(["status", "sequence_id", "taxonomy_id", "length", "taxonomy_ids_lca", "sequence_ids_lca", "protein_sequence", "taxonomic_lineage", "prokka_EC_number", "prokka_product", "swissprot_gene", "swissprot_EC_number", "swissprot_eggNOG", "swissprot_KO", "swissprot_Pfam", "swissprot_CAZy", "swissprot_TIGRFAMs"]).apply(lambda x: x.str.split(',').explode()).reset_index()
                tDF_explode_gene_annotations=tDF_explode[['taxonomic_lineage']].copy()
                
                DF_explode_gene_annotations=pd.concat([DF_explode_gene_annotations, tDF_explode_gene_annotations], ignore_index=True, axis=0)
                DF_explode_gene_annotations.drop_duplicates(inplace=True)
                
                tDF_explode['eggNOG'] =tDF_explode['swissprot_eggNOG'].str.split(",").str[0]
                tDF_explode_summary=tDF_explode.groupby(['eggNOG', "taxonomic_lineage"], as_index=False, dropna=False)['status'].count()
                tDF_explode_summary.rename(columns={'status':'number_of_reads'}, inplace=True)

                DF_explode_summary=DF_explode_summary.merge(tDF_explode_summary, on=['eggNOG', "taxonomic_lineage"], how='outer')
                DF_explode_summary["number_of_reads_x"].fillna(0, inplace=True)
                DF_explode_summary["number_of_reads_y"].fillna(0, inplace=True)
        
                DF_explode_summary['number_of_reads']=DF_explode_summary['number_of_reads_x']+DF_explode_summary['number_of_reads_y']
                DF_explode_summary.drop(columns=['number_of_reads_x', 'number_of_reads_y'], inplace=True)
        
        counter+=1
        


    DF_explode_summary.to_csv(out_name)

## Make count tables out of the functional annotation summary files produced above 
# pfam count table
print("Making functional count tables")

paths = glob.glob(os.path.join(pf_dir, "*"))
print("Pfam")
master_dict = {}
for txt in paths:

    (dir, file) = os.path.split(txt)
    sname = file.replace("_Pfam_results.csv", "")
 
    df = pd.read_csv(txt, index_col= 0)
    df = df.dropna(subset=['Pfam'])
    
    df1 = df.groupby(['Pfam'])['number_of_reads'].sum()
    dict1 = df1.to_dict()
    master_dict[sname] = dict1
    

pfam_table = pd.DataFrame.from_dict(master_dict, orient='index')
pfam_table = pfam_table.fillna(0)
out_name = os.path.join(func_dir, "pfam_count_table.txt")
pfam_table.to_csv(out_name, index = True, header=True, sep='\t')

# eggNOG count table
paths = glob.glob(os.path.join(og_dir, "*"))
print("eggNOG")
master_dict = {}

for txt in paths:

    (dir, file) = os.path.split(txt)
    sname = file.replace("_OG_results.csv", "")
 
    df = pd.read_csv(txt, index_col= 0)
    df = df.dropna(subset=['eggNOG'])
    
    df1 = df.groupby(['eggNOG'])['number_of_reads'].sum()
    dict1 = df1.to_dict()
    master_dict[sname] = dict1
    
og_table = pd.DataFrame.from_dict(master_dict, orient='index')
og_table = og_table.fillna(0)
out_name = os.path.join(func_dir, "OG_count_table.txt")
og_table.to_csv(out_name, index = True, header=True, sep='\t')

# KO count table
paths = glob.glob(os.path.join(ko_dir, "*"))
print("KO")
master_dict = {}

for txt in paths:

    (dir, file) = os.path.split(txt)
    sname = file.replace("_KO_results.csv", "")
 
    df = pd.read_csv(txt, index_col= 0)
    df = df.dropna(subset=['KO_ID'])
    
    df1 = df.groupby(['KO_ID'])['number_of_reads'].sum()
    dict1 = df1.to_dict()
    master_dict[sname] = dict1
    
ko_table = pd.DataFrame.from_dict(master_dict, orient='index')
ko_table = ko_table.fillna(0)
out_name = os.path.join(func_dir, "ko_count_table.txt")
ko_table.to_csv(out_name, index = True, header=True, sep='\t')

## make table of the taxonomic composition (Family level) of each functional category
# KO annotations
print("make taxonomic composition tables for each functional group")
wo_dir = os.path.join(abs_path, "function_taxonomy_composition/KO")
if not os.path.exists(wo_dir):
    os.makedirs(wo_dir)

paths = glob.glob(os.path.join(ko_dir, "*"))
for path in paths:
    
    (dir, file) = os.path.split(path)
    oname = file.replace("_annotated_KO_results.csv", "_KO_Family_composition.txt")
    outname = os.path.join(wo_dir, oname)
    
    if not os.path.exists(outname):
        
        df = pd.read_csv(path, index_col = 0)
        df[['Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species', 'NA']] = df['taxonomic_lineage'].str.split('; ', n = 8, expand=True)
        df.drop('NA', axis=1, inplace = True)
        df.drop('taxonomic_lineage', axis=1, inplace = True)

        df_summary= df.groupby(['KO_ID', "Family"], as_index=False, dropna=True)['number_of_reads'].sum()
        df_summary.sort_values(by=['number_of_reads'], ascending = False)

        ko_list = list(set(df_summary.KO_ID.to_list()))

        master_dict = defaultdict(list)
        for ko in ko_list:
            ko_dict = {}
            for tup in df_summary.itertuples():
                if tup.KO_ID is ko:
                    ko_dict[tup.Family] = tup.number_of_reads
    
            master_dict[ko] = ko_dict
        comp_dict = pd.DataFrame.from_dict(master_dict).fillna(0)
        comp_dict = comp_dict.div(comp_dict.sum(axis=0),axis=1)
        comp_dict = comp_dict.multiply(100)
        comp_dict.to_csv(outname, sep="\t", header = True)  
        
    else:
        continue
        
# Pfam annotations

wo_dir = os.path.join(abs_path, "function_taxonomy_composition/pfam")
if not os.path.exists(wo_dir):
    os.makedirs(wo_dir)
    
paths = glob.glob(os.path.join(pf_dir, "*"))
for path in paths:
    
    (dir, file) = os.path.split(path)
    oname = file.replace("_annotated_Pfam_results.csv", "_Pfam_Family_composition.txt")
    outname = os.path.join(wo_dir, oname)
    
    if not os.path.exists(outname):
        
        df = pd.read_csv(path, index_col = 0)
        df[['Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species', 'NA']] = df['taxonomic_lineage'].str.split('; ', n = 8, expand=True)
        df.drop('NA', axis=1, inplace = True)
        df.drop('taxonomic_lineage', axis=1, inplace = True)

        df_summary= df.groupby(['Pfam', "Family"], as_index=False, dropna=True)['number_of_reads'].sum()
        df_summary.sort_values(by=['number_of_reads'], ascending = False)

        pfam_list = list(set(df_summary.Pfam.to_list()))

        master_dict = defaultdict(list)
        for pfam in pfam_list:
            pfam_dict = {}
            for tup in df_summary.itertuples():
                if tup.Pfam is pfam:
                    pfam_dict[tup.Family] = tup.number_of_reads
    
            master_dict[pfam] = pfam_dict
        comp_dict = pd.DataFrame.from_dict(master_dict).fillna(0)
        comp_dict = comp_dict.div(comp_dict.sum(axis=0),axis=1)
        comp_dict = comp_dict.multiply(100)
        comp_dict.to_csv(outname, sep="\t", header = True)  
        
    else:
        continue
        


# eggNOG annotations

wo_dir = os.path.join(abs_path, "function_taxonomy_composition/eggNOG")
if not os.path.exists(wo_dir):
    os.makedirs(wo_dir)

paths = glob.glob(os.path.join(og_dir, "*"))
for path in paths:
    
    (dir, file) = os.path.split(path)
    oname = file.replace("_annotated_OG_results.csv", "_eggNOG_Family_composition.txt")
    outname = os.path.join(wo_dir, oname)
    
    if not os.path.exists(outname):
        
        df = pd.read_csv(path, index_col = 0)
        df[['Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species', 'NA']] = df['taxonomic_lineage'].str.split('; ', n = 8, expand=True)
        df.drop('NA', axis=1, inplace = True)
        df.drop('taxonomic_lineage', axis=1, inplace = True)

        df_summary= df.groupby(['eggNOG', "Family"], as_index=False, dropna=True)['number_of_reads'].sum()
        df_summary.sort_values(by=['number_of_reads'], ascending = False)

        og_list = list(set(df_summary.eggNOG.to_list()))

        master_dict = defaultdict(list)
        for og in og_list:
            og_dict = {}
            for tup in df_summary.itertuples():
                if tup.eggNOG is og:
                    og_dict[tup.Family] = tup.number_of_reads
    
            master_dict[og] = og_dict
        comp_dict = pd.DataFrame.from_dict(master_dict).fillna(0)
        comp_dict = comp_dict.div(comp_dict.sum(axis=0),axis=1)
        comp_dict = comp_dict.multiply(100)
        comp_dict.to_csv(outname, sep="\t", header = True)  
        
    else:
        continue
        
