import argparse
import pandas as pd
import numpy as np
import re
import platform
from Bio import SeqIO

# input
#   input_file_n - fasta-file with alignment
#   subtype_table_n - table with COMET results of subtyping sequences from input_file_n

# output
#   writes 2 alignments - the first one with all sequence; 
#       the sequences which were subtyped by COMET have new subtype in their name
#       the sequences which were not subtyped by COMET with enough support have NaN in addition to the original subtype in sequence name)
#   the second alignment contains only sequences which were not subtyped by COMET

def change_subtype(input_file_n, subtype_table_n):

    seqs = list(SeqIO.parse(input_file_n,'fasta')) # sequences from input file
    comet_table = pd.read_csv(subtype_table_n,"\t") # table with subtyping results from COMET
    
    # list for sequences with new names (with changed subtype)
    list_new_seqs = []
    # list for sequences which were not subtyped by COMET
    list_non = []
    
    clean_from = re.compile(r"[_\s,]")
    clean_to = "-"
    
    for index, row in comet_table.iterrows():
    
        for rec in seqs:
        
            if row['name'].startswith(rec.id[0:8]):
                recid = rec.id.split('_')
                subtype = recid[1]
                subtype_comet = re.sub(clean_from, clean_to, comet_table['subtype'][index])
                support = comet_table['bootstrap support'][index]

                if np.isnan(support):
                    rec.id = recid[0]+'_'+subtype+'-NaN_'+'_'.join(recid[2:])
                    rec.description = rec.id
                    list_new_seqs.append(rec)
                    list_non.append(rec)
                    break
                    
                else:
                    if subtype == subtype_comet:
                        list_new_seqs.append(rec)
                        break
                    else:
                        
                        rec.id = recid[0]+'_'+subtype_comet+'_'+'_'.join(recid[2:])
                        rec.description = rec.id
                        list_new_seqs.append(rec)
                        break
                seqs.remove(rec)
    print(list_new_seqs[0].id)

    output_subtyped = '.'.join(input_file_n.split('.'))[:-1] + "_comet_typed.fasta"
    output_not_subtyped = '.'.join(input_file_n.split('.'))[:-1] + "_comet_nottyped.fasta"

    print(output_subtyped)
    print(output_not_subtyped)
    
    SeqIO.write(list_new_seqs, open(output_subtyped, 'w'), 'fasta')
    SeqIO.write(list_non, open(output_not_subtyped, 'w'), 'fasta')


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    requiredNamed = parser.add_argument_group('required named arguments')
    requiredNamed.add_argument("--input_file", type=str,
                        help="name of fasta-file with alignment", required = True)
    requiredNamed.add_argument("--subtype_table", type=str,
                        help="name of table with results from COMET ", required = True)
    args = parser.parse_args()

    change_subtype(args.input_file, args.subtype_table)
        
#input_file_n = "D:\\MY_FILES\\DATA\\Lukashev\\HIV\\fasta\\aligned_curated\\fsu_foreign_pol_cds_2295-3300_950_2_aln_curated.fas"
#subtype_table = "D:\\MY_FILES\\DATA\\Lukashev\\HIV\\fasta\\aligned_curated\\2_comet.csv"