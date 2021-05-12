import argparse
import pandas as pd
from Bio import SeqIO


def fetch_foreign_blast(blast_out, fsu_fasta_n, foreign_fasta_n, output_f_name):

    blasto = pd.read_csv(blast_out,"\t") # blastoutput, format = 6

    current_seqid=  '' # seq id for current query seq
    cur_list_foreign_ids = [] # seq id of foreign sequences for current seq
    list_foreign_ids = [] # list of all foreign sequences

    for index, row in blasto.iterrows():

        if current_seqid != row.qseqid: #new query seq
            current_seqid = row.qseqid

            if cur_list_foreign_ids !=[]:

                f=0 # first foreign seq id to add to
                # cheking list of foreign seqs for previous query seq
                for i in range(len(cur_list_foreign_ids)): 
                    c = cur_list_foreign_ids[i].split('_')[-2] #country
                    y = cur_list_foreign_ids[i].split('_')[-1].split(':')[0] #collection year

                    #adds country and year of the first seq only if they are known
                    if (f==0) and (c!='NA') and (y!='NA'): 
                        countries = [cur_list_foreign_ids[i].split('_')[-2]]
                        years =  [cur_list_foreign_ids[i].split('_')[-1]]
                        list_foreign_ids.append(cur_list_foreign_ids[i])
                        f=1 
                        k=1
                        continue

                    if f==1:
                        #adds sseqid if country or collection year hasn't observed in previous sseqids
                        #allows unknown collection year
                        if ((c!='NA') and ((c not in countries) or (y not in years))):
                            list_foreign_ids.append(cur_list_foreign_ids[i])
                            countries.append(c)
                            if y!='NA':
                                years.append(y)
                            k+=1
                            if k>=7:
                                break

            #adds first cur_list_foreign_ids
            cur_list_foreign_ids = []
            if float(row.qcovs) > 95:
                cur_list_foreign_ids.append(row.sseqid)
        else:
            #adds sseqid to cur_list_foreign_ids
            if float(row.qcovs) > 95:
                cur_list_foreign_ids.append(row.sseqid)
    #deletes duplicates
    list_foreign_ids = set(list_foreign_ids)
    foreign_seq = SeqIO.parse(open(foreign_fasta_n),'fasta')
    list_foreign_recs = []
    for record in foreign_seq:
        if record.id in list_foreign_ids:
            list_foreign_recs.append(record)
    #print(len(list_foreign_recs))
    #print(len(list_foreign_ids))

    
    SeqIO.write(SeqIO.parse(open(fsu_fasta_n),'fasta'), open(output_f_name, 'w'), 'fasta')
    SeqIO.write(list_foreign_recs, open(output_f_name, 'a'), 'fasta')

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    
    requiredNamed = parser.add_argument_group('required named arguments')
    
    requiredNamed.add_argument("--q_fasta_file", type=str,
                        help="name of fasta-file with sequences used as query", required=True)
    requiredNamed.add_argument("--db_fasta_file", type=str,
                        help="name of fasta-file with sequences used to create the local database", required=True)
    requiredNamed.add_argument("--blast_output", type=str,
                        help="name of file with blast output (fmt 6)", required=True)
    requiredNamed.add_argument("--output_name", type=str,
                        help="name of output file", required=True)
    args = parser.parse_args()

    fetch_foreign_blast(args.blast_output, args.q_fasta_file, args.db_fasta_file, args.output_name)

#blast_out = "D:/FBB/kursa4_2/#5/blastn_output"
#foreign_fasta_n = "D:/FBB/kursa4_2/#5/foreign_align.fasta"
#fsu_fasta_n = "D:/FBB/kursa4_2/#5/USSR_align2300-3300.fasta"
#output_f_name = "D:/FBB/kursa4_2/#5/fsu_foreign_pol_cds_2300-3300_950_2.fasta"
#fetch_foreign_blast(blast_out, fsu_fasta_n, foreign_fasta_n, output_f_name)
    
    
#python fetch_foreign.py --q_fasta_file USSR_align2300-3300.fasta --db_fasta_file foreign_align.fasta --blast_output blastn_output --output_name fetch_output