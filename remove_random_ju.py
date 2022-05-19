import argparse
import matplotlib.pyplot as plt
import os
import random
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


def remove_random(fasta):
    '''
    Divides all sequences into groups by the first 5 characters in GenBank accession number
    Randomly removes k% sequences in groups which size exceed m. k and m are defined by user.
    Saves the resulting sequences in fasta format
    fasta - name of file with sequences in fasta_format
    Returns name of output file
    Input:
        fasta - path to file with sequences in fasta-format
    Output:
        fasta1 - path to output file
    '''

    #creates list of sequence names (the order of sequences is significant)
    
    records = list(SeqIO.parse(open(fasta), "fasta"))
    seq_names_ordered = []
    for rec in records:
        seq_names_ordered.append(rec.id)


    #dictionary with sequences
    record_dict = SeqIO.to_dict(records)

    #creating dictionary with names of groups (first five chars of GB ids) as keys and list of sequences' names in this
    #group as values
    GB_groups = dict()
    #this iteration includes ref_name
    for key in record_dict.keys():
            if key[:5] not in GB_groups:
                    GB_groups[key[:5]] = []
            GB_groups[key[:5]].append(key)

    #creating dictionary and list of number of seqs in each groups for hist
    GB_groups_num = dict()
    nums_for_hist = []
    for key in GB_groups.keys():
            GB_groups_num[key] = len(GB_groups[key])
            nums_for_hist.append(len(GB_groups[key]))

    #creating histogram
    #plt.figure(figsize=(15,10))
    plt.hist(nums_for_hist, bins = range(min(nums_for_hist),max(nums_for_hist),1))
    plt.title('Group size distribution', size = "25")
    plt.xlabel('Size of the group ', size = "20")
    plt.ylabel('Number of groups', size = "20")
    plt.show()

    print("Type in the maximal group size; then press Enter")
    max_size = int(input())
    print("Type in the percent of sequences to remove from the groups with size larger than adjusted; then press Enter")
    per_of_seq = float(input())


    # creating list with names of sequences in final sample
    final_sample_ids = []
    num = 0
    for value in GB_groups.values():
            if len(value) <= max_size:
                num = num + len(value)
                for i in value:
                        final_sample_ids.append(i)
            else:
                
                random_sample = random.sample(value, max(1,int(len(value)*per_of_seq/100)))
                num = num + (len(random_sample))
                for seq in random_sample:
                    final_sample_ids.append(seq)
    print('Number of sequences in resulting alignment {}'.format(num))

    #creating list with final sequences
    final_seq_list = []
    for name in seq_names_ordered:
        if name in final_sample_ids:
            #print(record)
            final_seq_list.append(record_dict[name])

    #writing sequences into file
    fasta1 =  os.path.splitext(fasta)[0] + "_random_" + str(max_size)+ '_' + str(per_of_seq) + ".fasta"
    SeqIO.write(final_seq_list, fasta1, "fasta")
    print('done')
    return(fasta1)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-input", "--input_file", type=str,
                        help="Input file in fasta format", required=True)

    args = parser.parse_args()

    remove_random(args.input_file)