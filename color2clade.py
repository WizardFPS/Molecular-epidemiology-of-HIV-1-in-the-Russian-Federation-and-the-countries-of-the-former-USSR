import argparse
import re
import platform
from Bio import SeqIO


#dictionary with correspondence between colors and clades
color_dict = {
'#999999': '06-cpx',
'#ff0033':'A1',
'#990000':'A1_strange',
'#0000ff':'B',
'#33ff33':'C',
'#cccc00':'G',
'#ff9933':'D',
'#cc00cc':'11-cpx',
'#00ffff':'03-AB',
'#999900':'18-cpx',
'#ffcccc':'63-02A1',
'#ccccff':'01-AE',
'#ff00ff':'F1',
'#009999':'02-AG',
#'#ffff33':'35-AD'
}

#input:
#   fasta_file_n - fast-file with alignment
#   tree_file_n - tree-file in nexus format with coloured leaves
#   output_dir - output directory
#output:
#   saves sequences from fasta_file_n to new alignments according color of each sequence in the tree
#   the names of alignments are defined in color_dict
#   new alignments will saved in output_dir
def color2clade(tree_file_n, fasta_file_n, output_dir):

    tree_file = open(tree_file_n) #file with tree
    tl_fl = 0 #taxlabels have began
    names_dict = {} #names_dict[seq_name] = colour
    for line in tree_file:
        if re.match('\ttaxlabels',line):
            tl_fl = 1
            continue
        if tl_fl == 1:
            if line == ';\n':
                tl_fl = 0
            else:
                try:
                    cur_color = re.search(r'#[0-9a-z]+', line).group()
                    
                except:
                    #print(line)
                    cur_color = '#000000'         
                if cur_color in color_dict.keys():
                        #print(cur_color)
                        seq_name = re.search(r"[A-Za-z0-9_\-\/]+",line).group().split('_')[0]
                        names_dict[seq_name] = color_dict[cur_color]
                
                        
    tree_file.close()
    #print(names_dict.keys())


    fasta_seq = SeqIO.parse(open(fasta_file_n),'fasta') #input fasta-file
    seq_dict = {}
    seqs_other = []
    for record in fasta_seq:
        #print(record)
        rec_name = record.id.split('_')[0]
        if rec_name in names_dict.keys():
            if names_dict[rec_name] not in seq_dict.keys():
                seq_dict[names_dict[rec_name]] = []
            seq_dict[names_dict[rec_name]].append(record)
        else:
            seqs_other.append(record)
    fasta_seq.close()


    
    for st_group in seq_dict:
        file_name = output_dir + st_group+'.fasta'
        SeqIO.write(seq_dict[st_group], open(file_name, "w"), "fasta")
    SeqIO.write(seqs_other, open(output_dir + 'other.fasta', "w"), "fasta")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    
    requiredNamed = parser.add_argument_group('required named arguments')
    
    requiredNamed.add_argument("--tree_file", type=str,
                        help="name of tree-file in nexus format", required = True)
    requiredNamed.add_argument("--fasta_file", type=str,
                        help="name of fasta-file with alignment", required = True)
    parser.add_argument("--output_dir", type=str,
                        help="output directory; \
                        the files will be saved in the directory of input fasta alignment if not defined")
    args = parser.parse_args()


    if not args.output_dir:
        print("The output files will be saved in the directory of input fasta alignment")
        if  platform.system() == 'Windows':
            output_dir = "\\".join(args.fasta_file.split("\\")[:-1]) + "\\"
        else:
            output_dir = "/".join(args.fasta_file.split("/")[:-1]) + "/"
        color2clade(args.tree_file, args.fasta_file, output_dir)
    else:
        color2clade(args.tree_file, args.fasta_file, args.output_dir)

#tree_file_n = 'D:/MY_FILES/DATA/Lukashev/HIV/trees/tree_yu_n2_sector.tree'
#fasta_file_n = 'D:/MY_FILES/DATA/Lukashev/HIV/fasta/final/fsu_foreign_pol_aln_man-st_ed.fasta'
#output_dir = 'D:/MY_FILES/DATA/Lukashev/HIV/fasta/final/clades/'