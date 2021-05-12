# названия файлов с уже отсортированными accessions 
get_from = ['Africa_2300-3300', 'Asia_2300-3300', 'Canada_2300-3300', 'Caribbean_2300-3300', 'CentralAmerica_2300-3300', 'Europe_2300-3300', 'MiddleEast_2300-3300', 'Oceania_2300-3300', 'Southamerica_2300-3300', 'USSR_2300-3300', 'UnitedStates_2300-3300']
# названия файлов с исходными выравниваниями
align_list = ['Africa.fasta', 'Asia.fasta', 'Canada.fasta', 'Caribbean.fasta', 'CentralAmerica.fasta', 'Europe.fasta', 'MiddleEast.fasta', 'Oceania.fasta', 'SouthAmerica.fasta', 'USSR.fasta', 'UnitedStates.fasta']

# функция создаёт список из accessions
def get_acc(file):  #
    output = []
    with open(file, 'r') as in_file:
        for line in in_file:
            output.append(line.split('_')[0])
    return(output)

# функция достаёт нужные последовательности из исходного файла
def get_seq(file, acc_lst):
    output_name = file.split('.')[0]+'_align2300-3300.fasta'
    with open(file, 'r') as in_file, open(output_name, 'w') as out_file:
        k = False
        for line in in_file:
            if k:
                out_file.write(line)
                k = False
            if line.split('_')[0][1:] in acc_lst:
                out_file.write(line)
                k = True
            else:
                pass
                
        return(out_file)           

# запускаем обе функции для каждой пары файлов (исходное выравнивание и отсортированные accessions)
for i in range(len(get_from)):
    acc_lst = get_acc(get_from[i])
    get_seq(align_list[i], acc_lst)
            