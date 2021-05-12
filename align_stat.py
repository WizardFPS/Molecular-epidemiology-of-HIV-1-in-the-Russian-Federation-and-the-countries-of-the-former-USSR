from collections import Counter
from matplotlib import pyplot as plt
import pandas as pd
import seaborn
seaborn.set()

names = ['01-AE', '02-AG', '03-AB', '06-cpx', '11-cpx', '63-02A1', 'A1', 'B', 'C', 'F1', 'G']

def get_data(name):
    align = name+'.fasta'
    country_lst = []
    with open(align, 'r') as al:
        for line in al:
            if line.startswith('>'):
                c_item = line[1:].split('_')[4]
                if ':' in c_item:
                    c_item = c_item.replace(':', '-')
                country_lst.append(c_item)
    dic_t = dict(Counter(country_lst))
    df = pd.DataFrame(dic_t.items(), columns = ['Country', 'N'])
    
    return(df)
   


def get_colors(df, name):
    tree = name + '_2300-3300'
    cd = {}
    cd = cd.fromkeys(df['Country'])
    with open(tree, 'r') as tr:
        for line in tr:
            for el in cd.keys():
                if el in line and not cd[el]:
                    try:
                        cd[el] ='#'+(line.split('#')[1])[:-2]
                    except:
                        cd[el] = '#000000'
    return(cd)

def lable_coord(n): #Координаты Total и FSU, размер шрифта числовых подписей столбцов
    if n >=1600:
        a, b = n, n-200
        f_s = 7
        return(a, b, f_s)
    if n >=900:
        a, b = n, n-80
        f_s = 12
        return(a, b, f_s)
    if n >= 200:
        a, b = n, n-20
        f_s = 12
        return(a, b, f_s)
    if n >= 50:
        a, b = n, n-10
        f_s = 12
        return(a, b, f_s)
    if n >= 10:
        a, b = n, n-3
        f_s = 12
        return(a, b, f_s)
    else:
        a, b = n, n-1
        f_s = 14
        return(a, b, f_s)
    
def count_FSU(name): #Подсчёт количества FSU последовательностей
    if name in ['03-AB','63-02A1', 'A1']:
        c = str(sum(colored_df[colored_df.Color != '#000000']['N']))
    else:
        c = str(sum(colored_df[colored_df.Color == '#0000ff']['N']))
    return(c)
for name in ['03-AB', '63-02A1', 'A1']:
    print(get_data(name)) 
    
  
for name in names:
    df = get_data(name)   
    df['Color'] = get_colors(df, name).values()

    colored_df = df[df.Color != '#000000'].reset_index(drop=True)
    
    colored_df = colored_df.append({'Country': 'Others', 'N': sum(df[df.Color == '#000000']['N']), 'Color': '#000000'}, ignore_index=True)
    print(colored_df)
    
    fig, ax = plt.subplots(figsize = (10,5))
    ax.bar(colored_df['Country'], colored_df['N'], width=0.4, color = colored_df['Color'])
    
    a, b, f_s = lable_coord(max(colored_df['N']))
    
    for index,data in enumerate(colored_df['N']):
        plt.text(x=index , y = data+1 , s=f"{data}" , fontdict=dict(fontsize=f_s))
    plt.tight_layout()
    plt.xticks(rotation=90)
    
    plt.text(len(colored_df), a, s = 'Total: ' + str(sum(colored_df['N'])))
    plt.text(len(colored_df), b, s = 'FSU: ' + count_FSU(name))

    plt.savefig(name+'_stat.png', dpi=700, format='png', bbox_inches='tight')