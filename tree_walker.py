import ete3
from ete3 import Tree
import pandas as pd

tree = Tree('test_small\\pastml\\test_named.nwk', format=1)
#print(tree)
dict_walk = {}
for node in tree.traverse("levelorder"):
    dict_walk[node.name] = list(node.children)




## Selecting one curtain state for each node (based on the highest probability)
df = pd.read_csv('test_small\\pastml\\marginal_probabilities.character_Country.model_F81.tab', sep='\t')
node_state = pd.DataFrame()
node_state['node'] = df['node']
df1 = df.drop('node', 1)
node_state['state'] = df1.idxmax(axis=1)

#print(node_state)
#print(dict_walk.keys())

to_rus = []
from_rus = []

for node_name in dict_walk:
    if len(dict_walk[node_name]) == 2:
        x=dict_walk[node_name][0].name
        y=dict_walk[node_name][1].name
        
        k = node_state[node_state['node']==node_name]['state'].values[0]
        l = node_state[node_state['node']==x]['state'].values[0]
        m = node_state[node_state['node']==y]['state'].values[0]
       
        #print(node_name, k, l, m)
       
        if 'RUS' in k:
            
            if not 'RUS' in l:
                #print(node_name+': '+k,x+': '+l)
                from_rus.append([node_name,x])
            if not 'RUS' in m:
                #print(node_name+': '+k,y+': '+m)
                from_rus.append([node_name,x])
        if not 'RUS' in k:
    
            if 'RUS' in l:
                #print(node_name+': '+k,x+': '+l)
                to_rus.append([node_name,x])
            if 'RUS' in m:
                #print(node_name+': '+k,y+': '+m)
                to_rus.append([node_name,x])


print(len(to_rus)+len(from_rus))
print(to_rus)
print(from_rus)