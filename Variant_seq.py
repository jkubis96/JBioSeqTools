import pandas as pd
import numpy as np

from tqdm import tqdm


enzymes = pd.read_excel('data/enzymes.xlsx')


for i in enzymes.index:
    enzymes['Sequence'][i] = ''.join(c.upper() for c in enzymes['Sequence'][i] if c.isalpha())



var_dic = {'N':['A', 'C', 'T', 'G'],'M':['A', 'C'],'R':['A', 'G'],'W':['A', 'T'],'S':['C', 'G'],'Y':['C', 'T'],'K':['G', 'T'],'V':['A', 'C', 'G'],'H':['A', 'C', 'T'],'D':['A', 'T', 'G'],'B':['C', 'T', 'G']}


sequence =  enzymes['Sequence'][i]

def silnia_i(n):
    s = 1
    for i in range(2, n+1):
        s *= i
    return s


def sub_var(sequence:str()):
    variation_num = []
    for element in sequence:
        if element in var_dic.keys():
            variation_num.append(len(var_dic[element]))
    
    variation_num = int(np.prod(variation_num))
    seq_list = [sequence]
    final_seq = []
    check = True

    while(check == True):
        for seq in seq_list:
            letter_list = []
            for letter in var_dic.keys():
                if letter in list(seq):
                    letter_list.append(letter)
            if len(letter_list) == 0:
                final_seq.append(seq)
            else:
                for nc in range(0, len(seq)):
                    if seq[nc] in var_dic.keys():
                        for j in var_dic[seq[nc]]:
                           tmp = list(seq)
                           tmp[nc] = j
                           tmp = "".join(tmp)
                           seq_list.append(tmp)
                           seq_list =  np.unique(seq_list).tolist()
            final_seq = np.unique(final_seq).tolist()
            if len(final_seq) == variation_num:
                check = False
                        
    return final_seq
    

        
           
new = sub_var(sequence)         
           


dictionary = {'restriction_place':[], 'sequence':[], "name":[]}
for n in tqdm(enzymes.index):
    tmp = sub_var(enzymes['Sequence'][n])
    for i in tmp:
        dictionary['restriction_place'].append(enzymes['Restriction_spot'][n])
        dictionary['sequence'].append(i)
        dictionary['name'].append(enzymes['Name'][n])
        
    
df = pd.DataFrame.from_dict(dictionary)

df.to_excel('restriction_enzymes.xlsx')
