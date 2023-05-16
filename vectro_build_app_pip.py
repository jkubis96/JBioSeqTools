# pip install JBioSeqTools==1.1.3

from JBioSeqTools import vector_build as vb
from JBioSeqTools import api as ap
import pandas as pd


# Czesc pierwsza
# zmienne inicjalne

project_name = 'test'

species = 'human'

vector_type = 'ssAAV'

#liczba transkryptów / genów 'n'

n = 1


project = vb.create_project(project_name)


# Czesc druga:

# Ładujemy projekt z czesci peirwszej i uzupelniamy o wszytskie elementy z czesci pierwszej zwiazanej z tabelka tranmskryptow


# ładujemy metadane
metadata = vb.load_metadata()


# ładujemy wszystkie informacje do elements


# informacje o iloci transkryptów i linkerów w formie ORFn - transcript ln - linker
# 1 x transcript - ['ORF1', 'linker1']
# 2 x transcripts - ['ORF1', 'linker1', 'ORF2']
# 3 x transcript - ['ORF1', 'linker1', 'ORF2']
# 4 x transcript - ['ORF1', 'linker1', 'ORF2', 'linker2', 'ORF3']


project['elements']['transcripts'] = ['ORF1', 'linker1']

# linkers

# sequence
# kiedy brak lub kiedy tylko jeden transcript
project['elements']['linkers']['linker1'] = 'GGGGGGG' or ''
# name
# kiedy brak lub kiedy tylko jeden transcript
project['elements']['linkers']['linker1_name'] = 'NAME' or ''


# fluorescence

# sequence
# jeżeli brak
project['elements']['fluorescence']['sequence'] = 'ATGGTGAGGGGGCCTGCGTAGC' or ''
# name
project['elements']['fluorescence']['name'] = 'NAME' or ''  # jeżeli brak

# fluorescence linker

# sequence
project['elements']['fluorescence']['linker'] = 'GGGGGG' or ''  # jeżeli brak
# name
project['elements']['fluorescence']['linker_name'] = 'NAME' or ''  # jeżeli brak


# promoter

# sequence
project['elements']['promoter']['linker'] = 'GAGGGAGTGGGATGGGATTGATGAG'  # musi być
# name
project['elements']['promoter']['linker_name'] = 'NAME'  # musi być


# ładujemy dane o transkryptach

# szykujemy tabelkę


# struktura w formie słownika
# transcripts_df = {'name': [], 'ORF': [], 'sequence':[]}

# 1 transcript 
# transcripts_df = {'name': ['SMN1'], 'ORF': ['ORF1'], 'sequence':['GGATGGATGCGGCTGAGCT']}

# 2 transcripts 
transcripts_df = {'name': ['SMN1', 'SMN2'], 'ORF': ['ORF1', 'ORF2'], 'sequence':['GGATGGATGCGGCTGAGCT','GGATGGATGCGGCTGAGCT']}

# analogicznie gdy 3,4,5 ... transkryptów, musza być podane w kolejnoci jak zdefiniował użytkownik 

project['transcripts']['sequences'] = pd.DataFrame(transcripts_df)

# decyzja użytkownika o optymalizacji lub jej braku

#optymalizacja + 
#n - liczba transkryptów

#zoptymalizowane
optimization_list = ['o']*n

#bez-optymalizacji
optimization_list = ['n']*n





#pedefined by user list of restriction enzymes
user_defined_enzymes = ['AluI', 'AvaII', 'Tsp45I'] 

#jeli brak to pusta lista
if len(user_defined_enzymes) == 0:
    user_defined_enzymes = []




###################################################################################

## KONIEC ŁADOWANIA DANYCH I PIERWSZEGO KROKU
#2 krok single-cell, na razie nie ma, będzie kompatybilne 

#3 krok analiza!


#usuwamy kodon stop z ostatniej transkrybowanej sekwencji
project = vb.check_stop(project, metadata['codons'])

#wzbogacamy sekwencje
project = vb.sequence_enrichment(project, metadata['codons'], species)


project = ap.add_chosen_sequence_variant(project, optimization_list)


project = vb.find_restriction_vector(project, metadata['restriction'])




#na podstawie predefiniowanej listy enzymów zostają naprawione miejsca restrykcyjne

# wybieranie miejsc restrykcyjnych na podstawie wyboru użytkownika
list_of_list = []
user_defined_enzymes = [e.upper() for e in user_defined_enzymes]
for trans in project['transcripts']['sequences'].index:
    if len(user_defined_enzymes) > 0:
        tmp = pd.DataFrame(project['transcripts']['sequences']['enzymes_df'][trans])
        tmp[0] = [t.upper() for t in tmp[0]]
        en = list(tmp[1][tmp[0].isin(user_defined_enzymes)])
        list_tmp = [item for sublist in en for item in sublist]
        list_of_list.append(list_tmp)
    else:
        list_of_list.append([])



project = ap.add_chosen_restriction(project, list_of_list)

project = vb.repair_restriction_vector(project, metadata['restriction'], metadata['codons'], species)


project = vb.vector_string(project, metadata['backbone'], vector_type)


project = vb.eval_vector(project, metadata['vectors'], vector_type)

title = 'ssAAV::id:1 \n Example'


#funkcja wytworzenia grafu dla wektora AAV
project, pl = vb.vector_plot_project(project, title)


#FRONTEND
#Do wywietlenia użytkownikowki
for i in  project['transcripts']['sequences'].index:
        print('-------------------------------------------------------------')
        print('name : ' + str( project['transcripts']['sequences']['ORF'][i] + ' -> ' +  project['transcripts']['sequences']['name'][i]))
        print('**************************************************************')
        print('Before optimization:')
        print('* GC % : ' + str( project['transcripts']['sequences']['vector_sequence_GC'][i]))
        print('* Mean codon frequence : ' + str( project['transcripts']['sequences']['vector_sequence_frequence'][i]))
        print('* Sequence : ' + str( project['transcripts']['sequences']['sequence'][i]))
        print('**************************************************************')
        print('After optimization and restriction repaired:')
        print('* GC % : ' + str( project['transcripts']['sequences']['optimized_vector_sequence_GC'][i]))
        print('* Mean codon frequence : ' + str( project['transcripts']['sequences']['optimized_vector_sequence_frequence'][i]))
        print('* Not removed restriction : ' + str( project['transcripts']['sequences']['not_repaired'][i]))
        print('* Sequence : ' + str( project['transcripts']['sequences']['vector_sequence'][i]))


#elementy wektora
project['vector']['elements']
# cala sekwencja wektora
project['vector']['fasta']

# graf wektora AAV
project['vector']['graph']
