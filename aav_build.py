import pandas as pd
import numpy as np
from tqdm import tqdm
from pandasgui import show


vectors = pd.read_excel('data/vectors.xlsx')


linkers = pd.read_excel('data/linkers.xlsx')


regulators = pd.read_excel('data/regulators.xlsx')


fluorescence_tag = pd.read_excel('data/fluorescence_tag.xlsx')


backbone = pd.read_excel('data/backbone.xlsx')


promotors = pd.read_excel('data/promotors.xlsx')

# codon_df = pd.read_excel('data/codon_frequence.xlsx')

# Prepare codon freq mini-db

# c = []
# g = []

# for n, codon in enumerate(codon_df['Triplet']):
#     c.append(int(codon_df['Triplet'][n].count('C')))
#     g.append(int(codon_df['Triplet'][n].count('G')))



# codon_df['GC_content'] = [x + y for x, y in zip(c, g)]

# codon_df.to_excel('codons.xlsx')

codons = pd.read_excel('data/codons.xlsx')

restriction = pd.read_excel('data/restriction_enzymes.xlsx')

def load_metadata(codons:str() = 'https://github.com/jkubis96/JBioSeqTools/blob/main/data/codons.xlsx?raw=true', vectors:str() = 'https://github.com/jkubis96/JBioSeqTools/blob/main/data/vectors.xlsx?raw=true', linkers:str() = 'https://github.com/jkubis96/JBioSeqTools/blob/main/data/linkers.xlsx?raw=true', regulators:str() = 'https://github.com/jkubis96/JBioSeqTools/blob/main/data/regulators.xlsx?raw=true', fluorescence_tag:str() = 'https://github.com/jkubis96/JBioSeqTools/blob/main/data/fluorescence_tag.xlsx?raw=true', backbone:str() = 'https://github.com/jkubis96/JBioSeqTools/blob/main/data/backbone.xlsx?raw=true', promotors:str() = 'https://github.com/jkubis96/JBioSeqTools/blob/main/data/promotors.xlsx?raw=true'):
    globals()['codons'] = pd.read_excel(codons)
    globals()['vectors'] = pd.read_excel(vectors)
    globals()['linkers'] = pd.read_excel(linkers)
    globals()['regulators'] = pd.read_excel(regulators)
    globals()['fluorescence_tag'] = pd.read_excel(fluorescence_tag)
    globals()['backbone'] = pd.read_excel(backbone)
    globals()['promotors'] = pd.read_excel(promotors)
    
load_metadata()    
    
transcripts = 'ATGAGAGGCGCTCGCGGCGCCCAAAAAAGGTGGTGGGATTTTCTCTGCGTTCTGCTCCTACTGCTTCGCGTCCAGACAGGCTCTTCTCAACCATCTGTGAGTCCAGGGGAACCGTCTCCACCATCCATCCATCCAGGAAAATCAGACTTAATAGTCCGCGTGGGCGACGAGATTAGGCTGTTATGCACTGATCCGGGCTTTGTCAAATGGACTTTTGAGATCCTGGATGAAACGAATGAGAATAAGCAGAATGAATGGATCACGGAAAAGGCAGAAGCCACCAACACCGGCAAATACACGTGCACCAACAAACACGGCTTAAGCAATTCCATTTATGTGTTTGTTAGAGATCCTGCCAAGCTTTTCCTTGTTGACCGCTCCTTGTATGGGAAAGAAGACAACGACACGCTGGTCCGCTGTCCTCTCACAGACCCAGAAGTGACCAATTATTCCCTCAAGGGGTGCCAGGGGAAGCCTCTTCCCAAGGACTTGAGGTTTATTCCTGACCCCAAGGCGGGCATCATGATCAAAAGTGTGAAACGCGCCTACCATCGGCTCTGTCTGCATTGTTCTGTGGACCAGGAGGGCAAGTCAGTGCTGTCGGAAAAATTCATCCTGAAAGTGAGGCCAGCCTTCAAAGCTGTGCCTGTTGTGTCTGTGTCCAAAGCAAGCTATCTTCTTAGGGAAGGGGAAGAATTCACAGTGACGTGCACAATAAAAGATGTGTCTAGTTCTGTGTACTCAACGTGGAAAAGAGAAAACAGTCAGACTAAACTACAGGAGAAATATAATAGCTGGCATCACGGTGACTTCAATTATGAACGTCAGGCAACGTTGACTATCAGTTCAGCGAGAGTTAATGATTCTGGAGTGTTCATGTGTTATGCCAATAATACTTTTGGATCAGCAAATGTCACAACAACCTTGGAAGTAGTAGATAAAGGATTCATTAATATCTTCCCCATGATAAACACTACAGTATTTGTAAACGATGGAGAAAATGTAGATTTGATTGTTGAATATGAAGCATTCCCCAAACCTGAACACCAGCAGTGGATCTATATGAACAGAACCTTCACTGATAAATGGGAAGATTATCCCAAGTCTGAGAATGAAAGTAATATCAGATACGTAAGTGAACTTCATCTAACGAGATTAAAAGGCACCGAAGGAGGCACTTACACATTCCTAGTGTCCAATTCTGACGTCAATGCTGCCATAGCATTTAATGTTTATGTGAATACAAAACCAGAAATCCTGACTTACGACAGGCTCGTGAATGGCATGCTCCAATGTGTGGCAGCAGGATTCCCAGAGCCCACAATAGATTGGTATTTTTGTCCAGGAACTGAGCAGAGATGCTCTGCTTCTGTACTGCCAGTGGATGTGCAGACACTAAACTCATCTGGGCCACCGTTTGGAAAGCTAGTGGTTCAGAGTTCTATAGATTCTAGTGCATTCAAGCACAATGGCACGGTTGAATGTAAGGCTTACAACGATGTGGGCAAGACTTCTGCCTATTTTAACTTTGCATTTAAAGGTAACAACAAAGAGCAAATCCATCCCCACACCCTGTTCACTCCTTTGCTGATTGGTTTCGTAATCGTAGCTGGCATGATGTGCATTATTGTGATGATTCTGACCTACAAATATTTACAGAAACCCATGTATGAAGTACAGTGGAAGGTTGTTGAGGAGATAAATGGAAACAATTATGTTTACATAGACCCAACACAACTTCCTTATGATCACAAATGGGAGTTTCCCAGAAACAGGCTGAGTTTTGGGAAAACCCTGGGTGCTGGAGCTTTCGGGAAGGTTGTTGAGGCAACTGCTTATGGCTTAATTAAGTCAGATGCGGCCATGACTGTCGCTGTAAAGATGCTCAAGCCGAGTGCCCATTTGACAGAACGGGAAGCCCTCATGTCTGAACTCAAAGTCCTGAGTTACCTTGGTAATCACATGAATATTGTGAATCTACTTGGAGCCTGCACCATTGGAGGGCCCACCCTGGTCATTACAGAATATTGTTGCTATGGTGATCTTTTGAATTTTTTGAGAAGAAAACGTGATTCATTTATTTGTTCAAAGCAGGAAGATCATGCAGAAGCTGCACTTTATAAGAATCTTCTGCATTCAAAGGAGTCTTCCTGCAGCGATAGTACTAATGAGTACATGGACATGAAACCTGGAGTTTCTTATGTTGTCCCAACCAAGGCCGACAAAAGGAGATCTGTGAGAATAGGCTCATACATAGAAAGAGATGTGACTCCCGCCATCATGGAGGATGACGAGTTGGCCCTAGACTTAGAAGACTTGCTGAGCTTTTCTTACCAGGTGGCAAAGGGCATGGCTTTCCTCGCCTCCAAGAATTGTATTCACAGAGACTTGGCAGCCAGAAATATCCTCCTTACTCATGGTCGGATCACAAAGATTTGTGATTTTGGTCTAGCCAGAGACATCAAGAATGATTCTAATTATGTGGTTAAAGGAAACGCTCGACTACCTGTGAAGTGGATGGCACCTGAAAGCATTTTCAACTGTGTATACACGTTTGAAAGTGACGTCTGGTCCTATGGGATTTTTCTTTGGGAGCTGTTCTCTTTAGGAAGCAGCCCCTATCCTGGAATGCCGGTCGATTCTAAGTTCTACAAGATGATCAAGGAAGGCTTCCGGATGCTCAGCCCTGAACACGCACCTGCTGAAATGTATGACATAATGAAGACTTGCTGGGATGCAGATCCCCTAAAAAGACCAACATTCAAGCAAATTGTTCAGCTAATTGAGAAGCAGATTTCAGAGAGCACCAATCATATTTACTCCAACTTAGCAAACTGCAGCCCCAACCGACAGAAGCCCGTGGTAGACCATTCTGTGCGGATCAATTCTGTCGGCAGCACCGCTTCCTCCTCCCAGCCTCTGCTTGTGCACGACGATGTCTGA'

n = 2

def load_transcripts(n):
    transcripts = {'name': [], 'ORF': [], 'sequence': []}
    for i in range(1,n+1):
        check = True
        check_name = True
        while (check == True or check_name == True):
            if str('ORF' + str(i)) not in globals() and check == True:
                globals()[str('ORF' + str(i))] = input('Writte sequence ' + str('ORF'+str(i)) + ': ').replace('\\n', '\n')
                globals()[str('ORF' + str(i))] = ''.join(c.upper() for c in eval(str('ORF' + str(i))) if c.isalpha())
            if str('ORF' + str(i) + '_gen') not in globals() and check_name == True:
                globals()[str('ORF' + str(i) + '_gen')] = input('Writte sequence name ' + str('ORF'+str(i)) + ': ')
                globals()[str('ORF' + str(i) + '_gen')] = eval(str('ORF' + str(i) + '_gen')).upper()

            test = eval(str('ORF'+str(i)))
            test = [eval(str('ORF'+str(i)))[y:y+3] for y in range(0, len(eval(str('ORF'+str(i)))), 3)]
            if ((len(test) == 0) or len(test[-1]) < 3):
                print("Wrong sequence " + str(i) + ". No three nucleotides repeat. Load right transcript")
                check = True
                del  globals()[str('ORF' + str(i))]
            else:
                check = False
            if (len(eval(str('ORF' + str(i) + '_gen'))) == 0):
                print("Wrong name.  Writte sequence name")
                check_name = True
                del  globals()[str('ORF' + str(i) + '_gen')]
            else:
                check_name = False
                
        transcripts['name'].append(eval(str('ORF' + str(i) + '_gen')).upper())
        transcripts['ORF'].append(str('ORF' + str(i)))
        transcripts['sequence'].append(''.join(c.upper() for c in eval(str('ORF' + str(i))) if c.isalpha()))
        del  globals()[str('ORF' + str(i) + '_gen')]
        del  globals()[str('ORF' + str(i))]
        
    transcripts = pd.DataFrame(transcripts)
    transcript_list = []
    for i in range(1,n+1):
        transcript_list.append(str('ORF'+str(i)))
        transcript_list.append(str('linker'+str(i)))
    
    transcript_list = transcript_list[0:len(transcript_list) - 1]
    
    return transcripts, transcript_list       
                
transcripts, transcript_list = load_transcripts(n)





def choose_promotor(promotors:pd.DataFrame()):
    if 'promotor' not in globals():
        for lin in promotors['id']:
            print('-------------------------------------------------------------')
            print('id : ' + str(lin))
            print('name : ' + str(promotors['name'][promotors['id'] == lin][lin-1]))
            print('specificity : ' + str(promotors['tissue'][promotors['id'] == lin][lin-1]))
            print('description : ' + str(promotors['description'][promotors['id'] == lin][lin-1]))
            print('role : ' + str(promotors['role'][promotors['id'] == lin][lin-1]))
            print('reference : ' + str(promotors['ref'][promotors['id'] == lin][lin-1]))
       
        check = True
        while (check == True):
            x = input('Write id for promotor: ')
            if (int(locals()['x']) > 0 and len(locals()['x']) > 0) and locals()['x'].isnumeric() and (int(locals()['x']) in range(0, len(promotors['role'])+1)):
    
                if x == str(0):
                    globals()[str('promotor_name')] = ''
                    globals()[str('promotor')] = ''
                else:
                    globals()[str('promotor_name')] = str(promotors['name'][promotors['id'] == eval(x)][eval(x)-1])
                    globals()[str('promotor')] = str(promotors['seq'][promotors['id'] == eval(x)][eval(x)-1])
                check = False

        
choose_promotor(promotors)



def choose_fluorescence(fluorescence_tag:pd.DataFrame(), linkers:pd.DataFrame()):
    check_f = True
    check_l = True
    while(check_f == True and check_l == True):
        if 'fluorescence' not in globals() and check_f == True:
            print('-------------------------------------------------------------')
            print('id : 0')
            print('Lack of fluorescence tag')
            for lin in fluorescence_tag['id']:
                print('-------------------------------------------------------------')
                print('id : ' + str(lin))
                print('name : ' + str(fluorescence_tag['name'][linkers['id'] == lin][lin-1]))
                print('description : ' + str(fluorescence_tag['description'][fluorescence_tag['id'] == lin][lin-1]))
                print('role : ' + str(fluorescence_tag['role'][fluorescence_tag['id'] == lin][lin-1]))
                print('reference : ' + str(fluorescence_tag['ref'][fluorescence_tag['id'] == lin][lin-1]))
    
            locals()['x'] = input('Write id for fluorescence tag: ')
            if (len(locals()['x']) > 0) and locals()['x'].isnumeric() and (int(locals()['x']) in range(0, len(fluorescence_tag['role'])+1) ):
                check_f = False
                if locals()['x'] == str(0):
                    globals()['fluorescence'] = ''
                    globals()['fluorescence_name'] = ''
                    globals()['fluoroscence_tag_linker_name'] = ''
                    globals()['fluoroscence_tag_linker'] = ''
                else:
                    globals()[str('fluorescence')] = str(fluorescence_tag['seq'][fluorescence_tag['id'] == int(locals()['x'])][int(locals()['x'])-1])
                    globals()[str('fluorescence_name')] = str(fluorescence_tag['name'][fluorescence_tag['id'] == int(locals()['x'])][int(locals()['x'])-1])

        if 'fluorescence' in globals():
            check_f = False
            
           
                
        if 'fluoroscence_tag_linker' not in globals() and check_l == True and check_f == False:
            print('-------------------------------------------------------------')
            print('id : 0')
            print('Lack of linker between last protein and fluorescence tag')
            for lin in linkers['id']:
                print('-------------------------------------------------------------')
                print('id : ' + str(lin))
                print('name : ' + str(linkers['name'][linkers['id'] == lin][lin-1]))
                print('description : ' + str(linkers['description'][linkers['id'] == lin][lin-1]))
                print('role : ' + str(linkers['role'][linkers['id'] == lin][lin-1]))
                
            
            locals()['l'] = input('Write id for linker: ')
            
            if (len(locals()['l']) > 0) and locals()['l'].isnumeric() and (int(locals()['l']) in range(0, len(linkers['role'])+1)):
                check_l = False
                if locals()['l'] == str(0):
                    globals()['fluoroscence_tag_linker_name'] = ''
                    globals()['fluoroscence_tag_linker'] = ''
                else:
                    globals()['fluoroscence_tag_linker_name'] = str(linkers['name'][linkers['id'] == int(locals()['l'])][int(locals()['l'])-1])
                    globals()['fluoroscence_tag_linker'] = str(linkers['seq'][linkers['id'] == int(locals()['l'])][int(locals()['l'])-1])

        if 'fluoroscence_tag_linker' in globals():
            check_f = False
     
        
    

    
choose_fluorescence(fluorescence_tag, linkers)



def choose_linkers(n:int(), linkers:pd.DataFrame()):
    if 'linker1' not in globals():
        print('-------------------------------------------------------------')
        print('id : 0')
        print('Lack of linker between proteins')
        for lin in linkers['id']:
            print('-------------------------------------------------------------')
            print('id : ' + str(lin))
            print('name : ' + str(linkers['name'][linkers['id'] == lin][lin-1]))
            print('description : ' + str(linkers['description'][linkers['id'] == lin][lin-1]))
            print('role : ' + str(linkers['role'][linkers['id'] == lin][lin-1]))
     
    
    for i in range(1,n):
        if str('linker' + str(i)) not in globals():
            check = True
            while (check == True):
                locals()['x'] = input('Write id for linker ' + str(i) + ': ')
                if (len(locals()['x']) > 0) and locals()['x'].isnumeric() and (int(locals()['x']) in range(0, len(linkers['role'])+1)):
                    if locals()['x'] == str(0):
                        globals()[str('linker'+str(i))] = ''
                        globals()[str('linker'+str(i) + '_name')] = ''

                    else:
                        globals()[str('linker'+str(i))] = str(linkers['seq'][linkers['id'] == int(locals()['x'])][int(locals()['x'])-1])
                        globals()[str('linker'+str(i) + '_name')] = str(linkers['name'][linkers['id'] == int(locals()['x'])][int(locals()['x'])-1])

                    check = False
         
        
choose_linkers(n, linkers)




def choose_enhancer(regulators:pd.DataFrame()):
    if 'enhancer' not in globals():
        print('-------------------------------------------------------------')
        print('id : 0')
        print('Lack of regulators')
        for lin in regulators['id']:
            print('-------------------------------------------------------------')
            print('id : ' + str(lin))
            print('name : ' + str(regulators['name'][regulators['id'] == lin][lin-1]))
            print('description : ' + str(regulators['description'][regulators['id'] == lin][lin-1]))
            print('role : ' + str(regulators['role'][regulators['id'] == lin][lin-1]))
            print('type : ' + str(regulators['type'][regulators['id'] == lin][lin-1]))
            print('reference : ' + str(regulators['ref'][regulators['id'] == lin][lin-1]))
       
        check = True
        while (check == True):
            x = input('Write id for regulator: ')
            if (len(locals()['x']) > 0) and locals()['x'].isnumeric() and (int(locals()['x']) in range(0, len(regulators['role'])+1)):
    
                if x == str(0):
                    globals()[str('enhancer_name')] = ''
                    globals()[str('enhancer')] = ''

                else:
                    globals()[str('enhancer_name')] = str(regulators['name'][regulators['id'] == eval(x)][eval(x)-1])
                    globals()[str('enhancer')] = str(regulators['seq'][regulators['id'] == eval(x)][eval(x)-1])

                check = False

        
choose_enhancer(regulators)




def ORF_STOP(transcripts:pd.DataFrame(), fluorescence:str(), codons:pd.DataFrame()):
    if len(transcripts['sequence']) > 1 and len(fluorescence) > 0:
        repaired = []
        for transcript in range(0,len(transcripts['sequence'])):
            test = [transcripts['sequence'][transcript][y:y+3] for y in range(0, len(transcripts['sequence'][transcript]), 3)]
            if test[-1] in list(codons['Triplet'][codons['Amino acid'] == '*']):
                repaired.append(transcripts['sequence'][transcript][0:len(transcripts['sequence'][transcript])-3])
            else:
                repaired.append(transcripts['sequence'][transcript])
    elif len(transcripts['sequence']) > 1 and len(fluorescence) == 0:
        repaired = []
        for transcript in range(0,len(transcripts['sequence'])):
            test = [transcripts['sequence'][transcript][y:y+3] for y in range(0, len(transcripts['sequence'][transcript]), 3)]
            if transcript < max(range(0,len(transcripts['sequence']))) and test[-1] in list(codons['Triplet'][codons['Amino acid'] == '*']):
                repaired.append(transcripts['sequence'][transcript][0:len(transcripts['sequence'][transcript])-3])
            else:
                repaired.append(transcripts['sequence'][transcript])
                
    transcripts['vector_sequence'] = repaired

    return transcripts


transcripts = ORF_STOP(transcripts, fluorescence, codons)


def codon_otymization(sequence:str(), codons:pd.DataFrame, species:str()):
    codons = codons[codons['Species'] == species]
    seq_codon = [sequence[y:y+3] for y in range(0, len(sequence), 3)]
    seq_codon_fr = [codons['Fraction'][codons['Triplet'] == seq][codons['Fraction'][codons['Triplet'] == seq].index[0]] for seq in seq_codon]
    seq_codon_fr = round(sum(seq_codon_fr) / len(seq_codon_fr),2)
    
    seq_codon_mean = ''.join(seq_codon).count('C')
    seq_codon_GC = (''.join(seq_codon).count('C') + ''.join(seq_codon).count('G')) / len(''.join(seq_codon)) * 100
    seq_aa = []
    for element in seq_codon:
        tmp = codons['Amino acid'][codons['Triplet'] == element]
        tmp = tmp.reset_index()
        seq_aa.append(tmp['Amino acid'][0])
        
    mean_GC = (len(sequence)-1)*58/100/len(sequence)*3
    
    seq_tmp = []
    
    for element in seq_aa:
        tmp = codons[codons['Amino acid'] == element].sort_values(['Fraction', 'GC_content'], ascending=[False, False])
        tmp = tmp.reset_index()
        seq_tmp.append(tmp['Triplet'][0])
        
    c = []
    g = []

    for n, codon in enumerate(seq_tmp):
        c.append(int(seq_tmp[n].count('C')))
        g.append(int(seq_tmp[n].count('G')))
    
    tmp2 = [x + y for x, y in zip(c, g)]
    df = np.array([seq_tmp, tmp2])
    seq_tmp_GC_1 = (''.join(seq_tmp).count('C') + ''.join(seq_tmp).count('G')) / len(''.join(seq_tmp)) * 100

    m = 1
    for i in tqdm(range(1, len(df[1]))):
        if m/(i) > mean_GC*1.05:
            tmp_np = df[0:2,i-1:i+1]
            aa_1 =  str(codons['Amino acid'][codons['Triplet'] == df[0,i-1]][codons['Amino acid'][codons['Triplet'] == df[0,i-1]].index[0]])
            aa_2 =  str(codons['Amino acid'][codons['Triplet'] == df[0,i]][codons['Amino acid'][codons['Triplet'] == df[0,i]].index[0]])
            tmp_1 = codons[codons['Amino acid'] == aa_1].sort_values(['Fraction', 'GC_content'], ascending=[False, False])
            tmp_1 = tmp_1.reset_index()
            fr1_up = tmp_1['Fraction'][0]
            tmp_1 = tmp_1[tmp_1['GC_content'] < int(df[1,i-1])]
            if len(tmp_1) > 0:
                tmp_1 = tmp_1.reset_index()
                fr1_down = tmp_1['Fraction'][0]
                diff1 = fr1_up - fr1_down
            else: 
                diff1 = 1000
            tmp_2 = codons[codons['Amino acid'] == aa_2].sort_values(['Fraction', 'GC_content'], ascending=[False, False])
            tmp_2 = tmp_2.reset_index()
            fr2_up = tmp_2['Fraction'][0]
            tmp_2 = tmp_2[tmp_2['GC_content'] < int(df[1,i])]
            if len(tmp_2) > 0:
                tmp_2 = tmp_2.reset_index()
                fr2_down = tmp_2['Fraction'][0]
                diff2 = fr2_up - fr2_down
            else: 
                diff2 = 1000


            if diff1 <= diff2 and diff1 != 1000:
               df[0,i-1] = tmp_1['Triplet'][0]
               df[1,i-1] = tmp_1['GC_content'][0]
               m += int(tmp_1['GC_content'][0])
            elif diff1 > diff2:
                df[0,i] = tmp_2['Triplet'][0]
                df[1,i] = tmp_2['GC_content'][0]
                m += int(tmp_2['GC_content'][0])
            elif diff1 == 1000 &  diff2 == 1000:
                next
        else:
            m += int(df[1,i])
                    
                  
    seq_tmp_GC_2 = (''.join(df[0]).count('C') + ''.join(df[0]).count('G')) / len(''.join(df[0])) * 100
    
    seq_aa_2 = []
    for element in df[0]:
        tmp = codons['Amino acid'][codons['Triplet'] == element]
        tmp = tmp.reset_index()
        seq_aa_2.append(tmp['Amino acid'][0])
        
    seq_codon_fr2 = [codons['Fraction'][codons['Triplet'] == seq][codons['Fraction'][codons['Triplet'] == seq].index[0]] for seq in df[0]]
    seq_codon_fr2 = round(sum(seq_codon_fr2) / len(seq_codon_fr2),2)
        
    df_final = {'status':[], 'sequence_na':[], 'sequence_aa':[], 'frequence':[], 'GC%': []}
    df_final['status'].append('not_optimized')
    df_final['status'].append('optimized')
    df_final['sequence_na'].append(''.join(seq_codon))
    df_final['sequence_na'].append(''.join(list(df[0])))
    df_final['sequence_aa'].append(''.join(seq_aa))
    df_final['sequence_aa'].append(''.join(seq_aa_2))
    df_final['frequence'].append(seq_codon_fr)
    df_final['frequence'].append(seq_codon_fr2)
    df_final['GC%'].append(seq_codon_GC)
    df_final['GC%'].append(seq_tmp_GC_2)
    
    df_final = pd.DataFrame(df_final)
    
    return df_final



def transcript_expression_enrichment(transcripts:pd.DataFrame(), codons:pd.DataFrame(), species:str()):
    transcripts['vector_sequence_GC'] = np.nan  
    transcripts['vector_sequence_frequence'] = np.nan  
    transcripts['optimized_vector_sequence'] = np.nan  
    transcripts['optimized_vector_sequence_GC'] = np.nan  
    transcripts['optimized_vector_sequence_frequence'] = np.nan  
    transcripts['sequence_aa'] = np.nan 
    for tn in range(0,len(transcripts['sequence'])):
            tmp = codon_otymization(transcripts['vector_sequence'][tn], codons, species)
            transcripts['vector_sequence_GC'][tn] = tmp['GC%'][0] 
            transcripts['vector_sequence_frequence'][tn] = tmp['frequence'][0] 
            transcripts['optimized_vector_sequence'] = tmp['sequence_na'][1] 
            transcripts['optimized_vector_sequence_GC'][tn] = tmp['GC%'][1] 
            transcripts['optimized_vector_sequence_frequence'][tn] = tmp['frequence'][1] 
            transcripts['sequence_aa'] = tmp['sequence_aa'][1]
            


    return transcripts

transcripts = transcript_expression_enrichment(transcripts, codons, 'human')



def choose_transcript_variant(transcripts:pd.DataFrame()):
    for i in transcripts.index:
        if str('ORF_sv' + str(i+1)) not in globals():
            print('-------------------------------------------------------------')
            print('name : ' + str(transcripts['ORF'][i] + ' -> ' + transcripts['name'][i]))
            print('**************************************************************')
            print('Before optimization:')
            print('* GC % : ' + str(transcripts['vector_sequence_GC'][i]))
            print('* Mean codon frequence : ' + str(transcripts['vector_sequence_frequence'][i]))
            print('**************************************************************')
            print('After optimization:')
            print('* GC % : ' + str(transcripts['optimized_vector_sequence_GC'][i]))
            print('* Mean codon frequence : ' + str(transcripts['optimized_vector_sequence_frequence'][i]))
            print('Choose sequence: optimized [o] or not optimized [n]')
    
            
            check = True
            while (check == True):
                locals()[str('ORF_sv' + str(i+1))] = input('Writte your choose [o/n]: ')
                if str('ORF_sv' + str(i+1)) in locals() and locals()[str('ORF_sv' + str(i+1))] == 'o' or str('ORF_sv' + str(i+1)) in locals() and locals()[str('ORF_sv' + str(i+1))] == 'n':
                    check = False
                    
                
        if str('ORF_sv' + str(i+1)) in locals() and locals()[str('ORF_sv' + str(i+1))] == 'o' :
            transcripts['vector_sequence'][i] = transcripts['optimized_vector_sequence'][i]
        if str('ORF_sv' + str(i+1)) in globals() and globals()[str('ORF_sv' + str(i+1))] == 'o' :
            transcripts['vector_sequence'][i] = transcripts['optimized_vector_sequence'][i]

        
    
    return transcripts       
                
transcripts_vector = choose_transcript_variant(transcripts)





def check_restriction(sequence:str(), restriction:pd.DataFrame()):

    enzyme_restriction = {'name':[], 'restriction_place':[], 'restriction_sequence':[], 'sequence':[], 'start':[], 'stop':[]}
    
    for r in tqdm(restriction.index):
        check = True
        if restriction['sequence'][r] in sequence:
            while(check == True):
                bmp = list(sequence)
                for n in range(0,len(restriction['sequence'][r])):
                    for j in range(n,len(bmp)-len(restriction['sequence'][r])):
                       lower = j
                       upper = j + len(restriction['sequence'][r])
                       if upper < len(bmp) and ''.join(bmp[lower:upper]) == restriction['sequence'][r]:
                            enzyme_restriction['name'].append(restriction['name'][r])
                            enzyme_restriction['restriction_sequence'].append(restriction['sequence'][r])
                            enzyme_restriction['restriction_place'].append(restriction['restriction_place'][r])
                            enzyme_restriction['sequence'].append(sequence)
                            enzyme_restriction['start'].append(lower)
                            enzyme_restriction['stop'].append(upper)
                            check = False

                               
    enzyme_restriction = pd.DataFrame.from_dict(enzyme_restriction)
    enzyme_restriction = enzyme_restriction.drop_duplicates()
    enzyme_restriction = enzyme_restriction.reset_index(drop=True)
    
    if len(enzyme_restriction) == 0:
        print('Any restriction places were found')
    
    return enzyme_restriction






def choose_restriction_to_remove(restriction_df:pd.DataFrame(), enzyme_list:list() = []):
    if len(restriction_df) != 0 and len(enzyme_list) == 0:
        for i in restriction_df.index:
            if str('ORF_sv' + str(i+1)) not in globals():
                print('-------------------------------------------------------------')
                print('id : ' + str(i))
                print('name : ' + restriction_df['name'][i])
                print('restriction_side : ' + restriction_df['restriction_place'][i])
    
        enzyme_list = []
        check = True
        enzyme_n = 1
        while (check == True):
            print('Provide enzyme id, if no restriction sites are relevant to your experiment or you have already provided all enzyme ids, write "x"')
            enzyme = input('Write enzyme '+str(enzyme_n) + ' id: ')
            if len(enzyme) != 0 and not enzyme.isalpha() and int(enzyme) in restriction_df.index:
                enzyme_n += 1
                enzyme_list.append(int(enzyme))
            elif len(enzyme) != 0 and enzyme.upper() == 'X':
                check = False
        
        enzyme_list = np.unique(enzyme_list)
    else:
        print('Any restriction places to choose')
        
    return np.asarray(enzyme_list)       




def repair_sequences(sequence:str(), codons:pd.DataFrame, restriction_df:pd.DataFrame(), restriction:pd.DataFrame(), enzyme_list:list(), species:str()):
    if len(restriction_df) != 0:
        not_repaired = []
        codons = codons[codons['Species'] == species]
        seq_codon = [sequence[y:y+3] for y in range(0, len(sequence), 3)]
        seq_codon_fr = [codons['Fraction'][codons['Triplet'] == seq][codons['Fraction'][codons['Triplet'] == seq].index[0]] for seq in seq_codon]
        seq_codon_fr = round(sum(seq_codon_fr) / len(seq_codon_fr),2)
        seq_codon_GC = (''.join(seq_codon).count('C') + ''.join(seq_codon).count('G')) / len(''.join(seq_codon)) * 100
        
        seq_aa = []
        for element in seq_codon:
            tmp = codons['Amino acid'][codons['Triplet'] == element]
            tmp = tmp.reset_index()
            seq_aa.append(tmp['Amino acid'][0])
        
        dic = {'seq':[], 'range':[], 'codon_n':[], 'aa':[]}
        n = 0
        for num, seq in enumerate(seq_codon):
            for i in range(0,3):
                dic['seq'].append(seq)
                dic['range'].append(n)
                dic['codon_n'].append(num)
                dic['aa'].append(seq_aa[num])
                n += 1
                
        dic = pd.DataFrame.from_dict(dic)
        
        print('\n Codon changing...')
        for eid in tqdm(enzyme_list):
            check = dic[['seq','codon_n']].drop_duplicates()
            check = ''.join(check['seq'])
            if restriction_df['sequence'][eid] in check:
                dic_tmp = dic[(dic['range'] >= restriction_df['start'][eid]) & (dic['range'] < restriction_df['stop'][eid])] 
                tmp = codons[codons['Triplet'].isin(np.unique(dic['seq']))]
                dictionary = {'seq':[], 'aa':[], 'triplet':[], 'freq':[], 'gc':[]}
                for i in np.unique(dic_tmp['seq']):
                    t = tmp[tmp['Triplet'] == i]
                    t = t.reset_index(drop = True)
                    t = tmp[tmp['Amino acid'] == t['Amino acid'][0]]
                    for n in t.index:
                        dictionary['seq'].append(i)
                        dictionary['aa'].append(t['Amino acid'][n])
                        dictionary['triplet'].append(t['Triplet'][n])
                        dictionary['freq'].append(t['Fraction'][n])
                        dictionary['gc'].append(t['GC_content'][n])
        
                
                dictionary = pd.DataFrame.from_dict(dictionary)
                dictionary = dictionary[~dictionary['triplet'].isin(dictionary['seq'])]
                dictionary = dictionary.sort_values(['freq', 'gc'], ascending=[False, False])
                dictionary = dictionary.reset_index()
        
        
                all_enzymes_variant = restriction[restriction['name'] == restriction_df['name'][eid]]
                
                
                seq_new = dic_tmp[['seq','codon_n']].drop_duplicates()
                seq_old = ''.join(seq_new['seq'])
                
                for d in dictionary.index:
                    seq_new = dic_tmp[['seq','codon_n']].drop_duplicates()
                    dictionary['seq'][d]
                    dictionary['triplet'][d]
                    seq_new['seq'][seq_new['seq'] == dictionary['seq'][d]] = dictionary['triplet'][d]
                    if ''.join(seq_new['seq']) in all_enzymes_variant['sequence']:
                        break
                    
                if seq_old == ''.join(seq_new['seq']):
                    print('It was impossieble to change restriction side for', + str(restriction_df['name'][eid]))
                    not_repaired.append(restriction_df['name'][eid])
                elif seq_old != ''.join(seq_new['seq']):
                    for new in seq_new.index:
                        dic['seq'][dic['codon_n'] == seq_new['codon_n'][new]]  = seq_new['seq'][new]
        
    
            
        final_sequence = dic[['seq','codon_n']].drop_duplicates()
        final_sequence = ''.join(final_sequence['seq'])
        
        if len(not_repaired) == 0:    
            print('\n Restriction place in sequence repaired...')
        else:
            print('\n Restriction place for:')
            for i in not_repaired:
                print('\n'+ str(i))
                
            print('\n were unable to optimize:')
            print('\n Rest of chosen restriction place in sequence repaired...')
    
    
        enzyme_restriction = {'name':[], 'restriction_place':[], 'restriction_sequence':[], 'sequence':[], 'start':[], 'stop':[]}
        
        print('\n Checking new restriction...')
        for r in tqdm(restriction.index):
            check = True
            if restriction['sequence'][r] in final_sequence:
                while(check == True):
                    bmp = list(final_sequence)
                    for n in range(0,len(restriction['sequence'][r])):
                        for j in range(n,len(bmp)-len(restriction['sequence'][r])):
                           lower = j
                           upper = j + len(restriction['sequence'][r])
                           if upper < len(bmp) and ''.join(bmp[lower:upper]) == restriction['sequence'][r]:
                                enzyme_restriction['name'].append(restriction['name'][r])
                                enzyme_restriction['restriction_sequence'].append(restriction['sequence'][r])
                                enzyme_restriction['restriction_place'].append(restriction['restriction_place'][r])
                                enzyme_restriction['sequence'].append(final_sequence)
                                enzyme_restriction['start'].append(lower)
                                enzyme_restriction['stop'].append(upper)
                                check = False
    
                                   
        enzyme_restriction = pd.DataFrame.from_dict(enzyme_restriction)
        enzyme_restriction = enzyme_restriction.drop_duplicates()
        enzyme_restriction = enzyme_restriction.reset_index(drop=True)
        enzyme_restriction = enzyme_restriction[~enzyme_restriction['name'].isin(restriction_df['name'])]
    
        if len(enzyme_restriction) == 0:
            print('\n Any new restriction places were created')
        else:
            print('\n New restriction places were created:')
            for name in enzyme_restriction['name']:
                print(name)
    
    else:
        enzyme_restriction = {'name':[], 'restriction_place':[], 'restriction_sequence':[], 'sequence':[], 'start':[], 'stop':[]}
        enzyme_restriction = pd.DataFrame.from_dict(enzyme_restriction)
        not_repaired = []
        final_sequence = sequence

    return final_sequence, not_repaired, enzyme_restriction


species = 'human'


def find_restriction_vector(transcripts_vector:pd.DataFrame(), restriction:pd.DataFrame()):
    transcripts_vector['enzymes_df'] = np.NaN
    for trans in transcripts_vector.index:
        transcripts_vector['enzymes_df'][trans] = check_restriction(str(transcripts_vector['vector_sequence'][trans]), restriction).to_records(index=True)
        
    return transcripts_vector
        
transcripts_vector = find_restriction_vector(transcripts_vector, restriction)


def choose_restriction_vector(transcripts_vector:pd.DataFrame(), restriction:pd.DataFrame()):
    transcripts_vector['enzymes'] = ''
    for trans in transcripts_vector.index:
        index = pd.DataFrame(transcripts_vector['enzymes_df'][trans])
        index.index = index['index']
        print("Choose enzymes for " + str(transcripts_vector['ORF']))
        transcripts_vector['enzymes'][trans] = choose_restriction_to_remove(index).tolist()
        
    return transcripts_vector
        
transcripts_vector = choose_restriction_vector(transcripts_vector, restriction)




def repair_restriction_vector(transcripts_vector:pd.DataFrame(), restriction:pd.DataFrame(), codons:pd.DataFrame(), species:str()):
    transcripts_vector['not_repaired'] = ''
    for trans in transcripts_vector.index:
        final_sequence, not_repaired, enzyme_restriction =  repair_sequences(transcripts_vector['vector_sequence'][trans], codons, transcripts_vector['enzymes_df'][trans], restriction, transcripts_vector['enzymes'][trans] , species)
        transcripts_vector['vector_sequence'][trans] = final_sequence
        transcripts_vector['not_repaired'][trans] = not_repaired
        transcripts_vector['enzymes_df'][trans] = enzyme_restriction.to_records(index=True)
        

        transcripts_vector['optimized_vector_sequence_GC'][trans] =  (transcripts_vector['vector_sequence'][trans].count('C') + transcripts_vector['vector_sequence'][trans].count('G')) / len(transcripts_vector['vector_sequence'][trans]) * 100
        
        seq_codon = [transcripts_vector['vector_sequence'][trans][y:y+3] for y in range(0, len(transcripts_vector['vector_sequence'][trans]), 3)]
        seq_codon = [codons['Fraction'][codons['Triplet'] == seq][codons['Fraction'][codons['Triplet'] == seq].index[0]] for seq in seq_codon]
        seq_codon = round(sum(seq_codon) / len(seq_codon),2)

        transcripts_vector['optimized_vector_sequence_frequence'][trans] = seq_codon
        
    return transcripts_vector
        
transcripts_vector = repair_restriction_vector(transcripts_vector, restriction, codons, species)



    
def vector_string(backbone:pd.DataFrame(), transcript_list:str(), vector_type:str()):
    backbone = backbone[backbone['vector_type'] == vector_type]
    vector1 = str(backbone['operators'][backbone['element'] == 'p1'][0])
    for i in transcript_list:
        vector1 = vector1 + ' + ' + str(i)
   
    vector1 = vector1 + str(backbone['operators'][backbone['element'] == 'p2'][1])
    
    return vector1


text_to_eval = vector_string(backbone, transcript_list, 'ssAAV')




def eval_vector(text_to_eval, vectors:pd.DataFrame(), transcripts_vector:pd.DataFrame(), vector_type:str()):
    vectors = vectors[vectors['vector_type'] == vector_type]
    for element, n in enumerate(vectors.index):
        locals()[str(vectors['component'][n])] = vectors['sequence'][n]
      
    for element, n in enumerate(transcripts.index):
        locals()[str(transcripts['ORF'][n])] = transcripts['vector_sequence'][n]
    
    
    elements = text_to_eval.split()
    tf =  [x != '+' for x in elements]
    elemensts = [i for indx,i in enumerate(elements) if tf[indx] == True]
    
    data_frame = {'element':[], 'sequence':[], 'start':[], 'end':[], 'length': []}
    
    start = 0
    for el in elemensts:
        if el not in globals() and el not in locals():
            print('\n')
            print('Variable -> ' + str(el) + ' <- was not found. Provide all required variables!!!')
            print('\n')
            return False, False
        else:
            data_frame['element'].append(str(el))
            data_frame['sequence'].append(eval(el))
            data_frame['start'].append(start + 1)
            start = start + int(len(eval(el)))
            data_frame['end'].append(start)
            data_frame['length'].append(len(eval(el)))

    
    fasta = eval(text_to_eval)
    data_frame = pd.DataFrame(data_frame)
    data_frame = data_frame[data_frame['length'] > 0]
    
    new_element = []
    for x in data_frame['element']:
        if 'break' in x:
            new_element.append('backbone_element')
        else:
            new_element.append(x)
            
    data_frame['element'] = new_element
    data_frame = data_frame.reset_index(drop=True)
    
    for n in data_frame.index: 
        if str(data_frame['element'][n]) + '_name' in globals():
            data_frame['element'][n] = data_frame['element'][n] + ' : ' + eval(str(data_frame['element'][n]) + '_name')
        elif str(data_frame['element'][n]) in list(transcripts['ORF']):
            data_frame['element'][n] = data_frame['element'][n] + ' : ' + str(transcripts['name'][transcripts['ORF'] == 'ORF1'][transcripts['name'][transcripts['ORF'] == 'ORF1'].index[0]])

    globals_list = ['regulators', 'regulators_name', 'promotor', 'promotor_name', 'enhancer', 'enhancer_name', 'fluorescence', 'fluorescence_tag', 'fluorescence_name', 'fluoroscence_tag_linker', 'fluoroscence_tag_linker_name', 'restriction', 'linker1', 'linker1_name', 'linker2', 'linker2_name', 'linker3', 'linker3_name', 'transcript_list', 'n', 'text_to_eval', 'transcripts_vector', 'backbone', 'codons', 'linkers', 'promotors', 'vectors', 'transcripts']        
    for g in globals_list:
        if g in globals():
            del globals()[g]
            
            
    return fasta, data_frame


vector, df_seq = eval_vector(text_to_eval, vectors, transcripts, 'ssAAV')


