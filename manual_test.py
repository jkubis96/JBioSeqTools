from jbst import seq_tools as st



# Test1

metadata = st.load_metadata() 
if metadata != None:
    print('\nTest1 passed!')
else:
    print('\nTest1 failed!')
    
    

# Test2
  
metadata_reduced = st.load_metadata(linkers = False, 
                                    loops = False, 
                                    regulators = False, 
                                    fluorescent_tag = False, 
                                    promoters = False, 
                                    polya = False, 
                                    marker = False, 
                                    utr5 = False, 
                                    utr3 = False) 


if metadata_reduced != None:
    print('\nTest2 passed!')
else:
    print('\nTest2 failed!')

# Test3

test_res1 = st.get_sequences_gene('SMN1', species = 'human', max_results = 20)
if test_res1 != None:
    print('\nTest3 passed!')
else:
    print('\nTest3 failed!')
    
# Test4

test_res2 = st.get_sequences_gene('SMN1', species = 'mouse', max_results = 20)
if test_res2 != None:
    print('\nTest4 passed!')
else:
    print('\nTest4 failed!')
    
# Test5
    
test_res3 = st.get_sequences_gene('SMN1', species = 'rat', max_results = 20)
if test_res3 != None:
    print('\nTest5 passed!')
else:
    print('\nTest5 failed!')
    
    
# Test6
    
test_res4 = st.get_sequences_gene('SMN1', species = 'both', max_results = 20)
if test_res4 != None:
    print('\nTest6 passed!')
else:
    print('\nTest6 failed!')
    
# Test7
    
test_res5 = st.get_sequences_gene('SMN1', species = 'both2', max_results = 20)
if test_res5 != None:
    print('\nTest7 passed!')
else:
    print('\nTest7 failed!')
    
# Test8
    
test_res6 = st.get_sequences_gene('SMN1', species = 'multi', max_results = 20)
if test_res6 != None:
    print('\nTest8 passed!')
else:
    print('\nTest8 failed!')
    
# Test9

test_res7 = st.get_sequences_accesion(['X81403.1', 'KJ890665.1'])
if test_res7 != None:
    print('\nTest9 passed!')
else:
    print('\nTest9 failed!')
    
# Test10

test_res8 = st.generate_fasta_string(test_res4)
if test_res8 != None:
    print('\nTest10 passed!')
else:
    print('\nTest10 failed!')
    
# Test11
    
try:
    st.write_fasta(test_res8, path = None, name = 'fasta_file_test')
    print('\nTest11 passed!')
except:
    print('\nTest11 failed!')

# Test12

test_res9 = st.MuscleMultipleSequenceAlignment(test_res8, output = None, gapopen = 10, gapextend = 0.5)
if test_res9 != None:
    print('\nTest12 passed!')
else:
    print('\nTest12 failed!')

# Test13

test_res10 = st.decode_alignments(test_res9)
if test_res10 != None:
    print('\nTest13 passed!')
else:
    print('\nTest13 failed!')

# Test14

try:
    st.write_aligments(test_res10, path = None, name = 'SMN1_aligments_file_test')
    print('\nTest14 passed!')
except:
    print('\nTest14 failed!')

# Test15

test_res11 = st.ExtractConsensuse(test_res9, refseq_sequences = None)
if test_res11 != None:
    print('\nTest15 passed!')
else:
    print('\nTest15 failed!')
    
# Test16  

test_res12 = st.DisplayAlignment(test_res9, color_scheme="Taylor", wrap_length=80, show_grid=True, show_consensus=True)
if test_res12 != None:
    print('\nTest16 passed!')
else:
    print('\nTest16 failed!')

# Test17

try:
    test_res12.savefig("aligments.svg")
    print('\nTest17 passed!')
except:
    print('\nTest17 failed!')

    

# Test18

# copy below sequence to the clipboard


#PAX3

#                             atg gcgatgagca gcggcggcag tggtggcggc gtcccggagc
#        61 aggaggattc cgtgctgttc cggcgcggca caggccagag cgatgattct gacatttggg
#       121 atgatacagc actgataaaa gcatatgata aagctgtggc ttcatttaag catgctctaa
#       181 agaatggtga catttgtgaa acttcgggta aaccaaaaac cacacctaaa agaaaacctg
#       241 ctaagaagaa taaaagccaa aagaagaata ctgcagcttc cttacaacag tggaaagttg
#       301 gggacaaatg ttctgccatt tggtcagaag acggttgcat ttacccagct accattgctt
#       361 caattgattt taagagagaa acctgtgttg tggtttacac tggatatgga aatagagagg
#       421 agcaaaatct gtccgatcta ctttccccaa tctgtgaagt agctaataat atagaacaaa
#       481 atgctcaaga gaatgaaaat gaaagccaag tttcaacaga tgaaagtgag aactccaggt
#       541 ctcctggaaa taaatcagat aacatcaagc ccaaatctgc tccatggaac tcttttctcc
#       601 ctccaccacc ccccatgcca gggccaagac tgggaccagg aaagccaggt ctaaaattca
#       661 atggcccacc accgccaccg ccaccaccac caccccactt actatcatgc tggctgcctc
#       721 catttccttc tggaccacca ataattcccc caccacctcc catatgtcca gattctcttg
#       781 atgatgctga tgctttggga agtatgttaa tttcatggta catgagtggc tatcatactg
#       841 gctattatat gggtttcaga caaaatcaaa aagaaggaag gtgctcacat tccttaaatt
#       901 aa
     
     
     
sequence = st.load_sequence()
if sequence != None:
    print('\nTest18 passed!')
else:
    print('\nTest18 failed!')

# Test19

sequence = st.clear_sequence(sequence)
if sequence != None:
    print('\nTest19 passed!')
else:
    print('\nTest19 failed!')

# Test20

dec = st.check_coding(sequence)
if dec in [True, False]:
    print('\nTest20 passed!')
else:
    print('\nTest20 failed!')

# Test21

dec2 = st.check_upac(sequence)
if dec2 in [True, False]:
    print('\nTest21 passed!')
else:
    print('\nTest21 failed!')

# Test21.5

sequence_r = st.reverse(sequence = sequence)
if sequence_r != None:
    print('\nTest21.5 passed!')
else:
    print('\nTest21.5 failed!')

# Test22

sequence_c = st.complement(sequence = sequence)
if sequence_c != None:
    print('\nTest22 passed!')
else:
    print('\nTest22 failed!')

# Test23

sequence_RNA = st.dna_to_rna(sequence = sequence, enrichment= False)
if sequence_RNA != None:
    print('\nTest23 passed!')
else:
    print('\nTest23 failed!')
    
# Test24
    
sequence_RNA2 = st.dna_to_rna(sequence = sequence, enrichment= True)
if sequence_RNA2 != None:
    print('\nTest24 passed!')
else:
    print('\nTest24 failed!')

# Test25

sequence_DNA = st.rna_to_dna(sequence_RNA)
if sequence_DNA != None:
    print('\nTest25 passed!')
else:
    print('\nTest25 failed!')

# Test26

prot = st.seuqence_to_protein(sequence_RNA, metadata)
if prot != None:
    print('\nTest26 passed!')
else:
    print('\nTest26 failed!')

# Test27

prot2 = st.seuqence_to_protein(sequence, metadata)
if prot2 != None:
    print('\nTest27 passed!')
else:
    print('\nTest27 failed!')

# Test28

pred1, dot1 = st.predict_structure(sequence_RNA, 
                  anty_sequence = '',
                  height=None, 
                  width=None, 
                  dis_alpha = 0.15, 
                  seq_force = 27, 
                  pair_force = 3, 
                  show_plot = True)

if pred1 != None:
    print('\nTest28 passed!')
else:
    print('\nTest28 failed!')

# Test29

try:
    pred1.savefig('smn1_structure.svg')
    print('\nTest29 passed!')
except:
    print('\nTest29 failed!')

# Test30

res_rnai =  st.FindRNAi(sequence, metadata_reduced, length = 23, n = 200, max_repeat_len = 3, max_off = 1, species = 'human', output = None, database_name = "refseq_select_rna",  evalue = 1e-3, outfmt =  5, word_size = 7, max_hsps = 20, reward = 1, penalty = -3, gapopen = 5, gapextend = 2, dust = "no", extension = 'xml')    
if len(res_rnai.index) > 0:
    print('\nTest30 passed!')
else:
    print('\nTest30 failed!')
    

# Test31

#examle loop (obtained from metadata['loops'])
loop_seq = 'CTCGAG'

res_rnai2 = st.loop_complementary_adjustment(res_rnai, loop_seq, min_length=3)
if len(res_rnai2.index) > 0:
    print('\nTest31 passed!')
else:
    print('\nTest31 failed!')


# Test32

# copy below sequence to the clipboard


# atgttaact agtggacttg tagtaagcaa
#       481 catgttctcc tatcatctcg cagccttggg actcatgccg tcattccaga tggaagggcg
#       541 aggtcgagtt aatcagctag gaggtgtttt cattaatgga cggccactgc ccaaccatat
#       601 acgactgaag attgtcgagc tagcggccca gggcgtccgt ccgtgcgtca tcagtagaca
#       661 gctgcgggtg tcacatggct gtgtcagtaa aatactccaa cgatatcaag aaaccggaag
#       721 tatccgacct ggggttattg gcggaagtaa accaagggtc gcaactccgg aagttgagaa
#       781 aaagatagaa caatacaaaa aagataatcc gggaattttc agttgggaga ttcgggatcg
#       841 gctgctgaag gaggggattt gtgaccgcag caccgtgcca agtgtgagct ccatcagtcg
#       901 agtattacgg agcaggttcc agaaatgtga ttctgatgac aatgacaatg acaatgacaa
#       961 tgaggacgac gatggcgatg acggcagtaa cagtagtgtg gcagacaggt ctgttaactt
#      1021 ctctgtcagc ggtctgctgt ccgacaataa aagcgacaaa agcgacaacg attccgattg
#      1081 tgaatcagag ccggggctat ctgtaaaacg gaagcaacgc cgcagtcgaa ctactttcac
#      1141 cgcggagcag ttggaggaac tggaaagagc ctttgaacga actcactatc cggatatata
#      1201 tacgcgagag gaattagcac aaagaacaaa gctaaccgag gcaagagtcc aagtatggtt
#      1261 tagtaaccga agagcgagat ggcggaaaca gatgggtagc aatcagctga cagccttgaa
#      1321 cagtatatta caagtgccac agggtatggg aacgccctct tatatgctgc acgagcctgg
#      1381 gtatccactc tcacataatg cagacaatct ttggcataga tcgtctatgg cccagtcatt
#      1441 acagtcattt ggtcagacaa taaaaccaga gaattcctac gccggtctta tggaaaacta
#      1501 tttatctcat tcatcacagc ttcatggtct tcctacacat agttcatccg atcccctctc
#      1561 atccacttgg tcatctcccg tgtccacttc cgttcctgcg ctaggataca cgccatctag
#      1621 tggccattac catcattact ctgatgtcac caaaagtact cttcattcat ataacgctca
#      1681 tattccttca gtcacaaaca tggagagatg ttcagttgat gacagtttgg ttgctttacg
#      1741 tatgaagtca cgtgagcatt ccgccgctct cagtttgatg caggtggcag acaacaaaat
#      1801 ggctacctca ttttga
     

sequence_test = st.load_sequence()
if sequence_test != None:
    print('\nTest32 passed!')
else:
    print('\nTest32 failed!')

# Test33

sequence_test = st.clear_sequence(sequence_test)
if sequence_test != None:
    print('\nTest33 passed!')
else:
    print('\nTest33 failed!')

# Test34

res_rnai3 = st.remove_specific_to_sequence(res_rnai2, sequence_test)
if len(res_rnai3.index) > 0:
    print('\nTest34 passed!')
else:
    print('\nTest34 failed!')

# Test35

# RNAi
#loop the same as above
loop_seq = 'CTCGAG'
# reducing the RNAi (ought to contain GC content between 30 and 60%)
res_rnai3 = res_rnai3[(res_rnai3['GC%'] >= 30) & (res_rnai3['GC%'] <= 60)].reset_index(drop = True)
if len(res_rnai3.index) > 0:
    print('\nTest35 passed!')
else:
    print('\nTest35 failed!')
    
# Test36

# structure of shRNA example 5' - senseRNAi - loop - antisenseRNAi (silencign sequence) - TT (tail) - 3'
shrna =  st.dna_to_rna(res_rnai3['RNAi_sense'][0] + loop_seq + res_rnai3['RNAi_seq'][0] + 'TT')
if shrna != None:
    print('\nTest36 passed!')
else:
    print('\nTest36 failed!')

# Test37

pred2, dot1 = st.predict_structure(res_rnai3['RNAi_sense'][0], 
                  anty_sequence = res_rnai3['RNAi_seq'][0],
                  height=None, 
                  width=None, 
                  dis_alpha = 0.35, 
                  seq_force = 27, 
                  pair_force = 8, 
                  show_plot = True)

if pred2 != None:
    print('\nTest37 passed!')
else:
    print('\nTest37 failed!')
    
    
    
# Test37.5

pred2, dot1 = st.predict_structure(shrna, 
                  anty_sequence = '',
                  height=None, 
                  width=None, 
                  dis_alpha = 0.35, 
                  seq_force = 27, 
                  pair_force = 8, 
                  show_plot = True)

if pred2 != None:
    print('\nTest37.5 passed!')
else:
    print('\nTest37.5 failed!')
    

# Test38

try:
    pred2.savefig('smn1_sh_structure.svg')
    print('\nTest38 passed!')
except:
    print('\nTest38 failed!')



# Test39


res_optimized1 = st.codon_otymization(sequence, metadata_reduced, species = 'human')
if len(res_optimized1.index) > 0:
    print('\nTest39 passed!')
else:
    print('\nTest39 failed!')
    
# Test40

res_optimized2 = st.codon_otymization(sequence, metadata_reduced, species = 'mouse')
if len(res_optimized2.index) > 0:
    print('\nTest40 passed!')
else:
    print('\nTest40 failed!')
    
# Test41

res_optimized3 = st.codon_otymization(sequence, metadata_reduced, species = 'rat')
if len(res_optimized3.index) > 0:
    print('\nTest41 passed!')
else:
    print('\nTest41 failed!')
    
# Test42

res_optimized4 = st.codon_otymization(sequence, metadata_reduced, species = 'both')
if len(res_optimized4.index) > 0:
    print('\nTest42 passed!')
else:
    print('\nTest42 failed!')
    
# Test43

res_optimized5 = st.codon_otymization(sequence, metadata_reduced, species = 'both2')
if len(res_optimized5.index) > 0:
    print('\nTest43 passed!')
else:
    print('\nTest43 failed!')
    
# Test44

res_optimized6 = st.codon_otymization(sequence, metadata_reduced, species = 'multi')
if len(res_optimized6.index) > 0:
    print('\nTest44 passed!')
else:
    print('\nTest44 failed!')

# Test45

res1, res2 = st.check_restriction(sequence, metadata_reduced)
if len(res1.index) > 0:
    print('\nTest45 passed!')
else:
    print('\nTest45 failed!')

# Test46

#provide 70 - SfcI and 5 - BbvI and accept 'x'

repaired_sequence0 = st.sequence_restriction_removal(sequence, metadata_reduced, restriction_places = [], species = 'human')
if len(repaired_sequence0.index) > 0:
    print('\nTest46 passed!')
else:
    print('\nTest46 failed!')


# Test47
   
repaired_sequence1 = st.sequence_restriction_removal(sequence, metadata_reduced, restriction_places = ['SfcI', 'BbvI'], species = 'human')
if len(repaired_sequence0.index) > 0:
    print('\nTest47 passed!')
else:
    print('\nTest47 failed!')
   
# Test48
    
if repaired_sequence0['sequence_na'][1] == repaired_sequence1['sequence_na'][1]:
    print('\nTest48 passed!')
else:
    print('\nTest48 failed!')   
        
# Test49

repaired_sequence2 = st.sequence_restriction_removal(sequence, metadata, restriction_places = ['AclI', 'AflIII'], species = 'mouse')
if len(repaired_sequence2.index) > 0:
    print('\nTest49 passed!')
else:
    print('\nTest49 failed!')
    
# Test50
    
repaired_sequence3 = st.sequence_restriction_removal(sequence, metadata, restriction_places = ['AclI', 'AflIII'], species = 'rat')
if len(repaired_sequence3.index) > 0:
    print('\nTest50 passed!')
else:
    print('\nTest50 failed!')
    
# Test51

repaired_sequence4 = st.sequence_restriction_removal(sequence, metadata, restriction_places = ['AclI', 'AflIII'], species = 'both')
if len(repaired_sequence4.index) > 0:
    print('\nTest51 passed!')
else:
    print('\nTest51 failed!')
    
# Test52
    
repaired_sequence5 = st.sequence_restriction_removal(sequence, metadata, restriction_places = ['AclI', 'AflIII'], species = 'both2')
if len(repaired_sequence5.index) > 0:
    print('\nTest52 passed!')
else:
    print('\nTest52 failed!')
    
# Test53

repaired_sequence6 = st.sequence_restriction_removal(sequence, metadata, restriction_places = ['AclI', 'AflIII'], species = 'multi')
if len(repaired_sequence6.index) > 0:
    print('\nTest53 passed!')
else:
    print('\nTest53 failed!')



# Test54

import pkg_resources
from jbst import vector_build as vb

fasta_r = vb.load_fasta(pkg_resources.resource_filename("jbst", "tests/fasta_vector_test.fasta"))
if fasta_r != None:
    print('\nTest54 passed!')
else:
    print('\nTest54 failed!')
    
# Test55

df_fasta = vb.decode_fasta_to_dataframe(fasta_r)
if len(df_fasta.index) > 0 or df_fasta != None:
    print('\nTest55 passed!')
else:
    print('\nTest55 failed!')
    
# Test56

df_fasta = vb.extract_fasta_info(df_fasta)
if len(df_fasta.index) > 0 or df_fasta != None:
    print('\nTest56 passed!')
else:
    print('\nTest56 failed!')
    
    
# Test57

gb_format = vb.get_genebank(df_fasta, 
                     name = 'viral_vector', 
                     definition = ' Synthetic viral plasmid vector')
    
if gb_format != None:
    print('\nTest57 passed!')
else:
    print('\nTest57 failed!')
    

    
# Test58

plot = vb.plot_vector(df_fasta, title = None, title_size = 20)
if plot != None:
    print('\nTest58 passed!')
else:
    print('\nTest58 failed!')
    
# Test59

try:
    print('\nTest59 passed!')
    plot.savefig('vector_from_fasta.svg')
except:
    print('\nTest59 failed!')





##############################################################################


from jbst import vector_build as vb


metadata_reduced = vb.load_metadata(linkers = False, 
                                    loops = False, 
                                    regulators = False, 
                                    fluorescent_tag = False, 
                                    promoters = False, 
                                    polya = False, 
                                    marker = False, 
                                    utr5 = False, 
                                    utr3 = False) 



# Expression vector:
    

input_dict = {
    
    'project_name':'test_expression',
    'vector_type':'ssAAV',
    'vector_function':'expression',
    'species':'human',
    'sequences':['ATGTTAACTAGTGGACTTGTAGTAAGCAACATGTTCTCCTATCATCTCGCAGCCTTGGGACTCATGCCGTCATTCCAGATGGAAGGGCGAGGTCGAGTTAATCAGCTAGGAGGTGTTTTCATTAATGGACGGCCACTGCCCAACCATATACGACTGAAGATTGTCGAGCTAGCGGCCCAGGGCGTCCGTCCGTGCGTCATCAGTAGACAGCTGCGGGTGTCACATGGCTGTGTCAGTAAAATACTCCAACGATATCAAGAAACCGGAAGTATCCGACCTGGGGTTATTGGCGGAAGTAAACCAAGGGTCGCAACTCCGGAAGTTGAGAAAAAGATAGAACAATACAAAAAAGATAATCCGGGAATTTTCAGTTGGGAGATTCGGGATCGGCTGCTGAAGGAGGGGATTTGTGACCGCAGCACCGTGCCAAGTGTGAGCTCCATCAGTCGAGTATTACGGAGCAGGTTCCAGAAATGTGATTCTGATGACAATGACAATGACAATGACAATGAGGACGACGATGGCGATGACGGCAGTAACAGTAGTGTGGCAGACAGGTCTGTTAACTTCTCTGTCAGCGGTCTGCTGTCCGACAATAAAAGCGACAAAAGCGACAACGATTCCGATTGTGAATCAGAGCCGGGGCTATCTGTAAAACGGAAGCAACGCCGCAGTCGAACTACTTTCACCGCGGAGCAGTTGGAGGAACTGGAAAGAGCCTTTGAACGAACTCACTATCCGGATATATATACGCGAGAGGAATTAGCACAAAGAACAAAGCTAACCGAGGCAAGAGTCCAAGTATGGTTTAGTAACCGAAGAGCGAGATGGCGGAAACAGATGGGTAGCAATCAGCTGACAGCCTTGAACAGTATATTACAAGTGCCACAGGGTATGGGAACGCCCTCTTATATGCTGCACGAGCCTGGGTATCCACTCTCACATAATGCAGACAATCTTTGGCATAGATCGTCTATGGCCCAGTCATTACAGTCATTTGGTCAGACAATAAAACCAGAGAATTCCTACGCCGGTCTTATGGAAAACTATTTATCTCATTCATCACAGCTTCATGGTCTTCCTACACATAGTTCATCCGATCCCCTCTCATCCACTTGGTCATCTCCCGTGTCCACTTCCGTTCCTGCGCTAGGATACACGCCATCTAGTGGCCATTACCATCATTACTCTGATGTCACCAAAAGTACTCTTCATTCATATAACGCTCATATTCCTTCAGTCACAAACATGGAGAGATGTTCAGTTGATGACAGTTTGGTTGCTTTACGTATGAAGTCACGTGAGCATTCCGCCGCTCTCAGTTTGATGCAGGTGGCAGACAACAAAATGGCTACCTCATTTTGA',
                 'ATGGCGATGAGCAGCGGCGGCAGTGGTGGCGGCGTCCCGGAGCAGGAGGATTCCGTGCTGTTCCGGCGCGGCACAGGCCAGAGCGATGATTCTGACATTTGGGATGATACAGCACTGATAAAAGCATATGATAAAGCTGTGGCTTCATTTAAGCATGCTCTAAAGAATGGTGACATTTGTGAAACTTCGGGTAAACCAAAAACCACACCTAAAAGAAAACCTGCTAAGAAGAATAAAAGCCAAAAGAAGAATACTGCAGCTTCCTTACAACAGTGGAAAGTTGGGGACAAATGTTCTGCCATTTGGTCAGAAGACGGTTGCATTTACCCAGCTACCATTGCTTCAATTGATTTTAAGAGAGAAACCTGTGTTGTGGTTTACACTGGATATGGAAATAGAGAGGAGCAAAATCTGTCCGATCTACTTTCCCCAATCTGTGAAGTAGCTAATAATATAGAACAAAATGCTCAAGAGAATGAAAATGAAAGCCAAGTTTCAACAGATGAAAGTGAGAACTCCAGGTCTCCTGGAAATAAATCAGATAACATCAAGCCCAAATCTGCTCCATGGAACTCTTTTCTCCCTCCACCACCCCCCATGCCAGGGCCAAGACTGGGACCAGGAAAGCCAGGTCTAAAATTCAATGGCCCACCACCGCCACCGCCACCACCACCACCCCACTTACTATCATGCTGGCTGCCTCCATTTCCTTCTGGACCACCAATAATTCCCCCACCACCTCCCATATGTCCAGATTCTCTTGATGATGCTGATGCTTTGGGAAGTATGTTAATTTCATGGTACATGAGTGGCTATCATACTGGCTATTATATGGGTTTCAGACAAAATCAAAAAGAAGGAAGGTGCTCACATTCCTTAAATTAA'],
    'sequences_names':['SMN1','PAX3'],
    'promoter_sequence':'GGGCTGGAAGCTACCTTTGACATCATTTCCTCTGCGAATGCATGTATAATTTCTACAGAACCTATTAGAAAGGATCACCCAGCCTCTGCTTTTGTACAACTTTCCCTTAAAAAACTGCCAATTCCACTGCTGTTTGGCCCAATAGTGAGAACTTTTTCCTGCTGCCTCTTGGTGCTTTTGCCTATGGCCCCTATTCTGCCTGCTGAAGACACTCTTGCCAGCATGGACTTAAACCCCTCCAGCTCTGACAATCCTCTTTCTCTTTTGTTTTACATGAAGGGTCTGGCAGCCAAAGCAATCACTCAAAGTTCAAACCTTATCATTTTTTGCTTTGTTCCTCTTGGCCTTGGTTTTGTACATCAGCTTTGAAAATACCATCCCAGGGTTAATGCTGGGGTTAATTTATAACTAAGAGTGCTCTAGTTTTGCAATACAGGACATGCTATAAAAATGGAAAGAT',
    'promoter_name':'TBG',
    'regulator_sequence':'CGATAATCAACCTCTGGATTACAAAATTTGTGAAAGATTGACTGGTATTCTTAACTATGTTGCTCCTTTTACGCTATGTGGATACGCTGCTTTAATGCCTTTGTATCATGCTATTGCTTCCCGTATGGCTTTCATTTTCTCCTCCTTGTATAAATCCTGGTTGCTGTCTCTTTATGAGGAGTTGTGGCCCGTTGTCAGGCAACGTGGCGTGGTGTGCACTGTGTTTGCTGACGCAACCCCCACTGGTTGGGGCATTGCCACCACCTGTCAGCTCCTTTCCGGGACTTTCGCTTTCCCCCTCCCTATTGCCACGGCGGAACTCATCGCCGCCTGCCTTGCCCGCTGCTGGACAGGGGCTCGGCTGTTGGGCACTGACAATTCCGTGGTGTTGTCGGGGAAGCTGACGTCCTTTCCATGGCTGCTCGCCTGTGTTGCCACCTGGATTCTGCGCGGGACGTCCTTCTGCTACGTCCCTTCGGCCCTCAATCCAGCGGACCTTCCTTCCCGCGGCCTGCTGCCGGCTCTGCGGCCTCTTCCGCGTCTTCGCCTTCGCCCTCAGACGAGTCGGATCTCCCTTTGGGCCGCCTCCCCGCATCGG',
    'regulator_name':'WPRE',
    'polya_sequence':'CAGACATGATAAGATACATTGATGAGTTTGGACAAACCACAACTAGAATGCAGTGAAAAAAATGCTTTATTTGTGAAATTTGTGATGCTATTGCTTTATTTGTAACCATTATAAGCTGCAATAAACAAGTTAACAACAACAATTGCATTCATTTTATGTTTCAGGTTCAGGGGGAGGTGTGGGAGGTTTTTTAAAGCAAGTAAAACCTCTACAAATGTGGTA',
    'polya_name':'SV40_late',
    'linkers_sequences':['GGAAGCGGAGAGGGCAGGGGAAGTCTTCTAACATGCGGGGACGTGGAGGAAAATCCCGGCCCC'],
    'linkers_names':['T2A'],
    'fluorescence_sequence':'ATGGTGAGCAAGGGCGAGGAGCTGTTCACCGGGGTGGTGCCCATCCTGGTCGAGCTGGACGGCGACGTAAACGGCCACAAGTTCAGCGTGTCCGGCGAGGGCGAGGGCGATGCCACCTACGGCAAGCTGACCCTGAAGTTCATCTGCACCACCGGCAAGCTGCCCGTGCCCTGGCCCACCCTCGTGACCACCCTGACCTACGGCGTGCAGTGCTTCAGCCGCTACCCCGACCACATGAAGCAGCACGACTTCTTCAAGTCCGCCATGCCCGAAGGCTACGTCCAGGAGCGCACCATCTTCTTCAAGGACGACGGCAACTACAAGACCCGCGCCGAGGTGAAGTTCGAGGGCGACACCCTGGTGAACCGCATCGAGCTGAAGGGCATCGACTTCAAGGAGGACGGCAACATCCTGGGGCACAAGCTGGAGTACAACTACAACAGCCACAACGTCTATATCATGGCCGACAAGCAGAAGAACGGCATCAAGGTGAACTTCAAGATCCGCCACAACATCGAGGACGGCAGCGTGCAGCTCGCCGACCACTACCAGCAGAACACCCCCATCGGCGACGGCCCCGTGCTGCTGCCCGACAACCACTACCTGAGCACCCAGTCCGCCCTGAGCAAAGACCCCAACGAGAAGCGCGATCACATGGTCCTGCTGGAGTTCGTGACCGCCGCCGGGATCACTCTCGGCATGGACGAGCTGTACAAGTAA',
    'fluorescence_name':'EGFP',
    'fluorescence_promoter_sequence':'CTCGACATTGATTATTGACTAGTTATTAATAGTAATCAATTACGGGGTCATTAGTTCATAGCCCATATATGGAGTTCCGCGTTACATAACTTACGGTAAATGGCCCGCCTGGCTGACCGCCCAACGACCCCCGCCCATTGACGTCAATAATGACGTATGTTCCCATAGTAACGCCAATAGGGACTTTCCATTGACGTCAATGGGTGGAGTATTTACGGTAAACTGCCCACTTGGCAGTACATCAAGTGTATCATATGCCAAGTACGCCCCCTATTGACGTCAATGACGGTAAATGGCCCGCCTGGCATTATGCCCAGTACATGACCTTATGGGACTTTCCTACTTGGCAGTACATCTACGTATTAGTCATCGCTATTACCATGGTCGAGGTGAGCCCCACGTTCTGCTTCACTCTCCCCATCTCCCCCCCCTCCCCACCCCCAATTTTGTATTTATTTATTTTTTAATTATTTTGTGCAGCGATGGGGGCGGGGGGGGGGGGGGGGCGCGCGCCAGGCGGGGCGGGGCGGGGCGAGGGGCGGGGCGGGGCGAGGCGGAGAGGTGCGGCGGCAGCCAATCAGAGCGGCGCGCTCCGAAAGTTTCCTTTTATGGCGAGGCGGCGGCGGCGGCGGCCCTATAAAAAGCGAAGCGCGCGGCGGGCGGGAGTCGCTGCGCGCTGCCTTCGCCCCGTGCCCCGCTCCGCCGCCGCCTCGCGCCGCCCGCCCCGGCTCTGACTGACCGCGTTACTCCCACAGGTGAGCGGGCGGGACGGCCCTTCTCCTCCGGGCTGTAATTAGCGCTTGGTTTAATGACGGCTTGTTTCTTTTCTGTGGCTGCGTGAAAGCCTTGAGGGGCTCCGGGAGGGCCCTTTGTGCGGGGGGAGCGGCTCGGGGGGTGCGTGCGTGTGTGTGTGCGTGGGGAGCGCCGCGTGCGGCTCCGCGCTGCCCGGCGGCTGTGAGCGCTGCGGGCGCGGCGCGGGGCTTTGTGCGCTCCGCAGTGTGCGCGAGGGGAGCGCGGCCGGGGGCGGTGCCCCGCGGTGCGGGGGGGGCTGCGAGGGGAACAAAGGCTGCGTGCGGGGTGTGTGCGTGGGGGGGTGAGCAGGGGGTGTGGGCGCGTCGGTCGGGCTGCAACCCCCCCTGCACCCCCCTCCCCGAGTTGCTGAGCACGGCCCGGCTTCGGGTGCGGGGCTCCGTACGGGGCGTGGCGCGGGGCTCGCCGTGCCGGGCGGGGGGTGGCGGCAGGTGGGGGTGCCGGGCGGGGCGGGGCCGCCTCGGGCCGGGGAGGGCTCGGGGGAGGGGCGCGGCGGCCCCCGGAGCGCCGGCGGCTGTCGAGGCGCGGCGAGCCGCAGCCATTGCCTTTTATGGTAATCGTGCGAGAGGGCGCAGGGACTTCCTTTGTCCCAAATCTGTGCGGAGCCGAAATCTGGGAGGCGCCGCCGCACCCCCTCTAGCGGGCGCGGGGCGAAGCGGTGCGGCGCCGGCAGGAAGGAAATGGGCGGGGAGGGCCTTCGTGCGTCGCCGCGCCGCCGTCCCCTTCTCCCTCTCCAGCCTCGGGGCTGTCCGCGGGGGGACGGCTGCCTTCGGGGGGGACGGGGCAGGGCGGGGTTCGGCTTCTGGCGTGTGACCGGCGGCTCTAGAGCCTCTGCTAACCATGTTCATGCCTTCTTCTTTTTCCTACAGCTCCTGGGCAACGTGCTGGTTATTGTGCTGTCTCATCATTTTGGCAAAGAATTG',
    'fluorescence_promoter_name':'CAG',
    'fluorescence_polya_sequence':'CTGTGCCTTCTAGTTGCCAGCCATCTGTTGTTTGCCCCTCCCCCGTGCCTTCCTTGACCCTGGAAGGTGCCACTCCCACTGTCCTTTCCTAATAAAATGAGGAAATTGCATCGCATTGTCTGAGTAGGTGTCATTCTATTCTGGGGGGTGGGGTGGGGCAGGACAGCAAGGGGGAGGATTGGGAAGAGAATAGCAGGCATGCTGGGGA',
    'fluorescence_polya_name':'bGH',
    'fluorescence_linker_sequence':'GGAAGCGGAGAGGGCAGGGGAAGTCTTCTAACATGCGGGGACGTGGAGGAAAATCCCGGCCCC',
    'fluorescence_linker_name':'T2A',
    'selection_marker_sequence':'ATGAGTATTCAACATTTCCGTGTCGCCCTTATTCCCTTTTTTGCGGCATTTTGCCTTCCTGTTTTTGCTCACCCAGAAACGCTGGTGAAAGTAAAAGATGCTGAAGATCAGTTGGGTGCACGAGTGGGTTACATCGAACTGGATCTCAACAGCGGTAAGATCCTTGAGAGTTTTCGCCCCGAAGAACGTTTTCCAATGATGAGCACTTTTAAAGTTCTGCTATGTGGCGCGGTATTATCCCGTATTGACGCCGGGCAAGAGCAACTCGGTCGCCGCATACACTATTCTCAGAATGACTTGGTTGAGTACTCACCAGTCACAGAAAAGCATCTTACGGATGGCATGACAGTAAGAGAATTATGCAGTGCTGCCATAACCATGAGTGATAACACTGCGGCCAACTTACTTCTGACAACGATCGGAGGACCGAAGGAGCTAACCGCTTTTTTGCACAACATGGGGGATCATGTAACTCGCCTTGATCGTTGGGAACCGGAGCTGAATGAAGCCATACCAAACGACGAGCGTGACACCACGATGCCTGTAGCAATGGCAACAACGTTGCGCAAACTATTAACTGGCGAACTACTTACTCTAGCTTCCCGGCAACAATTAATAGACTGGATGGAGGCGGATAAAGTTGCAGGACCACTTCTGCGCTCGGCCCTTCCGGCTGGCTGGTTTATTGCTGATAAATCTGGAGCCGGTGAGCGTGGAAGCCGCGGTATCATTGCAGCACTGGGGCCAGATGGTAAGCCCTCCCGTATCGTAGTTATCTACACGACGGGGAGTCAGGCAACTATGGATGAACGAAATAGACAGATCGCTGAGATAGGTGCCTCACTGATTAAGCATTGGTAA',
    'selection_marker_name':'Ampicillin',
    'restriction_list':['RsaI', 'MnlI', 'AciI', 'AluI', 'BmrI'],
    'optimize':True,
    'transcript_GC':58,
    'poly_len':7
}



project = vb.vector_create_on_dict(metadata_reduced, input_dict, show_plot=False)


#OUTPUT:
    
# Name of project
project['project']

# Graph of the designed vector
vector_plot = project['vector']['graph']
vector_plot.savefig('expression_vector.svg')

# Complete FASTA file of the designed vecotr
fasta = project['vector']['full_fasta']

fasta_sep = project['vector']['fasta']

# Complete GeneBank file of the designed vecotr
gene_bank = project['vector']['genebank']




## genes names
project['transcripts']['sequences']['name']

## proteins sequences
project['transcripts']['sequences']['sequence_aa']

## average codon frequency in the input sequence
project['transcripts']['sequences']['vector_sequence_frequence']

## GC% content in the input sequence
project['transcripts']['sequences']['vector_sequence_GC']

############################################################################

## average codon frequency in the output sequence
project['transcripts']['sequences']['optimized_vector_sequence_frequence']

## GC% content in the output sequence
project['transcripts']['sequences']['optimized_vector_sequence_GC']
    



# RNAi vecotr:
    

input_dict = {

    'project_name':'test_RNAi',
    'vector_type':'ssAAV',
    'vector_function':'rnai',
    'species':'human',
    'promoter_ncrna_name':'U6',
    'promoter_ncrna_sequence':'GAGGGCCTATTTCCCATGATTCCTTCATATTTGCATATACGATACAAGGCTGTTAGAGAGATAATTGGAATTAATTTGACTGTAAACACAAAGATATTAGTACAAAATACGTGACGTAGAAAGTAATAATTTCTTGGGTAGTTTGCAGTTTTAAAATTATGTTTTAAAATGGACTATCATATGCTTACCGTAACTTGAAAGTATTTCGATTTCTTGGCTTTATATATCTTGTGGAAAGGACGAAACACC',
    'rnai_sequence':'',
    'rnai_gene_name':'PAX3',
    'rnai_length':20,
    'overhang_3_prime':'UU',
    'loop_sequence':'TAGTGAAGCCACAGATGTAC',
    'sequences':['ATGTTAACTAGTGGACTTGTAGTAAGCAACATGTTCTCCTATCATCTCGCAGCCTTGGGACTCATGCCGTCATTCCAGATGGAAGGGCGAGGTCGAGTTAATCAGCTAGGAGGTGTTTTCATTAATGGACGGCCACTGCCCAACCATATACGACTGAAGATTGTCGAGCTAGCGGCCCAGGGCGTCCGTCCGTGCGTCATCAGTAGACAGCTGCGGGTGTCACATGGCTGTGTCAGTAAAATACTCCAACGATATCAAGAAACCGGAAGTATCCGACCTGGGGTTATTGGCGGAAGTAAACCAAGGGTCGCAACTCCGGAAGTTGAGAAAAAGATAGAACAATACAAAAAAGATAATCCGGGAATTTTCAGTTGGGAGATTCGGGATCGGCTGCTGAAGGAGGGGATTTGTGACCGCAGCACCGTGCCAAGTGTGAGCTCCATCAGTCGAGTATTACGGAGCAGGTTCCAGAAATGTGATTCTGATGACAATGACAATGACAATGACAATGAGGACGACGATGGCGATGACGGCAGTAACAGTAGTGTGGCAGACAGGTCTGTTAACTTCTCTGTCAGCGGTCTGCTGTCCGACAATAAAAGCGACAAAAGCGACAACGATTCCGATTGTGAATCAGAGCCGGGGCTATCTGTAAAACGGAAGCAACGCCGCAGTCGAACTACTTTCACCGCGGAGCAGTTGGAGGAACTGGAAAGAGCCTTTGAACGAACTCACTATCCGGATATATATACGCGAGAGGAATTAGCACAAAGAACAAAGCTAACCGAGGCAAGAGTCCAAGTATGGTTTAGTAACCGAAGAGCGAGATGGCGGAAACAGATGGGTAGCAATCAGCTGACAGCCTTGAACAGTATATTACAAGTGCCACAGGGTATGGGAACGCCCTCTTATATGCTGCACGAGCCTGGGTATCCACTCTCACATAATGCAGACAATCTTTGGCATAGATCGTCTATGGCCCAGTCATTACAGTCATTTGGTCAGACAATAAAACCAGAGAATTCCTACGCCGGTCTTATGGAAAACTATTTATCTCATTCATCACAGCTTCATGGTCTTCCTACACATAGTTCATCCGATCCCCTCTCATCCACTTGGTCATCTCCCGTGTCCACTTCCGTTCCTGCGCTAGGATACACGCCATCTAGTGGCCATTACCATCATTACTCTGATGTCACCAAAAGTACTCTTCATTCATATAACGCTCATATTCCTTCAGTCACAAACATGGAGAGATGTTCAGTTGATGACAGTTTGGTTGCTTTACGTATGAAGTCACGTGAGCATTCCGCCGCTCTCAGTTTGATGCAGGTGGCAGACAACAAAATGGCTACCTCATTTTGA'],
    'sequences_names':['SMN1'],
    'linkers_sequences':[''],
    'linkers_names':[''],
    'promoter_sequence':'GGGCTGGAAGCTACCTTTGACATCATTTCCTCTGCGAATGCATGTATAATTTCTACAGAACCTATTAGAAAGGATCACCCAGCCTCTGCTTTTGTACAACTTTCCCTTAAAAAACTGCCAATTCCACTGCTGTTTGGCCCAATAGTGAGAACTTTTTCCTGCTGCCTCTTGGTGCTTTTGCCTATGGCCCCTATTCTGCCTGCTGAAGACACTCTTGCCAGCATGGACTTAAACCCCTCCAGCTCTGACAATCCTCTTTCTCTTTTGTTTTACATGAAGGGTCTGGCAGCCAAAGCAATCACTCAAAGTTCAAACCTTATCATTTTTTGCTTTGTTCCTCTTGGCCTTGGTTTTGTACATCAGCTTTGAAAATACCATCCCAGGGTTAATGCTGGGGTTAATTTATAACTAAGAGTGCTCTAGTTTTGCAATACAGGACATGCTATAAAAATGGAAAGAT',
    'promoter_name':'TBG',
    'regulator_sequence':'CGATAATCAACCTCTGGATTACAAAATTTGTGAAAGATTGACTGGTATTCTTAACTATGTTGCTCCTTTTACGCTATGTGGATACGCTGCTTTAATGCCTTTGTATCATGCTATTGCTTCCCGTATGGCTTTCATTTTCTCCTCCTTGTATAAATCCTGGTTGCTGTCTCTTTATGAGGAGTTGTGGCCCGTTGTCAGGCAACGTGGCGTGGTGTGCACTGTGTTTGCTGACGCAACCCCCACTGGTTGGGGCATTGCCACCACCTGTCAGCTCCTTTCCGGGACTTTCGCTTTCCCCCTCCCTATTGCCACGGCGGAACTCATCGCCGCCTGCCTTGCCCGCTGCTGGACAGGGGCTCGGCTGTTGGGCACTGACAATTCCGTGGTGTTGTCGGGGAAGCTGACGTCCTTTCCATGGCTGCTCGCCTGTGTTGCCACCTGGATTCTGCGCGGGACGTCCTTCTGCTACGTCCCTTCGGCCCTCAATCCAGCGGACCTTCCTTCCCGCGGCCTGCTGCCGGCTCTGCGGCCTCTTCCGCGTCTTCGCCTTCGCCCTCAGACGAGTCGGATCTCCCTTTGGGCCGCCTCCCCGCATCGG',
    'regulator_name':'WPRE',
    'polya_sequence':'CAGACATGATAAGATACATTGATGAGTTTGGACAAACCACAACTAGAATGCAGTGAAAAAAATGCTTTATTTGTGAAATTTGTGATGCTATTGCTTTATTTGTAACCATTATAAGCTGCAATAAACAAGTTAACAACAACAATTGCATTCATTTTATGTTTCAGGTTCAGGGGGAGGTGTGGGAGGTTTTTTAAAGCAAGTAAAACCTCTACAAATGTGGTA',
    'polya_name':'SV40_late',
    'fluorescence_sequence':'ATGGTGAGCAAGGGCGAGGAGCTGTTCACCGGGGTGGTGCCCATCCTGGTCGAGCTGGACGGCGACGTAAACGGCCACAAGTTCAGCGTGTCCGGCGAGGGCGAGGGCGATGCCACCTACGGCAAGCTGACCCTGAAGTTCATCTGCACCACCGGCAAGCTGCCCGTGCCCTGGCCCACCCTCGTGACCACCCTGACCTACGGCGTGCAGTGCTTCAGCCGCTACCCCGACCACATGAAGCAGCACGACTTCTTCAAGTCCGCCATGCCCGAAGGCTACGTCCAGGAGCGCACCATCTTCTTCAAGGACGACGGCAACTACAAGACCCGCGCCGAGGTGAAGTTCGAGGGCGACACCCTGGTGAACCGCATCGAGCTGAAGGGCATCGACTTCAAGGAGGACGGCAACATCCTGGGGCACAAGCTGGAGTACAACTACAACAGCCACAACGTCTATATCATGGCCGACAAGCAGAAGAACGGCATCAAGGTGAACTTCAAGATCCGCCACAACATCGAGGACGGCAGCGTGCAGCTCGCCGACCACTACCAGCAGAACACCCCCATCGGCGACGGCCCCGTGCTGCTGCCCGACAACCACTACCTGAGCACCCAGTCCGCCCTGAGCAAAGACCCCAACGAGAAGCGCGATCACATGGTCCTGCTGGAGTTCGTGACCGCCGCCGGGATCACTCTCGGCATGGACGAGCTGTACAAGTAA',
    'fluorescence_name':'EGFP',
    'fluorescence_linker_sequence':'GGAAGCGGAGAGGGCAGGGGAAGTCTTCTAACATGCGGGGACGTGGAGGAAAATCCCGGCCCC',
    'fluorescence_linker_name':'T2A',
    'selection_marker_sequence':'ATGAGTATTCAACATTTCCGTGTCGCCCTTATTCCCTTTTTTGCGGCATTTTGCCTTCCTGTTTTTGCTCACCCAGAAACGCTGGTGAAAGTAAAAGATGCTGAAGATCAGTTGGGTGCACGAGTGGGTTACATCGAACTGGATCTCAACAGCGGTAAGATCCTTGAGAGTTTTCGCCCCGAAGAACGTTTTCCAATGATGAGCACTTTTAAAGTTCTGCTATGTGGCGCGGTATTATCCCGTATTGACGCCGGGCAAGAGCAACTCGGTCGCCGCATACACTATTCTCAGAATGACTTGGTTGAGTACTCACCAGTCACAGAAAAGCATCTTACGGATGGCATGACAGTAAGAGAATTATGCAGTGCTGCCATAACCATGAGTGATAACACTGCGGCCAACTTACTTCTGACAACGATCGGAGGACCGAAGGAGCTAACCGCTTTTTTGCACAACATGGGGGATCATGTAACTCGCCTTGATCGTTGGGAACCGGAGCTGAATGAAGCCATACCAAACGACGAGCGTGACACCACGATGCCTGTAGCAATGGCAACAACGTTGCGCAAACTATTAACTGGCGAACTACTTACTCTAGCTTCCCGGCAACAATTAATAGACTGGATGGAGGCGGATAAAGTTGCAGGACCACTTCTGCGCTCGGCCCTTCCGGCTGGCTGGTTTATTGCTGATAAATCTGGAGCCGGTGAGCGTGGAAGCCGCGGTATCATTGCAGCACTGGGGCCAGATGGTAAGCCCTCCCGTATCGTAGTTATCTACACGACGGGGAGTCAGGCAACTATGGATGAACGAAATAGACAGATCGCTGAGATAGGTGCCTCACTGATTAAGCATTGGTAA',
    'selection_marker_name':'Ampicillin',
    'restriction_list':['RsaI', 'MnlI', 'AciI', 'AluI', 'BmrI'],
    'optimize':True,
    'transcript_GC':58,
    'poly_len':7
}  


project = vb.vector_create_on_dict(metadata_reduced, input_dict, show_plot=False)



# Name of project
project['project']



# Graph of the designed vector
vector_plot = project['vector']['graph']
vector_plot.savefig('RNAi_vector.svg')

# Complete FASTA file of the designed vecotr
fasta = project['vector']['full_fasta']


# The FASTA file is divided into particular elements of the designed vector
fasta_sep = project['vector']['fasta']


# Complete GeneBank file of the designed vecotr
gene_bank = project['vector']['genebank']


# Top 1 designed RNAi in shRNA form information
# RNAi name
project['rnai']['name']

# RNAi sequence
project['rnai']['sequence']

# RNAi prediction structure
rnai_prediction = project['rnai']['figure']
rnai_prediction.savefig('rnai_predicted_structure.svg')


## *if occur user-defined sequences for additional expression in the plasmid vector
## genes names
project['transcripts']['sequences']['name']

## proteins sequences
project['transcripts']['sequences']['sequence_aa']

## average codon frequency in the input sequence
project['transcripts']['sequences']['vector_sequence_frequence']

## GC% content in the input sequence
project['transcripts']['sequences']['vector_sequence_GC']

############################################################################

## average codon frequency in the output sequence
project['transcripts']['sequences']['optimized_vector_sequence_frequence']

## GC% content in the output sequence
project['transcripts']['sequences']['optimized_vector_sequence_GC']



# transcription - mRNA:


input_dict = {
    
    'project_name':'test_invitro_transcription_mRNA',
    'vector_type':'transcription',
    'vector_function':'mrna',
    'species':'human',
    'sequences':['ATGTTAACTAGTGGACTTGTAGTAAGCAACATGTTCTCCTATCATCTCGCAGCCTTGGGACTCATGCCGTCATTCCAGATGGAAGGGCGAGGTCGAGTTAATCAGCTAGGAGGTGTTTTCATTAATGGACGGCCACTGCCCAACCATATACGACTGAAGATTGTCGAGCTAGCGGCCCAGGGCGTCCGTCCGTGCGTCATCAGTAGACAGCTGCGGGTGTCACATGGCTGTGTCAGTAAAATACTCCAACGATATCAAGAAACCGGAAGTATCCGACCTGGGGTTATTGGCGGAAGTAAACCAAGGGTCGCAACTCCGGAAGTTGAGAAAAAGATAGAACAATACAAAAAAGATAATCCGGGAATTTTCAGTTGGGAGATTCGGGATCGGCTGCTGAAGGAGGGGATTTGTGACCGCAGCACCGTGCCAAGTGTGAGCTCCATCAGTCGAGTATTACGGAGCAGGTTCCAGAAATGTGATTCTGATGACAATGACAATGACAATGACAATGAGGACGACGATGGCGATGACGGCAGTAACAGTAGTGTGGCAGACAGGTCTGTTAACTTCTCTGTCAGCGGTCTGCTGTCCGACAATAAAAGCGACAAAAGCGACAACGATTCCGATTGTGAATCAGAGCCGGGGCTATCTGTAAAACGGAAGCAACGCCGCAGTCGAACTACTTTCACCGCGGAGCAGTTGGAGGAACTGGAAAGAGCCTTTGAACGAACTCACTATCCGGATATATATACGCGAGAGGAATTAGCACAAAGAACAAAGCTAACCGAGGCAAGAGTCCAAGTATGGTTTAGTAACCGAAGAGCGAGATGGCGGAAACAGATGGGTAGCAATCAGCTGACAGCCTTGAACAGTATATTACAAGTGCCACAGGGTATGGGAACGCCCTCTTATATGCTGCACGAGCCTGGGTATCCACTCTCACATAATGCAGACAATCTTTGGCATAGATCGTCTATGGCCCAGTCATTACAGTCATTTGGTCAGACAATAAAACCAGAGAATTCCTACGCCGGTCTTATGGAAAACTATTTATCTCATTCATCACAGCTTCATGGTCTTCCTACACATAGTTCATCCGATCCCCTCTCATCCACTTGGTCATCTCCCGTGTCCACTTCCGTTCCTGCGCTAGGATACACGCCATCTAGTGGCCATTACCATCATTACTCTGATGTCACCAAAAGTACTCTTCATTCATATAACGCTCATATTCCTTCAGTCACAAACATGGAGAGATGTTCAGTTGATGACAGTTTGGTTGCTTTACGTATGAAGTCACGTGAGCATTCCGCCGCTCTCAGTTTGATGCAGGTGGCAGACAACAAAATGGCTACCTCATTTTGA'],
    'sequences_names':['SMN1'],
    'linkers_sequences':[''],
    'linkers_names':[''],
    'utr5_sequence':'GGGCTGGAAGCTACCTTTGACATCATTTCCTCTGCGAATGCATGTATAATTTCTACAGAACCTATTAGAAAGGATCACCCAGCCTCTGCTTTTGTACAACTTTCCCTTAAAAAACTGCCAATTCCACTGCTGTTTGGCCCAATAGTGAGAACTTTTTCCTGCTGCCTCTTGGTGCTTTTGCCTATGGCCCCTATTCTGCCTGCTGAAGACACTCTTGCCAGCATGGACTTAAACCCCTCCAGCTCTGACAATCCTCTTTCTCTTTTGTTTTACATGAAGGGTCTGGCAGCCAAAGCAATCACTCAAAGTTCAAACCTTATCATTTTTTGCTTTGTTCCTCTTGGCCTTGGTTTTGTACATCAGCTTTGAAAATACCATCCCAGGGTTAATGCTGGGGTTAATTTATAACTAAGAGTGCTCTAGTTTTGCAATACAGGACATGCTATAAAAATGGAAAGAT',
    'utr5_name':'SMN1',
    'utr3_sequence':'CAGACATGATAAGATACATTGATGAGTTTGGACAAACCACAACTAGAATGCAGTGAAAAAAATGCTTTATTTGTGAAATTTGTGATGCTATTGCTTTATTTGTAACCATTATAAGCTGCAATAAACAAGTTAACAACAACAATTGCATTCATTTTATGTTTCAGGTTCAGGGGGAGGTGTGGGAGGTTTTTTAAAGCAAGTAAAACCTCTACAAATGTGGTA',
    'utr3_name':'KIT',
    'polya_tail_x':50,
    'selection_marker_sequence':'ATGAGTATTCAACATTTCCGTGTCGCCCTTATTCCCTTTTTTGCGGCATTTTGCCTTCCTGTTTTTGCTCACCCAGAAACGCTGGTGAAAGTAAAAGATGCTGAAGATCAGTTGGGTGCACGAGTGGGTTACATCGAACTGGATCTCAACAGCGGTAAGATCCTTGAGAGTTTTCGCCCCGAAGAACGTTTTCCAATGATGAGCACTTTTAAAGTTCTGCTATGTGGCGCGGTATTATCCCGTATTGACGCCGGGCAAGAGCAACTCGGTCGCCGCATACACTATTCTCAGAATGACTTGGTTGAGTACTCACCAGTCACAGAAAAGCATCTTACGGATGGCATGACAGTAAGAGAATTATGCAGTGCTGCCATAACCATGAGTGATAACACTGCGGCCAACTTACTTCTGACAACGATCGGAGGACCGAAGGAGCTAACCGCTTTTTTGCACAACATGGGGGATCATGTAACTCGCCTTGATCGTTGGGAACCGGAGCTGAATGAAGCCATACCAAACGACGAGCGTGACACCACGATGCCTGTAGCAATGGCAACAACGTTGCGCAAACTATTAACTGGCGAACTACTTACTCTAGCTTCCCGGCAACAATTAATAGACTGGATGGAGGCGGATAAAGTTGCAGGACCACTTCTGCGCTCGGCCCTTCCGGCTGGCTGGTTTATTGCTGATAAATCTGGAGCCGGTGAGCGTGGAAGCCGCGGTATCATTGCAGCACTGGGGCCAGATGGTAAGCCCTCCCGTATCGTAGTTATCTACACGACGGGGAGTCAGGCAACTATGGATGAACGAAATAGACAGATCGCTGAGATAGGTGCCTCACTGATTAAGCATTGGTAA',
    'selection_marker_name':'Ampicillin',
    'restriction_list':['RsaI', 'MnlI', 'AciI', 'AluI', 'BmrI'],
    'optimize':True,
    'transcript_GC':58,
    'poly_len':7
}

project = vb.vector_create_on_dict(metadata_reduced, input_dict, show_plot=False)


##### Output:

# Name of project
project['project']


# Graph of the designed vector
vector_plot = project['vector']['graph']
vector_plot.savefig('transcription_mrna_vector.svg')


# Complete FASTA file of the designed vecotr
fasta = project['vector']['full_fasta']


# The FASTA file is divided into particular elements of the designed vector
fasta_sep = project['vector']['fasta']


# Complete GeneBank file of the designed vecotr
gene_bank = project['vector']['genebank']


## genes names
project['transcripts']['sequences']['name']

## proteins sequences
project['transcripts']['sequences']['sequence_aa']

## average codon frequency in the input sequence
project['transcripts']['sequences']['vector_sequence_frequence']

## GC% content in the input sequence
project['transcripts']['sequences']['vector_sequence_GC']

############################################################################

## average codon frequency in the output sequence
project['transcripts']['sequences']['optimized_vector_sequence_frequence']

## GC% content in the output sequence
project['transcripts']['sequences']['optimized_vector_sequence_GC']



# transcription - RNAi



##### Example dictionary:



input_dict = {

    'project_name':'test_invitro_transcription_RNAi',
    'vector_type':'transcription',
    'vector_function':'rnai',
    'species':'human',
    'rnai_sequence':'',
    'rnai_length':20,
    'overhang_3_prime':'UU',
    'rnai_gene_name':'KIT',
    'loop_sequence':'TAGTGAAGCCACAGATGTAC',
    'selection_marker_sequence':'ATGAGTATTCAACATTTCCGTGTCGCCCTTATTCCCTTTTTTGCGGCATTTTGCCTTCCTGTTTTTGCTCACCCAGAAACGCTGGTGAAAGTAAAAGATGCTGAAGATCAGTTGGGTGCACGAGTGGGTTACATCGAACTGGATCTCAACAGCGGTAAGATCCTTGAGAGTTTTCGCCCCGAAGAACGTTTTCCAATGATGAGCACTTTTAAAGTTCTGCTATGTGGCGCGGTATTATCCCGTATTGACGCCGGGCAAGAGCAACTCGGTCGCCGCATACACTATTCTCAGAATGACTTGGTTGAGTACTCACCAGTCACAGAAAAGCATCTTACGGATGGCATGACAGTAAGAGAATTATGCAGTGCTGCCATAACCATGAGTGATAACACTGCGGCCAACTTACTTCTGACAACGATCGGAGGACCGAAGGAGCTAACCGCTTTTTTGCACAACATGGGGGATCATGTAACTCGCCTTGATCGTTGGGAACCGGAGCTGAATGAAGCCATACCAAACGACGAGCGTGACACCACGATGCCTGTAGCAATGGCAACAACGTTGCGCAAACTATTAACTGGCGAACTACTTACTCTAGCTTCCCGGCAACAATTAATAGACTGGATGGAGGCGGATAAAGTTGCAGGACCACTTCTGCGCTCGGCCCTTCCGGCTGGCTGGTTTATTGCTGATAAATCTGGAGCCGGTGAGCGTGGAAGCCGCGGTATCATTGCAGCACTGGGGCCAGATGGTAAGCCCTCCCGTATCGTAGTTATCTACACGACGGGGAGTCAGGCAACTATGGATGAACGAAATAGACAGATCGCTGAGATAGGTGCCTCACTGATTAAGCATTGGTAA',
    'selection_marker_name':'Ampicillin'
}


project = vb.vector_create_on_dict(metadata_reduced, input_dict, show_plot=False)


# Name of project
project['project']



# Graph of the designed vector
vector_plot = project['vector']['graph']
vector_plot.savefig('transcription_rnai_vector.svg')



# Complete FASTA file of the designed vecotr
fasta = project['vector']['full_fasta']



# The FASTA file is divided into particular elements of the designed vector
fatsa_sep = project['vector']['fasta']


# Complete GeneBank file of the designed vecotr
gene_bank = project['vector']['genebank']



# Top 1 designed RNAi in shRNA form information
# RNAi name
project['rnai']['name']

# RNAi sequence
project['rnai']['sequence']

# RNAi prediction structure
rnai_prediction = project['rnai']['figure']
rnai_prediction.savefig('rnai_predicted_structure_transcription.svg')




# denovo - RNAi - siRNA



##### Example dictionary:




input_dict = {

    'project_name':'test_siRNA',
    'species':'human',
    'rnai_type':'sirna',
    'rnai_sequence':'',
    'rnai_length':20,
    'overhang_3_prime':'UU',
    'rnai_gene_name':'KIT',
    'loop_sequence':''
}

project = vb.create_sequence_from_dict(metadata_reduced, input_dict, show_plot=False)


# Name of project
project['project']



# Graph of the designed vector
structure_graph = project['rnai']['figure']


# Predicted RNAi data
data = project['rnai']['full_data']




# denovo - RNAi - shRNA



##### Example dictionary:




input_dict = {

    'project_name':'test_siRNA',
    'species':'human',
    'rnai_type':'sh',
    'rnai_sequence':'',
    'rnai_length':20,
    'overhang_3_prime':'UU',
    'rnai_gene_name':'KIT',
    'loop_sequence':'TAGTGAAGCCACAGATGTAC'
}

project = vb.create_sequence_from_dict(metadata_reduced, input_dict, show_plot=False)


# Name of project
project['project']



# Graph of the designed vector
structure_graph = project['rnai']['figure']


# Predicted RNAi data
data = project['rnai']['full_data']





# denovo - mRNA



##### Example dictionary:


input_dict = {

    'project_name':'test_mRNA',
    'species':'human',
    'sequence':'ATGTTAACTAGTGGACTTGTAGTAAGCAACATGTTCTCCTATCATCTCGCAGCCTTGGGACTCATGCCGTCATTCCAGATGGAAGGGCGAGGTCGAGTTAATCAGCTAGGAGGTGTTTTCATTAATGGACGGCCACTGCCCAACCATATACGACTGAAGATTGTCGAGCTAGCGGCCCAGGGCGTCCGTCCGTGCGTCATCAGTAGACAGCTGCGGGTGTCACATGGCTGTGTCAGTAAAATACTCCAACGATATCAAGAAACCGGAAGTATCCGACCTGGGGTTATTGGCGGAAGTAAACCAAGGGTCGCAACTCCGGAAGTTGAGAAAAAGATAGAACAATACAAAAAAGATAATCCGGGAATTTTCAGTTGGGAGATTCGGGATCGGCTGCTGAAGGAGGGGATTTGTGACCGCAGCACCGTGCCAAGTGTGAGCTCCATCAGTCGAGTATTACGGAGCAGGTTCCAGAAATGTGATTCTGATGACAATGACAATGACAATGACAATGAGGACGACGATGGCGATGACGGCAGTAACAGTAGTGTGGCAGACAGGTCTGTTAACTTCTCTGTCAGCGGTCTGCTGTCCGACAATAAAAGCGACAAAAGCGACAACGATTCCGATTGTGAATCAGAGCCGGGGCTATCTGTAAAACGGAAGCAACGCCGCAGTCGAACTACTTTCACCGCGGAGCAGTTGGAGGAACTGGAAAGAGCCTTTGAACGAACTCACTATCCGGATATATATACGCGAGAGGAATTAGCACAAAGAACAAAGCTAACCGAGGCAAGAGTCCAAGTATGGTTTAGTAACCGAAGAGCGAGATGGCGGAAACAGATGGGTAGCAATCAGCTGACAGCCTTGAACAGTATATTACAAGTGCCACAGGGTATGGGAACGCCCTCTTATATGCTGCACGAGCCTGGGTATCCACTCTCACATAATGCAGACAATCTTTGGCATAGATCGTCTATGGCCCAGTCATTACAGTCATTTGGTCAGACAATAAAACCAGAGAATTCCTACGCCGGTCTTATGGAAAACTATTTATCTCATTCATCACAGCTTCATGGTCTTCCTACACATAGTTCATCCGATCCCCTCTCATCCACTTGGTCATCTCCCGTGTCCACTTCCGTTCCTGCGCTAGGATACACGCCATCTAGTGGCCATTACCATCATTACTCTGATGTCACCAAAAGTACTCTTCATTCATATAACGCTCATATTCCTTCAGTCACAAACATGGAGAGATGTTCAGTTGATGACAGTTTGGTTGCTTTACGTATGAAGTCACGTGAGCATTCCGCCGCTCTCAGTTTGATGCAGGTGGCAGACAACAAAATGGCTACCTCATTTTGA',
    'sequence_name':'SMN1',
    'restriction_list':['RsaI', 'MnlI', 'AciI', 'AluI', 'BmrI'],
    'optimize':True,
    'transcript_GC':58,
    'poly_len':7
    
}



project = vb.create_sequence_from_dict(metadata_reduced, input_dict, show_plot=False)


# Name of project
project['project']





###############################################################################



# Itterative vector build from dict test!

from jbst.tests.jseq_tests import *


metadata_reduced = load_metadata(linkers = False, 
                                    loops = False, 
                                    regulators = False, 
                                    fluorescent_tag = False, 
                                    promoters = False, 
                                    polya = False, 
                                    marker = False, 
                                    utr5 = False, 
                                    utr3 = False) 

passed, errored = test_vectors(tests_dictonaries, tests_dictonaries_names, test_transcription_denovo_list, test_transcription_denovo_list_names, non, metadata, show_plot = False)




# if errored check using vector_create_on_dict() function change only 'input_dict' to dictionary name included in errored list:
    
project = vector_create_on_dict(metadata_reduced, input_dict_rnai_1, show_plot=False)

project['vector']['graph']

pd.DataFrame(project['rnai']['sequence'])



