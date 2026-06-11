import sys
import re
from Bio.Seq import Seq
import os
import pandas as pd

os.chdir(path=str(sys.argv[1]))    
file_list = list()
name_list = list()

#filtering reads with more than 200bp match and convert dna sequence to amino acid
for x in os.listdir():  
    if x.endswith("seq2.txt"):
        print("Running : "+x)
        try:
            with open(x, 'r') as f:
                P1P1p_list = list()
                for line in f:
                    line_list = line.split(',')
                    MD_tag = re.sub('[ATCG]',",",line_list[0])
                    MD_list = MD_tag.split(',')
                    MD_num = list()
                    #remove empty entries, and insert into new list
                    for i in range(0, len(MD_list)):
                        if MD_list[i] != '':
                            MD_num.append(int(MD_list[i]))
                            total_match = sum(MD_num)
                    if total_match > 250:
                        read = Seq(line_list[2])
                        motif = "GTCACATCTGACCGT"
                        codon_start = read.find(motif.upper())
                        aa_seq = read[codon_start+len(motif)-9:codon_start+len(motif)+189]
                        aa_site = aa_seq.translate()
                        if (str(aa_site)[0:3] == "SDR") and (str(aa_site)[-3:] == "NTI"):
                            P1P1p_list.append((str(aa_site[3:63])+","+str(aa_seq)))
                        #print(aa_site, codon_start+len(motif))
                        #print(line_list[0], line_list[2])  
                print("Total Usable Reads : "+str(len(P1P1p_list)))
                file_list.insert(len(file_list), P1P1p_list)
                name_list.insert(len(name_list), x)
        except OSError as oserr:
            print('OS error:', oserr) 

df = pd.DataFrame(file_list).transpose()
df.columns = name_list

#extract aa and dna sequence into another dataframe
for i in df.columns:
    aa_col_name = i[0:26]+"_aa"
    dna_col_name = i[0:26]+"_dna"
    df[[aa_col_name, dna_col_name]] = df[i].str.split(',', expand=True)
    df = df.drop(i, axis=1)

#check to make sure all files are saved into data frame
print(df.columns) 

#count aa frequencies and save into file aa_count.csv
all_aa = pd.DataFrame()
for i in df.columns:
    if i[-2:] == "aa":
        count_aa = df[i].value_counts(dropna=True)
        all_aa = pd.concat([all_aa, count_aa], axis =1)
all_aa = all_aa.fillna(0)
all_aa.to_csv("aa_count.csv")  
