from Bio import SeqIO
from re import finditer
from tqdm import tqdm

import re
#import sys
import pandas as pd
        

def printProgressBar(iteration, total, prefix = '', suffix = '', decimals = 1, length = 100, fill = '█'):
    percent = ("{0:." + str(decimals) + "f}").format(100 * (iteration / float(total)))
    filledLength = int(length * iteration // total)
    bar = fill * filledLength + '-' * (length - filledLength)
    print('\r%s |%s| %s%% %s' % (prefix, bar, percent, suffix), end = '\r')
    # Print New Line on Complete
    if iteration == total: 
        print()

def complementary(let):
    # A-T, C-G
    if let == 'A':
        return 'T'
    if let == 'T':
        return 'A'
    if let == 'C':
        return 'G'
    if let == 'G':
        return 'C'
    if let == 'N':
        return 'N'
    
def check_loc(new, prev):
    if new-prev>=1:
        return 1
    else:
        return 0
    
def check_12G(seq):
    if len([1 for x in seq.upper() if x=='G'])>=12:
        return 1
    else:
        return 0

def check_12C(seq):
    if len([1 for x in seq.upper() if x=='C'])>=12:
        return 1
    else:
        return 0
    
def check_ggg(seq):
    if 'GGG' in seq:
        cont = True
    else:
        cont = False
    return cont


def check_ccc(seq):
    if 'CCC' in seq:
        cont = True
    else:
        cont = False
    return cont


#CE
pattern1 = '(([G]{3,}\w{1,12}){3,})([G]{3,})'
#BULGE
pattern2 = '(([G][ATC]{0,1}[G][ATC]{0,1}[G][ATCG]{1,3}){3,})([G][ATC]{0,1}[G][ATC]{0,1}[G])'
#IRREGULAR
pattern3 = '(([G]{1,2}[ATC]{1,2}){7,})([G]{2})'

#CE
pattern4 = '(([C]{3,}\w{1,12}){3,})([C]{3,})'
#BULGE
pattern5 = '(([C][TAG]{0,1}[C][TAG]{0,1}[C][TACG]{1,3}){3,})([C][TAG]{0,1}[C][ATG]{0,1}[C])'
#IRREGULAR
pattern6 = '(([C]{1,2}[ATG]{1,2}){7,})([C]{2})'


region = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7',
            'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13',
            'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19',
            'chr20', 'chr21', 'chr22', 'chrX', 'chrY']


start_loc = []
end_loc = []
pattern = []
chrN = []
strand = []


seq_dict = {}

with open('../data/hg19.fa', mode='r') as handle:
    # process each record in .fa file if there's more than one
    for record in SeqIO.parse(handle, 'fasta'):
        identifier = record.id
        description = record.description
        # for any human genome assembly, there might be alternative chromosome contigs, 
        # but we only save the conventional ones from the "region" list
        if identifier in region:
            print(description)
            sequence = record.seq

            seq_dict[identifier] = str(sequence)


for ck in tqdm(region):
    #printProgressBar(i+1, len(region))
    #seq = seq_dict[ck]
    for p,t in zip([pattern1,pattern2,pattern3], ['classical_extended','bulge','agro']): #, 'agro100'
        prev = 0
        for m in finditer(p, seq_dict[ck], flags=re.IGNORECASE):
            if check_loc(m.start(), prev) and check_12G(m.group(0)):
                if t=='bulge':
                    if not check_ggg(m.group(0)):
                        continue
                start_loc.append(m.start())
                end_loc.append(m.end())
                pattern.append(t)
                chrN.append(ck)
                strand.append('+')
                prev = m.end()

    for p,t in zip([pattern4,pattern5,pattern6], ['classical_extended','bulge','agro']): #, 'agro100'
        prev = 0
        for m in finditer(p, seq_dict[ck], flags=re.IGNORECASE):
            if check_loc(m.start(), prev) and check_12C(m.group(0)):
                if t=='bulge':
                    if not check_ccc(m.group(0)):
                        continue
                start_loc.append(m.start())
                end_loc.append(m.end())
                pattern.append(t)
                chrN.append(ck)
                strand.append('-')
                prev = m.end()

pqs_id = []

for i in range(len(pattern)):
    if pattern[i]=='classical_extended':
        pqs_id.append(f'pqs{i+1}_ce')
    if pattern[i]=='bulge':
        pqs_id.append(f'pqs{i+1}_b')
    if pattern[i]=='agro':
        pqs_id.append(f'pqs{i+1}_ag')

with open('PQS_all.bed', 'w') as f:
    f.write('chr\tstart\tend\ttype\tn\tstrand\n')
    for st, en, p, c, s in zip(start_loc, end_loc, pqs_id, chrN, strand):
        f.write(f'{c}\t{st}\t{en}\t{p}\t.\t{s}\n')

g4type = {}
g4type['ce'] = 'Classical extended'
g4type['b'] = 'Bulge'
g4type['ag'] = 'Irregular'

df_pqs = pd.DataFrame({'id':pqs_id})
df_pqs['type'] = [g4type[x.split('_')[1]] for x in pqs_id]

N_all = len(df_pqs)
N_ce = len(df_pqs[df_pqs['type']=='Classical extended'])
N_b = len(df_pqs[df_pqs['type']=='Bulge'])
N_ag = len(df_pqs[df_pqs['type']=='Irregular'])

print('Total:',N_all)
print(f'Classical extended: {N_ce}, {N_ce*100/N_all}%')
print(f'Bulge: {N_b}, {N_b*100/N_all}%')
print(f'Irregular: {N_ag}, {N_ag*100/N_all}%')

