# -*- coding: utf-8 -*-
"""
Created on Thu Sep  5 22:53:39 2019

@author: samsung
"""
import os

# extract the protein id lines  
publicDir = 'E:/publicData/' 
org = 'Rat'
orgFull = 'RAT'

fas = open(os.path.join(publicDir, 'swiss_{}_201907.fasta'.format(org)), 'r')
protid = open(os.path.join(publicDir, 'swiss_{}_201907_ID.txt'.format(org)),'w')
for line in fas.readlines():
    if line.startswith('>'):
        protid.write(line)
fas.close()
protid.close()

# get the Accession and Symbol from ID lines
with open(os.path.join(publicDir, 'swiss_{}_201907_ID.txt'.format(org)),'r') as f:
    protid = [line.rstrip('\n') for line in f]

seq = ''.join(protid)
seq.count('|')
print('If the number of "|" is Twice of entries: {}'.format(seq.count('|') == 2*len(protid)))
# so it's OK to extract the Accession from between 2 "|"

print('If the number of "GN=" is Same of entries: {}'.format(seq.count('GN=') == len(protid)))
# so not every gene has the GN symbol

print('If the number of "OS=" is Same of entries: {}'.format(seq.count('OS=') == len(protid)))
print('If the number of "{}" is Same of entries: {}'.format(orgFull,seq.count('_{}'.format(orgFull)) == len(protid)))
# So it's ok to extract the protein name from between '_RAT' and ' OS='

# extract the info, and write in tab delimited txt
accession = []
symbol = []
protname = []
for i in protid:
    if 'GN=' in i:
        accession.append(i.split('|')[1])
        symbol.append(i.split('GN=')[1].split(' ')[0])
        protname.append(i.split('_{} '.format(orgFull))[1].split(' OS=')[0])
    else:
        accession.append(i.split('|')[1])
        symbol.append('Nosymbol')
        protname.append(i.split('_{} '.format(orgFull))[1].split(' OS=')[0])
        
# write to txt
fw = open(os.path.join(publicDir, 'AccSymPro{}.txt'.format(org)), 'w')
fw.write('\t'.join(['Accession','Symbol','ProteinName']))
fw.write('\n')    
for line in zip(accession, symbol, protname):
    fw.write('\t'.join(line))
    fw.write('\n')
fw.close()

fw = open(os.path.join(publicDir, 'Accession{}.txt'.format(org)), 'w')
fw.write('Accession')
fw.write('\n')    
for line in accession:
    fw.write(line)
    fw.write('\n')
fw.close()