# -*- coding: utf-8 -*-
"""
Created on Fri Sep  6 14:33:44 2019

@author: samsung
"""

import os
publicDir = 'E:/publicData/'
org = 'Hs'
orgFull = 'HUMAN'
outDir = publicDir + 'uniprotFunctions{}/'.format(org)

# extract the Accession and corresponding Function
txts = os.listdir(outDir)
accession = []
function = []
for root, dirs, files in os.walk(outDir):
    for i in files:
        accession.append(i.split('.')[0])
        f = open(os.path.join(root,i), 'r')
        function.append(f.readline())
        f.close()
print(len(accession) == len(function))

# 将所有下载好的accession-function, write to txt
fw = open(os.path.join(publicDir, 'AccFunction{}.txt'.format(org)), 'w')
fw.write('\t'.join(['Accession','Function']))
fw.write('\n')    
for line in zip(accession, function):
    fw.write('\t'.join(line))
    fw.write('\n')
fw.close()

# Merge with the AccSymProtname.txt
f1 = open(os.path.join(publicDir,'AccSymPro{}.txt'.format(org)),'r')
f2 = open(os.path.join(publicDir,'AccFunction{}.txt'.format(org)),'r')
flist1 = f1.readlines()[1:]
flist2 = f2.readlines()[1:]
f1.close()
f2.close()
# 定义各属性数据存储列表并存储数据
acc1 = []
sym = []
protname = []
acc2 = []
func = []
for i in flist1:
    acc1.append(i.split('\t')[0])
    sym.append(i.split('\t')[1])
    protname.append(i.split('\t')[2].strip())
for i in flist2:
    acc2.append(i.split('\t')[0])
    func.append(i.split('\t')[1].strip())
# merge
flist3 = []
for i in range(len(acc1)):
    if acc1[i] in acc2:
        x = acc2.index(acc1[i])
        s = '\t'.join([acc1[i], sym[i], protname[i], func[x]])
        s += '\n'
        flist3.append(s)
    else:
        pass
f3 = open(os.path.join(publicDir,'AccSymProFunc{}.txt'.format(org)), 'w')
f3.write('\t'.join(['Accession','Symbol','ProteinName', 'Function']))
f3.write('\n') 
f3.writelines(flist3)
f3.close()
flist3[0]
# count the number of accession without Gene Symbol or without Function summary
syms = []
funcs = []
for i in flist3:
    syms.append(i.split('\t')[1])
    funcs.append(i.split('\t')[-1].strip())
syms.count('Nosymbol')    # 149
funcs.count('Not found function summary.') # 4245


