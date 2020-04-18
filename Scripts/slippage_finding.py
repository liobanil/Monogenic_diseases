#These scripts creates a dictionary, with chromosome number as keys.
#Each with counts per biotype per range
#An addtional key 'counts' stores total number of counts per biotype.


#####################################################################
# Input: file with only the ranges by chromosome

import os
# load human genome annotation
human = {}
#For use modify the paths to the folder which cotains all cromosome files
for chrom in os.listdir('/home/jgonzale/Documentos/4_sem/humana/all_anot'): # for every chromosome file
    with open( os.path.join('/home/jgonzale/Documentos/4_sem/humana/all_anot', chrom),'r' ) as file:
        lines = { (int(i.split('\t')[3]), int(i.split('\t')[4]), i.split('\t')[2] ):[]  for i in file if not i.startswith('#')} # keep coordinates and range type
    human[int(chrom.split('.')[0][-2:])] = lines #{chromosome_number: [[start, end, type], ...] , ...}
    print('Chr',chrom.split('.')[0][-2:], 'ready...' )
#####################################################################


#####################################################################
# Input: files cpntaining chromosome number, start and end of range, disease ID (NCBI stable ID)
# Loade disease ranges
# Files_order:
# OMIM-Entry-Retrieval_singlegene3.txt --- 1 2 3 0
# OMIM-Entry-Retrieval_monogen3.txt --- 1 2 3 0
# MendelianRanges.txt --- 0 3 1 2
# KEGG_mart_export.txt --- 4 1 2 6
# DISGENET.txt --- 1 2 3 4
print('Enter your file name:')
file_name = input()
print('\nIs coma separated values 0/1:')
coma = int(input())
print('\nWrite your files column indexes in this order: chrom, start, end, diseaseID. separated by one space each:')
order= input().split(' ')
order= [int(i) for i in order]

print(order)
with open( os.path.join('/home/jgonzale/Documentos/4_sem/humana/enfermedades',file_name),'r') as file:
    if coma:
        disease = [ [int(i.split(',')[order[0]]), int(i.split(',')[order[1]]),  int(i.split(',')[order[2]]), i.split(',')[order[3]].strip('\n')] for i in file if not i.startswith('#') if i.split(',')[order[0]] in [str(chrom) for chrom in human.keys()] ]
    else:
        disease = [ [int(i.split('\t')[order[0]]), int(i.split('\t')[order[1]]),  int(i.split('\t')[order[2]]), i.split('\t')[order[3]].strip('\n')] for i in file if not i.startswith('#')  if i.split('\t')[order[0]] in [str(chrom) for chrom in human.keys()] ]

    ### [ [chrom, start, end, diseaseID], ...]

print(disease[:50])
#####################################################################


#####################################################################
#Obtains the 3 types of slippage (see article Methods)
# filter ranges, keep insiders or slipages to right and left 
for dis in disease:
    for rang in human[dis[0]]:
        if (dis[1] >= rang[0] and dis[2] <= rang[1]) or (dis[1] <= rang[0] and (dis[2] <= rang[1] and dis[2] >= rang[0])) or ( (dis[1] >= rang[0] and dis[1] <= rang[1]) and dis[2] >= rang[1] ):
            human[dis[0]][rang].append(dis[1:]) # { chromosome:{(start,end,type):[ [start,end,diseaseID], ... ] } , ...}
#####################################################################


#####################################################################
#Obtains all biotypes present on human genome

# get all interval types
types = {}
for chromosome in human:
    for interval in human[chromosome]:
        types[interval[2]] = 0
print(types)
#####################################################################


#####################################################################
# count for interval type

for chromosome in human: # for every chr
    human[chromosome]['counts'] = dict(types) # make a new key 'counts' and value interval types
    for interval in human[chromosome]: # for every interval in the chr
        if interval != 'counts': # if it is not the recently created key 'counts'
            human[chromosome]['counts'][interval[2]] += len(human[chromosome][interval]) # sum the number of diseases that land in that interval grouped by their interval type
#####################################################################


#####################################################################
#Write final results files

with open( os.path.join('/home/jgonzale/Documentos/4_sem/humana/monogen_res',  'RES_' + file_name ),'w') as out: 
    out.write('#this file contains the counts of monogenic disease intervals that land in genomic intervals. This file is named after the origianl Disease data file\n')
    out.write('#chr\ttype\tcount\n')
    for chromosome in human:
        print('\n\nCHROMOSOME ', chromosome )
        for typ in human[chromosome]['counts']:
            out.write( str(chromosome) + '\t' + typ + '\t' + str(human[chromosome]['counts'][typ]) + '\n')
            print(typ, human[chromosome]['counts'][typ])
#####################################################################


#####################################################################
#Visualization of results

# see results for specific interval types
for chromosome in human:
    print('\n\nCHROMOSOME ', chromosome )
    print('CDS: ', human[chromosome]['counts']['CDS'])
#####################################################################