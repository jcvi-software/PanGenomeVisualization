#--------------------------------#
# Reference FASTA File           #
#         Namita Dongre          #
#           7/27/2015            #
#--------------------------------#

import csv
from collections import OrderedDict


def main():
#---Reading Files
#---Centroid: Fasta for Genes & Cluster
#---Core.attfGI: Contig to Cluster
    reader_fasta = csv.reader(open('centroids.csv', 'rt', encoding="ascii"), delimiter=',')
    reader_core = csv.reader(open('Core_attfGI.csv', 'rt', encoding="ascii"), delimiter=',')
    reader_combined = csv.reader(open('combined_att_file.csv', 'rt', encoding="ascii"), delimiter=',')

#---Creating List of Necessary Cluster/Contig Pairs
#---Ignoring CL_INS for now and CONTEXT
    coreDict = OrderedDict()
    for row in reader_core:
        (contig, cluster, start, end) = (row[0], row[1], row[2], row[3])
        if not 'CONTEXT' in cluster:
            if not 'CL_INS' in cluster:
                cluster = cluster[cluster.index('_')+1::]
                coreDict[cluster] = (contig, start, end)
            else: coreDict[cluster] = (contig, start, end)


#---Create an ordered dictionary (same order as coreattList) of clusters
    key = ''
    values = ''
    centroidFileDict = dict()
    geneClusterDict = dict()
    for row in reader_fasta:
        row = [element for element in row if element != '']
        if '>' in row[0]:
            centroidFileDict[key] = values
            values = ''
            lines = row[0].split(' ')[0]
            key = lines[lines.index('_')+1::]
            geneClusterDict[lines[1]] = key
        else:
            values+=row[0]

    contigDict = OrderedDict([('1',''), ('2',''), ('3',''), ('4',''), ('Missing', '')])

    for clusterfGI in coreDict:
        (contig, start, end) = (coreDict[clusterfGI][0], coreDict[clusterfGI][1], coreDict[clusterfGI][2])
        if clusterfGI in centroidFileDict:
            contigDict[contig] += centroidFileDict[clusterfGI]
        else:
            contigDict[contig] += 'X' * (abs(int(end)-int(start)) + 1)


    for row in reader_combined:
        (gene, start, end) = (row[1], row[2], row[3])
        if not gene in geneClusterDict:
            contigDict['Missing'] += 'X' * (abs(int(end)-int(start)) + 1)

    referenceList = []
    for contig in contigDict:
        referenceList.append('>'+contig+'\n')
        referenceList.append(contigDict[contig]+'\n')

    with open("referenceFASTA.csv", "w") as f:
        f.writelines(referenceList)

        
main()
