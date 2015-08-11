#--------------------------------#
# Reference FASTA File           #
#         Namita Dongre          #
#           7/27/2015            #
#--------------------------------#

import csv, pickle
from collections import OrderedDict


def main():
#---Reading Files
#---Centroid: Fasta for Genes & Cluster
#---Core.attfGI: Contig to Cluster
    reader_fasta = csv.reader(open('centroids.csv', 'rt', encoding="ascii"), delimiter=',')
    reader_core = csv.reader(open('Core_attfGI.csv', 'rt', encoding="ascii"), delimiter=',')
    reader_combined = csv.reader(open('combined_att_file.csv', 'rt', encoding="ascii"), delimiter=',')

    # coreDictCreate(reader_core)
    # createClusterDicts(reader_fasta)

    coreDict = pickle.load(open('consensusCore.dict', 'rb'))
    centroidFileDict = pickle.load(open('centroidSequence.dict', 'rb'))
    geneClusterDict = pickle.load(open('geneToCluster.dict', 'rb'))

    contigDict = OrderedDict([('1',''), ('2',''), ('3',''), ('4',''), ('Missing', '')])

    (oldStart, oldEnd) = (0, 0)
    count = 0
    contigList = []
    contigList.append('>1' + '\n')
    contigStr = ''
    countPush = 0
    countinc = 1
    for clusterfGI in coreDict:
        (contig, start, end) = (coreDict[clusterfGI][0], int(coreDict[clusterfGI][1]), int(coreDict[clusterfGI][2]))
        if start > end:
            (start, end) = (end, start)
        if start == oldEnd:
            start+=1
        if clusterfGI in centroidFileDict:
            contigDict[contig] += centroidFileDict[clusterfGI]
        else:
            contigDict[contig] += 'X'*int(end-start+countPush)
            countPush+=1
        (oldStart, oldEnd) = (start, end)


    coordinateList = (9999999999999999999, -999999999999999999999999)
    for row in reader_combined:
        (gene, start, end) = (row[1], int(row[2]), int(row[2]))
        if start > end:
            (start, end) = (end, start)
        coordinateList = [min(start, coordinateList[0]), max(end, coordinateList[1])]

    contigDict['Missing'] += 'X' * modResult(coordinateList[0], coordinateList[1])

    referenceList = []
    for contig in contigDict:
        referenceList.append('>'+contig+'\n')
        referenceList.append(contigDict[contig]+'\n')

    with open("reference.fasta", "w") as f:
        f.writelines(referenceList)





################################################################################
def modResult(start, end):
    modResult = (abs(start - end))%3
    if modResult == 0:
        return (int(abs(start - end)/3))
    if modResult == 1:
        return (int(abs(start - end)/3) - 1)
    return (int(abs(start - end)/3) + 1)


#---Creating List of Necessary Cluster/Contig Pairs
#---Ignoring CL_INS for now and CONTEXT
def coreDictCreate(reader_core):
    (oldStart, oldEnd) = (0, 0)
    coreDict = OrderedDict()
    for row in reader_core:
        (contig, cluster, start, end, name, feature, aminoAcids) = (row[0], row[1], row[2], row[3], row[4], row[5], row[6])
        orientation = '+'
        (start, end) = (int(int(start)/3),int(int(end)/3))
        if start > end:
            (start, end, orientation) = (end, start, '-')
        if start == oldEnd:
            start+=1
        end = start+int(aminoAcids)
        if not 'CONTEXT' in cluster:
            if not 'CL_INS' in cluster:
                cluster = cluster[cluster.index('_')+1::]
                coreDict[cluster] = (contig, start, end)
            else: coreDict[cluster] = (contig, start, end)

    pickle.dump(coreDict, open('consensusCore.dict', 'wb'))


#---Create an ordered dictionary (same order as coreattList) of clusters
def createClusterDicts(reader_fasta):
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
    pickle.dump(centroidFileDict, open('centroidSequence.dict', 'wb'))
    pickle.dump(geneClusterDict, open('geneToCluster.dict', 'wb'))


main()
