#--------------------------------#
# Reference FASTA File           #
#         Namita Dongre          #
#           7/27/2015            #
#--------------------------------#

import csv, pickle, os
from collections import OrderedDict


def main():
    path, endPath = getFilePath()
    reader_core = csv.reader(open(path+'Core_attfGI.csv', 'rt', encoding="ascii"), delimiter=',')
    reader_combined = list(csv.reader(open(path+'combined_att_file.csv', 'rt', encoding="ascii"), delimiter=','))

    coreDict = pickle.load(open(path+'coreCluster.dict', 'rb'))
    centroidFileDict = pickle.load(open(path+'centroidSequence.dict', 'rb'))
    geneClusterDict = pickle.load(open(path+'geneToCluster.dict', 'rb'))

    contigDict = OrderedDict([('1',''), ('2',''), ('3',''), ('4',''), ('Missing', '')])

    (oldStart, oldEnd) = (0, 0)
    contigStr = ''
    countPush = 0
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

    with open(endPath+"reference.fasta", "w") as f:
        f.writelines(referenceList)


#------------------------------------------------------------------------------------------------------------------Methods
def getFilePath():
    path = input('Enter file path for input data:   ')
    endPath = input('Enter file path for output data:   ')
    if not(os.path.isdir(path)):
        print('Sorry the input path you specified does not exist. Make sure you end the path name with a backslash.')
        exit()
    if not(os.path.isdir(endPath)):
        print('Sorry the output path you specified does not exist. Make sure you end the path name with a backslash.')
        exit()
    if not(os.path.isfile(path+'fGI_stats.csv') and os.path.isfile(path+'Core_attfGI.csv')):
        print('Sorry the path you specified does not contain the required files fGI_stats.csv and Core_attfGI.csv')
        exit()
    return path, endPath


def modResult(start, end):
    modResult = (abs(start - end))%3
    if modResult == 0:
        return (int(abs(start - end)/3))
    if modResult == 1:
        return (int(abs(start - end)/3) - 1)
    return (int(abs(start - end)/3) + 1)

def createCoreClusterDict(reader_core):
    coreDict = dict()
    for row in reader_core:
        (cluster, start, end, contig, proteinName) = (row[1], int(int(row[2])/3), int(int(row[3])/3), row[0], row[4])
        if start > end:
            (start, end) = (end, start)
        if not 'INS' in cluster and not 'CONTEXT' in cluster:
            clusterID = cluster[cluster.index('_')+1::]
            coreDict[clusterID] = [contig, start, end]
    return coreDict

main()
