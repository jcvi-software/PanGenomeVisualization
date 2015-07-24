import csv
from collections import OrderedDict


def main():
#---Reading Files
#---Centroid: Fasta for Genes & Cluster
#---Core.attfGI: Contig to Cluster
    reader_fasta = csv.reader(open('centroids.csv', 'rt', encoding="ascii"), delimiter=',')
    reader_coreatt = csv.reader(open('Core_attfGI.csv', 'rt', encoding="ascii"), delimiter=',')

#---Creating List of Necessary Cluster/Contig Pairs
#---Ignoring CL_INS for now and CONTEXT
    coreattList = []
    cluster2contig = dict()
    for row in reader_coreatt:
        (contigID, clusterID) = (int(row[0]), row[1])
        if not(':' in clusterID) and clusterID.count('_') == 1: 
            clusterID = clusterID[clusterID.index('_')+1::]
            coreattList.append([clusterID, contigID])
            cluster2contig[clusterID]=contigID

##    for item in coreattList:
##        print(item)

#---Create an ordered dictionary (same order as coreattList) of clusters
    centroidDict = OrderedDict()

    for item in coreattList:
        centroidDict[item[0]] = 0

#---Retrieve data from centroid        
    fastaList = []
    for row in reader_fasta:
        row = [element for element in row if element != ''][0]
        row = row.split(' ')
        if row[0][0] == '>':
            key = row[0][row[0].index('_')+1::]
            fastaList.append('>'  + key)
        else: fastaList.append(row[0])

#---Match each cluster to sequence 
    key = fastaList[0][1::]
    values = ''
    for index in range(1,len(fastaList)):
        if '>' in fastaList[index]:
            centroidDict[key] = values
            key = fastaList[index][1::]
            values = ''
        else:
            values+=fastaList[index]


##    for key in centroidDict:
##        print(key +': '+  centroidDict[key])

#---Dictionary of concateated sequences for each contig
    contigDict = OrderedDict([(1, ''), (2, ''), (3, ''), (4, ''), (5, '')])
    for clusterID in centroidDict:
        if clusterID in cluster2contig:
            contig = cluster2contig[clusterID]
            contigDict[contig] += centroidDict[clusterID]


#---Convert dictionary of contig:sequence into list to save into file
    referenceList = []
    for key in contigDict:
        referenceList.append('>' + str(key))
        referenceList.append(contigDict[key])

    with open("referenceFASTA.csv", "w") as f:
        f.writelines(referenceList)

    

        
        
        
        
main()
