#-------------------------#
#   Namita Dongre         #
#   8/7/2015              #
#   Genome GFF Tracks     #
#-------------------------#

import csv, pickle

def main():
    reader_combined = list(csv.reader(open('combined_att_file.csv', 'rt', encoding="ascii"), delimiter=','))
    reader_panoct = list(csv.reader(open('panoct.csv', 'rt', encoding="ascii"), delimiter=','))
    reader_fgi = list(csv.reader(open('fGI_stats.csv', 'rt', encoding="ascii"), delimiter=','))
    reader_core = list(csv.reader(open('Core_attfGI.csv', 'rt', encoding="ascii"), delimiter=','))
    genomeListing = open('db.txt', 'r')
    genomeListing = list(genomeListing)

    genome = 'E00001'
    panoctList = createClusterGenomeList(reader_panoct)
    coreDict, fgiDict, maxLenInsertion = createCoreClusterDict(reader_core)

    for genome in genomeListing:
       if '\n' in genome:
           genome = genome[0:-1]
       genomeDict = createfgiInsertDict(reader_fgi, genome)
       referenceList = createfGIFeatures(genomeDict, coreDict, fgiDict, genome, maxLenInsertion)
       writeFile(genome, referenceList)
       print(genome + ' File Created')


    #-------------Iterate through the fGI_stats file
    #-----------------Create a dictionary {Unique fGI Inserts: (Start Cluster, End Cluster, Feature)}
def createfgiInsertDict(reader_fgi, genome):
    genomeDict=dict()
    for row in reader_fgi:
        (insert, numberOfGenomes, genomes) = (row[0], row[1], row[2])
        if genome in genomes:
            if ':' not in insert:
                startCluster = insert[insert.rindex('E')+1:-1]
                endCluster = 'X'
                feature = 'Break'
            elif ('STOP_CORE' in insert) or ('U_CORE' in insert):
                startCluster = insert[insert.index('E')+1:insert.index(':')-1]
                endCluster = insert[insert.rindex('E')+1:-1]
                feature = insert[insert.rindex(':')+1:insert.rindex('E')+1]
            else:
                startCluster = insert[insert.index('E')+1:insert.index(':')-1]
                endCluster = 'X'
                feature = 'Break'
            genomeDict[insert] = [startCluster, endCluster, feature]
    return genomeDict

    #-------------Iterate through the Core.attfGI file
    #-----------------Create a dictionary {Core Cluster: (start coordinate, end coordinate, contig, proteinName)}
    #-----------------Create a dictionary {fGI: (contig, start coordinate, end coordinate)}

def createCoreClusterDict(reader_core):
    coreDict = dict()
    fgiDict = dict()
    maxLenInsertion = []
    for row in reader_core:
        (cluster, start, end, contig, proteinName) = (row[1], int(int(row[2])/3), int(int(row[3])/3), row[0], row[4])
        if start > end:
            (start, end) = (end, start)
        if not 'INS' in cluster and not 'CONTEXT' in cluster:
            clusterID = cluster[cluster.index('_')+1::]
            coreDict[clusterID] = [start, end, contig, proteinName]
        elif not 'CONTEXT' in cluster:
            coordString = str(start) + str(end)
            fgiDict[coordString] = (contig, cluster)
            maxLenInsertion.append(abs(end-start))
    return coreDict, fgiDict, max(maxLenInsertion)



    #-------------Iterate through dictionary of fGI inserts and through clusters within each fGI insert
    #-----------------Create a list of lists for GFF files of fGI clusters
    #-----------------Clusters missing in the Core file will be placed in another contig for now since no coordinates available with respect to Core
    #-----------------Cluster inserts on fGIs (not START and End core clusters will be denoted as a subfeature)
    #-----------------Some fGIs exist across multiple contigs, these fGIs will appear extremely long on JBrowse
def createfGIFeatures(genomeDict, coreDict, fgiDict, genome, maxLenInsertion):
    referenceList = []
    count = 1
    oldEnd = 0
    for insert in genomeDict:
        (startCluster, endCluster, feature) = (genomeDict[insert][0], genomeDict[insert][1], genomeDict[insert][2])
        (start, startInsert, startContig, proteinName) = (coreDict[startCluster][0], coreDict[startCluster][1], coreDict[startCluster][2], coreDict[startCluster][3])
        (end, endInsert, endContig, proteinName) = ('', '', '', '')
        if feature != 'Break':
            (endInsert, end, endContig, proteinName) = (coreDict[endCluster][0], coreDict[endCluster][1], coreDict[endCluster][2], coreDict[endCluster][3])
        else: end = startInsert

        orientation = '+'
        if start > end:
            (start, end, startInsert, endInsert, orientation) = (end, start, endInsert, startInsert, '-')

        annotation = 'Note=' + insert +';Contig' + startContig +';Genome=' + genome +';Name=' + feature + ' fGI Insertion'
        coordString = str(startInsert) + str(endInsert)
        if feature == 'STOP_CORE' or feature == 'U_CORE':
            print(feature)
            # if coordString in fgiDict:
            #     referenceList.append([startContig, 'PanOCT', feature, startInsert, endInsert, '.', orientation, '.', annotation])
            # if (endInsert - startInsert) <= maxLenInsertion:
            #     referenceList.append([startContig, 'PanOCT', feature, startInsert, endInsert, '.', orientation, '.', annotation])
    return referenceList

#---------List of cluster, gene, gene, gene, etc.
def createClusterGenomeList(reader_panoct):
    panoctList = list()
    for row in reader_panoct:
        row = [element for element in row if element != '']
        panoctList.append(row)
    return panoctList


def writeFile(genome, referenceList):
    with open(genome+'.gff', 'w') as f:
        writer = csv.writer(f, delimiter='\t', quoting=csv.QUOTE_NONE)
        writer.writerows(referenceList)


main()
