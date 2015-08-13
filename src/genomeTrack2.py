#-------------------------#
#   Namita Dongre         #
#   8/7/2015              #
#   Genome GFF Tracks     #
#-------------------------#

import csv, pickle, os

def main():
    path, endPath = getFilePath()
    reader_fgi = list(csv.reader(open(path+'fGI_stats.csv', 'rt', encoding="ascii"), delimiter=','))
    reader_core = list(csv.reader(open(path+'Core_attfGI.csv', 'rt', encoding="ascii"), delimiter=','))
    genomeListing = list(open(path+'db.txt', 'r'))

    genomeClusterDict = pickle.load(open(path+'genomeCluster.dict', 'rb'))
    genomeLocusDict = pickle.load(open(path+'genomeLocus.dict', 'rb'))
    coreDict, fgiDict = createCoreClusterDict(reader_core)

    genomeIdDict = {}
    index= 3
    for genome in genomeListing:
        if '\n' in genome:
            genome = genome[0:-1]
        genomeIdDict[genome] = index
        index+=1

    for genome in genomeIdDict:
       genomeDict = createfgiInsertDict(reader_fgi, genome)
       referenceList = createfGIFeatures(genomeDict, coreDict, fgiDict, genomeClusterDict, genomeLocusDict, genome, genomeIdDict[genome])
       writeFile(endPath, genome, referenceList)

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
    return coreDict, fgiDict

def createfGIFeatures(genomeDict, coreDict, fgiDict, genomeClusterDict, genomeLocusDict, genome,idCode):
    referenceList = []
    count = 1
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

        insertClusterList = clusterInsertSplit(insert)
        geneStr = ''
        bpLen = 0
        if len(insertClusterList) != 0:
            geneStr+= ';Genes='
        elif feature == 'Break':
            bpLen = end - start
        else: bpLen = endInsert-startInsert
        for cluster in insertClusterList:
            cluster = cluster[0:-1]
            (name, centroidGenome, centroidLocus, locus) = (genomeClusterDict[cluster][0], genomeClusterDict[cluster][1], genomeClusterDict[cluster][2], genomeClusterDict[cluster][idCode])
            (genome, genomeContig, start, end, length)  = (genomeLocusDict[locus][0], genomeLocusDict[locus][1], genomeLocusDict[locus][2], genomeLocusDict[locus][3], genomeLocusDict[locus][4])
            geneStr += cluster+'-'+str(length)+' bp, '
            bpLen += length

        if feature == 'STOP_CORE' or feature == 'U_CORE' or feature=='Break':
            referenceList.append([startContig, 'PanOCT', feature, startInsert, startInsert+bpLen, '.', orientation, '.', annotation+geneStr])
    return referenceList

def clusterInsertSplit(insert):
    return insert.split(':')[1:-1]


def writeFile(outPath, genome, referenceList):
    with open(outPath + genome+'.gff', 'w') as f:
        writer = csv.writer(f, delimiter='\t', quoting=csv.QUOTE_NONE)
        writer.writerows(referenceList)
    print(genome+'.gff File Created.')


main()
