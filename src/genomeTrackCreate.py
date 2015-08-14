#-------------------------#
#   Namita Dongre         #
#   8/13/2015             #
#   Genome GFF Tracks     #
#-------------------------#

import csv, pickle, os

def main():
    # path = 'ExampleData/input_files/'
    # endPath = 'ExampleData/test/'
    path, endPath = getFilePath()
    reader_fgi = list(csv.reader(open(path+'fGI_stats.csv', 'rt', encoding="ascii"), delimiter=','))
    reader_core = list(csv.reader(open(path+'Core_attfGI.csv', 'rt', encoding="ascii"), delimiter=','))
    genomeListing = list(open(path+'db.txt', 'r'))

    genomeClusterDict = pickle.load(open(path+'genomeCluster.dict', 'rb'))
    genomeLocusDict = pickle.load(open(path+'genomeLocus.dict', 'rb'))
    coreDict, fgiDict = createCoreClusterDict(reader_core)

    #genome = 'E00002'
    genomeIdDict = {}
    index= 3
    for genome2 in genomeListing:
        if '\n' in genome2:
            genome2 = genome2[0:-1]
        genomeIdDict[genome2] = index
        index+=1

    for genome in genomeIdDict:
       genomeDict = createfgiInsertDict(reader_fgi, genome)
       referenceList = createfGIFeatures(genomeDict, coreDict, fgiDict, genomeClusterDict, genomeLocusDict, genome, genomeIdDict[genome])
       writeFile(endPath, genome, referenceList)

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
    clusterList = []
    count = 1
    for insert in genomeDict:
        (startCluster, endCluster, feature) = (genomeDict[insert][0], genomeDict[insert][1], genomeDict[insert][2])
        (start, startInsert, startContig, proteinName) = (coreDict[startCluster][0], coreDict[startCluster][1], coreDict[startCluster][2], coreDict[startCluster][3])
        (end, endInsert, endContig, proteinName) = ('', '', '', '')
        if feature != 'Break':
            (endInsert, end, endContig, proteinName) = (coreDict[endCluster][0], coreDict[endCluster][1], coreDict[endCluster][2], coreDict[endCluster][3])
        else: end = startInsert

        orientation = '+'
        if feature != 'Break':
            temp = start
            coordList= [start, end, startInsert, endInsert]
            coordList.sort()
            (start, startInsert, endInsert, end) = (coordList[0], coordList[1], coordList[2], coordList[3])
            if temp == end:
                orientation = '-'

        annotation = 'ID=' + feature+ str(count)+ ';Note=' + insert +';Contig=' + startContig +';Genome=' + genome +';Name=' + feature + ' Insertion'

        insertClusterList, startEnd = clusterInsertSplit(insert, feature)
        geneStr = ''
        bpLen = len(insertClusterList)

        startCoord = startInsert
        if bpLen != 0:
            geneStr+= ';Genes='

        totalLen = 0
        increment = 0
        for cluster in insertClusterList:
            cluster = cluster[0:-1]
            locus = genomeClusterDict[cluster][idCode]
            totalLen +=  genomeLocusDict[locus][4]
            geneStr += cluster+', '

        if bpLen != 0:
            increment = int(round(((totalLen/startCoord)/bpLen),-1))

        if feature=='Break':
            referenceList.append([startContig, 'PanOCT', feature, start, startInsert+totalLen, '.', orientation, '.', annotation+geneStr])
            referenceList.append([startContig, 'PanOCT', feature+'_Cluster', start, startInsert, '.', orientation, '.','Name=CL_'+startEnd[0] + ';Parent=' + feature+ str(count) + ';Note=Core Cluster; Contig=' + startContig + ';Genome=' + genome + ';Index=1'])
        else:
            referenceList.append([startContig, 'PanOCT', feature, start, end, '.', orientation, '.', annotation+geneStr])
            referenceList.append([startContig, 'PanOCT', feature+'_Cluster', start, startInsert, '.', orientation, '.', 'Name=CL_'+startEnd[0] + ';Parent=' + feature+ str(count) + ';Note=Core Cluster; Contig=' + startContig + ';Genome=' + genome + ';Index=1'])
            referenceList.append([startContig, 'PanOCT', feature+'_Cluster', endInsert, end, '.', orientation, '.', 'Name=CL_'+startEnd[1] + ';Parent=' + feature+ str(count) + ';Note=Core Cluster; Contig=' + startContig + ';Genome=' + genome +';Index=1'])
        count+=1

#-----------------Code for subfeatures
        # if bpLen != 0:
        #     referenceList.append([startContig, 'PanOCT', feature+'_Cluster', start, startInsert, '.', orientation, '.', 'Name=CL_'+startEnd[0] + ';Parent=' + feature+ str(count) + ';Note=Core Cluster; Contig=' + startContig + ';Genome=' + genome + ';Index=1'])
        #     if feature=='Break' and bpLen != 0:
        #         referenceList.append([startContig, 'PanOCT', feature, start, startInsert + totalLen, '.', orientation, '.', annotation+geneStr])
        #     elif feature == 'U_CORE' or feature == 'STOP_CORE' and bpLen != 0:
        #         referenceList.append([startContig, 'PanOCT', feature, startInsert, endInsert, '.', orientation, '.', annotation+geneStr])
        #         referenceList.append([startContig, 'PanOCT', feature+'_Cluster', endInsert, end, '.', orientation, '.', 'Name=CL_'+startEnd[1] + ';Parent=' + feature+ str(count) + ';Note=Core Cluster; Contig=' + startContig + ';Genome=' + genome +';Index=1'])

            # for cluster in insertClusterList:
            #     cluster = cluster[0:-1]
            #     (name, centroidGenome, centroidLocus, locus) = (genomeClusterDict[cluster][0], genomeClusterDict[cluster][1], genomeClusterDict[cluster][2], genomeClusterDict[cluster][idCode])
            #     (genome, genomeContig, start, end, length)  = (genomeLocusDict[locus][0], genomeLocusDict[locus][1], genomeLocusDict[locus][2], genomeLocusDict[locus][3], genomeLocusDict[locus][4])
            #     clusterAnnote = 'Name=' + name + ';Parent=' + feature+ str(count) +';Cluster=' + str(cluster) + ';Locus_ID=' + locus + ';Contig=' + startContig +';Genome=' + genome +';Length=' + str(length)
            #     clusterList.append([startContig, 'PanOCT', 'NonCoreCluster', startCoord, startCoord+length, '.', orientation, '.', clusterAnnote])
            #     startCoord += (length + increment)

            #referenceList.extend(clusterList)

        # else:
        #     if feature == 'Break':
        #         referenceList.append([startContig, 'PanOCT', feature, start, startInsert, '.', orientation, '.', 'Name=Break_Cluster;Contig=' + startContig + ';Genome=' + genome])
        #     else:
        #         referenceList.append([startContig, 'PanOCT', feature+'_Cluster', start, startInsert, '.', orientation, '.', 'Name=CL_'+startEnd[0] + ';Parent=' + feature+ str(count) + ';Note=Core Cluster; Contig=' + startContig + ';Genome=' + genome + ';Index=1'])
        #         referenceList.append([startContig, 'PanOCT', feature, start, end, '.', orientation, '.', annotation+geneStr])
        #         referenceList.append([startContig, 'PanOCT', feature+'_Cluster', endInsert, end, '.', orientation, '.', 'Name=CL_'+startEnd[1] + ';Parent=' + feature+ str(count) + ';Note=Core Cluster; Contig=' + startContig + ';Genome=' + genome +';Index=1'])


    return referenceList

def clusterInsertSplit(insert, feature):
    start = insert[insert.index('E')+1:-1]
    end = ''
    if feature != 'Break':
        start = insert[insert.index('E')+1:insert.index(':')]
        end = insert[insert.rindex('E')+1::]
        startEnd = [start[0:-1], end[0:-1]]
    elif feature == 'Break' and ':' in insert:
        start = insert[insert.index('E')+1:insert.index(':')]
        startEnd = [start, end]
    else: startEnd = [start, end]
    return insert.split(':')[1:-1], startEnd


def writeFile(outPath, genome, referenceList):
    with open(outPath + genome+'.gff', 'w') as f:
        writer = csv.writer(f, delimiter='\t', quoting=csv.QUOTE_NONE)
        writer.writerows(referenceList)
    print(genome+'.gff File Created.')


main()
