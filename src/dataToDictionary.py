#------------------------------------------------#
#   Namita Dongre                                #
#   8/13/2015                                    #
#   Dictionaries for Pangenome Visualization     #
#------------------------------------------------#

import csv, pickle, os
from collections import OrderedDict

def main():
    path, endPath = getFilePath()
    reader_combined = list(csv.reader(open(path+'combined_att_file.csv', 'rt', encoding="ascii"), delimiter=','))
    reader_all = list(csv.reader(open(path+'all_clusters_members.csv', 'rt', encoding="ascii"), delimiter=','))
    reader_fasta = csv.reader(open(path+'centroids.csv', 'rt', encoding="ascii"), delimiter=',')
    reader_core = csv.reader(open(path+'Core_attfGI.csv', 'rt', encoding="ascii"), delimiter=',')

    genomeClusterDict = createGenomeClusterDict(reader_all)
    genomeLocusDict = createGenomeLocusDict(reader_combined)
    coreDict = coreDictCreate(reader_core)
    centroidFileDict, geneClusterDict = createClusterDicts(reader_fasta)

    pickle.dump(genomeClusterDict, open(endPath+'genomeCluster.dict', 'wb'))
    pickle.dump(genomeLocusDict, open(endPath+'genomeLocus.dict', 'wb'))
    pickle.dump(centroidFileDict, open(endPath+'centroidSequence.dict', 'wb'))
    pickle.dump(geneClusterDict, open(endPath+'geneToCluster.dict', 'wb'))
    pickle.dump(coreDict, open(endPath+'coreCluster.dict', 'wb'))

#------------------------------------------------------------------------------------------------------------------Methods
def getFilePath():
    path = input('Enter file path for input data:   ')
    endPath = input('Enter file path for output data. We suggest saving in the same location as input data:   ')
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

def createGenomeClusterDict(reader_all):
    genomeClusterDict = dict()
    for row in reader_all:
        (clusterID, name, centroidGenome, centroidLocus, locus1, locus2, locus3, locus4, locus5, locus31, locus32, locus33, locus34, locus35) = (row[0], row[1], row[2], row[3], row[6], row[7], row[8], row[9], row[10], row[11], row[12], row[13], row[14], row[15])
        genomeClusterDict[clusterID] = [name, centroidGenome, centroidLocus, locus1, locus2, locus3, locus4, locus5, locus31, locus32, locus33, locus34, locus35]
    return genomeClusterDict

def createGenomeLocusDict(reader_combined):
    genomeLocusDict = dict()
    for row in reader_combined:
        (genomeContig, geneLocus, start, end, name, genome) = (row[0], row[1], int(row[2]), int(row[3]), row[4], row[5])
        if start > end:
            (start, end) = (end, start)
        genomeLocusDict[geneLocus] = [genome, genomeContig, start, end, int(round((end-start)/3, -1)), name]
    return genomeLocusDict


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

    return coreDict


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
    return centroidFileDict, geneClusterDict



main()
