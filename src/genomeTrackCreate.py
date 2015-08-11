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
    coreDict, clusterContigDict = createCoreClusterDict(reader_core)
    # genomeDict = createfgiInsertDict(reader_fgi, genome)
    # geneDict = createGeneDict(panoctList, reader_combined, genome)
    # referenceList = createfGIFeatures(genomeDict, coreDict, genome)
    # referenceList.extend(createGeneFeatures(clusterContigDict, geneDict, genome))
    # writeFile(genome, referenceList)

    for genome in genomeListing:
       if '\n' in genome:
           genome = genome[0:-1]
       genomeDict = createfgiInsertDict(reader_fgi, genome)
       #geneDict = createGeneDict(panoctList, reader_combined, genome)
       referenceList = createfGIFeatures(genomeDict, coreDict, genome)
       #referenceList.extend(createGeneFeatures(clusterContigDict, geneDict, genome))
       #referenceList = createGeneFeatures(clusterContigDict, geneDict, genome)
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
    #-----------------Create a dictionary {Cluster: (Contig, start, end)}
def createCoreClusterDict(reader_core):
    coreDict = dict()
    clusterContigDict = dict()
    for row in reader_core:
        (cluster, start, end, contig, proteinName) = (row[1], int(int(row[2])/3), int(int(row[3])/3), row[0], row[4])
        if not 'INS' in cluster and not 'CONTEXT' in cluster:
            clusterID = cluster[cluster.index('_')+1::]
            coreDict[clusterID] = [start, end, contig, proteinName]
            clusterContigDict[clusterID] = (contig, start, end)
    return coreDict, clusterContigDict



    #-------------Iterate through dictionary of fGI inserts and through clusters within each fGI insert
    #-----------------Create a list of lists for GFF files of fGI clusters
    #-----------------Clusters missing in the Core file will be placed in another contig for now since no coordinates available with respect to Core
    #-----------------Cluster inserts on fGIs (not START and End core clusters will be denoted as a subfeature)
    #-----------------Some fGIs exist across multiple contigs, these fGIs will appear extremely long on JBrowse
def createfGIFeatures(genomeDict, coreDict, genome):
    referenceList = []
    count = 1
    oldEnd = 0
    for insert in genomeDict:
        (startCluster, endCluster, feature) = (genomeDict[insert][0], genomeDict[insert][1], genomeDict[insert][2])
        (start, startInsert, contig, endContig) = (min(coreDict[startCluster][0], coreDict[startCluster][1]), max(coreDict[startCluster][0], coreDict[startCluster][1]), coreDict[startCluster][2], '')
        orientation = '+'
        if feature != 'Break':
            (endInsert, end, endContig) = (min(coreDict[endCluster][0], coreDict[endCluster][1]), max(coreDict[endCluster][0], coreDict[endCluster][1]), coreDict[endCluster][2])
        else:
            end = startInsert
        if start > end:
            if feature == 'U_CORE':
                (startInsert, endInsert, start, end, orientation) = (end, start, endInsert, startInsert, '-')
            else: (start, end, startInsert, endInsert, orientation) = (end, start, endInsert, startInsert, '-')

        if start == oldEnd:
            start+=1
        annotation = 'ID='+feature+str(count)+';Name='+feature+';Note='+insert + ';Genome=' + genome
        if contig != endContig and endContig != '':
            annotation += ';Start_Contig=' + contig + ';End_Contig=' + endContig
            print(True)
        else: annotation += ';Contig=' + contig
        referenceList.append([contig, 'PanOCT', feature, start, end, '.', orientation, '.', annotation])

        if feature != 'Break':
            insert2 = insert[insert.index(':')+1:insert.rindex(':')]
            if len(insert2) != 0:
                if feature == 'U_CORE':
                    referenceList.append([contig, 'PanOCT', 'fgi_u', startInsert, endInsert, '.', orientation, '.', 'Parent='+feature+str(count)+';Index=1;Note='+insert2])
                else:
                    referenceList.append([contig, 'PanOCT', 'fgi_s', startInsert, endInsert, '.', orientation, '.', 'Parent='+feature+str(count)+';Index=1;Note='+insert2])
            else:
                referenceList.append([contig, 'PanOCT', 'fgi_s', startInsert, endInsert, '.', orientation, '.', 'Parent='+feature+str(count)+';Index=1'])
        else: referenceList.append([contig, 'PanOCT', 'breakSub', start, end, '.', orientation, '.', 'Parent='+feature+str(count)+';Index=1'])
        oldEnd = end
        count+=1

    return referenceList

#---------List of cluster, gene, gene, gene, etc.
def createClusterGenomeList(reader_panoct):
    panoctList = list()
    for row in reader_panoct:
        row = [element for element in row if element != '']
        panoctList.append(row)
    return panoctList


#---------Dictionary of {geneID:(start, end, cluster)} for genome; from Combined.att_file
def createGeneDict(panoctList, reader_combined, genome):
    geneDict = dict()
    count = 0
    for row in reader_combined:
        (reference, geneID, start, stop, info, genomeID) = (row[0], row[1], int(int(row[2])/3), int(int(row[3])/3), row[4], row[5])
        if genomeID == genome:
            count+=1
            geneDict[geneID] = [start, stop, reference, info]
        else: continue
    for gene in geneDict:
        for row in panoctList:
            if gene in row:
                geneDict[gene].append(row[0])
        if len(geneDict[gene]) == 4:
            geneDict[gene].append('None')
    return geneDict

#---------Creating Gene features
def createGeneFeatures(clusterContigDict, geneDict, genome):
    referenceList = []
    oldEnd = 0
    for gene in geneDict:
        (startGenome, endGenome, reference, info, cluster) = (geneDict[gene][0], geneDict[gene][1], geneDict[gene][2], geneDict[gene][3], geneDict[gene][4])
        orientation = orientation2 = '+'
        if int(startGenome) > int(endGenome):
            (startGenome, endGenome, orientation2) = (endGenome, startGenome, '-')
        if startGenome == oldEnd:
            startGenome+=1
        annotation = 'Name= ' + info + '  |    Cluster ' + cluster + ';Gene_ID= ' + gene + ';Cluster_ID= CL_' + cluster + ';Genome_Contig=' + reference +'; Genome_Location=' + str(startGenome) + '..' + str(endGenome) + ' (' + orientation2 + ' strand)'
        if cluster in clusterContigDict:
            (contig, start, end) = (clusterContigDict[cluster][0], clusterContigDict[cluster][1], clusterContigDict[cluster][2])
            if int(start) > int(end):
                (start, end, orientation) = (end, start, '-')
            referenceList.append([contig, 'PanOCT', 'Gene', start, end, '.', orientation, '.', annotation])
        else:
            referenceList.append(['Missing', 'PanOCT', 'Gene', startGenome, endGenome, '.', orientation, '.', annotation])
        oldEnd = endGenome
    return referenceList

def writeFile(genome, referenceList):
    with open(genome+'.gff', 'w') as f:
        writer = csv.writer(f, delimiter='\t', quoting=csv.QUOTE_NONE)
        writer.writerows(referenceList)


main()
