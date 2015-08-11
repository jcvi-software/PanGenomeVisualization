import csv

def main():
    reader_combined = csv.reader(open('combined_att_file.csv', 'rt', encoding="ascii"), delimiter=',')
    reader_panoct = csv.reader(open('panoct.csv', 'rt', encoding="ascii"), delimiter=',')
    reader_fgi = csv.reader(open('fGI_stats.csv', 'rt', encoding="ascii"), delimiter=',')
    reader_core = csv.reader(open('Core_attfGI.csv', 'rt', encoding="ascii"), delimiter=',')

    referenceList = []

    #-------------Iterate through the fGI_stats file
    #-----------------Create a dictionary {Unique fGI Inserts: (Start Cluster, End Cluster, Feature)}
    genome = 'E00002'
    genomeDict=dict() 
    for row in reader_fgi:
        (insert, numberOfGenomes, genomes) = (row[0], row[1], row[2])
        if genome in genomes:
            if ':' not in insert:
                startCluster = insert[insert.rindex('E')+1:-1]
                endCluster = 'X'
                feature = 'Break'
            elif 'STOP_CORE' in insert or 'U_CORE' in insert:
                startCluster = insert[insert.index('E')+1:insert.index(':')-1]
                endCluster = insert[insert.rindex('E')+1:-1]
                feature = insert[insert.rindex(':')+1:insert.rindex('E')+1]
            else:
                print(insert)
                startCluster = insert[insert.index('E')+1:insert.index(':')-1]
                endCluster = 'X'
                feature = 'Break'
            genomeDict[insert] = [startCluster, endCluster, feature]

    #-------------Iterate through the Core.attfGI file
    #-----------------Create a dictionary {Core Cluster: (start coordinate, end coordinate, contig, proteinName)}
    coreDict = dict()
    for row in reader_core:
        (cluster, start, end, contig, proteinName) = (row[1], row[2], row[3], row[0], row[4])
        if not 'INS' in cluster and not 'CONTEXT' in cluster:
            clusterID = cluster[cluster.index('_')+1::]
            coreDict[clusterID] = [start, end, contig, proteinName]
    
    #-------------Iterate through dictionary of fGI inserts and through clusters within each fGI insert
    #-----------------Create a list of lists for GFF files of fGI clusters
    #-----------------Clusters missing in the Core file will be placed in another contig for now since no coordinates available with respect to Core
    #-----------------Cluster inserts on fGIs (not START and End core clusters will be denoted as a subfeature)
    #-----------------Some fGIs exist across multiple contigs, these fGIs will appear extremely long on JBrowse
    count = 1
    for insert in genomeDict:
        (startCluster, endCluster, feature) = (genomeDict[insert][0], genomeDict[insert][1], genomeDict[insert][2])
        (start, contig, endContig) = (coreDict[startCluster][0], coreDict[startCluster][2], '')
        orientation = '+'
        print(feature)
        if feature != 'Break':
            end = coreDict[endCluster][1]
            endContig = coreDict[endCluster][2]
        else: end = coreDict[startCluster][1]
        if int(start) > int(end):
            (start, end, orientation) = (end, start, '-')
        annotation = 'ID='+feature+str(count)+';Name='+feature+';Note='+insert + ';Genome=' + genome
        if contig != endContig and endContig != '':
            annotation += ';Start_Contig=' + contig + ';End_Contig=' + endContig
        else: annotation += 'Contig=' + contig
        referenceList.append([contig, 'PanOCT', feature, start, end, '.', orientation, '.', annotation])

        (startCluster, endCluster, feature) = (genomeDict[insert][0], genomeDict[insert][1], genomeDict[insert][2])
        (start, contig) = (coreDict[startCluster][1], coreDict[startCluster][2])
        orientation = '+'
        if endCluster != 'X':
            end = coreDict[endCluster][0]
            if int(start) > int(end):
                (start, end, orientation) = (end, start, '-')
            referenceList.append([contig, 'PanOCT', 'fgi', start, end, '.', orientation, '.', 'Parent='+feature+str(count)+';Index=1;Note='+insert[insert.index(':')+1:insert.rindex(':')]])
        count+=1


#---------Dictionary of {geneID:(start, end)} for genome; from Combined.att_file
    genomeDict = dict()
    for row in reader_combined:
        (reference, geneID, start, stop, info, genomeID) = (row[0], row[1], row[2], row[3], row[4], row[5])
        if genome == genomeID:
            genomeDict[geneID] = [start, stop, reference, info]
    
#---------List of cluster, gene, gene, gene, etc.
    panoctList = list()
    for row in reader_panoct:
        row = [element for element in row if element != '']
        panoctList.append(row)

#---------Appending cluster to values of genomeDict...>{geneID:(start, end, cluster)}
    for gene in genomeDict:
        for row in panoctList:
            if gene in row:
                genomeDict[gene].append(row[0])
        if len(genomeDict[gene]) == 4:
            genomeDict[gene].append('None')

#---------Dictionary of {Cluster: (Contig, start, end)}
    clusterContigDict = dict()
    for row in reader_core:
        key = row[1][row[1].index('_')+1::]
        clusterContigDict[key] = (row[0], row[2], row[3])


#---------Creating Gene features        
    referenceList = list()
    for gene in genomeDict:
        #print(gene, len(genomeDict[gene]), genomeDict[gene])
        (startGenome, endGenome, reference, info, cluster) = (genomeDict[gene][0], genomeDict[gene][1], genomeDict[gene][2], genomeDict[gene][3], genomeDict[gene][4])
        orientation = orientation2 = '+'
        if int(startGenome) > int(endGenome):
            (startGenome, endGenome, orientation2) = (endGenome, startGenome, '-')

        annotation = 'Name= ' + info + '  |    Cluster ' + cluster + ';Gene_id= ' + gene + ';Cluster ID= CL_' + cluster + ';Genome Contig=' + reference +'; Genome Location=' + startGenome + '..' + endGenome + ' (' + orientation2 + ' strand)'

        if cluster in clusterContigDict:
            (contig, start, end) = (clusterContigDict[cluster][0], clusterContigDict[cluster][1], clusterContigDict[cluster][2])
            if int(start) > int(end):
                (start, end, orientation) = (end, start, '-')
            referenceList.append([contig, 'PanOCT', 'Gene', start, end, '.', orientation, '.', annotation])

        else: referenceList.append(['Other', 'PanOCT', 'Gene', startGenome, endGenome, '.', orientation, '.', annotation])


    with open( "e00002GFF1.csv", "w") as f:
        writer = csv.writer(f, delimiter='\t')
        writer.writerows(referenceList)
        

main()
