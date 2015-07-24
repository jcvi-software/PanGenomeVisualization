import csv

def main():
    reader_combined = csv.reader(open('combined_att_file.csv', 'rt', encoding="ascii"), delimiter=',')
    reader_panoct = csv.reader(open('panoct.csv', 'rt', encoding="ascii"), delimiter=',')
    reader_core = csv.reader(open('Core_attfGI.csv', 'rt', encoding="ascii"), delimiter=',')
    reader_genome = csv.reader(open('e00001GFF.csv', 'rt', encoding="ascii"), delimiter='\t')

#---------Dictionary of {geneID:(start, end)} for genome; from Combined.att_file
    genomeID = 'E00001'
    genomeDict = dict()
    for row in reader_combined:
        (reference, geneID, start, stop, info, genome) = (row[0], row[1], row[2], row[3], row[4], row[5])
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
    gffList = list()
    for gene in genomeDict:
        #print(gene, len(genomeDict[gene]), genomeDict[gene])
        (startGenome, endGenome, reference, info, cluster) = (genomeDict[gene][0], genomeDict[gene][1], genomeDict[gene][2], genomeDict[gene][3], genomeDict[gene][4])
        orientation = orientation2 = '+'
        if int(startGenome) > int(endGenome):
            (startGenome, endGenome, orientation2) = (endGenome, startGenome, '-')

        annotation = 'Name= ' + info + '  |    Cluster ' + cluster + '; Gene_id= ' + gene + '; Description= CL_' + cluster + '; Genome_contig=' + reference +'; Location_Genome=' + startGenome + '..' + endGenome + ' (' + orientation2 + ' strand)'

        if cluster in clusterContigDict:
            (contig, start, end) = (clusterContigDict[cluster][0], clusterContigDict[cluster][1], clusterContigDict[cluster][2])
            if int(start) > int(end):
                (start, end, orientation) = (end, start, '-')
            gffList.append([contig, 'PanOCT', 'Gene', start, end, '.', orientation, '.', annotation])

        else: gffList.append(['Other', 'PanOCT', 'Gene', startGenome, endGenome, '.', orientation, '.', annotation])

#---------
    for row in reader_genome:
        if len(row) != 0 and row[2] != 'CL':
            gffList.append(row)

    with open( "e00001GFF1.csv", "w") as f:
        writer = csv.writer(f, delimiter='\t')
        writer.writerows(gffList)
        

main()
