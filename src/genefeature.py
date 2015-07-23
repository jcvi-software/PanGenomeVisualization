import csv

def main():
    #---E00001
    reader_Core_attfGI = csv.reader(open('Core_attfGI.csv', 'rt', encoding="ascii"), delimiter=',')
    reader_combined = csv.reader(open('combined_att_file.csv', 'rt', encoding="ascii"), delimiter=',')
    reader_panoct = csv.reader(open('panoct.csv', 'rt', encoding="ascii"), delimiter=',')

    print(True)
    genomeID = 'E00001'
    e00001GeneDict = dict()
    for row in reader_combined:
        (reference, geneID, start, stop, info, genome) = (row[0], row[1], row[2], row[3], row[4], row[5])
        if genome == genomeID:
            e00001GeneDict[geneID] = [start, stop]

    panoctList = list()
    for row in reader_panoct:
        row = [element for element in row if element != '']
        panoctList.append(row)

    clusterDict = []
    for gene in genomeDict:
        for row in panoctList:
            if gene in row:
                clusterDict[gene].append(row[0])
                print(clusterDict[gene])

##    gffList = list()
##    for gene in e0001GeneDict:
##        (start, end) = e0001GeneDict[gene]
##
##
##    with open( "e00001GFF.csv", "w") as f:
##        writer = csv.writer(f, delimiter='\t')
##        writer.writerows(e00001List)
        

main()
