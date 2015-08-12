#-------------------------#
#   Namita Dongre         #
#   8/7/2015              #
#   Genome GFF Tracks     #
#-------------------------#

import csv, pickle

def main():
    reader_combined = list(csv.reader(open('combined_att_file.csv', 'rt', encoding="ascii"), delimiter=','))
    reader_all = list(csv.reader(open('all_clusters_members.csv', 'rt', encoding="ascii"), delimiter=','))

    genomeClusterDict = createGenomeClusterDict(reader_all)
    genomeLocusDict = createGenomeLocusDict(reader_combined)

    # print(genomeLocusDict['SR69_RS17940'])
    #
    # for row in genomeClusterDict:
    #     if genomeClusterDict[row][12] == 'SR69_RS17940':
    #         print(genomeClusterDict[row])

    pickle.dump(genomeClusterDict, open('genomeCluster.dict', 'wb'))
    pickle.dump(genomeLocusDict, open('genomeLocus.dict', 'wb'))


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




main()
