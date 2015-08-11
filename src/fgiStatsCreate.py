#--------------------------------#
# Genome fGI Insert Data for GFF #
#         Namita Dongre          #
#           7/24/2015            #
#--------------------------------#


import csv


def main():
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
    count = 1
    for insert in genomeDict:
        (startCluster, endCluster, feature) = (genomeDict[insert][0], genomeDict[insert][1], genomeDict[insert][2])
        (start, contig) = (coreDict[startCluster][0], coreDict[startCluster][2])
        orientation = '+'
        print(feature)
        if feature != 'Break':
            end = coreDict[endCluster][1]
        else: end = coreDict[startCluster][1]
        if int(start) > int(end):
            (start, end, orientation) = (end, start, '-')
        referenceList.append([contig, 'PanOCT', feature, start, end, '.', orientation, '.', 'ID='+feature+str(count)+';Name='+feature+';Note='+insert])

        (startCluster, endCluster, feature) = (genomeDict[insert][0], genomeDict[insert][1], genomeDict[insert][2])
        (start, contig) = (coreDict[startCluster][1], coreDict[startCluster][2])
        orientation = '+'
        if endCluster != 'X':
            end = coreDict[endCluster][0]
            if int(start) > int(end):
                (start, end, orientation) = (end, start, '-')
            referenceList.append([contig, 'PanOCT', 'fgi_CL', start, end, '.', orientation, '.', 'Parent='+feature+str(count)+';Index=1;Note='+insert[insert.index(':')+1:insert.rindex(':')]])
        count+=1


    with open("e00002GFF.csv", "w") as f:
        writer = csv.writer(f, delimiter='\t')
        writer.writerows(referenceList)

#--------------Method to split up fGI cluster strings into list of cluster IDs + orientation
def clusterSplit(stng):
    clusters = []
    stng = stng[stng.index('E')+1::]
    while ':' in stng:
        clusters.append(stng[0:stng.index(':')])
        stng = stng[stng.index(':')+1::]
    if 'E' in stng:
        clusters.append(stng[stng.index('E')+1:len(stng)])
    else: clusters.append(stng[0:len(stng)])
    return clusters

main()

