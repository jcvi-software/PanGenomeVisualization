#------------------------------#
#   Gene GFF Program           #
#        Namita Dongre         #
#------------------------------#

import csv


def main():
    #---Consensus
    reader_fgi = csv.reader(open('fGI_stats.csv', 'rt', encoding="ascii"), delimiter=',')
    reader_core = csv.reader(open('Core_attfGI.csv', 'rt', encoding="ascii"), delimiter=',')
    referenceList = []

    #-----------------Created a dictionary of {Insertion ID: (cluster1_direction, cluster2_direction)}
    #-----------------------In some cases, only 1 cluster present
    #-----------------------Direction denoted by 5' or 3'
    consensusDict = dict()
    for row in reader_fgi:
        if 'CL_INS' in row[0]:
            consensusDict[row[0]] = (row[1], row[2])

    #-----------------Created a list of GFF data for fGI inserts and Clusters of consensus
    (oldStart, oldEnd) = (0, 0)
    for row in reader_core:
        (contig, insertionID, start, end, name, feature, aminoAcids) = (row[0], row[1], row[2], row[3], row[4], row[5], row[6])
        orientation = '+'
        (start, end) = (int(int(start)/3),int(int(end)/3))
        if start > end:
            (start, end, orientation) = (end, start, '-')
        if start == oldEnd:
            start+=1
        end = start+int(aminoAcids)-1
        if 'CL_INS' in insertionID:
            (cluster1, cluster2) = consensusDict[insertionID]
            note = cluster1 +' '
            if cluster2 != '':
                note += cluster2
            annotation = 'Name=' + insertionID +'; Note=' + note + '; fGI_info=' + name;
            referenceList.append([contig, 'PanOCT', 'fGI_INS', start, end, '.', orientation, '.', annotation])
        else:
            annotation = 'Name= ' + insertionID[insertionID.index('_')+1::] + '; Note= ' + name
            referenceList.append([contig, 'PanOCT', feature, start, end, '.', orientation, '.', annotation])
        (oldStart, oldEnd) = (start, end)

    with open("consensus.gff", "w") as f:
        writer = csv.writer(f, delimiter='\t')
        writer.writerows(referenceList)

main()
