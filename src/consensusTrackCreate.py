#----------------------------------------- ----#
#        Creating the Consensus Track          #
#              Namita Dongre                   #
#----------------------------------------------#

import csv, os


def main():
    #---Consensus
    path, endPath = getFilePath()
    reader_fgi = csv.reader(open(path+'fGI_stats.csv', 'rt', encoding="ascii"), delimiter=',')
    reader_core = csv.reader(open(path+'Core_attfGI.csv', 'rt', encoding="ascii"), delimiter=',')
    referenceList = []

    consensusDict = dict()
    for row in reader_fgi:
        if 'CL_INS' in row[0]:
            consensusDict[row[0]] = (row[1], row[2])

    referenceList = createConsensusFile(consensusDict, reader_core)
    writeFile(endPath, referenceList)

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


def createConsensusFile(consensusDict, reader_core):
    (oldStart, oldEnd) = (0, 0)
    referenceList = []
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
    return referenceList

def writeFile(outpath, referenceList):
    with open(outpath +'consensus.gff', 'w') as f:
        writer = csv.writer(f, delimiter='\t', quoting=csv.QUOTE_NONE)
        writer.writerows(referenceList)
    print('consensus.gff File Created')

#-------------------------------------------------------------------------------------------------------------------Main
main()
