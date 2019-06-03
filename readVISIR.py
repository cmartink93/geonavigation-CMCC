import csv
import sys



def readVISIRgraphFiles(edgesFileName, weightsFileName, coordsFileName):
    '''Reads VISIR input files and creates lists of all relevant information
    edgesFileName: file containing edges information
    weightsFileName: file containing weight information
    coordsFileName: file containing coordinates of nodes
    
    returns:
    startList: each row of the edges file has a start and end node of each edge. The startList has all of the start nodes.
    endList: all end nodes of edge file
    wList: weights for each edge. Goes in same order as edge file edges. For example, the first row of the start/end lists corresponds to the first row of the wList.
    lonList: same concept as edge file and the start/end lists. lonList is a list of all node longitudes.
    latList: list of all node latitudes.
    '''

    weights = open(weightsFileName, 'r')
    edges   = open(edgesFileName,   'r')
    coords  = open(coordsFileName,  'r')

    wList     = makeWeightDict(csv.reader(weights))
    startList = makeList(csv.reader(edges), 0, 0) 
    edges.seek(0) 
    endList = makeList(csv.reader(edges), 1, 0)
    lonList = makeList(csv.reader(coords), 0, 1)
    coords.seek(0) 
    latList = makeList(csv.reader(coords), 1, 1)

    
    return startList, endList, wList, lonList, latList


def makeList(reader, column, varType):
    '''Uses input VISIR files to create lists of data.
    reader: input reader that reads one of the VISIR csv input files.
    column: for a file that has two columns (i.e. edge or coord files), specifies which column to make a list of
    varType: indicates whether to convert input data from a string into either an integer (0) or float (1)
    
    returns: list of input data from csv file
    '''
    listName = []
    for row in reader:
        if varType == 0:
            w = int(row[column])
        else:
            w = float(row[column])
        listName.append(w)
    return listName

def makeWeightDict(reader):
    '''Makes a list of weights using the VISIR input file
    reader: input reader that reads the weight VISIR csv input file.
    
    returns: a list of weight lists. Each element of the returned list, d, is a list of all edge weights at a given time step.
    '''
    d = []
    for row in reader:
        l = []
        for col in row:
            w = float(col)
            l.append(w)
        d.append(l)
    return d

def makeDict(fileName):
    '''Makes a dictionary from an input csv file, specifically the parmsDict.csv or nameDict.csv files
    fileName: name of input csv file
    
    returns: dictionary of all flags. key = flag name, value = flag itself
    '''
    fileDict = {}
    f = open(fileName, 'rU')
    for row in csv.reader(f, dialect='excel'):
        values = []
        values.append(row[0])
        values.append(row[2])
        values.append(row[3])
        fileDict[row[1]] = values
    return fileDict
