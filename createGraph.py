import digraph as d
import readVISIR as ri
import readGMSH as rg
import math as m
import intersection as inter
import calcWeights as cwts
import edgelists as ed
import sys
import time
import nx_pylab as nxp
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap as bmap
from scipy import interpolate
import numpy as np




def makeAGraph(files, parmsDict, nameDict):
    '''Creates a graph based on input data and user specifications.
    files: List of locations/file names for retrieving input data
    parmsDict: csv file dictionary of directory information/flags
    nameDict: csv file dictionary of user specifications/flags
    
    returns: graph with all necessary properties
    '''
    if parmsDict['VISIRbenchmark_flag'][0] == '1':
        t0 = time.time() 
        edgeFile = files[0]
        weightFile = files[1]
        coordsFile = files[3]        
        graph = makeVISIRGraph(edgeFile, weightFile, coordsFile, nameDict)
        tf = time.time()-t0
        print "time to make graph: " + str(tf)

            
    if parmsDict['VISIRbenchmark_flag'][0] == '2':
        mesh = files[0]       
        t0 = time.time()
        graph = GMSHGraph(mesh,nameDict,files)
        tf = time.time()-t0
        print "time to make graph: " + str(tf)
        ed.write_edgelist(graph, files[4])
    
    return graph  



def GMSHGraph(msh,nameDict,files):
    '''Creates a graph based on GMSH mesh. Adds relevant graph properties.
    msh: GMSH mesh file to be used to create graph (.msh file)
    nameDict: csv file dictionary of user specifications/flags
    files: list of directories/file names for file retrieval
    
    returns: graph    
    '''
    mesh = open(msh, 'r')
    graph = d.DiGraph()
    nodeID = rg.getNodes(mesh)
    mesh.seek(0) 
    start, end, coast = rg.getEdges(mesh)
    mesh.close()
    keys = nodeID.keys()
    for i in range(len(keys)):
        graph.add_node(keys[i])
    for i in range(len(start)-1):
        graph.add_edge(start[i], end[i])
    
    gTot = gmshExtraEdges(graph, nameDict)
    
    nProps = gmshNodeProps(gTot, nodeID)
    
    '''Code for checking for intersections and removing nodes'''
    #t0i = time.time()
    #eList = inter.intCheck(gTot, gCoast)
    #tfi = time.time() - t0i
    #print "time to check for intersections: " + str(tfi)
    #print "new number of edges: " + str(len(gTot.edges()))
    
    '''Code for interpolating current magnitudes and creating dictionaries that can then be used to calculate weights based on current magnitudes'''
    #flowArrays = cwts.getFlowArrays(nameDict)
    #vesselSpeedMax = float(nameDict['vessel_speed_max'][0])
    #iModule = flowArrays[4]
    #if vesselSpeedMax < iModule.max():   
        #print "ERROR: Maximum vessel speed does not exceed maximum current speed within domain"
        #sys.exit()    
    #t0funct = time.time()
    #fuDict, fvDict = cwts.interpField(nProps, flowArrays, files)
    #tffunct = time.time() - t0funct
    #print "function time: " + str(tffunct)
   
    '''Adds edge properties. Right now current interpolation is not properly functioning so it will not be an input, but once interpolation is up and running all relevant information will need to be an input for adding graph properties '''
    gmshEdgeProps(gTot, nodeID, nameDict)#, flowArrays, fuDict, fvDict, vesselSpeedMax) 
   
   
   
    return gTot
 
    
def gmshNodeProps(gTot, nodeID):
    '''Adds longitude, latitude, and use properties to the nodes in the graph. Lon/lat are coordinates in degrees, 'use' tracks whether or not the node has been visited by the shortest path algorithm.
    gTot: input graph
    nodeID: dictionary of nodes in mesh. key = node ID, value = tuple, (x_coord, y_coord, z_coord)
    
    returns: list of properties. Each element is a list. Element 0 is a list of longitudes, element 1 is a list of latitudes, and element 2 is a list of node ID numbers
    '''
    nodeProps = []
    nodeLons = []
    nodeLats = []
    nodeNames = []
    for node in gTot.nodes():
        gTot.node[node]['lon'] = cwts.calcLon(nodeID[node][0], nodeID[node][1])
        gTot.node[node]['lat'] = cwts.calcLat(nodeID[node][2])
        gTot.node[node]['use'] = -1
        nodeLons.append(gTot.node[node]['lon'])
        nodeLats.append(gTot.node[node]['lat'])
        nodeNames.append(node)
    nodeProps.append(nodeLons)
    nodeProps.append(nodeLats)
    nodeProps.append(nodeNames)
    return nodeProps
        
def gmshEdgeProps(gTot, nodeID, nameDict):#, flowArrays, fuDict, fvDict, vesselSpeedMax): #have calcEdgeWeight return COG value?
    '''Adds edge properties to graph. The properties are: weight, COG, and dist (great circle distance between two nodes of an edge)
    gTot: graph to add properties to
    nodeID: dictionary of nodes in mesh. key = node ID, value = tuple, (x_coord, y_coord, z_coord)
    nameDict: csv file dictionary of user specifications/flags    
    '''
    
    '''Code for when interpolation works'''
        #for edge in gTot.edges():
            #gTot.edge[edge[0]][edge[1]]['weight'] = cwts.calcEdgeWeight(gTot, edge[0], edge[1], vesselSpeedMax, fuDict, fvDict)
            #print gTot.edge[edge[0]][edge[1]]['weight']
            #gTot.edge[edge[0]][edge[1]]['cog'] = cwts.calcEdgeBearing(gTot, edge[0], edge[1])
            #gTot.edge[edge[0]][edge[1]]['dist'] = cwts.calcDist(gTot, edge[0], edge[1])
        
              
    
    '''Temporary code until the interpolation works (uses distance as weight). Can be deleted when interpolation works and above code should be used instead. However, below code could also be used to handle different kinds of routes. In other words, the above block of code could be for the optimal route and the below code could be for the geodetic route.'''
    for edge in gTot.edges():
        gTot.edge[edge[0]][edge[1]]['weight'] = cwts.calcDist(gTot, edge[0], edge[1])
        gTot.edge[edge[0]][edge[1]]['cog'] = cwts.calcEdgeBearing(gTot, edge[0], edge[1])
        gTot.edge[edge[0]][edge[1]]['dist'] = cwts.calcDist(gTot, edge[0], edge[1])
   

def gmshExtraEdges(graph, nameDict):
    '''Adds extra edges to graph based on user-specified degree.
    graph: graph to add edges to
    nameDict: csv file dictionary with user specifications/flags.
    
    returns: graph with added edges
    '''
    print "# g nodes: " + str(len(graph.nodes()))
    print "# g edges: " + str(len(graph.edges()))
    gTot = d.DiGraph()
    i=0
    nProcessed = 0
    sp={}
    t0 = time.time()
    deg = int(nameDict['con_deg'][0])
    for node in graph.nodes():
        nhbrs = []
        
        for node1st in graph.neighbors(node):
            nhbrs.append(node1st) #1st order
            if deg > 1:
                for node2nd in graph.neighbors(node1st):
                    nhbrs.append(node2nd) #2nd order
                    if deg > 2:
                        for node3rd in graph.neighbors(node2nd):
                            nhbrs.append(node3rd) #3rd order 
                            if deg > 3:
                                for node4th in graph.neighbors(node3rd):
                                    nhbrs.append(node4th) #4th order
                                    if deg > 4:
                                        for node5th in graph.neighbors(node4th):
                                            nhbrs.append(node5th) #5th order 
                                            if deg > 5:
                                                for node6th in graph.neighbors(node5th):
                                                    nhbrs.append(node6th) #6th order                             
      
        for nde in nhbrs:
            if node != nde: # GM: comment out
                gTot.add_edge(node, nde) 
    
        i=i+1 # GM check if i++1
        nProcessed = nProcessed + 1 # GM check if nProcessed++1


    print "gTot nodes: " + str(len(gTot.nodes()))
    print "gTot edges: " + str(len(gTot.edges()))
    print "Nodes processed: " + str(nProcessed)
    tf = time.time() - t0
    print "time to go through " + str(nProcessed) + " nodes: " + str(tf)
 
    return gTot
        



def readCoords(extremesFileName):
    '''Extracts coordinates from a file
    extremesFileName: file containing coordinates of start/target nodes
    
    returns: list of start/target node coordinates in sequential order
    '''
    extremes = open(extremesFileName, 'rU')
    l = []
    for row in extremes:
        l.append(row)
    return l


def coordsToExtremes(g, coords, startOrTarget):
    '''Find the node in the graph that is closest to the user's desired start or target coordinates.
    g: graph to search
    coords: coordinate list (in the form [start_lat,start_lon,target_lat,target_lon]) of desired start or target
    startOrTarget: flag to specify whether method should look at first two elements in list (user's desired start coords) or last two elements in list (user's desired end coords)
    
    returns: number ID of node that best matches user's specifications
    '''
    if startOrTarget == 'START':
        inLat = float(coords[0])
        inLon = float(coords[1])
    else:
        inLat = float(coords[2])
        inLon = float(coords[3])  
    bestDiff = 100000
    for n in g.nodes():
        lat = g.node[n]['lat']
        lon = g.node[n]['lon']
        diff = cwts.calcDistPreproc(inLat, inLon, lat, lon)
        if diff <= bestDiff:
            bestDiff = diff
            best = n
    return best



def VISIRGraph(startList, endList, wList, lonList, latList, nameDict): #TO DO: make it add interpolation properties (similar to what gmshEdgeProps() does)
    '''Creates graph using input data from VISIR. Adds relevant properties
    startList: list of start nodes for each edge
    endList: list of end nodes for each edges
    wList: list of weights for each edge
    lonList: list of longitudes for each node
    latList: list of latitudes for each node
    nameDict: csv file dictionary of user specifications/flags
    
    returns: graph
    '''
   
    g = d.DiGraph()
    for i in range(len(startList)):
        if nameDict['is_tdep'][0] == '0':
            if i < len(lonList):
                g.add_node(i+1, lon = lonList[i], lat = latList[i], use = -1)
                g.add_edge(startList[i], endList[i], weight = wList[i][0])# to generalize for multiple keys (assuming new key file structure is the same as for the weights), simply add: keyname = keyDict[i]
                #dictLoop(g, i+1, startList[i], endList[i], wDict)
            else:
                g.add_edge(startList[i], endList[i], weight = wList[i][0])# to generalize for multiple keys (assuming new key file structure is the same as for the weights), simply add: keyname = keyDict[i]
                #dictLoop(g, i+1, startList[i], endList[i], wDict)
        if nameDict['is_tdep'][1] == '1':
            if i < len(lonList):
                g.add_node(i+1, lon = lonList[i], lat = latList[i], use = -1)
                g.add_edge(startList[i], endList[i], weight = wList[i])# to generalize for multiple keys (assuming new key file structure is the same as for the weights), simply add: keyname = keyDict[i]
                #dictLoop(g, i+1, startList[i], endList[i], wDict)
            else:
                g.add_edge(startList[i], endList[i], weight = wList[i])# to generalize for multiple keys (assuming new key file structure is the same as for the weights), simply add: keyname = keyDict[i]
                #dictLoop(g, i+1, startList[i], endList[i], wDict) 
    for edge in g.edges():
        g.edge[edge[0]][edge[1]]['cog'] = cwts.calcEdgeBearing(g, edge[0], edge[1])
        g.edge[edge[0]][edge[1]]['dist'] = cwts.calcDist(g, edge[0], edge[1])
    return g

    
def makeVISIRGraph(edgeFile, weightFile, coordsFile, nameDict):
    '''Retrieves VISIR input data and feeds it to graph method to create a graph
    edgeFile: file containing a list of all edges
    weightFile: file containing a list of weights for each edge
    coordsFile: file containing a list of (lat,lon) coordinates for each node
    nameDict: csv file dictionary of user specifications/flags
    
    returns: graph
    '''
    startList, endList, wList, lonList, latList = ri.readVISIRgraphFiles(edgeFile, weightFile, coordsFile)
    g = VISIRGraph(startList, endList, wList, lonList, latList, nameDict)    
    return g


    
def getExtreme(extremesList, nodeType):
    '''Will retrieve start or end node ID. Enter nodeType = 0 for start and nodeType = 1 for target node.
    extremesList: list of start and end node IDs
    nodeType: numerical flag specifying where to look in extremes list for node ID
    
    returns: integer node ID of the start or end extreme
    '''
    
    if nodeType == 0:
        n = int(extremesList[0].strip("\n"))
    else:
        n = int(extremesList[1].strip("\n"))
    return n



#NOTE: the following methods are only used/useful for visualizing the graph mid-creation, for example to see which edges were removed.

def makeLabels(inPath):
    '''Creates label dictionary
    inPath: list of nodes in path
    
    returns: dictionary of labels to be used for plotting. Key=node ID, value= string of node ID    
    '''
    labels={}
    for node in inPath:
        labels[node] = str(node)
    return labels

def mPos(graph):
    '''Creates position dictionary of nodes to plot
    graph: graph to plot
    
    returns: position dictionary. key= node ID, value= tuple, (longitude, latitude)
    '''
    posDict = {}
    for node in graph.nodes():
        posDict[node] = (graph.node[node]['lon'], graph.node[node]['lat'])
    return posDict   

def plotTemp(graph,gRem,gCoast):
    '''Plots graph mid-processing to see if certain functions, like edge removal, are working properly.
    graph: main graph to evaluate
    gRem: graph of all removed edges
    gCoast: graph of coastal nodes    
    '''
    fig = plt.figure(figsize=(8.0, 5.0))
    ax = fig.add_subplot(111)

    plt.xlabel('E Long. [deg]', fontsize=20)
    plt.ylabel('N Lat. [deg]', fontsize=20) 
    #plt.grid(True, which='both')
    
    fig.suptitle('Shortest Path', fontsize=20)

    for tick in ax.xaxis.get_major_ticks():
        tick.label.set_fontsize(15)
    for tick in ax.yaxis.get_major_ticks():
        tick.label.set_fontsize(15)     
    pDict = mPos(graph)
    rpDict = mPos(gRem)
    cpDict = mPos(gCoast)
    #m = bmap(llcrnrlon = -95, llcrnrlat = 10, urcrnrlon =-80, urcrnrlat = 30, resolution = 'h', projection = 'merc',lon_0=-90, lat_0=20)
    #m.drawcoastlines()
    nxp.draw_networkx_edges(graph, pos=pDict, edgelist=graph.edges(), edge_color='black', arrows=False)
    nxp.draw_networkx_edges(gRem, pos=rpDict, edgelist=gRem.edges(), edge_color='mediumvioletred', arrows=False, width=3)
    nxp.draw_networkx_nodes(gCoast, pos=cpDict, nodelist=gCoast.nodes(), node_color='limegreen', node_size=25)
    #nxp.draw_networkx_nodes(graph, pos=pDict, nodelist=graph.nodes(), node_color='black', node_size=10)
    #nxp.draw_networkx_labels(graph, pos=pDict, labels=makeLabels(graph.nodes()), fontsize=11)
    #nxp.draw_networkx_labels(gCoast, pos=cpDict, labels=makeLabels(gCoast.nodes()), fontsize=11)
    
    plt.show()
    
def remGraph(eList, nodeID):
    '''Creates a graph of all removed edges
    eList: list of removed edges, where each element is a tuple edge, (node1, node2)
    nodeID: dictionary of nodes in mesh. key = node ID, value = tuple, (x_coord, y_coord, z_coord)
    
    returns: graph of all removed edges    
    '''
    gRem = d.DiGraph()
    for edge in eList:
        gRem.add_edge(*edge)
    gmshNodeProps(gRem,nodeID)
    return gRem

def coastGraph(coastList, nodeID):
    '''Creates a graph of all coastal nodes
    coastList: list of all coastal nodes
    nodeID: dictionary of nodes in mesh. key = node ID, value = tuple, (x_coord, y_coord, z_coord)
    
    returns: graph of all coastal nodes
    '''
    gCoast = d.DiGraph()
    for i in range(len(coastList)):
        gCoast.add_node(coastList[i], lon = calcLon(nodeID[coastList[i]][0], nodeID[coastList[i]][1]), lat = calcLat(nodeID[coastList[i]][2]))
    return gCoast

def coastGraphGEODAS(shoreline):
    '''Creates coastal graph based on shoreline data from GEODAS
    shoreline: shoreline database (.a00 file)
    
    returns: graph of all coastal nodes    
    '''
    gCoast = d.DiGraph()
    coastFile = open(shoreline, 'r')
    nodeID = 0
    for line in coastFile.readlines():
        gCoast.add_node(nodeID, lon = float(line.split()[0]), lat = float(line.split()[1]))
        nodeID = nodeID + 1
    coastFile.close()
    return gCoast    
