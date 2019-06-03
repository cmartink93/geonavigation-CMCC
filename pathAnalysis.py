import os, sys
import time
import createGraph as cg
import digraph as d
import readVISIR as ri
import writeInterface as wi
import weighted as wgt
import pathInfo as pi
import edgelists as ed
import coords as c
import readGMSH as rg
import calcWeights as cwts



def findPath(files, outLog, outPath, parmsDict, nameDict, graph):
    '''Calculates and writes out shortest path sequence/cost based on specified input data.
    files: list of directories/file names for file retrieval
    outLog: name/directory of file that path statistics will be written to
    outPath: name/directory of file that path nodes will be written be to
    parmsDict: csv file dictionary with directory specifications/flags
    nameDict: csv file dictionary with user specifications/flags
    graph: graph to be used to find shortest path
    
    returns:
    path: list of all nodes in shortest path
    ex: list of extremes in the form [startID, targetID]
    '''
    if parmsDict['VISIRbenchmark_flag'][0] == '1':
        extremeFile = files[2]
        exStr = cg.readCoords(extremeFile)
        startID = cg.getExtreme(exStr, 0)
        targetID = cg.getExtreme(exStr, 1)  
      

        
    if parmsDict['VISIRbenchmark_flag'][0] == '2':
        extremeFile = files[1] 
        coordList = cg.readCoords(extremeFile)
    
        startID = cg.coordsToExtremes(graph, coordList, 'START')
        targetID = cg.coordsToExtremes(graph, coordList, 'TARGET')
        
    ex = []    
    ex.append(startID)
    ex.append(targetID)        
    print startID
    print targetID
    print graph.node[startID]['lat'], graph.node[startID]['lon']
    print graph.node[targetID]['lat'], graph.node[targetID]['lon']
    
   
    
    #optRoutines = ['0','2','4','5']
    #outLogNames = ['dijkstra_ss','dijkstra_sm','bellman_ford','astar']
    #dirName = str(parmsDict['route_name'][0])
    
    #for i in range(4):
        #print 'working on: ' + outLogNames[i]
        
        #out_folder = '../out/outGmsh/' + dirName + '/'
        #outLFilename = out_folder + outLogNames[i] + "_LOG.txt"
        #outPFilename = out_folder + outLogNames[i] + "_PATH.txt"
        #nameDict['opt_routine'][0] = optRoutines[i]
        
        #outLog = open(outLFilename, 'w')
        #outPath = open(outPFilename, 'w')         
        
    t0 = time.time()
    cost,path,tf = pi.selectOpt(nameDict, graph, startID, targetID)
    tfp = time.time()-t0
    nVis = pi.numVisited(graph)
    print "time to find path: " + str(tfp)
    print "great circle route for this path should be: " + str(cwts.calcDist(graph, startID,targetID))
    print "calculated cost is: " + str(cost)
    wi.printPath(path, outPath)
    wi.printLog(path, startID, targetID, cost, tf, nVis, outLog, graph)
    #outLog.close()
    #outPath.close()        
        
    
    return path, ex
   



def runOnce(files, parmsDict, nameDict):
    '''Runs through the entire shortest path procedure,from
    dijkstra to reporting the results.
    files: list of directories/file names for file retrieval
    parmsDict: csv file dictionary with directory specifications/flags
    nameDict: csv file dictionary with user specifications/flags
    
    returns: 
    graph: final graph used
    path: list of all nodes in shortest path
    ex: list of extremes in the form [startID, targetID]
    '''

    
    if parmsDict['VISIRbenchmark_flag'][0] == '1':
        outLFilename = files[4]
        outPFilename = files[5]
        loadGraphName = files[8]
        #fix loading graph for VISIR graphs
        
    if parmsDict['VISIRbenchmark_flag'][0] == '2':
        outLFilename = files[2]
        outPFilename = files[3]
        loadGraphName = files[5]
        mshName = files[0]
        
    
        
    outLog = open(outLFilename, 'w')
    outPath = open(outPFilename, 'w')   
    
    if nameDict['need_previous_graph'][0] == '0':
        graph = cg.makeAGraph(files, parmsDict, nameDict)
        
        
    if nameDict['need_previous_graph'][0] == '1':
        t0 = time.time()
        graph = ed.read_edgelist(loadGraphName, create_using=d.DiGraph(), nodetype=int, edgetype=float)
        mesh = open(mshName, 'r')
        nodeID = rg.getNodes(mesh)
        mesh.close()
        cg.gmshNodeProps(graph, nodeID)
        tf = time.time() - t0
        print "time to read in and add properties to graph: " + str(tf)
       
          
    
    t0 = time.time()
    path, ex = findPath(files, outLog, outPath, parmsDict, nameDict, graph)
    tf = time.time() - t0
    
    outLog.write("TOTAL JOB COMPUTATION TIME: " + str(tf) + "\n")    

    outLog.close()
    outPath.close()
    
    return graph, path, ex
    


