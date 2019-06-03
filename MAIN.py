import sys, csv, os
import readVISIR as ri
import pathAnalysis as pa
import setStructure as ss
import drawPath as dp
import course_plot as cplot
import edgelists as ed
import time







def MAIN():    
#if __name__ == '__main__':
    '''This program takes input data, can create a graph (or load a graph, if specified), add graph properties, and calculate the shortest path from a start to target node using a user-specified algorithm. It will then record all information, provide a graphical representation of the path, and provide a COG plot for reference. See README.txt for details about directory structure and necessary inputs. To run, type "python MAIN()" at the command line.'''
   
    print "Running..."  
    
    parmsDict = ri.makeDict('../in/parmsDict.csv')

    
    ##VISIR
    if parmsDict['VISIRbenchmark_flag'][0] == '1':
        nFile = open('../in/inVisir/' + str(parmsDict['parms_name'][0]), 'r')

        for num in nFile:
            
            gdtFiles = ss.makeFileList(parmsDict, num, 0)
            statFiles = ss.makeFileList(parmsDict, num, 1)
            
            nameDict = ri.makeDict(statFiles[6])
            t0 = time.time()
            print "Finding shortest path..."
            gGdt, pathGdt, exGdt = pa.runOnce(gdtFiles, parmsDict, nameDict)
            gStat, pathStat, exStat = pa.runOnce(statFiles, parmsDict, nameDict)
            tf = time.time() - t0
            print "total runOnce time: " + str(tf)         
            
            t0 = time.time()
            print "Rendering..."
            graphGdt = dp.drawMain(gGdt, gdtFiles[9], pathGdt, nameDict, exGdt)
            graphStat = dp.drawMain(gStat, statFiles[9], pathStat, nameDict, exStat)
            tf = time.time() - t0
            print "total drawMain time: " + str(tf)
            
            t0 = time.time()
            cplot.cogPlotsMain(pathGdt, graphGdt, nameDict, gdtFiles[10])
            cplot.cogPlotsMain(pathStat, graphStat, nameDict, statFiles[10])
            tf = time.time() - t0
            print "total cogPlots time: " + str(tf)            
            
        nFile.close() 
    
    ##GMSH
    if parmsDict['VISIRbenchmark_flag'][0] == '2':
        
        files = ss.makeFileList(parmsDict)
        nameDict = nameDict = ri.makeDict(files[9])
        
        t0 = time.time()
        print "Finding shortest path..."        
        gTot, path, ex = pa.runOnce(files, parmsDict, nameDict)
        tf = time.time() - t0
        print "total runOnce time: " + str(tf) 
        
        #sys.exit()
        
        t0 = time.time()
        print "Rendering..."
        usedGraph = dp.drawMain(gTot, files[6], path, nameDict, ex)
        tf = time.time() - t0
        print "total drawMain time: " + str(tf)
    
        
        t0 = time.time()
        cplot.cogPlotsMain(path, usedGraph, nameDict, files[7])
        tf = time.time() - t0
        print "total cogPlots time: " + str(tf) 
        
        
    
    print "Run complete."
           



MAIN()