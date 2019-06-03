import sys
import math as m

#step one: reads as if non-directional
#step two: TODO: make directional, improve efficiency

def getNodes(msh):
    '''Retrieves all nodes in mesh file
    msh: mesh file (.msh)
    
    returns: dictionary of nodes in mesh. key= node ID, value = tuple, (x_coord, y_coord, z_coord)
    '''
    nodeID = {}
    for line in msh.readlines()[5:]: #fix range to only look at $Nodes section, also fix so it finds the first node automatically
        if len(line.split()) != 4:
            continue
        else:
            if line.split()[0] == '$EndNodes':
                break 
            nodeID[int(line.split()[0])] = (float(line.split()[1]), float(line.split()[2]), float(line.split()[3].strip('\n')))
    return nodeID

def getEdges(msh): #right now set to create undirected edges
    '''Retrieves edges from mesh file
    msh: mesh file (.msh)
    
    returns: 
    start: list of start nodes in edges
    end: list of end nodes in edges. corresponds with start list, i.e. start[0] and end[0] are the start and end nodes of an edge
    coast: list of coastal nodes, not necessarily useful with new meshes. Should be evaluated.
    '''
    start = []
    end = []
    coast = []
    for line in msh.readlines(): #fix range to only look at $Elements section
       
        if len(line.split()) < 6:
            continue
        #if line.split()[1] == '1': #for old, uniform meshes
            #start.append(int(line.split()[5]))
            #end.append(int(line.split()[6].strip('\n')))
            
            #start.append(int(line.split()[6].strip('\n')))
            #end.append(int(line.split()[5]))
            
        elif line.split()[1] == '2':
            start.append(int(line.split()[5]))
            end.append(int(line.split()[6]))
            
            start.append(int(line.split()[6]))
            end.append(int(line.split()[5]))
            
            start.append(int(line.split()[5]))
            end.append(int(line.split()[7].strip('\n')))
            
            start.append(int(line.split()[7].strip('\n')))
            end.append(int(line.split()[5])) 
            
            start.append(int(line.split()[6]))
            end.append(int(line.split()[7].strip('\n')))  
            
            start.append(int(line.split()[7].strip('\n')))
            end.append(int(line.split()[6]))
            
        #elif line.split()[1] == '15': #for old, uniform meshes
            #if line.split()[5] != '1':
                #if line.split()[5] != '2':
                    #coast.append(int(line.split()[5]))
        elif line.split()[1] == '1': 
            if int(line.split()[5]) not in coast:
                coast.append(int(line.split()[5]))  
            if int(line.split()[6]) not in coast:
                coast.append(int(line.split()[6]))             
            
    return start, end, coast


 
 

