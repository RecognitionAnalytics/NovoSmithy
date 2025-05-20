import numpy as np
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from sklearn.manifold import MDS
from scipy.ndimage import fourier_gaussian
from mpl_toolkits.axes_grid1 import make_axes_locatable
import networkx as nx

sns.set_style("whitegrid")
sns.set_context("paper")
#############################################################################################################################################
#   plotting functions    


     

def PlotKMCNetworks(activeAminos,successDwellTimes, dwellTimes, passes)    :
    #create three 3d  subplots with the x y and z views
    fig, ax = plt.subplots(2, 2, figsize =(6,6), subplot_kw={'projection': '3d'})
    ax = np.ravel(ax)
    
    PlotColorCoords(activeAminos,np.log( successDwellTimes+1), 0, ax[0], 'Success Dwell Times')
    PlotColorCoords(activeAminos,np.log(dwellTimes+1), 1, ax[1], 'Dwell Times')

    PlotColorCoords(activeAminos,np.log(passes+1), 2, ax[2], 'Passes')
    
    plt.tight_layout()
    plt.show()
    
def PlotKMC( electronTimes, diffusions,title):
    

    _, ax = plt.subplots(1, 3, figsize =(9,3))
    bins = np.linspace(0, np.min([np.max(diffusions),100]) , 50)
    ax[0].hist(diffusions,bins=bins )
    ax[0].set_xlabel('Average Diffusion nm^2/ns')
    ax[0].set_ylabel('Count')
    ax[0].set_title('Diffusion Distribution')
    
    ax[1].hist(np.log10(electronTimes), bins=50)
    ax[1].set_xlabel('Transit Times log10(ns)')
    ax[1].set_ylabel('Count')
    ax[1].set_title('Transit Times')
    
    
    #ACS Omega paper page G  using this to give a rough estimate of the current
    currents_sum = np.log10(6.24e-18/(np.array(electronTimes)*1e-9)+1e-25)+9
    ax[2].hist(currents_sum, bins=50)
    ax[2].set_xlabel('Current log10(nA)')
    ax[2].set_ylabel('Count')
    ax[2].set_title('Single Electron Currents')
    plt.suptitle(title)
    plt.tight_layout()
    plt.show()     
   

def plotDensity(struct,fig, ax, title, units,extent,redoxPoints, activeAminos):            

        c=ax.imshow(struct , extent=extent, cmap='bwr')
        divider = make_axes_locatable(ax)
        cax = divider.append_axes('right', size='5%', pad=0.05)
        cb=fig.colorbar(c, cax=cax, orientation='vertical')
        cb.set_label(units)
        
        ax.set_title(title)
        ax.set_xlabel('Warped nm')    
        ax.set_ylabel('Warped nm')    
        
        redoxPoints = np.array(redoxPoints)
        #normalize the metric
        redoxPoints = (redoxPoints-np.min(redoxPoints))/(np.max(redoxPoints)-np.min(redoxPoints))
        #include the start and end points 
        cc=0
        for atom in activeAminos:
            if atom['startInject'] or atom['endEject']:
                ax.plot( atom['embedding'][0],atom['embedding'][1],'go',markersize=10)
            else:
                ax.plot( atom['embedding'][0],atom['embedding'][1], '.',color=(0,redoxPoints[cc],0),markersize=5)        
            cc+=1
            
 

def CreateProteinManifold(atom_COM, verbose=False):
    """ All atoms will be used to create a 2D map of the protein using the euclidean distance between the atoms as the metric
        embedding will be added to the atom_COM dictionary
    """
    newCoords =np.array( [acid['centerOfMass'] for acid in atom_COM])
    if verbose:
        print('Potentials mapped')
    embedding = MDS(n_components=2, normalized_stress='auto')
    X_transformed =embedding.fit_transform ( newCoords )
    
    for i in range(len(atom_COM)):
        atom_COM[i]['embedding'] = X_transformed[i]
    #convert newCoords to transformed array
    if verbose:
        print('Embedding Found')   
        
    
def PlotPotentialMap(atom_COM):
     
    potentials = [acid['potenial_EV'] for acid in atom_COM]
    totalPotentials = [acid['totalPotential'] for acid in atom_COM]
    print('Total pot (eV)', np.max(totalPotentials), np.min(totalPotentials))
    
    charges = [acid['charge'] for acid in atom_COM]
    print('charges (e)', np.max(charges), np.min(charges))
    redoxPoints = [acid['redoxEnergy_EV'] for acid in atom_COM]    
    vibrationPoints = [acid['std'] for acid in atom_COM]
   
    X_transformed = np.array( [acid['embedding'] for acid in atom_COM]    )

    #get max and min of transformed coordinates
    maxX,maxY = np.max(X_transformed[:,0]) ,np.max(X_transformed[:,1]) 
    minX,minY = np.min(X_transformed[:,0]), np.min(X_transformed[:,1])
    xdim=100
    ydim=100
    lX=(xdim-1)/( maxX-minX)
    lY= (ydim-1)/( maxY-minY)
    apbsStruct = np.zeros((xdim,ydim),dtype=np.float32)
    chargePots = np.zeros((xdim,ydim),dtype=np.float32)
    redoxStruct = np.zeros((xdim,ydim),dtype=np.float32)+np.max(redoxPoints)
    chargeStruct = np.zeros((xdim,ydim),dtype=np.float32)
    vibrateStruct = np.zeros((xdim,ydim),dtype=np.float32)
    
    
    for i in range(len(X_transformed)):
        x = int((X_transformed[i,0]-minX)*lX)
        y = int((X_transformed[i,1]-minY)*lY)
        apbsStruct[y,x] =  potentials[i]
        chargePots[y,x] =  totalPotentials[i] 
        if redoxPoints[i] >0:
            redoxStruct[y,x] =  redoxPoints[i]
        if charges[i] != 0:
            chargeStruct[y,x] =  charges[i]
        if vibrationPoints[i] !=0:
            vibrateStruct[y,x] = np.log10( vibrationPoints[i] )
        #add in field potential
        
    #plt.imshow(microStruct , origin='lower')     
    print('Smoothing Energy Map')
    sigma2=2
    apbsStruct = np.fft.ifftn(fourier_gaussian(np.fft.fftn(apbsStruct), sigma=sigma2)).real 
    apbsStruct = np.max(potentials)/np.max(apbsStruct)*apbsStruct
    
    chargePots = np.fft.ifftn(fourier_gaussian(np.fft.fftn(chargePots), sigma=sigma2)).real 
    chargePots = np.max(totalPotentials)/np.max(chargePots)*chargePots
    
    redoxStruct = np.fft.ifftn(fourier_gaussian(np.fft.fftn(redoxStruct), sigma=2)).real 
    min = np.min( [x for x in  redoxPoints if x>0])
    eRange = np.max(redoxPoints)-min
    redoxStruct = eRange* (redoxStruct- np.min(redoxStruct))/( np.max(redoxStruct)-min) +min
    
    chargeStruct = np.fft.ifftn(fourier_gaussian(np.fft.fftn(chargeStruct), sigma=1)).real 
    chargeStruct = np.max(charges)/np.max(chargeStruct)*chargeStruct
    
    vibrateStruct = np.fft.ifftn(fourier_gaussian(np.fft.fftn(vibrateStruct), sigma=1)).real
    vibrateStruct = np.max(np.log10( vibrationPoints))/np.max(vibrateStruct)*vibrateStruct
    
    fig, ax = plt.subplots(2,3,figsize=(13.5,8))
    ax=np.ravel(ax)

    def plotDensitySub(ax,metric,title, colorbar):
        c=ax.imshow(metric , extent=[minX,maxX,maxY,minY], cmap='bwr')
        divider = make_axes_locatable(ax)
        cax = divider.append_axes('right', size='5%', pad=0.05)
        cb=fig.colorbar(c, cax=cax, orientation='vertical')
        cb.set_label(colorbar)
        ax.set_title(title)
        ax.set_xlabel('Warped nm')    
        ax.set_ylabel('Warped nm')    
        #include the start and end points 
        for atom in atom_COM:
            if atom['startInject'] or atom['endEject']:
                ax.plot( atom['embedding'][0],atom['embedding'][1],'kx',markersize=10)
        
    plotDensitySub(ax[5],vibrateStruct,'Vibration (std)', 'log$_{10}$(nm)')
    plotDensitySub(ax[4],chargePots+redoxStruct,'Energy Landscape', '(EV)')
    plotDensitySub(ax[3],chargePots,'Partial Charge Voltages', 'Potential (EV)')
    plotDensitySub(ax[2],apbsStruct,'Adaptive Poisson-Boltzmann Voltages', 'Potential (EV)')
    plotDensitySub(ax[1],redoxStruct,'Redox Potential', 'Potential (EV)')
    plotDensitySub(ax[0],chargeStruct,'Partial Charges', 'Charge (e)')
     
    plt.tight_layout()
    plt.show()
    
def PlotColorCoords(activeAminos,metric, colorIndex, ax, title, coordOrder=(1,0,2), invertColors=False):
    """Plot the selected amino acids with a color based on the metric

    Keyword arguments:
    activeAminos -- the list of amino acids to plot (list of dictionaries)
    metric -- the metric to use for the color (list)
    colorIndex -- the index of the RGB color to use (0,1,2)
    ax -- the matplotlib axis to plot on
    title -- the title of the plot
    coordOrder -- the order of the projection coordinates (1,0,2)
    invertColors -- invert the color scale (False)
    """  
    cc=0
    
    #normalize the metric
    minMetric = np.min(metric)
    maxMetric = np.max(metric)
    distMetric = maxMetric-minMetric 
    if distMetric ==0:
        distMetric =1
    
    
    for acid in activeAminos:
        newCoords=(acid['centerOfMass'])
        b= (metric[cc]-minMetric)/distMetric
        if invertColors:
            b=1-b
        color = np.zeros(3)
        color[colorIndex]=b
        ax.plot (newCoords[ coordOrder[0]], newCoords[coordOrder[1]], newCoords[ coordOrder[2]],'o',color=color, alpha=0.3)
        
        cc+=1
     
    ax.set_title(title)   
    ax.set_xlabel('(nm)')
    ax.set_ylabel('(nm)')
    ax.set_zlabel('(nm)') 
    
def PlotBiColorCoords(activeAminos,metric,  ax, title, coordOrder=(1,0,2) ):
    """Plot the selected amino acids with a divergent color based on the metric

    Keyword arguments:
    activeAminos -- the list of amino acids to plot (list of dictionaries)
    metric -- the metric to use for the color (list)
    ax -- the matplotlib axis to plot on
    title -- the title of the plot
    coordOrder -- the order of the projection coordinates (1,0,2)
    """  
    
    cc=0
    if len(metric)==0:
        return
    minMetric = np.min(metric)
    if minMetric ==0:
        minMetric =1
    maxMetric = np.max(metric)
    if maxMetric ==0:
        maxMetric =1
    
    for acid in activeAminos:
        newCoords=(acid['centerOfMass'])
         
        if metric[cc]>0:
            color=(0,0,metric[cc]/maxMetric)
        else :
            color=(0,metric[cc]/minMetric,0)
            
        ax.plot (newCoords[ coordOrder[0]], newCoords[coordOrder[1]], newCoords[ coordOrder[2]],'o',color=color, alpha=0.3)
        
        cc+=1
     
    ax.set_title(title)   
    ax.set_xlabel('(nm)')
    ax.set_ylabel('(nm)')
    ax.set_zlabel('(nm)') 
    
    
def PlotPDBProjections(activeAminos,atom_COM,injectionAminos,exitAminos):
    #create three 3d  subplots with the x y and z views
    fig, ax = plt.subplots(2, 2, figsize =(8,8), subplot_kw={'projection': '3d'})
    ax = np.ravel(ax)

    redoxEnergiesPlot = [acid['redoxEnergy_EV'] for acid in activeAminos]
    PlotColorCoords(activeAminos, redoxEnergiesPlot, 2, ax[0], 'Standard Potential (EV)', coordOrder=(1,0,2), invertColors=True)
    
    stdPlot =np.log10( np.array([acid['std'] for acid in activeAminos]) )
    PlotColorCoords(activeAminos, stdPlot, 2, ax[1], 'Vibrational Motion', coordOrder=(1,0,2), invertColors=False)

    metric = [acid['potenial_EV'] for acid in atom_COM]
    newCoords = []
    cc=0
    minMetric = np.min(metric)
    maxMetric = np.max(metric)
    if minMetric != maxMetric:
        for acid in atom_COM:
            newCoords=(acid['centerOfMass'])
                
            if metric[cc]>0:
                a=1-metric[cc]/maxMetric
                color=(metric[cc]/maxMetric,0,0)
            else :
                a=1-metric[cc]/minMetric
                color=(0,0,metric[cc]/minMetric)
                
            ax[2].plot (newCoords[ 1], newCoords[0], newCoords[ 2],'.',color=color, alpha=a)
            
            cc+=1
        ax[2].set_title('Potential  ')
    else :
        aminos =list(set([acid['amino'] for acid in atom_COM]))
        for amino in aminos:
            newCoords=[]
            for acid in atom_COM:
                if acid['amino']==amino:
                    newCoords.append (acid['centerOfMass'])
            newCoords = np.array(newCoords)
         
            ax[2].plot (newCoords[ :,1], newCoords[:,0], newCoords[:, 2],'.', label = amino)
        ax[2].set_title('Amino Acid  ')
        ax[2].legend()
    
    charges=[acid['charge'] for acid in atom_COM]
    chargeList = []
    chargedAminos = []
    for i in range(len(charges)):
        if charges[i]!=0:
            chargeList.append(charges[i])
            chargedAminos.append(atom_COM[i])
    
    PlotBiColorCoords(chargedAminos,chargeList,  ax[3], 'Partial Charges', coordOrder=(1,0,2) )
   
    for amino in injectionAminos:
        startCoords = activeAminos[amino['aminoIndex']]['centerOfMass']
        for sax in ax:
            sax.plot (startCoords [ 1], startCoords[0], startCoords[ 2],'o',color='r' )
        
    for amino in exitAminos:
        endCoords = activeAminos[amino['aminoIndex']]['centerOfMass']
        for sax in ax:
            sax.plot (endCoords [ 1], endCoords[0], endCoords[ 2],'o',color='r' )
                
        
    plt.tight_layout()
    plt.suptitle('PDB 3D structure of the redox cofactors')
    plt.show()   
    
  
    
def PlotGraphRates(reorg, gammas, distanceRates, energyRates, transferRates, voltageRates):
    #yes it is a lot of graphs
    fig, ax = plt.subplots(2, 2, figsize=(10, 10))
    ax = np.ravel(ax)
    v, b = np.histogram(gammas, bins=100)
    # make bar plot of histogram
    ax[0].bar(b[:-1], np.log10(v+1), width=b[1]-b[0],
            color='b', align='edge', alpha=0.5)
    ax[0].set_xlabel('Gamma (A$^{-1}$)')
    ax[0].set_ylabel('log10(Count)')
    ax[0].set_title('Gamma Distribution')

    rates = np.array(distanceRates)
    r=rates[:,0]
    ax[1].plot(rates[:, 0], np.log10(rates[:, 1])-9, '.', label='Dutton')
    ax[1].plot(rates[:, 0], np.log10(rates[:, 2])-9, '.', label='Vibrate')
    ax[1].plot(rates[:, 0], np.log10(rates[:, 3])-9, '.', label='Published')
    
    ax[1].set_title('Distance Prefactor')
    ax[1].set_ylabel('Transfer Rate log$_{10}$(ns$^{-1}$)')
    ax[1].set_xlabel('Edge-Edge distance (A)')
    #ax[1].set_xlim([2, 28])
    #ax[1].set_ylim([-4, 14])
    ax[1].legend()

    rates = np.array(energyRates)
    ax[2].plot(rates[:, 0], np.log10( (rates[:, 1]*1e-9)), '.', label='Forward')
    ax[2].plot(rates[:, 0], np.log10( (rates[:, 2]*1e-9)), '.', label='Backward')
    ax[2].set_title("ΔF$\u2021$")
    ax[2].set_ylabel('Transfer Rate log$_{10}$(ns$^{-1}$)')
    ax[2].set_xlabel('ΔG$_0$ + eΔV (eV)')
    ax[2].legend()

    rates = np.array(transferRates)
   
    
    ax[3].plot(rates[:, 0], np.log10(rates[:, 1]*1e-9), '.', label='Dutton')
    ax[3].plot(rates[:, 0], np.log10(rates[:, 2]*1e-9), '.', label='Vibrate')
    ax[3].plot(rates[:, 0], np.log10(rates[:, 3]*1e-9), '.', label='Published')
    
    dG=.85- 1.08
    ax[3].plot(r, 15-.6*r-3.1*(reorg)**2/reorg-9,'-',alpha=.25,label='Dutton Y-Y') 
    ax[3].plot(r, 15-.6*r-3.1*(dG+reorg)**2/reorg-9,'-',alpha=.25,label='Dutton Y-W') 
    ax[3].plot(r, 15-.6*r-3.1*(reorg-dG)**2/reorg-9,'-',alpha=.25,label='Dutton W-Y') 
    
    ax[3].set_title('K$_{Hopping}$ for distance, redox energies, and voltages')
    ax[3].set_ylabel('Transfer Rate log$_{10}$(ns$^{-1}$)')
    ax[3].set_xlabel('Edge-Edge distance (A)')
    ax[3].legend()

    plt.show()
    
def PlotVoltageRates(voltageRates):    
    fig,ax = plt.subplots(1,2,figsize=(10,5))

    lvoltageRates=( np.array(voltageRates))
    v,b = np.histogram(np.log10(lvoltageRates[:,1]), bins=100)
    ax[0].plot(b[:-1], np.log10(v), '.', label='Forward')
    v,b = np.histogram(np.log10(lvoltageRates[:,2]), bins=100)
    ax[0].plot(b[:-1], np.log10(v), '.', label='Backward')
    ax[0].set_xlabel('Voltage Rate')
    ax[0].set_ylabel('log10(Count)')

    ax[1].plot(lvoltageRates[:,3], np.log10(lvoltageRates[:,4]), '.', label='Forward')
    ax[1].plot(lvoltageRates[:,3], np.log10(lvoltageRates[:,5]), '.', label='Back')
    ax[1].set_xlabel('Distance (A)')
    ax[1].set_ylabel('log10(Transfer Rate)')

    plt.suptitle('Changes in Transfer Rate due to Voltage Differences')
    plt.legend()
    plt.show()          
    
    

    
def PlotPaths(activeAminos,shortestPaths,ax,injectionAminos,exitAminos):
    visited = np.zeros(len(activeAminos))
    lastTime =0 
    for path in shortestPaths:
        if lastTime!=path['time']:
            lastTime = path['time']
            for node in path['path']:
                if node <len(visited):
                    visited[node] += 1
    
    
     
    PlotColorCoords(activeAminos, visited, 2, ax, 'Paths', coordOrder=(1,0,2), invertColors=False)
    
    cc=0
    lastTime =0 
    for i  in  range(0,len( shortestPaths ),1+int(np.floor(len( shortestPaths )/100))):
        path=shortestPaths[ i]
        lastTime = path['time']
        newCoords = []
        for index in path['path']:
            if index in activeAminos:
                coords=(activeAminos[index]['centerOfMass'])
                newCoords.append(coords)
        #plot out the paths
        if len(newCoords)>0:
            newCoords = np.array(newCoords)
            ax.plot (newCoords[:, 1], newCoords[:, 0], newCoords[:, 2],'-',color=(cc/100.0,0,1-cc/100.0),alpha= (1-cc/100.0)*.5+.25 )
        cc+=1
        if cc>99:
            break
        
    for i  in  range(3 ):
        path=shortestPaths[ i]
        lastTime = path['time']
        newCoords = []
        for index in path['path']:
            if index in activeAminos:
                coords=(activeAminos[index]['centerOfMass'])
                newCoords.append(coords)
        #plot out the paths
        if len(newCoords)>0:
            newCoords = np.array(newCoords)
            ax.plot (newCoords[:, 1], newCoords[:, 0], newCoords[:, 2],'-',color='g' ) 

    for amino in injectionAminos:
        startCoords = activeAminos[amino['aminoIndex']]['centerOfMass']
        ax.plot (startCoords [ 1], startCoords[0], startCoords[ 2],'o',color='r' )
        
    for amino in exitAminos:
        endCoords = activeAminos[amino['aminoIndex']]['centerOfMass']
        ax.plot (endCoords [ 1], endCoords[0], endCoords[ 2],'o',color='r' )
            

   
def PlotCriticalNode(G,shortest_Paths, axis,label):
    #collect all the lengths and determine what the common longest jump is
    totalLength = []
    distances  = []
    minMaxDist = []
    for path in shortest_Paths:   
        maxDist = 0
        
        dists = [] 
        for i in range(len( path['path'])-1):
            dist = G[path['path'][i]][path['path'][i+1]]['time'] 
            dists.append(dist)
            if dist > maxDist:
                maxDist = dist
        minMaxDist.append(maxDist)
        distances.extend(dists)
        totalLength.append(np.sum(dists))
                
    v,b,_=axis[0].hist(np.log10( totalLength), bins=10, label=label)
    mostProb =10** b[np.argmax(v)]
    axis[0].set_title('Total Travel Time')
    axis[0].set_xlabel('Length (ns)')
    axis[0].set_ylabel('Count')
    axis[0].legend()


    axis[1].hist(np.log10(distances), bins=20, label=label)
    
    axis[1].set_title('Inidividual Jump Time')
    axis[1].set_xlabel('Time log$_{10}$(ns)')
    axis[1].set_ylabel('Count')
    axis[1].legend()

    axis[2].hist(np.log10(minMaxDist), bins=20, label=label)
    axis[2].set_title('Max jump time in each path')
    axis[2].set_xlabel('Time log$_{10}$(ns)')
    axis[2].set_ylabel('Count')
    axis[2].legend()

    print('Network : ', label)
    print('Total paths tested: ', len(shortest_Paths))
    print(f'The average travel time is { mostProb:.2e} ns')
    print(f'The average rate is {1/(1e-9* mostProb):.2e} 1/s')
    #print(f'All paths must make at least one jump >= { np.min(minMaxDist):.2e} ns jump giving a max current of { 1.602e-19/np.min(minMaxDist):.2e} nA (assuming 1 electron in the molecule at a time)')
    print('\n\n')
         
 
 
  
        
def PlotDistancePaths(G_connected, shortestPaths_connected, weight):
    #collect all the lengths and determine what the common longest jump is
    totalLength = []
    distances  = []
    minMaxDist = []
    
    maxPathStartIndex=0
    maxPathEndIndex=0
    for path in shortestPaths_connected:   
        maxDist = 0
        
        dists = [] 
        for i in range(len( path['path'])-1):
            dist = G_connected[path['path'][i]][path['path'][i+1]][weight]  
            dists.append(dist)
            if dist > maxDist:
                maxDist = dist
                maxPathEndIndex = path['path'][i+1]
                maxPathStartIndex = path['path'][i]
                
        minMaxDist.append(maxDist)
        distances.extend(dists)
        totalLength.append(np.sum(dists))
            
    print(G_connected.nodes[maxPathStartIndex]['activeAmino'], G_connected.nodes[maxPathEndIndex]['activeAmino'])
    _,ax = plt.subplots(1,3,figsize=(8,3))
            
    ax[0].hist(totalLength, bins=50)
    ax[0].set_title('Total Path length')
    ax[0].set_xlabel('Length (nm)')
    ax[0].set_ylabel('Count')


    ax[1].hist(distances, bins=20)
    ax[1].set_title('Length of jumps')
    ax[1].set_xlabel('length (nm)')
    ax[1].set_ylabel('Count')

    ax[2].hist(minMaxDist, bins=20)
    ax[2].set_title('Max jump in each path')
    ax[2].set_xlabel('length (nm)')
    ax[2].set_ylabel('Count')


    print('Total paths tested: ', len(shortestPaths_connected))
    print(f'The average tunnel distance is { np.mean(distances):.2f} nm')
    print(f'All paths must make at least one jump >= { np.min(minMaxDist):.2f} nm jump')
    plt.tight_layout( h_pad=2)
    plt.show()
     
    
def DrawGraph(G, layout, subPlot, title,labelLongest, injectionAminos, exitAminos):
    plt.subplot(subPlot)
    pos=layout(G)
    nx.draw_networkx_nodes(G, pos, node_size=20)
    nx.draw_networkx_edges(G, pos, alpha=0.4)
    
     
    nx.draw_networkx_nodes(G, pos, nodelist=[amino['aminoIndex'] for amino in injectionAminos] , node_color='r', node_size=30)
    nx.draw_networkx_nodes(G, pos, nodelist=[amino['aminoIndex'] for amino in exitAminos], node_color='g', node_size=30)
    
    if labelLongest:
        #find the three longest edges and label them
        edgeList = list(G.edges(data='time'))
        edgeList.sort(key=lambda x: x[2], reverse=True)
        
        edges ={}
        for edge in edgeList:
            edges[(edge[0],edge[1])]=f"{edge[2]:.2e}"
        
        nx.draw_networkx_edge_labels(
            G, pos,
            edge_labels=edges,
            font_color='red'
        )

    plt.title(title)    
    
def DrawLayout(Gs, layout,   title, injectionAminos, exitAminos,labelLongest=False):
    plt.figure(figsize=(5*len(Gs),5))
    for i in range(len(Gs)):
        DrawGraph(Gs[i][1], layout, 100+len(Gs)*10+ i+1, Gs[i][0], labelLongest, injectionAminos, exitAminos)
    plt.suptitle(title)
    plt.show()     
    
  
def PlotNetworkMetrics(G,activeAminos):
    clusters = list( nx.clustering(G ).values())
    degree = G.degree() 
    degree = [degree[node] for node in G.nodes()]

    centrality =list(  nx.eigenvector_centrality(G, weight='time',tol=.1).values())
    #degCent = nx.degree_centrality(G_vibrate)
    betCent =list( nx.betweenness_centrality(G, weight='time', normalized=True, endpoints=True).values())

    popularity = np.zeros(len(activeAminos))
    for clique in nx.find_cliques(G):
        for node in clique:
            popularity[node] += 1

    #create three 3d  subplots with the x y and z views
    fig, ax = plt.subplots(2, 2, figsize =(8,8), subplot_kw={'projection': '3d'})
    ax = np.ravel(ax)

    PlotColorCoords(activeAminos, degree, 2, ax[0], 'Degree', coordOrder=(1,0,2))
    PlotColorCoords(activeAminos, centrality, 1, ax[1], 'Centrality', coordOrder=(1,0,2))
    PlotColorCoords(activeAminos, betCent, 0, ax[2], 'Choke Points', coordOrder=(1,0,2))
    PlotColorCoords(activeAminos, popularity,0 , ax[3], 'Cliques', coordOrder=(1,0,2))
    
    plt.show()    
    return popularity 


def CompareWithPublished(activeAminos,G_min,transferRates):

    rates =[x.split("\t") for x in """43	48	8.481	8.481	4563.516	0.51
36	55	1.702	1.702	4661.247	0.55
35	54	0.023	24.311	975.706	0.59
42	54	6.165	6.165	13.142	0.65
43	55	3.502	3.502	559.671	0.69
54	55	4.293	4.293	269.791	0.86
55	76	1.606	1.606	26.667	0.94
48	55	4.176	4.176	44.744	1
55	69	9.752	0.009	4.441	1.01
36	43	3.072	3.072	1.42	1.02
35	55	0.015	15.964	41.042	1.03
42	55	5.081	5.081	11.552	1.03
55	88	0.845	0.845	1.888	1.04
43	77	4.153	4.153	28.054	1.05
42	43	6.31	6.31	327.05	1.06
36	69	11.08	0.01	0.692	1.06
43	70	6.422	6.422	10.354	1.06
48	77	5.376	5.376	2.989	1.06
42	48	8.476	8.476	45.609	1.07
35	36	0.024	25.621	325.759	1.08
35	42	0.007	7.848	0.441	1.09
55	70	0.918	0.918	0.447	1.09
48	76	6.703	6.703	8.679	1.1
43	76	3.963	3.963	5.817	1.14
36	54	6.589	6.589	26.219	1.15
36	70	7.127	7.127	11.684	1.16
43	54	3.125	3.125	2.025	1.25
43	89	3.208	3.208	10.422	1.25
55	89	0.089	0.089	0.33	1.25
36	42	3.183	3.183	0.628	1.27
48	54	1.095	1.095	0.18	1.28
36	88	1.104	1.104	0.237	1.37
55	77	2.653	2.653	1.108	1.39
43	88	0.469	0.469	0.099	1.42
48	70	0.104	0.104	0.011	1.42
35	43	0.0004	0.327	0.012	1.61
""".split('\n') if x!='']

    nodeDict = {}
    cc=0
    for aa in activeAminos:
        nodeDict[aa['index']] = cc
        cc+=1
        
    transferRates = np.array(transferRates)
    pubRates = np.array( [[float(x[-1])*10, float(x[-3]),float(x[-4]), int(x[0]), int(x[1])] for x in rates])
    _,ax=plt.subplots(1,3,figsize=(12.5,4))

    plt.suptitle('Transfer Rates for min model compared with published rates')
    ax[0].semilogy (pubRates[:,0],pubRates[:,1],'.r', label='Published Forward')
    ax[0].semilogy(pubRates[:,0],pubRates[:,2],'.r', label='Published Back')

    calcRates=[]
    compares =""
    for pubrate in pubRates:
        
        fNode = nodeDict[int( pubrate[-1])]
        tNode = nodeDict[int( pubrate[-2])]
        forwardRate = G_min[ fNode ] [ tNode]['rate']
        backRate = G_min[ tNode ] [ fNode]['rate']
        brate =np.max( [ forwardRate ,backRate])
        drate= np.max([     pubrate[1],pubrate[2]])
        
        compares += f'{int( pubrate[-1])}\t{int( pubrate[-2])}\t{brate}\t{drate}\t{np.abs(brate-drate) }\n'
        calcRates .append( [ pubrate[0],brate ,forwardRate,backRate, pubrate[1],pubrate[2]])
    calcRates = np.array(calcRates)
    #print(compares)
    ax[0].semilogy(calcRates[:,0],calcRates[:,2],'.k', label='CG Forward')
    ax[0].semilogy(calcRates[:,0],calcRates[:,3],'.k', label='CG Back')

    ax[0].set_xlabel('Distance (A)')
    ax[0].set_ylabel('Rate (1/ns)')
    ax[0].legend()
    
    # ax[1].plot(calcRates[:,2],calcRates[:,-2],'.', label='Pub. Forward')
    # ax[1].plot(calcRates[:,3],calcRates[:,-1],'.', label='Pub. Back')

    # ax[1].plot([.002,20 ],[.002,20],  label='Equal')

    # ax[1].set_title('Linear Rates')
    # ax[1].set_xlabel('Coarse Grain Rate (1/ns)')
    # ax[1].set_ylabel('Published Rate (1/ns)')

    ax[1].loglog(calcRates[:,2],calcRates[:,-2], '.', label='Forward')
    ax[1].loglog(calcRates[:,3],calcRates[:,-1],'.',  label='Back')
    ax[1].loglog([.002,20 ],[.002,20],  label='Equal')

    ax[1].set_title('Log Log Rates')
    ax[1].set_xlabel('Coarse Grain Rate (1/ns)')
    ax[1].set_ylabel('Published Rate (1/ns)')
    ax[1].legend()


    ax[2].set_title('Residuals of rate comparison')
    ax[2].hist((calcRates[:,-2]-calcRates[:,2]),bins=20, alpha=.3, label='Forward')
    ax[2].hist((calcRates[:,-1]-calcRates[:,3]),bins=20, alpha=.3, label='Back')
    ax[2].set_xlabel('Residual |Pub-CG|(1/ns)')
    ax[2].set_ylabel('Count')

    plt.legend()
    plt.show()

    print(f"MAE forward:{np.mean(np.abs((calcRates[:,-2]-calcRates[:,2]))):.2e} 1/ns")
    print(f"MAE back:{np.mean(np.abs((calcRates[:,-1]-calcRates[:,3]))):.2e} 1/ns")
    print(f"MAE all:{np.mean( np.concatenate([np.abs(calcRates[:,-2]-calcRates[:,2]) ,np.abs((calcRates[:,-1]-calcRates[:,3]))])):.2e} 1/ns")    
    
    
def PlotComparison():
    attemptFrequency = 2*np.pi/hbar_eV/np.sqrt( 4*np.pi*reorgE_EV*kbT_eV  )   # attempt prefactor 2 pi/ hbar
    beta = 1/(kbT_eV)  # 1/kT in 1/eV

    rates =[x.split("\t") for x in """35	36	0.024	25.621	325.759	1.08	4.51
    35	42	0.007	7.848	0.441	1.09	9.097
    35	43	0	0.327	0.012	1.61	11.542
    35	54	0.023	24.311	975.706	0.59	3.437
    35	55	0.015	15.964	41.042	1.03	5.45
    36	42	3.183	3.183	0.628	1.27	9.529
    36	43	3.072	3.072	1.42	1.02	8.49
    36	54	6.589	6.589	26.219	1.15	6.674
    36	55	1.702	1.702	4661.247	0.55	2.421
    42	43	6.31	6.31	327.05	1.06	4.85
    42	48	8.476	8.476	45.609	1.07	6.683
    42	54	6.165	6.165	13.142	0.65	7.124
    42	55	5.081	5.081	11.552	1.03	7.069
    43	48	8.481	8.481	4563.516	0.51	3.419
    43	54	3.125	3.125	2.025	1.25	8.177
    43	55	3.502	3.502	559.671	0.69	4.136
    48	54	1.095	1.095	0.18	1.28	10.181
    48	55	4.176	4.176	44.744	1	5.98
    54	55	4.293	4.293	269.791	0.86	4.737
    36	69	11.08	0.01	0.692	1.06	8.594
    36	70	7.127	7.127	11.684	1.16	7.384
    36	88	1.104	1.104	0.237	1.37	9.799
    43	70	6.422	6.422	10.354	1.06	7.35
    43	76	3.963	3.963	5.817	1.14	7.441
    43	77	4.153	4.153	28.054	1.05	6.308
    43	88	0.469	0.469	0.099	1.42	10.328
    43	89	3.208	3.208	10.422	1.25	6.903
    48	70	0.104	0.104	0.011	1.42	12.525
    48	76	6.703	6.703	8.679	1.1	7.543
    48	77	5.376	5.376	2.989	1.06	8.233
    55	69	9.752	0.009	4.441	1.01	6.901
    55	70	0.918	0.918	0.447	1.09	9.022
    55	76	1.606	1.606	26.667	0.94	6.001
    55	77	2.653	2.653	1.108	1.39	8.64
    55	88	0.845	0.845	1.888	1.04	7.775
    55	89	0.089	0.089	0.33	1.25	8.734
    """.split('\n') if x!='']
    #Donor	Acceptor	Forward	Backward	g	R_COM	r_Min
    nodeDict = {}
    cc=0
    for aa in activeAminos:
        nodeDict[aa['index']] = cc
        cc+=1
    rmin=[]
    rateForward =[]
    rateBack= []
    aindex=[]
    aindex2=[]
    erf=[]
    erb=[]

    def CalcRate(r,reorgE_EV):
        vibrate = np.mean([donor['std'], acceptor['std']]) 
        
        gammaV= -1.4
        distance_rate_vibrate =attemptFrequency * np.exp(gammaV * r)
        preVibrate = 2.7/7* np.exp(gammaV/2 * r)
        distance_rate_vibrate=preVibrate**2* np.exp(gammaV*gammaV*vibrate*vibrate/2)  
        preG= (.01e-9)* np.pi*preVibrate**2/hbar_eV *np.exp( 3/2*gammaV*gammaV*vibrate*vibrate) 
        
        dE = acceptor['redoxEnergy_EV']-  donor['redoxEnergy_EV'] 

        dF_F = (reorgE_EV + dE)**2/(4*reorgE_EV)
        
        g_F=preG/ np.sqrt(reorgE_EV* dF_F)
        
        energy_rate_forward = np.exp(-beta * dF_F)
        
        #get the backard rate for returning to acceptor
        dF_B = (reorgE_EV - dE )**2/(4*reorgE_EV)
        g_B=preG/ np.sqrt(reorgE_EV* dF_B)
        
        energy_rate_back = np.exp(-beta * dF_B)
        
        
        fR=distance_rate_vibrate * energy_rate_forward/(1+g_F)*1e9 *.05
        bR= distance_rate_vibrate * energy_rate_back/(1+g_B)*1e9 *.05
        return fR,bR
    erf2=[]
    erb2=[]
    for x in rates:
        #Donor	Acceptor	Forward	Backward	g	R_COM	r_Min
        index,index2,forward,back,g,r_com,r_min= int(x[0]), int(x[1]), float(x[2]), float(x[3]),float(x[4]), float(x[5]),float(x[6])
        
        acceptor = activeAminos[nodeDict[index]]    
        donor = activeAminos[nodeDict[index2]]
        
        dE = acceptor['redoxEnergy_EV']-  donor['redoxEnergy_EV'] 
        
        if dE==0:
            rmin.append(r_min)
            rateForward.append(forward)
            rateBack.append(back)
            aindex.append(nodeDict[index])
            aindex2.append(nodeDict[index2])
            
            fR,bR = CalcRate(r_min, .2)
            erf.append(fR)
            erb.append(bR)
            
            fR,bR = CalcRate(r_min, .5)
            erf2.append(fR)
            erb2.append(bR)
        
    rmin=np.array(rmin)    
    plt.semilogy(rmin, rateForward,'o', label='Published Rates')
    #plt.plot(rmin, rateBack,'o')

    plt.plot(rmin, erf,'.', label='Calculated Reorg. .2V')
    #plt.plot(rmin, erb,'.')

    plt.plot(rmin, erf2,'.', label='Calculated Reorg. .5V')
    #plt.plot(rmin, erb2,'.')


    dG=.85- 1.08
    #lt.plot(rmin,10**( 15-.6*rmin-3.1*(reorgE_EV)**2/reorgE_EV-9),'-',alpha=.25,label='Dutton Y-Y') 
    #plt.plot(rmin,10**( 15-.6*rmin-3.1*(dG+reorgE_EV)**2/reorgE_EV-9),'-',alpha=.25,label='Dutton Y-W') 
    #plt.plot(rmin,10**( 15-.6*rmin-3.1*(reorgE_EV-dG)**2/reorgE_EV-9),'-',alpha=.25,label='Dutton W-Y') 
    plt.ylim(0,20)
    plt.ylabel('Rate (1/ns)')     
    plt.xlabel('Edge-Edge Distance (nm)')
    plt.legend()
    plt.show()    
    
    
    
def PlotTransitTimes(axTime, axRate,axCurrent, shortestPaths, label):
    times = np.array([x['time'] for x in shortestPaths if x['time'] >0 ])
    occupancy = np.array([x['occupancy'] for x in shortestPaths if x['time'] >0 ])
    
    times_sum = np.log10(  times  )

    #ACS Omega paper page G
    rates_sum  = np.log10(1.0/(times *1e-9) )
    v,b,_=axTime.hist(times_sum  ,   label = label )
    mostProb =10** b[np.argmax(v)]
    
    axTime.set_xlabel('log10(Transit Time(ns))')
    axTime.legend()
    
    axRate.hist(rates_sum ,     label = label)
    axRate.set_xlabel('Log10(K$_{ET}$ ($s^{-1}$))')
    axRate.legend()  
    
    current = np.log10( np.abs( 1.602e-19*occupancy/(times*1e-9)))+9
    axCurrent.hist(current ,     label = label)
    axCurrent.set_xlabel('Current Log10(nA)')
    axCurrent.legend()
    
    
    print('Network : ', label)
    print('Total paths tested: ', len(shortestPaths))
    print(f'The average travel time is { mostProb:.2e} ns')
    print(f'The molecule K(ET) is {1/(1e-9* mostProb):.2e} 1/s')
    #print(f'All paths must make at least one jump >= { np.min(minMaxDist):.2e} ns jump giving a max current of { 1.602e-19/np.min(minMaxDist):.2e} nA (assuming 1 electron in the molecule at a time)')
    print('\n\n')
        
                
def PlotShortedPaths(models):
    fig = plt.figure(figsize =(16.5,4))

    axPath = fig.add_subplot(1, 4, 1, projection='3d')
    axRate=  fig.add_subplot(1, 4, 2)
    axTime= fig.add_subplot(1, 4, 3)
    axCurrent= fig.add_subplot(1, 4, 4)
    #KMC_Plotting. PlotPaths(activeAminos,shortestPaths_static,axPath,injectionAminos,exitAminos)
    for model in models:
        PlotTransitTimes(axTime, axRate, axCurrent, model['ShortestPaths'], model['name'])

    axTime.set_title('Total Molecule Transit Time')
    axRate.set_title('Total Molecule Rate')

    plt.tight_layout(h_pad=4)
    plt.show()        
    

def PlotPathMetrics(G, shortestPaths, title): 

    nodeAminos = {}
    activeAminos = []
    for node in G.nodes():
        if 'activeAmino' in  G.nodes[node]:
            nodeAminos[node] = {'outRate':1/G.nodes[node]['outrate'], "aa": G.nodes[node]['activeAmino'], 'visited':0,
                                  'embedding': G.nodes[node]['activeAmino']['embedding'], 'redoxEnergy_EV': G.nodes[node]['activeAmino']['redoxEnergy_EV']}
            activeAminos.append(G.nodes[node]['activeAmino'])
                
    lastTime =0 
    for path in shortestPaths:
        if lastTime!=path['time']:
            lastTime = path['time']
            for node in path['path']:
                if node in nodeAminos:
                    nodeAminos[node]['visited'] += 1
                
    
    redoxPoints = [nodeAminos[ node]['redoxEnergy_EV'] for node in nodeAminos]  
    minRP = np.min(redoxPoints)
    maxRP = np.max(redoxPoints)

    X_transformed = np.array( [nodeAminos[ node]['embedding'] for node in nodeAminos]    )

    #get max and min of transformed coordinates
    maxX,maxY = np.max(X_transformed[:,0]) ,np.max(X_transformed[:,1]) 
    minX,minY = np.min(X_transformed[:,0]), np.min(X_transformed[:,1])
    xdim=100
    ydim=100
    lX=(xdim-1)/( maxX-minX)
    lY= (ydim-1)/( maxY-minY)
    visitStruct = np.zeros((xdim,ydim),dtype=np.float32)
    dwellStruct = np.zeros((xdim,ydim),dtype=np.float32) 
    
    visited = []
    outTimes = []   
    for node in nodeAminos:
        x,y= nodeAminos[node]['embedding']
        visit = nodeAminos[node]['visited']
        dwell = nodeAminos[node]['outRate']
        x = int((x-minX)*lX)
        y = int((y-minY)*lY)
        if visit>0:
            visitStruct[y,x] =  ( visit)**.5
            visited.append(visit)
        if dwell>0:
            dwellStruct[y,x] = np.log  (dwell)
            outTimes.append(dwell)
        
        
    print('Smoothing Energy Map')
    sigma2=3
    visitStruct = np.fft.ifftn(fourier_gaussian(np.fft.fftn(visitStruct), sigma=sigma2)).real 
    visitStruct = visitStruct/np.max(visitStruct)
    
    dwellStruct = np.fft.ifftn(fourier_gaussian(np.fft.fftn(dwellStruct), sigma=3)).real 
    dwellStruct = dwellStruct /np.max(dwellStruct)

    

    fig, ax = plt.subplots(1,2,figsize=(8,4))
    ax=np.ravel(ax)
    plt.suptitle(title)
    
                
    plotDensity(visitStruct,fig, ax[0], "Visited", 'Count$^{1/2}$', [minX,maxX,maxY,minY],redoxPoints,activeAminos)
    plotDensity(dwellStruct,fig, ax[1], "Rest Time", 'Log(ns)',[minX,maxX,maxY,minY],redoxPoints,activeAminos)   
            
    for i  in  range(2 ):
            path=shortestPaths[ i]
            lastTime = path['time']
            newCoords = []
            for index in path['path']:
                if 'activeAmino' in  G.nodes[index]:
                    coords=G.nodes [index]['activeAmino']['embedding']
                    newCoords.append(coords)
            #plot out the paths
            newCoords = np.array(newCoords)
            ax[0].plot (newCoords[:, 0], newCoords[:, 1] ,'-',color='g' ,alpha=.5)          
            ax[1].plot (newCoords[:, 0], newCoords[:, 1] ,'-',color='g'  ,alpha=.5)          
            
    plt.tight_layout()
    plt.show()                   