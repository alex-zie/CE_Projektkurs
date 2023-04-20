import numpy as np
import matplotlib.pyplot as plt

E = 1e4     # E-Modul
A = 0.25   # Querschnittsfläche der Balken

nodes = []
bars = []

height = 10
numSegments = 10
hs = height/numSegments # Höhe jedes Segments

for i in range(numSegments):
    nodes.append([0, 0, i*hs])         
    nodes.append([hs, 0, i*hs])     
    nodes.append([0, hs, i*hs])       
    nodes.append([hs, hs, i*hs])  


# x- und y-Richtung
for i in range(numSegments):
    bars.append([4*i, 4*i+1])
    bars.append([4*i+2, 4*i+3])
    bars.append([4*i, 4*i+2])
    bars.append([4*i+1, 4*i+3])

# z-Richtung
for i in range(numSegments-1):
    bars.append([4*i, 4*i+4])
    bars.append([4*i+1, 4*i+5])
    bars.append([4*i+2, 4*i+6])
    bars.append([4*i+3, 4*i+7])

# Kreuzstreben
for i in range(numSegments-1):
    bars.append([4*i, 4*i+5])
    bars.append([4*i+1, 4*i+4])
    bars.append([4*i+2, 4*i+7])
    bars.append([4*i+3, 4*i+6])
    bars.append([4*i+1, 4*i+7])
    bars.append([4*i+3, 4*i+5])
    bars.append([4*i+2, 4*i+4])
    bars.append([4*i, 4*i+6])
    


nodes = np.array(nodes).astype(float)
bars = np.array(bars)

# Äußere Kräfte
P = np.zeros_like(nodes)
#P[7,2] = 50

# Lager Verschiebung
Ur = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0] # 4 Festlager = 4*3 blockierte Freiheitsgrade

# Freiheitsgrade (1 = beweglich, 0 = fest) # evtl. booleans benutzen?
DOFCON = np.ones_like(nodes).astype(int)
# Festlager
DOFCON[0,:] = 0
DOFCON[1,:] = 0
DOFCON[2,:] = 0
DOFCON[3,:] = 0

#%% Truss structural analysis
def TrussAnalysis():
    NN = len(nodes)
    NE = len(bars)
    DOF = 3 # weil wir uns in 3D befinden
    NDOF = DOF*NN
    
    # Geometrie
    d = nodes[bars[:,1],:] - nodes[bars[:,0],:] # Endknoten - Anfangsknoten
    L = np.sqrt((d**2).sum(axis=1))             # Länge der Balken (Euklidische Norm)
    orientations = d.T/L                        # Richtungsvektoren der Balken
    trans = np.concatenate((-orientations.T,orientations.T), axis=1) # Transformationsvektor

    # Steifigkeitsmatrix
    K = np.zeros([NDOF,NDOF])
    for k in range(NE):
        aux  = 3*bars[k,:]
        index = np.r_[aux[0]:aux[0]+3,aux[1]:aux[1]+3]
        ES = np.dot(trans[k][np.newaxis].T*E*A, trans[k][np.newaxis])/L[k] # lokale Steifigkeiten
        K[np.ix_(index,index)] = K[np.ix_(index,index)] + ES
    
    freeDOF = DOFCON.flatten().nonzero()[0]
    supportDOF = (DOFCON.flatten() == 0).nonzero()[0]
    Kff = K[np.ix_(freeDOF,freeDOF)]
    Kfr = K[np.ix_(freeDOF,supportDOF)]
    Krf = Kfr.T
    Krr = K[np.ix_(supportDOF,supportDOF)] # für die Lagerkräfte
    Pf = P.flatten()[freeDOF]
    #Uf = np.linalg.solve(Kff,Pf) 
    Uf = np.linalg.lstsq(Kff,Pf)[0] # Deformation an jedem Freiheitsgrad # least squares damit auch überbestimmte Systeme fkt.
    U = DOFCON.astype(float).flatten()
    U[freeDOF] = Uf
    U[supportDOF] = Ur
    U = U.reshape(NN,DOF)
    u = np.concatenate((U[bars[:,0]],U[bars[:,1]]),axis=1)
    N = E*A/L[:]*(trans[:]*u[:]).sum(axis=1) # externe Kräfte
    R = (Krf[:]*Uf).sum(axis=1) + (Krr[:]*Ur).sum(axis=1) # Reaktionskräfte
    R = R.reshape(4,DOF) # 4 ist die Anzahl der Lager
    return np.array(N), np.array(R), U

def Plot(nodes,c,lt,lw,lg):
    plt.subplot(projection='3d')
    plt.gca().set_aspect('equal')
    #plt.gca(projection='3d')
    for i in range(len(bars)):
        xi, xf = nodes[bars[i,0],0], nodes[bars[i,1],0]
        yi, yf = nodes[bars[i,0],1], nodes[bars[i,1],1]
        zi, zf = nodes[bars[i,0],2], nodes[bars[i,1],2]
        line, = plt.plot([xi, xf],[yi, yf],[zi, zf],color=c,linestyle=lt,linewidth=lw)
    line.set_label(lg)
    plt.legend(prop={'size': 7})

#%% Result
N, R, U = TrussAnalysis()
print('Axial Forces (positive = tension, negative = compression)')
print(N[np.newaxis].T)
print('Reaction Forces (positive = upward, negative = downward')
print(R)
print('Deformation at nodes')
print(U)
Plot(nodes,'gray','--',1,'Undeformed')
scale = 5
Dnodes = U*scale + nodes
Plot(Dnodes,'red','-',2,'Deformed')
plt.show()
#plt.savefig('fig-1.png', dpi=300)
















