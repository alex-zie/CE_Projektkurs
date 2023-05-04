import numpy as np
import matplotlib.pyplot as plt

E = 1e4     # E-Modul
A = 0.111   # Querschnittsfläche der Balken
length = 10
height = 10
ls = 1
hs = 1

# nodes = []
# bars = []

# nodes.append([-37.5,0,200])
# nodes.append([37.5,0,200])
# nodes.append([-37.5,37.5,100])
# nodes.append([37.5,37.5,100])
# nodes.append([37.5,-37.5,100])
# nodes.append([-37.5,-37.5,100])
# nodes.append([-100,100,0])
# nodes.append([100,100,0])
# nodes.append([100,-100,0])
# nodes.append([-100,-100,0])

# bars.append([0,1])
# bars.append([3,0])
# bars.append([2,1])
# bars.append([4,0])
# bars.append([5,1])
# bars.append([3,1])
# bars.append([4,1])
# bars.append([2,0])
# bars.append([5,0])
# bars.append([5,2])
# bars.append([4,3])
# bars.append([2,3])
# bars.append([5,4])
# bars.append([9,2])
# bars.append([6,5])
# bars.append([8,3])
# bars.append([7,4])
# bars.append([6,3])
# bars.append([7,2])
# bars.append([9,4])
# bars.append([8,5])
# bars.append([9,5])
# bars.append([6,2])
# bars.append([7,3])
# bars.append([8,4])

nodes = []
bars = []

# for 10m high tower with segement size hs we need at least
# height / hs segements. For now i would just ignore the last segement
# size != hs
numSegmentsA = np.ceil(length / ls).astype('int')
numSegmentsT = np.ceil(height / hs).astype('int')
for i in range(numSegmentsT):
    nodes.append([0, 0, i * hs])  # Left Top
    nodes.append([hs, 0, i * hs])  # Right Top
    nodes.append([0, hs, i * hs])  # Left Bottom
    nodes.append([hs, hs, i * hs])  # Right Bottom
offset = len(nodes)
# Erstelle die nodes des Auslegers in positive x Richtung
for i in range(1, numSegmentsA + 1):
    nodes.append([hs + i * ls, 0, (numSegmentsT - 1) * hs])  # Left Top
    nodes.append([hs + i * ls, hs, (numSegmentsT - 1) * hs])  # Right Top
    nodes.append([hs + i * ls, 0, (numSegmentsT - 2) * hs])  # Left Bottom
    nodes.append([hs + i * ls, hs, (numSegmentsT - 2) * hs])  # Right Bottom

# Turm
# x- und y-Richtung (LT für Left Top usw.)
for i in range(numSegmentsT):
    bars.append([4 * i, 4 * i + 1])  # LT -> RT
    bars.append([4 * i + 2, 4 * i + 3])  # LB -> RB
    bars.append([4 * i, 4 * i + 2])  # LT -> LB
    bars.append([4 * i + 1, 4 * i + 3])  # RT -> RB

# z-Richtung
for i in range(numSegmentsT - 1):
    bars.append([4 * i, 4 * i + 4])  # LT
    bars.append([4 * i + 1, 4 * i + 5])  # RT
    bars.append([4 * i + 2, 4 * i + 6])  # LB
    bars.append([4 * i + 3, 4 * i + 7])  # RB

# Kreuzstreben (+1 jeweils für nächste Stufe)
for i in range(numSegmentsT - 1):
    bars.append([4 * i, 4 * i + 5])  # LT -> RT+1
    # bars.append([4 * i + 1, 4 * i + 4])  # RT -> RT+1
    bars.append([4 * i + 3, 4 * i + 6])  # RB -> LB+1
    # bars.append([4 * i + 2, 4 * i + 7])  # LB -> RB+1
    bars.append([4 * i + 1, 4 * i + 7])  # RT -> RB+1
    # bars.append([4 * i + 3, 4 * i + 5])  # RB -> RT+1
    bars.append([4 * i + 2, 4 * i + 4])  # LB -> LT+1
    # bars.append([4 * i, 4 * i + 6])  # LT -> LB+1

# Ausleger
bars.append([offset - 4, offset])  # LB -> LT
bars.append([offset - 1, offset + 1])  # RB -> RT
bars.append([offset - 8, offset + 2])  # LB-1 -> LB
bars.append([offset - 5, offset + 3])  # RB-1 -> RB
# x- und y-Richtung (LT für Left Top usw.)
for i in range(numSegmentsA):
    bars.append([4 * i + offset, 4 * i + offset + 1])  # LT -> RT
    bars.append([4 * i + 2 + offset, 4 * i + 3 + offset])  # LB -> RB
    bars.append([4 * i + offset, 4 * i + 2 + offset])  # LT -> LB
    bars.append([4 * i + 1 + offset, 4 * i + 3 + offset])  # RT -> RB

# z-Richtung
for i in range(numSegmentsA - 1):
    bars.append([4 * i + offset, 4 * i + 4 + offset])  # LT
    bars.append([4 * i + 1 + offset, 4 * i + 5 + offset])  # RT
    bars.append([4 * i + 2 + offset, 4 * i + 6 + offset])  # LB
    bars.append([4 * i + 3 + offset, 4 * i + 7 + offset])  # RB

# convert python list to np.array
nodes = np.array(nodes).astype(float)
bars = np.array(bars)

# Lager
supports = np.ones_like(nodes).astype(int)

#Äußere Kräfte
P = np.zeros_like(nodes)
P[-1,2] = -1e9
P[-2,2] = -1e9
# P[-3,2] = -1e9
# P[-4,2] = -1e9


#Support Displacement
Ur = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]

#Freiheitsgrade (1 = beweglich, 0 = fest)
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
    #Uf = np.linalg.solve(Kff,Pf) # Deformation an jedem Freiheitsgrad
    Uf = np.linalg.lstsq(Kff,Pf)[0] # Deformation an jedem Freiheitsgrad
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
    for i in range(len(bars)):
        xi, xf = nodes[bars[i,0],0], nodes[bars[i,1],0]
        yi, yf = nodes[bars[i,0],1], nodes[bars[i,1],1]
        zi, zf = nodes[bars[i,0],2], nodes[bars[i,1],2]
        line, = plt.plot([xi, xf],[yi, yf],[zi, zf],color=c,linestyle=lt,linewidth=lw)
    line.set_label(lg)
    plt.legend(prop={'size': 14})

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
















