from crane import crane
from fem import FEM
import numpy as np
import matplotlib.pyplot as plt

if __name__ == "__main__":
    h_b = 5e-2  # Höhe des Querschnitts der Balken in m
    b_b = 5e-2  # Breite des Querschnitts der Balken in m
    E = 210e9  # E-Modul in Pa
    A = h_b * b_b  # Querschnittsfläche der Balken in m^2
    rho = 7850  # Dichte in kg/m^3
    g = 9.81  # m/s^2
    myCrane = crane(10, 10, 1, 0.5, A, rho, E)
    nodes = myCrane.nodes
    bars = myCrane.bars
    # P = np.zeros_like(nodes)
    # P[16, 0] = 1
    # P[17, 0] = 1
    # P[18, 0] = 1
    # P[19, 0] = 1

    # # Lager Verschiebung
    # Ur = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]  # 4 Festlager = 4*3 blockierte Freiheitsgrade

    # # Freiheitsgrade (1 = beweglich, 0 = fest) # evtl. booleans benutzen?
    # DOFCON = np.ones_like(nodes).astype(int)
    # # Festlager
    # DOFCON[0, :] = 0
    # DOFCON[1, :] = 0
    # DOFCON[2, :] = 0
    # DOFCON[3, :] = 0
    fem = FEM(myCrane)
    N, R, U = fem.TrussAnalysis()
    print('\nAxial Forces (positive = tension, negative = compression)')
    # Anschaulichkeit
    print(N[np.newaxis].T)
    print('\nReaction Forces (positive = upward, negative = downward')
    print(R)
    print('\nDeformation at nodes')
    print(U)
    fem.Plot(nodes, bars, 'gray', '--', 1, 'Undeformed')
    scale = 1
    # Berechne die neue Position der Knoten
    Dnodes = U * scale + nodes
    fem.Plot(Dnodes, bars, 'red', '-', 2, 'Deformed')
    fem.plotPoint(Dnodes[-1])
    fem.plotPoint(Dnodes[-2])
    plt.show()
    # plt.savefig('fig-1.png', dpi=300)
