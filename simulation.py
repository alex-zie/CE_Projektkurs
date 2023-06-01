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

    myCrane = crane(1, 7, 7, 2, 2, 0.0225, rho, E)

    nodes = myCrane.nodes
    bars = myCrane.bars

    points = [] # indices of points where an external force is applied
    # for i in range(-1, -4, -1):
    #     myCrane.addExternalForce(i, 0, 0, -500e3)
    #     points.append(i)
    fem = FEM(myCrane, False)
    #fem.addWind(2000, 0, -1)

    # Veranschauung
    print('\nAxial Forces (positive = tension, negative = compression)')
    print(fem.N)
    # print('\nReaction Forces (positive = upward, negative = downward')
    # print(fem.R)
    # print(fem.getTension())
    # print('\nDeformation at nodes')
    # print("At node", np.where(U == U.max())[0], "is U:", U[-1:-5:-1])
    # print('\nDeformation')
    # print(fem.U)
    fem.Plot(nodes, bars, 'gray', '--', 1, 'Undeformed')
    scale = 1
    # Berechne die neue Position der Knoten
    Dnodes = fem.U * scale + nodes
    fem.Plot(Dnodes, bars, 'red', '-', 2, 'Deformed')
    # fem.Plot(Dnodes, bars[myCrane.x_side], 'green', '-', 2, 'Selected bars')
    # fem.Plot(Dnodes, bars[myCrane.y_side], 'yellow', '-', 2, 'Selected bars')
    
    for i in points:
        fem.plotPoint(Dnodes[i])

    plt.show()
    # plt.savefig('fig-1.png', dpi=300)
