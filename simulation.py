from crane import crane
from fem import FEM
import numpy as np
import matplotlib.pyplot as plt

if __name__ == "__main__":
    h_b = 10e-2  # Höhe des Querschnitts der Balken in m
    b_b = 10e-2  # Breite des Querschnitts der Balken in m
    E = 210e9  # E-Modul in Pa
    A = h_b * b_b  # Querschnittsfläche der Balken in m^2
    rho = 7850  # Dichte in kg/m^3

    myCrane = crane(1, 10, 5, 1, 1, h_b * b_b, rho, E)

    nodes = myCrane.nodes
    bars = myCrane.bars

    points = []  # indices of points where an external force is applied
    #for i in range(-1, -5, -1):
    #    myCrane.addExternalForce(i, 0, 0, -500e3 / 4)
    #    points.append(i)
    fem = FEM(myCrane, True)
    #fem.addWind(28, 0, -1)
    #print(len(fem.N[np.where(fem.N < 0)[0]]))
    #print(len(fem.F_krit()[np.where(fem.N < 0)[0]] > fem.N[np.where(fem.N < 0)[0]]))
    #print(fem.check_bending_force())
    # Veranschauung
    print('\nAxial Forces (positive = tension, negative = compression)')
    # print(fem.N)
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
    #colors = fem.paintBars(bars)
    Dnodes = fem.U * scale + nodes
    fem.Plot(Dnodes, bars, 'red', '-', 2, 'Deformed')
    plot_sides_x = True
    plot_sides_y = True
    if plot_sides_x:
        fem.Plot(Dnodes, bars[myCrane.x_negative_side], 'yellow', '-', 2, 'Selected bars')
        fem.Plot(Dnodes, bars[myCrane.x_positive_side], 'green', '-', 2, 'Selected bars')
        pass
    if plot_sides_y:
        fem.Plot(Dnodes, bars[myCrane.y_negative_side], 'yellow', '-', 2, 'Selected bars')
        fem.Plot(Dnodes, bars[myCrane.y_positive_side], 'green', '-', 2, 'Selected bars')
        pass


    # for i in points:
    #     fem.plotPoint(Dnodes[i])
    plt.show()
    # plt.savefig('fig-1.png', dpi=300)
