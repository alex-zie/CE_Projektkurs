from crane import crane_1
from crane import crane_2_1
from crane import crane_2_2
from fem import FEM
import numpy as np
import matplotlib.pyplot as plt

if __name__ == "__main__":
    tower_height = 10 # Maximalhöhe des Krans [m]
    jib_length = 10 # Maximallänge des Auslegers [m]
    length_segments = 1 # Länge der Segmente (nicht der Stäbe) [m]
    h_b = 10e-2  # Höhe des Querschnitts der Balken [m]
    b_b = 10e-2  # Breite des Querschnitts der Balken [m]
    A = h_b * b_b  # Querschnittsfläche der Balken in m^2
    E = 210e9  # E-Modul in Pa
    rho = 7850  # Dichte in kg/m^3
    load = 0 # angebrachte Last in kg

    # myCrane = crane_1(tower_height, jib_length, length_segments, A, rho, E)
    # myCrane = crane_2_1(10, 10, 1, A, rho, E)
    myCrane = crane_2_2(10, 10, 1, A, rho, E)

    nodes = myCrane.nodes
    bars = myCrane.bars

    # Gewicht
    for i in myCrane.tip_nodes:
       myCrane.addExternalForce(i, 0, 0, -load/len(myCrane.tip_nodes))

    # Gegengewicht
    for i in myCrane.counterweight_nodes:
       myCrane.addExternalForce(i, 0, 0, load/len(myCrane.counterweight_nodes))
    
    fem = FEM(myCrane, True)

    #fem.addWind(28, 0, -1) #TODO maybe pass arrays to wind
    
    # Visualisierung
    fem.display(tension=True)

    # Windangriffsfläche hervorheben
    # plot_sides_x = False
    # plot_sides_y = False
    # if plot_sides_x:
    #     fem.Plot(Dnodes, bars[myCrane.x_negative_side], 'yellow', '-', 2, 'neg. x')
    #     fem.Plot(Dnodes, bars[myCrane.x_positive_side], 'orange', '-', 2, 'pos. x')
    #     pass
    # if plot_sides_y:
    #     fem.Plot(Dnodes, bars[myCrane.y_negative_side], 'green', '-', 2, 'neg. y')
    #     fem.Plot(Dnodes, bars[myCrane.y_positive_side], 'cyan', '-', 2, 'pos. y')
    #     pass

    
