from fem import FEM
from truss import Truss
import numpy as np
import matplotlib.pyplot as plt


def originalExample():
    global myTruss

    nodes.append([-37.5, 0, 200])
    nodes.append([37.5, 0, 200])
    nodes.append([-37.5, 37.5, 100])
    nodes.append([37.5, 37.5, 100])
    nodes.append([37.5, -37.5, 100])
    nodes.append([-37.5, -37.5, 100])
    nodes.append([-100, 100, 0])
    nodes.append([100, 100, 0])
    nodes.append([100, -100, 0])
    nodes.append([-100, -100, 0])

    bars.append([0, 1])
    bars.append([3, 0])
    bars.append([2, 1])
    bars.append([4, 0])
    bars.append([5, 1])
    bars.append([3, 1])
    bars.append([4, 1])
    bars.append([2, 0])
    bars.append([5, 0])
    bars.append([5, 2])
    bars.append([4, 3])
    bars.append([2, 3])
    bars.append([5, 4])
    bars.append([9, 2])
    bars.append([6, 5])
    bars.append([8, 3])
    bars.append([7, 4])
    bars.append([6, 3])
    bars.append([7, 2])
    bars.append([9, 4])
    bars.append([8, 5])
    bars.append([9, 5])
    bars.append([6, 2])
    bars.append([7, 3])
    bars.append([8, 4])

    myTruss = Truss(nodes, bars, 0.111, 7850, 1e4)
    myTruss.addSupport(6, 0, 0, 0)
    myTruss.addSupport(7, 0, 0, 0)
    myTruss.addSupport(8, 0, 0, 0)
    myTruss.addSupport(9, 0, 0, 0)

    myTruss.addExternalForce(0, 0, 0, -10)
    myTruss.addExternalForce(1, 0, 0, -10)


def tetrahedron():
    global myTruss

    nodes.append([0, 0, -3])
    nodes.append([20, 0, -3])
    nodes.append([10, 9, -3])
    nodes.append([10, 6, 3])

    bars.append([0, 1])
    bars.append([1, 2])
    bars.append([0, 2])
    bars.append([0, 3])
    bars.append([1, 3])
    bars.append([2, 3])

    myTruss = Truss(nodes, bars, 0.111, 7850, 1e4)
    myTruss.addSupport(0, 1, 1, 0)
    myTruss.addSupport(1, 1, 1, 0)
    myTruss.addSupport(2, 1, 1, 0)

    myTruss.addExternalForce(3, 0, 0, -100)


def triangle():
    global myTruss

    nodes.append([0, 0, 0])
    nodes.append([0, 1, 0])
    nodes.append([1, 0, 0])

    bars.append([0, 1])
    bars.append([1, 2])
    bars.append([0, 2])

    myTruss = Truss(nodes, bars, 0.111, 7850, 1e4)
    myTruss.addSupport(0, 0, 0, 0)
    # myTruss.addSupport(1, 0, 0, 0)
    # myTruss.addSupport(2, 0, 0, 0)

    myTruss.addExternalForce(2, 0, 0, -10)


def bridge():
    global myTruss

    nodes.append([0, 0, 0])
    nodes.append([1, 0, 0])
    # nodes.append([2, 0, 0])
    # nodes.append([3, 0, 0])
    # nodes.append([4, 0, 0])

    bars.append([0, 1])
    # bars.append([1,2])
    # bars.append([2,3])
    # bars.append([3,4])

    myTruss = Truss(nodes, bars, 0.111, 7850, 1e4)
    myTruss.addSupport(0, 0, 0, 0)

    myTruss.addExternalForce(1, -10, 0, 0)  # funktioniert nicht


nodes = []
bars = []
myTruss = None

bridge()

fem = FEM(myTruss)

N, R, U = fem.TrussAnalysis()
print(myTruss.F)
print('\nAxial Forces (positive = tension, negative = compression)')
# Anschaulichkeit
print(N[np.newaxis].T)
print('\nReaction Forces (positive = upward, negative = downward')
print(R)
print('\nDeformation at nodes')
print(U)

fem.Plot(np.array(nodes).astype(float), np.array(bars), 'gray', '--', 1, 'Undeformed')
scale = 5
# Berechne die neue Position der Knoten
Dnodes = U * scale + nodes
fem.Plot(np.array(Dnodes), np.array(bars), 'red', '-', 2, 'Deformed')
plt.show()
