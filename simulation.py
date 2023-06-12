from crane import crane_1
from crane import crane_2_1
from crane import crane_2_2
from fem import FEM
import numpy as np
import matplotlib.pyplot as plt

if __name__ == "__main__":
   tower_height = 10 # maximum height of crane  [m]
   jib_length = 10 # maximum lengths of jib [m]
   length_segments = 1 # length of segments (not bars) [m]
   h_b = 10e-2  # heights of bars profile  [m]
   b_b = 10e-2  # width of bars profile  [m]
   A = h_b * b_b  # area of bars profile  [m^2]
   E = 210e9  # E-module [Pa]
   rho = 7850  # density [kg/m^3]
   load = 500e3 # attached weight [kg]

   # myCrane = crane_1(tower_height, jib_length, length_segments, A, rho, E)
   # myCrane = crane_2_1(10, 10, 1, A, rho, E)
   myCrane = crane_2_2(10, 10, 1, A, rho, E)

   nodes = myCrane.nodes
   bars = myCrane.bars

   # weight
   for i in myCrane.tip_nodes:
      myCrane.addExternalForce(i, 0, 0, -load/len(myCrane.tip_nodes))

   # counter weight
   for i in myCrane.counterweight_nodes:
      myCrane.addExternalForce(i, 0, 0, load/len(myCrane.counterweight_nodes))
   
   fem = FEM(myCrane, True)

   # fem.addWind(28, 0, -1) #TODO maybe pass arrays to wind
   
   # visualization
   fem.display(tension=True)

   # highlight target surface of the wind 
   # fem.plot(nodes, bars, 'gray', '--', 1)
   # fem.plot(nodes, bars[myCrane.x_positive_side], 'green', '-', 2, 'pos. x')
   # fem.plot(nodes, bars[myCrane.x_negative_side], 'green', '-', 2, 'neg. x')
   # fem.plot(nodes, bars[myCrane.y_positive_side], 'green', '-', 2, 'pos. y')
   # fem.plot(nodes, bars[myCrane.y_negative_side], 'green', '-', 2, 'neg. y')

   plt.show()