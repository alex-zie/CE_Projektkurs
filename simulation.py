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
   h_b = 5e-2  # heights of bar profile  [m]
   b_b = 5e-2  # width of bar profile  [m]
   A = h_b * b_b  # area of bar profile  [m^2]
   E = 210e9  # E-module [Pa]
   rho = 7850  # density [kg/m^3]
   load = 500e3 # attached weight [N]

   #myCrane = crane_1(tower_height, jib_length, length_segments, A, rho, E)
   myCrane = crane_2_1(10, 10, 1, A, rho, E)
   #myCrane = crane_2_2(10, 10, 1, A, rho, E)

   nodes = myCrane.nodes
   bars = myCrane.bars

   # weight
   for i in myCrane.tip_nodes:
      myCrane.addExternalForce(i, 0, 0, -load/len(myCrane.tip_nodes))

   # counter weight
   for i in myCrane.counterweight_nodes:
      myCrane.addExternalForce(i, 0, 0, -3*load/len(myCrane.counterweight_nodes))
   
   fem = FEM(myCrane, own_weight=True)
   fem.display(scale=1, tension=True)
   #fem.optimize_crossections(625e-4, 200e6)
   #fem.homogenize_tensions(625e-4, 200e6)
   
   # fem.addWind(28, 1, 1)
   
   # visualization
   fem.display(scale=1, tension=True)

   # highlight critical bars
   # plt.subplot(projection='3d')
   # fem.plot(nodes, bars[fem.getTension() > 200e6], 'lightskyblue', '-', 2)          
   # plt.show()