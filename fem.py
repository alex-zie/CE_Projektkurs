import numpy as np
from truss import Truss
import matplotlib.pyplot as plt


class FEM:

    def __init__(self, truss: Truss, own_weight: bool):
        self.truss = truss
        self.NN = len(truss.nodes)
        self.NE = len(truss.bars)
        self.DOF = 3  # because we are in 3D 
        self.NDOF = self.DOF * self.NN  # total number of degrees of freedom
        self.own_weight = own_weight
        self.firstCreate = True  # so that weights don't get added multiple times
        self.TrussAnalysis()
        self.firstCreate = False

    @property
    def N(self):
        """Return axial forces"""
        return self.__dict__["N"]

    @property
    def R(self):
        """Return reactional forces"""
        return self.__dict__["R"]

    @property
    def U(self):
        """Return deformations"""
        return self.__dict__["U"]

    def F_krit(self):
        """Return critical bending force """
        F_k = np.pi ** 2 * self.truss.E * self.truss.I / (self.truss.lengths ** 2)
        F_k[40:44] = 1.43 ** 2 * F_k[0:4]
        return F_k

    def check_bending_force(self):
        b = self.N[np.where(self.N < 0)[0]] < self.F_krit()[np.where(self.N < 0)[0]]
        return len(np.where(~b)[0]) == 0

    def TrussAnalysis(self, p=False):
        """
        returns: axial forces, reactional forces, displacements
        """
        if p and self.firstCreate:
            message = "Simulating truss with " + str(self.NN) + " nodes and " + str(self.NE) + " bars"
            message += " considering gravity." if self.own_weight else "."
            print(message)

        E = self.truss.E
        A = self.truss.A
        L = self.truss.lengths
        trans = np.concatenate((-self.truss.orientations.T, self.truss.orientations.T),
                               axis=1)  # transformation vector local -> global
        K = self.computeStiffnessMatrix(E, A, L, trans)
        # print("Determinant K:", np.linalg.det(K))
        freeDOF = self.truss.supports.flatten().nonzero()[0]  # check, which nodes are FG > 0 
        supportDOF = (self.truss.supports.flatten() == 0).nonzero()[0]  # Knoten mit Lagern
        Kff = K[np.ix_(freeDOF, freeDOF)]  # totally movable nodes 
        Kfr = K[np.ix_(freeDOF, supportDOF)]  # partially movable nodes  # are we shuurree?
        Krf = Kfr.T
        # Krr = K[np.ix_(supportDOF, supportDOF)]  # for the support forces

        # weight
        if (self.firstCreate and self.own_weight):
            weights = np.zeros_like(self.truss.F)
            weights[:, 2] = -self.computeWeight()
            self.truss.addExternalForces(weights)

        F = self.truss.F.flatten()[freeDOF]  # matrix of forces, like K, with not-null entries, as definied before
        Uf = np.linalg.lstsq(Kff, F, rcond=None)[0]
        U = self.truss.supports.astype(float).flatten()
        U[freeDOF] = Uf
        U[supportDOF] = self.truss.Ur
        U = U.reshape(self.NN, self.DOF)
        u = np.concatenate((U[self.truss.bars[:, 0]], U[self.truss.bars[:, 1]]),axis=1)  # displacement vector for each element
        N = E * A / L[:] * (trans[:] * u[:]).sum(axis=1)  # internal forces
        R = (Krf[:] * Uf).sum(axis=1)  # + (Krr[:] * self.truss.Ur).sum(axis=1)  # reactional forces
        # update N, R, U
        self.__dict__["N"] = N
        self.__dict__["R"] = R
        self.__dict__["U"] = U

    def computeStiffnessMatrix(self, E, A, L, trans):
        K = np.zeros([self.NDOF, self.NDOF])
        for k in range(self.NE):
            aux = self.DOF * self.truss.bars[k, :]
            # Indices of elements of effected nodes for position of sum
            index = np.r_[aux[0]:aux[0] + self.DOF, aux[1]:aux[1] + self.DOF]
            # local stiffness, np.newaxis for 0 lines(4x4)
            ES = np.dot(trans[k][np.newaxis].T * E * A, trans[k][np.newaxis]) / L[k]
            # global stiffness: sum of single stiffnesses, position !
            K[np.ix_(index, index)] = K[np.ix_(index, index)] + ES
        return K

    def computeWeight(self):
        bars = self.truss.bars
        nOutgoingBars = np.zeros_like(self.truss.nodes[:, 0])

        # compute for each node (index) the number of outgoing bars
        for bar in bars:
            nOutgoingBars[bar] = nOutgoingBars[bar] + 1

        masses = self.truss.mass
        weights = np.zeros_like(self.truss.nodes[:, 0])
        # add a weight to the node equalling the mass of a connected bar
        for i in range(len(bars)):
            weights[bars[i]] = weights[bars[i]] + masses[i]

        return 9.81 * weights / nOutgoingBars

    def addWind(self, speed, axis, dir):
        """
        speed in m/s
        axis: 0 (x-Axis) or 1 (y-Axis)
        dir: 1 (positive direction) or -1 (negative direction)
        """
        if not (dir == -1 or dir == 1):
            raise Exception("dir has to be either -1 or 1")
        if not (axis == 0 or axis == 1):
            raise Exception("axis can only take values 0 (x), 1 (y)")

        message = "Simulating wind with " + str(speed) + " m/s (" + str(round(speed * 3.6, 2)) + " km/h) blowing in "
        message += "negative " if dir == -1 else "positive "
        message += "x" if axis == 0 else "y"
        message += "-direction."
        print(message)

        if axis == 0:
            if dir == -1:
                bars = self.truss.bars[self.truss.x_negative_side]
                lengths = self.truss.lengths[self.truss.x_negative_side]
            else:
                bars = self.truss.bars[self.truss.x_positive_side]
                lengths = self.truss.lengths[self.truss.x_positive_side]
        else:
            if dir == -1:
                bars = self.truss.bars[self.truss.y_negative_side]
                lengths = self.truss.lengths[self.truss.y_negative_side]
            else:
                bars = self.truss.bars[self.truss.y_positive_side]
                lengths = self.truss.lengths[self.truss.y_positive_side]

        nOutgoingBars = np.zeros_like(self.truss.nodes[:, 0])
        # compute for each node (index) the number of outgoing bars
        for bar in bars:
            nOutgoingBars[bar] = nOutgoingBars[bar] + 1

        # height of the bars
        area_bars = self.truss.A ** 0.5 * lengths

        # add wind force to node for each bar that is connected to it
        # wind pressure: 0.5 * rho * v^2 [m/s] * A [m^2]
        forces = np.zeros_like(self.truss.nodes[:, 0])
        for i in range(len(bars)):
            forces[bars[i]] = forces[bars[i]] + 0.5 * area_bars[
                i] * 1.2 * speed ** 2  # 1.2 is the air density at sea level

        w_forces = np.zeros_like(self.truss.F)
        nOutgoingBars[nOutgoingBars == 0] = 1  # prevent division by zero
        w_forces[:, axis] = dir * forces / nOutgoingBars
        self.truss.addExternalForces(w_forces)
        self.TrussAnalysis()

    def getTension(self):
        return self.N / self.truss.A

    def reset(self):
        """Removes all external forces including gravity and wind"""
        self.firstCreate = True
        self.truss.reset()
        self.TrussAnalysis()

    def map_value_to_color(self, value):
        if value < 0 or value > 1:
            raise Exception("The input value must lay between 0 and one!")

        # define the color scale
        color_scale = [
            (0, (0, 255, 0)),      # green
            (0.25, (255, 255, 0)),  # yellow
            (0.5, (255, 165, 0)),  # orange
            (0.75, (255, 0, 0)),    # red
            (1, (128, 0, 128))     # purple
        ]
        
        # find the appropriate color range for the value
        for i in range(len(color_scale) - 1):
            if value <= color_scale[i + 1][0]:
                break
        
        # interpolate between the colors
        start_value, start_color = color_scale[i]
        end_value, end_color = color_scale[i + 1]
        ratio = (value - start_value) / (end_value - start_value)
        color = (
            int(start_color[0] + ratio * (end_color[0] - start_color[0])),
            int(start_color[1] + ratio * (end_color[1] - start_color[1])),
            int(start_color[2] + ratio * (end_color[2] - start_color[2]))
        )
        
        # convert the color to HEX format
        hex_color = '#{:02x}{:02x}{:02x}'.format(*color)
        
        return hex_color
    
    def getColorMap(self):
        tension = self.getTension()
        min = np.min(tension)
        max = np.max(tension)
        color = []
        for bar in self.truss.bars:
            force = (tension[bar[0]] + tension[bar[1]]) / 2
            value = (force - min) / (max - min) # normalise value 
            color.append(self.map_value_to_color(value))
        return color

    def display(self, scale=1, external_forces=True, tension=False):
        minTension = np.min(self.getTension())/1e6 # minimum tension [MPa]
        maxTension = np.max(self.getTension())/1e6 # maximum tension [MPa]

        # plot crane without deformations
        self.plot(self.truss.nodes, self.truss.bars, 'gray', '--', 1)

        # plot crane with deformations 
        dnodes = self.U * scale + self.truss.nodes # displacement of nodes

        if not tension:
            self.plot(dnodes, self.truss.bars, 'red', '-', 2)
        else:
            # visualize normal forces with colors 
            colors = self.getColorMap()
            self.plot(dnodes, self.truss.bars, colors, '-', 2)
            # legend
            color_scale = [
                (0, 255, 0),    # green
                (255, 255, 0),  # yellow
                (255, 165, 0),  # orange
                (255, 0, 0),    # red
                (128, 0, 128)]  # purple
        
            lines = [plt.Line2D([0], [0], color='#{:02x}{:02x}{:02x}'.format(*color_scale[0]), lw=2),
                     plt.Line2D([0], [0], color='#{:02x}{:02x}{:02x}'.format(*color_scale[1]), lw=2),
                     plt.Line2D([0], [0], color='#{:02x}{:02x}{:02x}'.format(*color_scale[2]), lw=2),
                     plt.Line2D([0], [0], color='#{:02x}{:02x}{:02x}'.format(*color_scale[3]), lw=2),
                     plt.Line2D([0], [0], color='#{:02x}{:02x}{:02x}'.format(*color_scale[4]), lw=2)]
            
            labels = [str(int(minTension))+" MPa",
                      str(int(minTension + 0.25*(maxTension-minTension)))+" MPa",
                      str(int(minTension + 0.5*(maxTension-minTension)))+" MPa",
                      str(int(minTension + 0.75*(maxTension-minTension)))+" MPa",
                      str(int(maxTension))+" MPa"]

            plt.legend(lines, labels, title="Spannungen")

        if external_forces:
            for i in self.truss.force_points:
                self.plotPoint(dnodes[i])

        plt.suptitle(self.truss)
        plt.title("spannung: ["+str(int(minTension))+", "+str(int(maxTension))+"] MPa", fontsize=10)
        # save plot
        # plt.savefig('figures/fig1', dpi=600)
        plt.show()


    def plot(self, nodes, bars, color, lt, lw, lg=None):
        if isinstance(color, str):
            color = [color] * len(bars)
        plt.subplot(projection='3d')
        plt.gca().set_aspect('equal')
        # plt.gca(projection='3d')
        for i in range(len(bars)):
            # each start and end coordinate
            xi, xf = nodes[bars[i, 0], 0], nodes[bars[i, 1], 0]
            yi, yf = nodes[bars[i, 0], 1], nodes[bars[i, 1], 1]
            zi, zf = nodes[bars[i, 0], 2], nodes[bars[i, 1], 2]
            # plot elements 
            line, = plt.plot([xi, xf], [yi, yf], [zi, zf], color=color[i], linestyle=lt, linewidth=lw)
        if lg is not None:
            line.set_label(lg)
            plt.legend(prop={'size': 7})
        

    def plotNode(self, node):
        plt.subplot(projection='3d')
        plt.plot(self.truss.nodes[node, 0], self.truss.nodes[node, 1], self.truss.nodes[node, 2], 'bo')

    def plotPoint(self, point):
        plt.subplot(projection='3d')
        plt.plot(point[0], point[1], point[2], 'bo')
