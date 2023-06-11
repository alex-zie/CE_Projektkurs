import numpy as np
from truss import Truss
import matplotlib.pyplot as plt


class FEM:

    def __init__(self, truss: Truss, own_weight: bool):
        self.truss = truss
        self.NN = len(truss.nodes)
        self.NE = len(truss.bars)
        self.DOF = 3  # weil wir uns in 3D befinden
        self.NDOF = self.DOF * self.NN  # Gesamtanzahl der Freihetsgrade
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
                               axis=1)  # Transformationsvektor lokal -> global
        K = self.computeStiffnessMatrix(E, A, L, trans)
        # print("Determinant K:", np.linalg.det(K))
        freeDOF = self.truss.supports.flatten().nonzero()[0]  # Prüfe, welche Knoten FG > 0 haben
        supportDOF = (self.truss.supports.flatten() == 0).nonzero()[0]  # Knoten mit Lagern
        Kff = K[np.ix_(freeDOF, freeDOF)]  # Vollkommen Bewegliche knoten
        Kfr = K[np.ix_(freeDOF, supportDOF)]  # Teilweise bewegliche Knoten # sicher?
        Krf = Kfr.T
        # Krr = K[np.ix_(supportDOF, supportDOF)]  # für die Lagerkräfte

        # Weight
        if (self.firstCreate and self.own_weight):
            weights = np.zeros_like(self.truss.F)
            weights[:, 2] = -self.computeWeight()
            self.truss.addExternalForces(weights)

        F = self.truss.F.flatten()[freeDOF]  # Kraftmatrix passend zu K mit nicht null Einträgen, wie oben definiert
        Uf = np.linalg.lstsq(Kff, F, rcond=None)[0]
        U = self.truss.supports.astype(float).flatten()
        U[freeDOF] = Uf
        U[supportDOF] = self.truss.Ur
        U = U.reshape(self.NN, self.DOF)
        u = np.concatenate((U[self.truss.bars[:, 0]], U[self.truss.bars[:, 1]]),
                           axis=1)  # Verschiebungsvektor für die einzelnen Elemente
        N = E * A / L[:] * (trans[:] * u[:]).sum(axis=1)  # interne Kräfte
        R = (Krf[:] * Uf).sum(axis=1)  # + (Krr[:] * self.truss.Ur).sum(axis=1)  # Reaktionskräfte
        # update N, R, U
        self.__dict__["N"] = N
        self.__dict__["R"] = R
        self.__dict__["U"] = U

    def computeStiffnessMatrix(self, E, A, L, trans):
        K = np.zeros([self.NDOF, self.NDOF])
        for k in range(self.NE):
            aux = self.DOF * self.truss.bars[k, :]
            # Indices der vom Element betroffenen Knoten für Position der Summierung
            index = np.r_[aux[0]:aux[0] + self.DOF, aux[1]:aux[1] + self.DOF]
            # lokale Steifigkeiten, np.newaxis für 0 Zeile(4x4)
            ES = np.dot(trans[k][np.newaxis].T * E * A, trans[k][np.newaxis]) / L[k]
            # Globale Steifigkeiten durch Summierung der Einzelsteifigkeiten, Position !
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

        # Height of the bars
        area_bars = self.truss.A ** 0.5 * lengths

        # add wind force to node for each bar that is connected to it
        # Wind pressure: 0.5 * rho * v^2 [m/s] * A [m^2]
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

        # Define the color scale
        color_scale = [
            (0, (0, 255, 0)),      # Green
            (0.2, (255, 255, 0)),  # Yellow
            (0.4, (255, 165, 0)),  # Orange
            (0.6, (255, 0, 0)),    # Red
            (1, (128, 0, 128))     # Purple
        ]
        
        # Find the appropriate color range for the value
        for i in range(len(color_scale) - 1):
            if value <= color_scale[i + 1][0]:
                break
        
        # Interpolate between the colors
        start_value, start_color = color_scale[i]
        end_value, end_color = color_scale[i + 1]
        ratio = (value - start_value) / (end_value - start_value)
        color = (
            int(start_color[0] + ratio * (end_color[0] - start_color[0])),
            int(start_color[1] + ratio * (end_color[1] - start_color[1])),
            int(start_color[2] + ratio * (end_color[2] - start_color[2]))
        )
        
        # Convert the color to HEX format
        hex_color = '#{:02x}{:02x}{:02x}'.format(*color)
        
        return hex_color
    
    def getColorMap(self):
        min = np.min(self.N)
        max = np.max(self.N)
        color = []
        for bar in self.truss.bars:
            force = (self.N[bar[0]] + self.N[bar[1]]) / 2
            value = (force - min) / (max - min) # Wert normalisieren
            color.append(self.map_value_to_color(value))
        return color

    def display(self, scale=1, external_forces=True, tension=False):
        # Zeichne undeformierten Kran
        self.plot(self.truss.nodes, self.truss.bars, 'gray', '--', 1)

        # Zeichne deformierten Kran
        dnodes = self.U * scale + self.truss.nodes # Verschiebung der Knoten

        if not tension:
            self.plot(dnodes, self.truss.bars, 'red', '-', 2)
        else:
            # Normalkräfte farblich visualisieren
            colors = self.getColorMap()
            self.plot(dnodes, self.truss.bars, colors, '-', 2)

        if external_forces:
            for i in self.truss.force_points:
                self.plotPoint(dnodes[i])

        plt.suptitle(self.truss)
        plt.title("forces: ["+str(int(np.min(self.N)/1000))+", "+str(int(np.max(self.N)/1000))+"] kN", fontsize=10)
        # Graphik speichern
        plt.savefig('figures/fig1', dpi=600)
        plt.show()
        
    


    def plot(self, nodes, bars, color, lt, lw, lg=None):
        if isinstance(color, str):
            color = [color] * len(bars)
        plt.subplot(projection='3d')
        plt.gca().set_aspect('equal')
        # plt.gca(projection='3d')
        for i in range(len(bars)):
            # Jeweilige Start und Endkoordiante
            xi, xf = nodes[bars[i, 0], 0], nodes[bars[i, 1], 0]
            yi, yf = nodes[bars[i, 0], 1], nodes[bars[i, 1], 1]
            zi, zf = nodes[bars[i, 0], 2], nodes[bars[i, 1], 2]
            # Plotte die Elemente
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
