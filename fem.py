import numpy as np
from truss import Truss
import matplotlib.pyplot as plt


class FEM:

    def __init__(self, truss: Truss, own_weight: bool):
        """
        truss : Truss
            the Truss that is to be simulated
        own_weight : bool
            whether or not its own weight should be accounted for
        """
        self.truss = truss
        self.num_nodes = len(truss.nodes)
        self.num_bars = len(truss.bars)
        self.DOF = 3  # because we are in 3D 
        self.NDOF = self.DOF * self.num_nodes  # total number of degrees of freedom
        self.own_weight = own_weight
        self.wind = False # True if wind is being simulated
        self.wind_dir = -1 # Diection of wind (0 to 3, -1 if no wind)
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

    def F_crit(self):
        """Return critical bending force """
        F_k = np.pi ** 2 * self.truss.E * self.truss.I / (self.truss.lengths ** 2)
        F_k[40:44] = 1.43 ** 2 * F_k[0:4]
        return F_k

    def check_bending_force(self):
        b = self.N[np.where(self.N < 0)[0]] < self.F_crit()[np.where(self.N < 0)[0]]
        return len(np.where(~b)[0]) == 0

    def TrussAnalysis(self, p=False):
        """
        Performs a statical analysis using the finite element method.
        Modifies the attributes N, R, U
        """
        if p and self.firstCreate:
            message = "Simulating truss with " + str(self.num_nodes) + " nodes and " + str(self.num_bars) + " bars"
            message += " considering gravity." if self.own_weight else "."
            print(message)

        E = self.truss.E
        A = self.truss.A
        L = self.truss.lengths
        trans = np.concatenate((-self.truss.orientations.T, self.truss.orientations.T),
                               axis=1)  # transformation vector local -> global
        K = self.computeStiffnessMatrix(E, A, L, trans)
        freeDOF = self.truss.supports.flatten().nonzero()[0]  # check, which nodes are FG > 0 
        supportDOF = (self.truss.supports.flatten() == 0).nonzero()[0]  # nodes with supports
        Kff = K[np.ix_(freeDOF, freeDOF)]  # split K appropriately
        Kfr = K[np.ix_(freeDOF, supportDOF)]
        Krf = Kfr.T

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
        U = U.reshape(self.num_nodes, self.DOF)
        u = np.concatenate((U[self.truss.bars[:, 0]], U[self.truss.bars[:, 1]]),axis=1)  # displacement vector for each element
        N = E * A / L[:] * (trans[:] * u[:]).sum(axis=1)  # internal forces
        R = (Krf[:] * Uf).sum(axis=1) # reactional forces
        # update N, R, U
        self.__dict__["N"] = N
        self.__dict__["R"] = R
        self.__dict__["U"] = U

    def computeStiffnessMatrix(self, E, A, L, trans):
        if isinstance(A, float) or isinstance(A, int):
            A = A*np.ones(self.num_bars)
        K = np.zeros([self.NDOF, self.NDOF])
        for k in range(self.num_bars):
            aux = self.DOF * self.truss.bars[k, :]
            # Indices of elements of effected nodes for position of sum
            index = np.r_[aux[0]:aux[0] + self.DOF, aux[1]:aux[1] + self.DOF]
            # local stiffness, np.newaxis for 0 lines(4x4)
            ES = np.dot(trans[k][np.newaxis].T * E * A[k], trans[k][np.newaxis]) / L[k]
            # global stiffness: sum of single stiffnesses, position !
            K[np.ix_(index, index)] = K[np.ix_(index, index)] + ES
        return K

    def computeWeight(self):
        """
        Returns a force vector that assigns a force to each node by averaging the weight of the connected bars
        """
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
        
        A = self.truss.A*np.ones(self.num_bars)
        
        self.wind = True

        message = "Simulating wind with " + str(speed) + " m/s (" + str(round(speed * 3.6, 2)) + " km/h) blowing in "
        message += "negative " if dir == -1 else "positive "
        message += "x" if axis == 0 else "y"
        message += "-direction."
        print(message)

        if axis == 0:
            if dir == -1:
                bars = self.truss.bars[self.truss.x_positive_side]
                lengths = self.truss.lengths[self.truss.x_positive_side]
                A = A[self.truss.x_positive_side]
                self.wind_dir = 0
            else:
                bars = self.truss.bars[self.truss.x_negative_side]
                lengths = self.truss.lengths[self.truss.x_negative_side]
                A = A[self.truss.x_negative_side]
                self.wind_dir = 1
        else:
            if dir == -1:
                bars = self.truss.bars[self.truss.y_positive_side]
                lengths = self.truss.lengths[self.truss.y_positive_side]
                A = A[self.truss.y_positive_side]
                self.wind_dir = 2
            else:
                bars = self.truss.bars[self.truss.y_negative_side]
                lengths = self.truss.lengths[self.truss.y_negative_side]
                A = A[self.truss.y_negative_side]
                self.wind_dir = 3

        nOutgoingBars = np.zeros_like(self.truss.nodes[:, 0])
        # compute for each node (index) the number of outgoing bars
        for bar in bars:
            nOutgoingBars[bar] = nOutgoingBars[bar] + 1

        # height of the bars
        area_bars = np.multiply(A ** 0.5, lengths)

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

    def optimize_crossections(self, max_A, crit_tension):
        """
        Minimizes the crossections in the interval [A, max_A] based on a critical tension

        returns: Array of crossections
        """
        abs_max_tension = np.max([np.abs(np.min(self.getTension())), np.abs(np.max(self.getTension()))])

        if isinstance(self.truss.A, float) or isinstance(self.truss.A, int):
            min_A = self.truss.A
            self.truss.A = np.zeros(len(self.truss.bars))
        else:
            min_A = np.max(self.truss.A)

        while abs_max_tension > crit_tension:
            critical_bars = np.abs(self.getTension()) > 0.999*crit_tension # indices where tension is exceeded # 0.995 to make it slightly smaller so that increased weight force is accounted for
            self.truss.A[critical_bars] = np.abs(self.N[critical_bars]) / (0.999*crit_tension) # increase crossections here
            self.truss.A[self.truss.A < min_A] = min_A # keep crossections whithin allowed range
            self.truss.A[self.truss.A > max_A] = max_A
            self.truss._computeMass() # update mass
            
            self.TrussAnalysis() # recalculate forces
            abs_max_tension = np.max([np.abs(np.min(self.getTension())), np.abs(np.max(self.getTension()))])       
        
        return self.truss.A

    def homogenize_tensions(self, max_A, crit_tension):
        """
        Changes crossections such that tension is equal

        returns: Array of crossections
        """
        abs_max_tension = np.max([np.abs(np.min(self.getTension())), np.abs(np.max(self.getTension()))])

        if isinstance(self.truss.A, float) or isinstance(self.truss.A, int):
            min_A = self.truss.A
            self.truss.A = np.zeros(len(self.truss.bars))
        else:
            min_A = np.max(self.truss.A)

        while abs_max_tension > crit_tension:
            self.truss.A = np.abs(self.N) / (0.995*crit_tension) # increase crossections here
            self.truss.A[self.truss.A < min_A] = min_A # keep crossections whithin allowed range
            self.truss.A[self.truss.A > max_A] = max_A

            self.TrussAnalysis() # recalculate forces
            abs_max_tension = np.max([np.abs(np.min(self.getTension())), np.abs(np.max(self.getTension()))])
            self.truss._computeMass() # update mass
                
        return self.truss.A


    def reset(self):
        """
        Removes all external forces including gravity and wind
        """
        self.firstCreate = True
        self.wind = False
        self.wind_dir = -1
        self.truss.reset()
        self.TrussAnalysis()

    def map_value_to_color(self, value, color_map):
        """
        Returns a hexadecimal color number as a string for a value between -1 and 1 based
        on the given color map
        """
        if value < -1 or value > 1:
            raise Exception("The input value must lay between -1 and 1!")
        
        # find the appropriate color range for the value
        for i in range(len(color_map) - 1):
            if value < color_map[i + 1][0]:
                break
        
        # interpolate between the colors
        start_value, start_color = color_map[i]
        end_value, end_color = color_map[i + 1]
        ratio = (value - start_value) / (end_value - start_value)
        color = (
            int(start_color[0] + ratio * (end_color[0] - start_color[0])),
            int(start_color[1] + ratio * (end_color[1] - start_color[1])),
            int(start_color[2] + ratio * (end_color[2] - start_color[2]))
        )
        
        # convert the color to HEX format
        hex_color = '#{:02x}{:02x}{:02x}'.format(*color)
        
        return hex_color
    
    def getColorMap(self, values, min, max):
        """
        values : array-like
        min :
            minimum value to be assigned a color
        max :
            maximum value to be assigned a color

        Returns:
            * a color array that maps a color to each element in values
            * the color map that was used
        """

        # define the color map
        color_map = [
            (-1, (128, 0, 128)),    # purple
            (-2/3, (0, 0, 255)),   # blue
            (-1/3, (0, 255, 255)), # cyan
            (0, (0, 255, 0)),       # green
            (1/3, (255, 255, 0)),  # yellow
            (2/3, (255, 165, 0)),  # orange
            (1, (200, 0, 64))       # dark red
        ]

        abs_max = np.max([np.abs(min), np.abs(max)]) # so that 0 is later mapped to zero
        color = []
        for value in values:
            color.append(self.map_value_to_color(value/abs_max, color_map)) # normalize value to fit interval [-1, 1]
        return color, color_map

    def display(self, scale=1, external_forces=True, tension=False, wind=False):
        """
        scale : float
            the scale in which the deformations should be displayed
        external_forces : bool
            display points where single external are being applied
        tension : bool
            highlight the tension in each bar using a color map
        wind : bool
            highlight the surface at which wind is attacking
        """
        if tension and wind:
            print("It is not recommended to have both tension and wind exposed surfaces displayed!")

        tensions = self.getTension()
        minTension = np.min(tensions)
        maxTension = np.max(tensions)

        if np.abs(minTension) >= 1e7 or np.abs(maxTension) >= 1e7:
            minTension = minTension/1e6  
            maxTension = maxTension/1e6
            tensions = tensions/1e6
            unit = 'MPa'
        else: # both are < 1e7
            minTension = minTension/1e3  
            maxTension = maxTension/1e3
            tensions = tensions/1e3
            unit = 'KPa'

        # plot crane without deformations
        self.plot(self.truss.nodes, self.truss.bars, 'gray', '--', 1)

        # plot crane with deformations 
        dnodes = self.U * scale + self.truss.nodes # displacement of nodes

        if not tension:
            self.plot(dnodes, self.truss.bars, 'red', '-', 2)
        else:
            # visualize tensions with colors 
            colors, color_map = self.getColorMap(tensions, minTension, maxTension)
            self.plot(dnodes, self.truss.bars, colors, '-', 2)
            
            # legend
            color_map = np.array(color_map)
            abs_max = np.max([np.abs(minTension), np.abs(maxTension)])
            legend_values = []
            legend_values.append(minTension) # start at minimum tension
            for value in color_map[:, 0]:
                # it could be that there are colors defined for nonexisting tensions
                # we don't want those
                if value*abs_max <= minTension:
                    continue
                legend_values.append(value*abs_max) # fill up till maximum is reached

            lines = []
            legend_colors, _ = self.getColorMap(legend_values, minTension, maxTension)
            for color in legend_colors:
                lines.append(plt.Line2D([0], [0], color=color, lw=2))
            lines.reverse()
            
            labels = []
            for value in legend_values:
                labels.append(str(int(value))+" "+unit)
            labels.reverse()

            plt.legend(lines, labels, title="Spannungen")

        if external_forces:
            for i in self.truss.force_points:
                self.plotPoint(dnodes[i])

        if wind:
            if not self.wind:
                print("There is no wind.")
            else:
                if self.wind_dir == 0:
                    self.plot(dnodes, self.truss.bars[self.truss.x_positive_side], 'lightskyblue', '-', 2, 'wind exposed area')
                elif self.wind_dir == 1:
                    self.plot(dnodes, self.truss.bars[self.truss.x_negative_side], 'lightskyblue', '-', 2, 'wind exposed area')
                elif self.wind_dir == 2:
                    self.plot(dnodes, self.truss.bars[self.truss.y_positive_side], 'lightskyblue', '-', 2, 'wind exposed area')
                elif self.wind_dir == 3:
                    self.plot(dnodes, self.truss.bars[self.truss.y_negative_side], 'lightskyblue', '-', 2, 'wind exposed area')

        plt.suptitle(self.truss)
        plt.title("Höhe: "+str(self.truss.height)+" m\n"
                + "Länge: "+str(self.truss.length)+" m\n"
                + "Segmentlänge: "+str(self.truss.ls)+" m\n"
                + "Masse: "+str(int(np.sum(self.truss.mass)))+" kg", fontsize=8, y = 0.95, loc="right")
        # save plot
        # plt.savefig('figures/fig1', dpi=600)
        plt.show()


    def plot(self, nodes, bars, color, lt, lw, lg=None):
        """
        nodes : number array
        bars : number array
        color : string or array of string
        lt : string
            line style
        lw : number
            linewidth
        """
        if isinstance(color, str):
            color = [color] * len(bars)
        plt.subplot(projection='3d')
        plt.gca().set_aspect('equal')
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
        """
        Plots the point at index node
        """
        plt.subplot(projection='3d')
        plt.plot(self.truss.nodes[node, 0], self.truss.nodes[node, 1], self.truss.nodes[node, 2], 'bo')

    def plotPoint(self, point):
        """
        Plots a 3D point
        """
        plt.subplot(projection='3d')
        plt.plot(point[0], point[1], point[2], 'bo')
