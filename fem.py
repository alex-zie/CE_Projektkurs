import numpy as np
from truss import Truss
import matplotlib.pyplot as plt


class FEM:
    """
    class which is able to perfom a finite element simulation on the given truss
    """
    def __init__(self, truss: Truss, own_weight: bool):
        """
        :param truss: truss to perfom analysis on
        :param own_weight: boolean whether to apply own weight or not

        Performs a FEM based analysis of the given truss.
        The analysis is automatically performed when creating an instance of this class.
        """
        self.truss = truss
        self.num_nodes = len(truss.nodes)
        self.num_bars = len(truss.bars)
        self.DOF = 3  # because we are in 3D 
        self.NDOF = self.DOF * self.num_nodes  # total number of degrees of freedom
        self.own_weight = own_weight
        self.wind = False
        self.wind_dir = -1
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
        """Returns critical bending forces for each element as a matrix shaped like self.bars. """
        F_k = np.pi ** 2 * self.truss.E * self.truss.I / (self.truss.lengths ** 2)
        F_k[40:44] = 1.43 ** 2 * F_k[0:4]
        return F_k

    def check_bending_force(self):
        """Checks if all elements avoid bending under the applied forces. True if no element excedes
        the critical bending force, False otherwise"""
        b = self.N[np.where(self.N < 0)[0]] < self.F_krit()[np.where(self.N < 0)[0]]
        return len(np.where(~b)[0]) == 0

    def TrussAnalysis(self, p=False):
        """
        :param p: Determines if you want a printed message showing the characteristics of your model

        Updates the axial forces, reaction forces and deformation of self.truss by performing a full FEM based
        truss analysis.
        """
        # Printing characteristics of the model
        if p and self.firstCreate:
            message = "Simulating truss with " + str(self.num_nodes) + " nodes and " + str(self.num_bars) + " bars"
            message += " considering gravity." if self.own_weight else "."
            print(message)

        # E-Module and area of each bar
        E = self.truss.E
        A = self.truss.A

        # lengths of all bars (see also: truss.computeLengths() for details)
        L = self.truss.lengths

        # Transformation Matrix local coordinate system -> global coordinate system
        trans = np.concatenate((-self.truss.orientations.T, self.truss.orientations.T),
                               axis=1)
        # Global stiffness matrix
        K = self.computeStiffnessMatrix(E, A, L, trans)

        # Find indices of nodes which degrees of freedom (DOF) are greater than 0
        freeDOF = self.truss.supports.flatten().nonzero()[0]
        # Find indices of nodes which DOF are equal to 0
        supportDOF = (self.truss.supports.flatten() == 0).nonzero()[0]

        # Stiffness matrix entries with full DOF
        Kff = K[np.ix_(freeDOF, freeDOF)]

        # Stiffness matrix entries with some DOF
        Kfr = K[np.ix_(freeDOF, supportDOF)]
        Krf = Kfr.T

        # Stiffness matrix entries with no DOF
        Krr = K[np.ix_(supportDOF, supportDOF)]  # for the support forces

        # Apply weight to the truss
        if self.firstCreate and self.own_weight:
            weights = np.zeros_like(self.truss.F)
            weights[:, 2] = -self.computeWeight()
            self.truss.addExternalForces(weights)

        # Reshape applied forces at entries with free DOF
        P = self.truss.F.flatten()[freeDOF]  # matrix of forces, like K, with not-null entries, as definied before

        # Solve the fundamental System F = K U to get the displacement at the nodes with free DOF
        Uf = np.linalg.lstsq(Kff, P, rcond=None)[0]

        # Create deformation matrix
        U = self.truss.supports.astype(float).flatten()
        # Fill the displacement at the nodes with free DOF
        U[freeDOF] = Uf
        # Fill the displacement at the nodes with no DOF (with zeros)
        U[supportDOF] = self.truss.Ur

        # Necessary reshape operation
        U = U.reshape(self.num_nodes, self.DOF)
        # displacement vector for each bar/element -> [deformation of start node, deformation of end node]
        u = np.concatenate((U[self.truss.bars[:, 0]], U[self.truss.bars[:, 1]]),
                           axis=1)

        # Calculate axial forces for each element
        N = E * A / L[:] * (trans[:] * u[:]).sum(axis=1)
        # Calculate reaction forces for the fixed noddes
        R = (Krf[:] * Uf).sum(axis=1) + (Krr[:] * self.truss.Ur).sum(axis=1)

        # update N, R, U with resulting forces and deformations
        self.__dict__["N"] = N
        self.__dict__["R"] = R
        self.__dict__["U"] = U

    def computeStiffnessMatrix(self, E, A, L, trans):
        """
        :param E: E-Module of the bars
        :param A: Area of the bars
        :param L: Length of each bar, shaped like self.bars
        :param trans: Transformation matrix, Trasforming each bar from local to global coordinate system

        Returns the global stiffnes matrix
        """
        if isinstance(A, float) or isinstance(A, int):
            A = A*np.ones(self.num_bars)

        # Create stiffness matrix
        K = np.zeros([self.NDOF, self.NDOF])

        for k in range(self.num_bars):
            # auxiliary variable
            aux = self.DOF * self.truss.bars[k]

            # Indices of the elements degrees of freedom at the end and start of the node in the stiffness matrix
            index = np.r_[aux[0]:aux[0] + self.DOF, aux[1]:aux[1] + self.DOF]

            # stiffness of the element
            ES = np.dot(trans[k][np.newaxis].T * E * A[k], trans[k][np.newaxis]) / L[k]

            # global stiffness: sum of single stiffness matrices in the correct order given by index
            K[np.ix_(index, index)] += ES

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

        Applies a force to wind exposed nodes based on the area of the bars that are connected to it.
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
        """
        Returns the tension in each bar
        """
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
    
    # # not needed because it does the same as optimize_tensions()
    # def homogenize_tensions(self, max_A, crit_tension):
    #     """
    #     Changes crossections such that tension is equal
    #     returns: Array of crossections
    #     """
    #     abs_max_tension = np.max([np.abs(np.min(self.getTension())), np.abs(np.max(self.getTension()))])
    #     if isinstance(self.truss.A, float) or isinstance(self.truss.A, int):
    #         min_A = self.truss.A
    #         self.truss.A = np.zeros(len(self.truss.bars))
    #     else:
    #         min_A = np.max(self.truss.A)
    #     while abs_max_tension > crit_tension:
    #         self.truss.A = np.abs(self.N) / (0.995*crit_tension) # increase crossections here
    #         self.truss.A[self.truss.A < min_A] = min_A # keep crossections whithin allowed range
    #         self.truss.A[self.truss.A > max_A] = max_A
    #         self.TrussAnalysis() # recalculate forces
    #         abs_max_tension = np.max([np.abs(np.min(self.getTension())), np.abs(np.max(self.getTension()))])
    #         self.truss._computeMass() # update mass
                
    #     return self.truss.A

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
    
    def display(self, scale=1, external_forces=True, tension=False, wind=False, saveFig=None):
        """
        scale : float
            the scale in which the deformations should be displayed
        external_forces : bool
            display points where single external are being applied
        tension : bool
            highlight the tension in each bar using a color map
        wind : bool
            highlight the surface at which wind is attacking
        saveFig : string
            name of the figure that should be saved
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
                + "Segmentlänge: "+str(np.round(self.truss.ls, 2))+" m\n"
                + "Masse: "+str(int(np.sum(self.truss.mass)))+" kg", fontsize=8, y = 0.95, loc="right")
        
        if saveFig is not None:
            plt.savefig('figures/'+str(saveFig)+'.pdf', dpi=600)

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
        plt.gca().view_init(elev=19, azim=-57, roll=0)
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
