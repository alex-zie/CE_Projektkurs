# here we commemorate the old cubical crane (version 0) which broke down upon the slightest wind breeze

def tower(self, nodes, bars):
        # Turm erstellen
        # for 10m high tower with segement size ls we need at least
        # height / ls segements. For now i would just ignore the last segement
        # size != ls

        # Nodes des Turms
        for i in range(self.nST):
            nodes.append([0, 0, i * self.ls])  # Left Top
            nodes.append([self.ls, 0, i * self.ls])  # Right Top
            nodes.append([0, self.ls, i * self.ls])  # Left Bottom
            nodes.append([self.ls, self.ls, i * self.ls])  # Right Bottom

        # Bars des Turms
        # x- und y-Richtung (LT für Left Top usw.)
        for i in range(self.nST):
            bars.append([4 * i, 4 * i + 1])  # LT -> RT
            self.y_side.append(len(bars) - 1)
            bars.append([4 * i + 2, 4 * i + 3])  # LB -> RB
            bars.append([4 * i, 4 * i + 2])  # LT -> LB
            bars.append([4 * i + 1, 4 * i + 3])  # RT -> RB

        # z-Richtung
        for i in range(self.nST - 1):
            bars.append([4 * i, 4 * i + 4])  # LT
            self.y_side.append(len(bars) - 1)
            bars.append([4 * i + 1, 4 * i + 5])  # RT
            self.y_side.append(len(bars) - 1)
            bars.append([4 * i + 2, 4 * i + 6])  # LB
            bars.append([4 * i + 3, 4 * i + 7])  # RB

        # Kreuzstreben (+1 jeweils für näclste Stufe)
        for i in range(self.nST - 1):
            bars.append([4 * i, 4 * i + 5])  # LT -> RT+1
            self.y_side.append(len(bars) - 1)
            # bars.append([4 * i + 1, 4 * i + 4])  # RT -> RT+1
            bars.append([4 * i + 3, 4 * i + 6])  # RB -> LB+1
            # bars.append([4 * i + 2, 4 * i + 7])  # LB -> RB+1
            bars.append([4 * i + 1, 4 * i + 7])  # RT -> RB+1
            # bars.append([4 * i + 3, 4 * i + 5])  # RB -> RT+1
            bars.append([4 * i + 2, 4 * i + 4])  # LB -> LT+1
            # bars.append([4 * i, 4 * i + 6])  # LT -> LB+1

def gegenausleger(self, nodes, bars, offsetT):
    # Gegenausleger erstellen

    # Nodes des Gegenauslegers in negative x Richtung, x Koordinaten mit +1, weil der Turm auf x= 0 bix x=1 steht
    for i in range(1, self.nSGA + 1):
        nodes.append([-(self.ls + i * self.ls) + 1, 0, (self.nST - 1) * self.ls])  # Left Top
        nodes.append([-(self.ls + i * self.ls) + 1, self.ls, (self.nST - 1) * self.ls])  # Right Top
        nodes.append([-(self.ls + i * self.ls) + 1, 0, (self.nST - 2) * self.ls])  # Left Bottom
        nodes.append([-(self.ls + i * self.ls) + 1, self.ls, (self.nST - 2) * self.ls])  # Right Bottom

    # Bars des Gegenausleger
    bars.append([offsetT - 2, offsetT + 1])  # hinten oben
    bars.append([offsetT - 4, offsetT])  # vorne oben
    bars.append([offsetT - 6, offsetT + 3])  # hinten unten
    bars.append([offsetT - 8, offsetT + 2])  # vorne unten

    # x- und y-Richtung (LT für Left Top usw.)
    for i in range(self.nSGA):
        bars.append([4 * i + offsetT, 4 * i + offsetT + 1])  # LT -> RT ?
        bars.append([4 * i + 2 + offsetT, 4 * i + 3 + offsetT])  # LB -> RB
        bars.append([4 * i + offsetT, 4 * i + 2 + offsetT])  # LT -> LB
        bars.append([4 * i + 1 + offsetT, 4 * i + 3 + offsetT])  # RT -> RB

    # z-Richtung
    for i in range(self.nSGA - 1):
        bars.append([4 * i + offsetT, 4 * i + 4 + offsetT])  # LT
        bars.append([4 * i + 1 + offsetT, 4 * i + 5 + offsetT])  # RT
        bars.append([4 * i + 2 + offsetT, 4 * i + 6 + offsetT])  # LB
        bars.append([4 * i + 3 + offsetT, 4 * i + 7 + offsetT])  # RB

    # Kreuzstreben vorne
    bars.append([offsetT - 8, offsetT])
    for i in range(self.nSGA // 2):
        bars.append([offsetT + 8 * i, offsetT + 6 + 8 * i])
        if self.nSGA % 2 != 0 and i == self.nSGA // 2:  # TODO checken, ob es allgemein funktioniert
            break
        bars.append([offsetT + 6 + 8 * i, offsetT + 8 + 8 * i])
    # Kreusstreben hinten
    bars.append([offsetT - 6, offsetT + 1])
    for i in range(self.nSGA // 2):
        bars.append([offsetT + 1 + 8 * i, offsetT + 7 + 8 * i])
        if self.nSGA % 2 != 0 and i == self.nSA // 2 - 1:
            break
        bars.append([offsetT + 7 + 8 * i, offsetT + 9 + 8 * i])

def ausleger(self, nodes, bars, offsetT, offsetTG):
    # Nodes des Ausleger in positive x Richtung
    for i in range(1, self.nSA + 1):
        nodes.append([self.ls + i * self.ls, 0, (self.nST - 1) * self.ls])  # Left Top
        nodes.append([self.ls + i * self.ls, self.ls, (self.nST - 1) * self.ls])  # Right Top
        nodes.append([self.ls + i * self.ls, 0, (self.nST - 2) * self.ls])  # Left Bottom
        nodes.append([self.ls + i * self.ls, self.ls, (self.nST - 2) * self.ls])  # Right Bottom

    # Bars des Ausleger
    bars.append([offsetT - 4, offsetTG])  # LB -> LT
    bars.append([offsetT - 1, offsetTG + 1])  # RB -> RT
    bars.append([offsetT - 8, offsetTG + 2])  # LB-1 -> LB
    bars.append([offsetT - 5, offsetTG + 3])  # RB-1 -> RB
    # x- und y-Richtung (LT für Left Top usw.)
    for i in range(self.nSA):
        bars.append([4 * i + offsetTG, 4 * i + offsetTG + 1])  # LT -> RT
        bars.append([4 * i + 2 + offsetTG, 4 * i + 3 + offsetTG])  # LB -> RB
        bars.append([4 * i + offsetTG, 4 * i + 2 + offsetTG])  # LT -> LB
        bars.append([4 * i + 1 + offsetTG, 4 * i + 3 + offsetTG])  # RT -> RB

    # z-Richtung
    for i in range(self.nSA - 1):
        bars.append([4 * i + offsetTG, 4 * i + 4 + offsetTG])  # LT
        bars.append([4 * i + 1 + offsetTG, 4 * i + 5 + offsetTG])  # RT
        bars.append([4 * i + 2 + offsetTG, 4 * i + 6 + offsetTG])  # LB
        bars.append([4 * i + 3 + offsetTG, 4 * i + 7 + offsetTG])  # RB

    # Kreuzstreben vorne
    bars.append([offsetT - 7, offsetTG])
    for i in range(self.nSA // 2):
        bars.append([offsetTG + 8 * i, offsetTG + 6 + 8 * i])
        if self.nSGA % 2 != 0 and i == self.nSA // 2 - 1:
            break
        bars.append([offsetTG + 6 + 8 * i, offsetTG + 8 + 8 * i])

    # Kreusstreben hinten
    bars.append([offsetT - 5, offsetTG + 1])
    for i in range(self.nSA // 2):
        bars.append([offsetTG + 1 + 8 * i, offsetTG + 7 + 8 * i])
        if self.nSGA % 2 != 0 and i == self.nSA / 2 - 1:
            break
        bars.append([offsetTG + 7 + 8 * i, offsetTG + 9 + 8 * i])