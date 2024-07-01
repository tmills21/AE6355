import math

class vehicle:
    def __init__(self, mass, diameter, noseRadius, deltac):
        self.mass = mass # kg
        self.diamter = diameter # m
        self.nodeRadius = noseRadius # m
        self.deltac = math.radians(deltac) # degrees to radians
        self.LDrat = 0 # purely ballistic TODO change

        self.ballisticCoeff = self.computeBallisticCoeff()

    def computeBallisticCoeff(self):
        # TODO fix this
        return 200