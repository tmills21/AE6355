import math

class vehicle:
    def __init__(self, mass, diameter, noseRadius, deltac, sigma):
        self.mass = mass # kg
        self.diamter = diameter # m
        self.nodeRadius = noseRadius # m
        self.deltac = math.radians(deltac) # degrees to radians
        self.LDrat = 0 # purely ballistic TODO change
        self.sigma = math.radians(sigma)

        self.ballisticCoeff = self.computeBallisticCoeff()

    def computeBallisticCoeff(self):
        # TODO fix this
        return 200
    
if __name__ == "__main__":
    pass