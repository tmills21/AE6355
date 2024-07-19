import numpy as np

class vehicle:
    def __init__(self, mass, coneDiameter, noseRadius, deltac, sigma):
        self.mass = mass # kg
        self.coneDiamter = coneDiameter # m
        self.coneRadius = self.coneDiamter / 2 # m
        self.nodeRadius = noseRadius # m
        self.deltac = np.radians(deltac) # degrees to radians
        self.sigma = np.radians(sigma) # degrees to radians

        self.CL = 0 # assume purely ballistic
        self.LDrat = 0 # assume purely ballistic
        self.angleOfAttack = 0 # assume purely ballistic

        self.cpmax = 2.0 # assume Newtonian flow

        self.Aref = self.computeAref()
        self.CN = self.computeCN()
        self.CA = self.computeCA()
        self.CD = self.computeCD()
        self.ballisticCoeff = self.computeBallisticCoeff()

    def computeAref(self):
        return np.pi * max(self.nodeRadius, self.coneRadius)**2
    
    def computeCN(self):
        
        # assuming singlular truncated cone blunt body
        term1 = 1 - ( self.nodeRadius / self.coneRadius )**2 * np.cos(self.deltac) * np.cos(self.deltac)
        term2 = np.cos(self.deltac) * np.cos(self.deltac) * np.sin(self.angleOfAttack) * np.cos(self.angleOfAttack)
        return self.cpmax * term1 * term2

    def computeCA(self):
        
        # assuming singlular truncated cone blunt body
        term1 = 0.5 * ( 1 - (np.sin(self.deltac))**4 ) * ( self.nodeRadius / self.coneRadius )**2
        term2 = np.sin(self.deltac) * np.sin(self.deltac) * np.cos(self.angleOfAttack) * np.cos(self.angleOfAttack) + \
                0.5 * np.cos(self.deltac) * np.cos(self.deltac) * np.sin(self.angleOfAttack) * np.sin(self.angleOfAttack)
        term3 = 1 - ( self.nodeRadius / self.coneRadius )**2 * np.cos(self.deltac) * np.cos(self.deltac)
        return self.cpmax * ( term1 + term2 * term3 )
    
    def computeCD(self):
        return self.CN * np.sin(self.angleOfAttack) + self.CA * np.cos(self.angleOfAttack)

    def computeBallisticCoeff(self):
        return self.mass / ( self.CD * self.Aref )
    
if __name__ == "__main__":
    stardust = vehicle(46, 0.8128, 0.2202, 60, 0)