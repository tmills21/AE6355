from vehicleGeometry import vehicle
from RK4sim import RK4withLayeredAtmosphere

import numpy as np
import math

class deorbit:
    def __init__(self, entryVehicle, RK4, ecc1, hp1, nmax, vatm):
        self.entryVehicle = entryVehicle # vehicle object
        self.RK4 = RK4 # RK4withLayeredAtmosphere object
        self.ecc1 = ecc1
        self.hp1 = hp1 # meters
        self.nmax = - abs(nmax) # g's
        self.vatm = vatm # m/s
        self.ue = self.takeDimensionOut(self.vatm)

        self.a1 = self.computea(self.ecc1, self.hp1) # meters
        self.p1 = self.computep(self.a1, self.ecc1) # meters
        self.alpha1 = self.computealpha(self.a1)

        # base values using thetaD = 180 degrees
        self.thetaD = math.radians(180) # degrees
        self.rD = self.computerD(self.p1, self.ecc1) # meters
        self.relambda = self.computelambda(self.rD)

        # get a ue that is 80% of the ue returned from equation (76)
        if self.vatm == str():
            self.computeEntryVels(0.99)  

    def computea(self, ecc, hp):
        return ( self.RK4.Re + hp ) / ( 1 - ecc )
    
    def computep(self, a, ecc):
        return a * ( 1 - ecc**2 )
    
    def computerD(self, p, ecc):
        return p / ( 1 + ecc * math.cos(self.thetaD) )
    
    def computelambda(self, rD):
        return rD / (self.RK4.Re + self.RK4.hatm)
    
    def computealpha(self, a):
        return a / (self.RK4.Re + self.RK4.hatm)
    
    def computeEntryVels(self, scale):
        self.ue = scale * self.computeMinUe()
        self.vatm = self.addDimensionBack(self.ue)
    
    def computeMinUe(self):
        numer = 2 * self.relambda * ( self.relambda - 1 ) * ( 1 - self.ecc1**2 ) * self.alpha1**2
        denom = self.relambda * ( 1 - self.ecc1**2 ) * self.alpha1**2 - 2 * self.alpha1 + self.relambda
        return math.sqrt( numer / denom )
    
    def addDimensionBack(self, val):
        scale = math.sqrt( self.RK4.mu / self.RK4.Re )
        return val * scale
    
    def takeDimensionOut(self, val):
        ans = 0
        if val != str():
            scale = math.sqrt( self.RK4.mu / self.RK4.Re )
            ans = self.vatm / scale
        return ans
    
    def computeu1(self):
        numer = 2 * self.alpha1 - self.relambda
        return math.sqrt(numer / self.alpha1)
    
    def computeu2(self):
        discrim = self.ue**2 + 2 * ( 1 - self.relambda )
        return math.sqrt(discrim)
    
    def computeDeltau(self):
        u1 = self.computeu1()
        u2 = self.computeu2()
        return math.sqrt(u1**2 + u2**2 - 2* u1 * u2 * math.cos(0))
    
    def computeOvershoot(self):
        vel = self.vatm # m/s
        lower_gammae = -0.5 # degrees
        upper_gammae = -50 # degrees

        gammaDif = 100
        prevGammae = lower_gammae

        while (gammaDif > 10**-4):
            model = RK4withLayeredAtmosphere(self.entryVehicle, vel, prevGammae, self.RK4.hatm)
            capure = model.runSim(abortOnSkip = True)[0]
            if capure:
                upper_gammae = prevGammae
            else:
                lower_gammae = prevGammae
            nextGammae = ( upper_gammae + lower_gammae ) / 2
            gammaDif = abs(prevGammae - nextGammae)
            prevGammae = nextGammae
            
        print('overshoot found, gamma = ' + str(prevGammae)) 

        # return deorbit position, delta v, ve, and gammae
        return prevGammae
    
    def AllenEggersGammaFromNmax(self, vel):

        # equation (93)
        return math.asin( ( 53.3 * self.RK4.H * self.nmax) / vel**2 ) # radians

    def computeUndershoot(self):
        vel = self.vatm # m/s
        startgamma = math.degrees(self.AllenEggersGammaFromNmax(vel))
        gammaincrement = -0.1

        model = RK4withLayeredAtmosphere(self.entryVehicle, vel, startgamma, self.RK4.hatm)
        ach = model.runSim(abortOnNmax = True, nmaxLim = self.nmax)[0]

        otherBoundFound = False
        if ach:
            while not otherBoundFound:
                startgamma += gammaincrement
                model = RK4withLayeredAtmosphere(self.entryVehicle, vel, startgamma, self.RK4.hatm)
                otherBoundFound = model.runSim(abortOnNmax = True, nmaxLim = self.nmax)[0]

            good_gammae = startgamma
            bad_gammae = startgamma - gammaincrement

        else:
            while not otherBoundFound:
                startgamma -= gammaincrement
                model = RK4withLayeredAtmosphere(self.entryVehicle, vel, startgamma, self.RK4.hatm)
                otherBoundFound = model.runSim(abortOnNmax = True, nmaxLim = self.nmax)[0]

            good_gammae = startgamma
            bad_gammae = startgamma + gammaincrement

        gammaDif = 100
        prevGammae = ( good_gammae + bad_gammae ) / 2

        while (gammaDif > 10**-4):
            model = RK4withLayeredAtmosphere(self.entryVehicle, vel, prevGammae, self.RK4.hatm)
            ach = model.runSim(abortOnNmax = True, nmaxLim = self.nmax)[0]
            if ach:
                good_gammae = prevGammae
            else:
                bad_gammae = prevGammae
            nextGammae = ( good_gammae + bad_gammae ) / 2
            gammaDif = abs(prevGammae - nextGammae)
            prevGammae = nextGammae
            
        print('undershoot found, gamma = ' + str(prevGammae)) 

        # return deorbit position, delta v, ve, and gammae
        return prevGammae
    
    def gammaeForMinDv(self):
        numer = self.relambda * ( 1 - self.ecc1 **2) * ( self.ue**2 + 2 * ( 1 - self.relambda ) )
        denom = 2 * self.alpha1 - self.relambda
        arccosArg = self.alpha1 / self.ue * math.sqrt(numer / denom)
        if arccosArg > 1 or arccosArg < -1:
            return None
        
        gamrad = math.acos( arccosArg )
        gamdeg = - math.degrees(gamrad)
        print('min dv found, gamma = ' + str(gamdeg)) 

        return gamdeg
    
    def validCorridor(self, min, max, toTest):
        return (min <= toTest <= max)
    
    def computeCorridor(self):
        valid = False
        overshootGamma = self.computeOvershoot()
        undershootGamma = self.computeUndershoot()

        minDVgamma = self.gammaeForMinDv()
        if minDVgamma is not None:
            valid = self.validCorridor(undershootGamma, overshootGamma, minDVgamma)
            if valid:
                deltau = self.computeDeltau()
                deltav = self.addDimensionBack(deltau)
                return [self.vatm, minDVgamma, self.rD, deltav, undershootGamma, overshootGamma]
                # return [undershootGamma, minDVgamma, overshootGamma, float(math.degrees(self.thetaD)), float(self.vatm)]

        thetaDoptions = np.linspace(190, 270, 20)

        # loop over thetaD values
        if not valid:
            for thetaD in thetaDoptions:
                print(thetaD, 0.99)
                self.thetaD = math.radians(thetaD) # degrees
                self.rD = self.computerD(self.p1, self.ecc1) # meters
                self.relambda = self.computelambda(self.rD)
                self.computeEntryVels(0.99) 

                overshootGamma = self.computeOvershoot()
                undershootGamma = self.computeUndershoot()

                minDVgamma = self.gammaeForMinDv()
                if minDVgamma is not None:
                    valid = self.validCorridor(undershootGamma, overshootGamma, minDVgamma)
                if valid:
                    deltau = self.computeDeltau()
                    deltav = self.addDimensionBack(deltau)
                    # return [undershootGamma, minDVgamma, overshootGamma, float(thetaD), float(self.vatm)]
                    return [self.vatm, minDVgamma, self.rD, deltav, undershootGamma, overshootGamma]
        
        velocityScaleOptions = np.linspace(0.98, 0.9, 9)

        # if still not valid, loop over different ue's too
        if not valid:
            for thetaD in thetaDoptions:
                for scaling in velocityScaleOptions:
                    print(thetaD, scaling)
                    self.thetaD = math.radians(thetaD) # degrees
                    self.rD = self.computerD(self.p1, self.ecc1) # meters
                    self.relambda = self.computelambda(self.rD)
                    self.computeEntryVels(scaling) 

                    overshootGamma = self.computeOvershoot()
                    undershootGamma = self.computeUndershoot()

                    minDVgamma = self.gammaeForMinDv()
                    if minDVgamma is not None:
                        valid = self.validCorridor(undershootGamma, overshootGamma, minDVgamma)

                    if valid:
                        deltau = self.computeDeltau()
                        deltav = self.addDimensionBack(deltau)
                        return [self.vatm, minDVgamma, self.rD, deltav, undershootGamma, overshootGamma]
                        # return [undershootGamma, minDVgamma, overshootGamma, float(thetaD), float(self.vatm)]

        return [0, 0, 0, 0, 0]





if __name__ == "__main__":
    stardust = vehicle(46, 0.8128, 0.2202, 60, 0)
    RK = RK4withLayeredAtmosphere(stardust, 11067.63087, -10, 120000)

    orb = deorbit(stardust, RK, 0.1, 400000, 30, '')
    # print(orb.validCorridor(-7.713486509983248, -3.451846122741699, -7.670675341282269))
    print(orb.computeCorridor())

