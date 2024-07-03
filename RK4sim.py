import numpy as np
import matplotlib.pyplot as plt

from vehicleGeometry import vehicle

class RK4withLayeredAtmosphere:
    def __init__(self, vehicleObject, v0, gamma0, h0):
        self.vehicleObject = vehicleObject

        self.v0 = v0 # m/s
        self.gamma0 = np.radians(gamma0) # degrees to radians
        self.h0 = h0 # m

        self.timeStep = 0.1 # s
        self.H = 7250 # m for Earth
        self.Re = 6378130 # m
        self.mu = 398600441800000 # m^3/s^2
        self.g0 = 9.806 # m/s^2
        self.g0p = 9.806 # m^2/s^2-m'
        self.Rbar = 8314.32 # N-m/kmol-K
        self.R = 287 # J/kg-K
        self.MW0 = 28.9644 # kg/kmol

    def dVdt(self, vel, gam, alt):
        val = - ( self.density(alt) * vel ** 2 ) / ( 2 * self.vehicleObject.ballisticCoeff ) - self.gravity(alt) * np.sin(gam)
        return val

    def dGamdt(self, vel, gam, alt):
        term1 = (vel * np.cos(gam)) / (self.Re + alt) 
        term2 = (self.density(alt) * vel) / (2 * self.vehicleObject.ballisticCoeff) * self.vehicleObject.LDrat * np.cos(self.vehicleObject.sigma) 
        term3 = - (self.gravity(alt) * np.cos(gam)) / vel
        return term1 + term2 + term3

    def dhdt(self, vel, gam):
        return vel * np.sin(gam)

    def density(self, alt):

        # # exponential atmosphere
        # rho = 1.225 * np.exp(-alt / self.H)
        # return rho

        # convert altitude from m to km
        altkm = alt / 10**3

        # TODO check tops of intervals for 0-86
        # TODO return an error if bad input?
        layers = [(0, 11.0102), (11.0102, 20.0631), (20.0631, 32.1619), (32.1619, 47.3501), (47.3501, 51.4125), 
                  (51.4125, 71.802), (71.802, 86), (86, 100), (100, 110), (110, 120), (120, 150)]
        def getLayer(layers, altkm):
            for i, (bottom, top) in enumerate(layers):
                if bottom <= abs(altkm) <= top:
                    return i
                
        layer = getLayer(layers, altkm)
        Li = [-6.5, 0, 1, 2.8, 0, -2.8, -2, 1.6481, 5, 10, 20] # K/km' or K/km
        Tmi = [288.15, 216.65, 216.65, 228.65, 270.65, 270.65, 214.65, 186.945, 210.65, 260.65, 360.65] # K
        
        # TODO check layer 8
        Pi = [101325, 22631.95, 5474.79, 868.01, 110.9, 66.94, 3.956, 0.34418, 0.029073, 0.006801, 0.002247] # N/m^2

        li = Li[layer]
        tmi = Tmi[layer]
        pi = Pi[layer]

        layerFloor = np.floor(layers[layer][0])

        # 1976 US Standard Atmosphere Model
        if layer <= 6:

            # convert altitude from m to m'
            altmprime = ( self.Re * alt ) / ( self.Re + alt )

            # convert altitude from m' to km'
            altkmprime = altmprime / 10**3

            tm = tmi + li * ( altkmprime - layerFloor ) 

            if li:
                base = tmi / tm
                exponent = ( self.g0p * self.MW0 ) / ( self.Rbar * 10**-3 * li )
                p = pi * base**exponent

            else:
                exponent = ( - self.g0p * self.MW0 * ( altmprime - layerFloor * 10**3 ) ) / ( self.Rbar * tmi )
                p = pi * np.exp(exponent)

            rho = ( p * self.MW0 ) / ( self.Rbar * tm )

        # 1962 Standard Atmosphere
        else:
            b = 3.31 * 10**-7 # 1/m

            # # pressure
            # power = - ( ( self.g0 / ( self.R * li * 10**-3 ) ) * ( 1 + b * ( tmi / ( li * 10**-3 ) - layerFloor * 10**3 ) ) )
            # exponential = ( ( self.g0 * b ) / ( self.R * li * 10**-3 ) * ( alt - layerFloor * 10**3 ) )
            # tm = tmi + li * ( altkm - layerFloor )
            # p = pi * ( tm / tmi )**power * np.exp(exponential)

            T = [186.946, 210.02, 257, 349.49]
            MWi = [28.9644, 28.88, 28.56, 28.08]
            t = T[layer - 7]
            mwi = MWi[layer - 7]

            # TODO check if this is right
            rhoi = ( pi * mwi ) / ( self.Rbar * t )

            # TODO is this R or Rbar
            power = - ( ( self.g0 / ( self.R * li * 10**-3 ) ) * ( ( ( self.R * li * 10**-3 ) / self.g0 ) + 1 + b * ( tmi / ( li * 10**-3 ) - layerFloor * 10**3 ) ) )
            exponential = ( ( self.g0 * b ) / ( self.R * li * 10**-3 ) * ( alt - layerFloor * 10**3 ) )
            tm = tmi + li * ( altkm - layerFloor )
            rho = rhoi * ( tm / tmi )**power * np.exp(exponential)

        return rho

    def gravity(self, alt):

        # inverse square gravity law
        return self.g0 * ( self.Re / ( self.Re + alt ) )**2
    
    def runSim(self):
        v = self.v0
        gam = self.gamma0
        h = self.h0

        vHistory = np.array([v])
        gamHistory = np.array([gam])
        hHistory = np.array([h])

        while ( h > 0 ):
            k1 = self.dVdt(v, gam, h)
            L1 = self.dGamdt(v, gam, h)
            m1 = self.dhdt(v, gam)

            k2 = self.dVdt(v + self.timeStep / 2 * k1, gam + self.timeStep / 2 * L1, h + self.timeStep / 2 * m1)
            L2 = self.dGamdt(v + self.timeStep / 2 * k1, gam + self.timeStep / 2 * L1, h + self.timeStep / 2 * m1)
            m2 = self.dhdt(v + self.timeStep / 2 * k1, gam + self.timeStep / 2 * L1)

            k3 = self.dVdt(v + self.timeStep / 2 * k2, gam + self.timeStep / 2 * L2, h + self.timeStep / 2 * m2)
            L3 = self.dGamdt(v + self.timeStep / 2 * k2, gam + self.timeStep / 2 * L2, h + self.timeStep / 2 * m2)
            m3 = self.dhdt(v + self.timeStep / 2 * k2, gam + self.timeStep / 2 * L2)

            k4 = self.dVdt(v + self.timeStep * k3, gam + self.timeStep * L3, h + self.timeStep * m3)
            L4 = self.dGamdt(v + self.timeStep * k3, gam + self.timeStep * L3, h + self.timeStep * m3)
            m4 = self.dhdt(v + self.timeStep * k3, gam + self.timeStep * L3)

            T4A = k1 + 2 * k2 + 2 * k3 + k4
            T4B = L1 + 2 * L2 + 2 * L3 + L4
            T4C = m1 + 2 * m2 + 2 * m3 + m4

            v += self.timeStep / 6 * T4A
            gam += self.timeStep / 6 * T4B
            h += self.timeStep / 6 * T4C

            vHistory = np.append(vHistory, v)
            gamHistory = np.append(gamHistory, gam)
            hHistory = np.append(hHistory, h)

        plt.plot(vHistory, hHistory)
        plt.xlabel('Velocity (m/s)')
        plt.ylabel('Altitude (m)')
        plt.show()

        return 'done'

    
if __name__ == "__main__":
    stardust = vehicle(46, 0.8128, 0.2202, 60, 0)
    test = RK4withLayeredAtmosphere(stardust, 11067.63087, -10, 120000)

    altitude = np.arange(0, 150000, 1)
    density = np.zeros_like(altitude)

    for i, alt in enumerate(altitude):
        density[i] = test.density(alt)

    plt.plot(density, altitude)
    plt.xlabel('density')
    plt.ylabel('Altitude (m)')
    plt.show()

    # print(test.density(49000))