import numpy as np
import math
import matplotlib.pyplot as plt

from vehicleGeometry import vehicle

class RK4Planar:
    def __init__(self, vehicleObject, v0, gamma0, h0):
        self.vehicleObject = vehicleObject

        self.v0 = v0 # m/s
        self.gamma0 = np.radians(gamma0) # degrees to radians
        self.h0 = h0 # m

        self.timeStep = 0.1 # s

        # for Earth
        self.H = 7250 # m
        self.Re = 6378130 # m
        self.hatm = 120000 # m
        self.mu = 398600441800000 # m^3/s^2
        self.g0 = 9.806 # m/s^2
        self.g0p = 9.806 # m^2/s^2-m'
        self.Rbar = 8314.32 # N-m/kmol-K
        self.R = 287 # J/kg-K
        self.MW0 = 28.9644 # kg/kmol
        self.K = 1.74153 * 10**-4

    def dVdt(self, vel, gam, alt):
        val = - ( self.density(alt) * vel ** 2 ) / ( 2 * self.vehicleObject.ballisticCoeff ) - self.gravity(alt) * np.sin(gam)
        return val

    def dGamdt(self, vel, gam, alt):
        term1 = ( vel * np.cos(gam) ) / ( self.Re + alt ) 
        term2 = ( self.density(alt) * vel ) / ( 2 * self.vehicleObject.ballisticCoeff ) * self.vehicleObject.LDrat * np.cos(self.vehicleObject.sigma) 
        term3 = - ( self.gravity(alt) * np.cos(gam) ) / vel
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
    
    def runSim(self, abortOnSkip = False, abortOnNmax = False, nmaxLim = 0, title = ''):
        v = self.v0
        gam = self.gamma0
        h = self.h0

        vHistory = np.array([v])
        gamHistory = np.array([gam])
        hHistory = np.array([h])
        decelHistory = np.array([])

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

            if abortOnSkip:
                if gam > 0:
                    return [False]

            vHistory = np.append(vHistory, v)
            gamHistory = np.append(gamHistory, gam)
            hHistory = np.append(hHistory, h)

            decel = ( vHistory[-1] - vHistory[-2] ) / self.timeStep
            decelGs = decel / self.g0
            decelHistory = np.append(decelHistory, decelGs)

            if abortOnNmax and (decelGs < nmaxLim):
                return [False]

        if abortOnSkip or abortOnNmax:
            return [True]

        else:
            fig = plt.figure()
            ax = fig.add_subplot(1, 1, 1)
            ax.plot(vHistory, hHistory)
            ax.set_xlabel('Velocity (m/s)')
            ax.set_ylabel('Altitude (m)')
            ax.set_title(title)
            # plt.show()

            return [vHistory, gamHistory, hHistory, decelHistory, fig]
        
    def SuttonGravesHeat(self, vHistory, hHistory):
        length = len(vHistory)
        qstag = [0] * length
    
        # equation (128)
        for i in range(length):
            qstag[i] = self.K * math.sqrt( self.density(hHistory[i]) / self.vehicleObject.noseRadius ) * vHistory[i]**3 / 10000; # W/cm^2
        
        tHistory = [0 + i * self.timeStep for i in range(length)]
        fig = plt.figure()
        ax = fig.add_subplot(1, 1, 1)
        ax.plot(tHistory, qstag)
        ax.set_xlabel('Time (s)')
        ax.set_ylabel('Convective Stagnation-Point Heat Rate ' + r'$\text{(W/cm}^2 \text{)}$')
        ax.set_title('Aerodynamic Heat Rate')
        
        return [qstag, fig]
    
    def integratedHeat(self, qstag):
        length = len(qstag)
        qtotal = [0] * length

        for i in range(length-1):
            area = (self.timeStep)/2 * (qstag[i+1] + qstag[i]) # J/cm^2
            qtotal[i+1] = qtotal[i] + area
    
        tHistory = [0 + i * self.timeStep for i in range(length)]
        fig = plt.figure()
        ax = fig.add_subplot(1, 1, 1)
        ax.plot(tHistory, qtotal)
        ax.set_xlabel('Time (s)')
        ax.set_ylabel('Stagnation-Point Integrated-Heat Rate ' + r'$\text{(J/cm}^2 \text{)}$')
        ax.set_title('Total Heat Rate')
        
        return [qtotal, fig]


    
if __name__ == "__main__":
    stardust = vehicle(46, 0.8128, 0.2202, 60, 0, 0, '')
    test = RK4Planar(stardust, 8977.441267279422, -8.270241244214445, 120000)
    simulation_result = test.runSim(abortOnNmax = True, nmaxLim = -30)

    # altitude = np.arange(3000, 8000, 1)
    # density = np.zeros_like(altitude)

    # for i, alt in enumerate(altitude):
    #     density[i] = test.density(alt)

    # plt.plot(density, altitude)
    # plt.xlabel('density')
    # plt.ylabel('Altitude (m)')
    # plt.show()

    # print(test.density(49000))