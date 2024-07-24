class TPS:
    def __init__(self, beta, vel, gam):
        self.beta = beta
        self.vel = vel
        self.gam = gam
        self.targetTemp = 523 # Kelvin

    def TPSequation(self, vel, gam, beta, thick):
        term1 = ( gam - (-6.5) ) / 3.5
        term2 = ( vel - 10.45) / 2.55
        term3 = ( beta - 60 ) / 40
        term4 = ( thick - 6 ) / 4 

        blTemp = 353.71725424 + 60.407666667 * term1 + 25.948222222 * term2+ 112.06988889 * term3 + \
            -170.5173333 * term4 + term1 * ( term2 * -80.552375 ) + term1 * ( term3 * 2.053875 ) + \
            term2 * ( term3 * 21.202625 ) + term1 * ( term4 * -39.463625 ) + term2 * ( term4 * -42.249875 ) + \
            term3 * ( term4 * 7.366375 ) + term1 * ( term1 * 48.46320339 ) + term2 * ( term2 * 31.12820339) + \
            term3 * ( term3 * 29.82320339 ) + term4 * ( term4 * 158.84820339 )
        
        return blTemp
    
    def computeThickness(self):
        start = {}
        lower = 2
        upper = 10

        start[self.TPSequation(self.vel, self.gam, self.beta, lower)] = lower
        start[self.TPSequation(self.vel, self.gam, self.beta, upper)] = upper
        
        lower = start[min(start.keys())]
        upper = start[max(start.keys())]

        diff = 100
        tol = 10**-4

        while (diff > tol):
            next = ( upper + lower ) / 2
            val = self.TPSequation(self.vel, self.gam, self.beta, next)
            diff = abs(self.targetTemp - val)

            if val > self.targetTemp:
                upper = next
            else:
                lower = next

        return next

if __name__ == "__main__":
    test = TPS(80, 8.5, -5)
    print(test.computeThickness())