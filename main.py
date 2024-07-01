from vehicleGeometry import vehicle
from RK4sim import RK4withLayeredAtmosphere

stardust = vehicle(46, 0.8128, 0.2202, 60, 0)
test = RK4withLayeredAtmosphere(stardust, 11000, -10, 120000)
print(test.dVdt(10000, -1, 110000))