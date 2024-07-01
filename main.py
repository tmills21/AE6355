from vehicleGeometry import vehicle
from RK4sim import RK4withLayeredAtmosphere

stardust = vehicle(46, 0.8128, 0.2202, 60, 0)
test = RK4withLayeredAtmosphere(stardust, 11067.63087, -10, 120000)
print(test.runSim())