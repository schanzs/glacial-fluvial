backgroundU = 0.5/1000 # in m/yr
D0 = 200/1000 # in m
a = 0.02/1000
erosion_depth_threshold = 0.2
glacial_discharge_sw = True
glacial_sed_supply_sw = True
isostacy_sw = True
dt = 5
# glacier timestep needs to be 0.01 to capture annual ELA variability and slow growth
dt_g = 0.01 
dt_i = 10.
dx = 2000
nodes = 100
initial_slope = 0.02
initial_sed_depth = 0.1