
""" Variables """
backgroundU = 0.01/1000          # rock uplift rate, in m/yr
D0 = 0.1                        # initial median grain size, in meters

""" Switches to change model type """
glacial_discharge_sw = True     # Boolean switch to turn glacial meltwater input to fluvial system on or off
glacial_sed_supply_sw = True    # Boolean switch to turn glacially eroded bedload supply to fluvial system on or off
isostacy_sw = True             # Boolean switch to turn isostatic calculations on or off
sediment_transport = True      # Boolean switch to turn sediment transport calculations on/off

""" Time """
dt = 1                          # model timestep, years
dt_g = 0.01                     # model timestep for glacial erosion. Through trial and error, found 0.01 is good to capture annual ELA variability, years
dt_i = 10.                      # model timestep to calculate isostatic adjustments, years

""" Model set up """
dx = 2000                       # distance between model nodes, meters. Found 2000 keeps sediment transport stable for all grain sizes modelled
nodes = 100                     # number of nodes in the model
initial_slope = 0.02            # initial slope of model system. Overriden in block uplift scenario
initial_sed_depth = 0.1         # initial depth of sediment, meters
a = 0.02/1000                   # attrition rate, Sternberg's Law
erosion_depth_threshold = 0.2   # thickness of sediment, m, that prevents bedrock erosion in cover effect