# RanGES: River and Glacial Erosion with Sediment

### Overview
This 1-D finite difference model simulates erosion by a mixed alluviul bedrock river and an overriding glacier. The two systems feedback in the form of 1) valley widening and sediment storage; 2) glacial supply to sediment loads; 3) meltwater; and 4) increased bedload size. The model couples two previously distinct erosion regimes, and is intended for geomorphologists and geoscientists interested in the development of topography during glacial-interglacial cycles.

### Installation instructions
To use this software, download the .py and .sh files. The code was developed using Python3, and can be run as a Python script in your preferred method. Inputs are user-defined when driver.py is called. The bash file is set up to run multiple simulations by iterating through input lists. 

### Example usage
To run the script with multiple simulations, update the desired simulations in the bash script, then run the bash script by calling your preferred method, for example:
$ bash run_glacial_fluvial.sh. 

To run the script with only one simulation (in this example, uplift = 1 mm/yr, initial grain size = 5 cm, attrition rate is 0.0001 m/m, and glaciation is on), run the python file by calling: 
$ python driver.py 0.001 0.05 0.0001 True

### Contributing and bug reports
Submit any bug reports or contributions to our issue tracker (https://github.com/schanzs/glacial-fluvial/issues). For substantial changes, please fork this repository and submit a pull request for review.

### Functionality
The model works through:
1) creating an initial fluvial landscape in part 1 of driver.py. This uses the inputs to simulate fluvial erosion over 5 My, starting with a 3000 m high plateau. Outputs a set of text files that can be used as input in further runs and an image of the 5 My long profile. Will automatically continue to step 2 unless otherwise edited. 
2) applying glacial and fluvial erosion to your landscape in part 2 of driver.py. The default runs for 800,000 years with 100,000 glacial cycles. This will output channel geometry and erosion every 1000 years, saved within the script as a dictionary and exported as a single Pandas pickle file.

Driver.py organizes the model, but erosion components are defined in glacial_fluvial_functions.py, which defines an object class used in driver.py.

Further details on the erosion components are in Schanz and Yanites, in submission at Journal of Geophysical Research - Earth Surface. 


