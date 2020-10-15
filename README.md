# glacial-fluvial
## clever model names? long profile erosion by glaciers, streams and sediment (e, g, st, se?)

the following sections correspond to review criteria from JOSS.

### Statement of need
This 1-D finite difference model simulates erosion of a mixed alluviul bedrock river by an overriding glacier and stream power, with feedbacks between the two erosional regimes in the form of valley widening and sediment storage; glacial supply to sediment loads; and meltwater. The model couples two previously distinct erosion regimes, and is intended for geomorphologists and geoscientists interested in the development of topography during glacial-interglacial cycles.

### Installation instructions
To use this software, download the .py files. The code was developed using Python3, and can be run as a Python script in your preferred method. Inputs are specified in inputs.py. If you do not have an initial topography, use spin_up_run.py to create one. Finally, analysis_run.py is run to apply glaciation onto your input topography.

### Example usage


### API documentation

### Contributing and bug reports
Submit any bug reports or contributions to our issue tracker (https://github.com/schanzs/glacial-fluvial/issues). For substantial changes, please fork this repository and submit a pull request for review.

### Functionality & tests



The model works through:
1) altering inputs in inputs.py that include model set up in nodes, dx, dt, uplift rate, grain size, etc as well as switches to turn on/off sediment supply, glacial-fluvial feedbacks, and isostatic adjustments.
2) creating an initial fluvial landscape in spin_up_run.py. This uses the inputs to simulate fluvial erosion over 5 My, starting with a 3000 m high plateau. Outputs a set of text files that can be used as input in the glacial function - or you can create your own set of text files to use as input to the next step.
3) applying glacial and fluvial erosion to your landscape in analysis_run.py. The default runs for 800,000 years with 100,000 glacial cycles. This will output channel geometry and erosion every 1000 years as text files, saved in a new subdirectory.

The erosion components are explained below. Briefly, fluvial bedrock erosion occurs via shear stress with a threshold sediment depth to stop erosion. Sediment is transported with a modified Meyer-Peter Mueller equation. Glacial bedrock erosion is driven by quarrying and is a power law relation to the basal sliding velocity.

## Fluvial erosion.
Bedrock erosion is proportional to shear stress (Howard & Kerby, 1983; Whipple and Tucker, 1999) as:

E=k_b τ_b^a

where a is equal to 1, corresponding to n = 0.7 in the stream power equation, k_b is a constant, and basal shear stress, τ_b, is the product of water density, gravitational acceleration, water depth, and channel slope.

A sediment depth threshold is used to imitate a sediment cover effect, and can be manually set in inputs.py.

## Sediment transport.
The transport capacity of sediment is set using a modified Meyer-Peter and Muller (1948) equation from Wong and Parker (2006) in which the transport capacity is a function of basal shear stress, critical shear stress, grain diameter, and grain density:

Q_t=3.97ρ_s W[τ_b/(ρ_s-ρ)gD- τ_c^* ]^(3/2) D^(3/2) √((ρ_s-ρ)/ρ g)

Grain size decreases downstream following Sternberg's (1875) law:

D=D_0 e^(-a_s x)

Sediment mass is conserved through Exner's equation (Exner, 1920; 1925).

## Glacial erosion.
The glacial erosion component is based off MacGregor et al (2000) with quarrying-law updates from Iverson (2012) that sets erosion as dependent on sliding velocity:

E_q=C_1 U_s^0.46   	


