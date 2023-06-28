# Bearing Capacity and Settlement of Shallow Foundations
Produces spread footing design chart using elastic theory method from Bowle's Foundation Analysis and Design, 5th ed.

This code is presently only applicable to computing settlement at the centerpoint of the foundation

Settlement is estimated from the following equation:

Settlement:           ΔH = 𝑞 ⋅ B'⋅ (1-ν²)/Es ⋅ m ⋅ 𝐼𝑠 ⋅ 𝐼𝑓


Instructions for use:
In "design_chart_imperial.py":
Specify the foundation embedment depth, "D", and shape ratio (B/L) as "Shape" for the shallow foundation.

Next, specify the backfill and foundation parameters:

gamma_backfill = unit weight of structural backfill in pcf
gamma_soil = unit weight of the foundation soil in pcf
Z = depth to hard layer in feet
phi = foundation soil friction angle in degrees
c = foundation soil cohesion in psf
nu = Poisson's ratio of foundation soil.

If you want to read in a file with layer information, set modulus_file=True.
The csv file format should contain a row for each layer and have columns named "layer_number", "top", "bottom", "modulus"
The first layer "top" is at the base of the foundation.

If modulus_file=False, the constant Es value will be used.
