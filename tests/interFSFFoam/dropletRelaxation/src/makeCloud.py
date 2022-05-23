import os
import math

d_ = 0.000025
step_ = 1.2*d_
dx_ = 0.02
dz_ = 0.000125

cloudFile = open('cloudPositions', 'w')

cloudFile.write("\
FoamFile\n\
{\n\
    format      ascii;\n\
    class       vectorField;\n\
    object      cloudPositions;\n\
}\n\
\n\
(\n\
")

for j in range(int(dz_/step_) - 1):
    for i in range(int(dx_/step_) - 1):
            cloudFile.write("(" + str(step_*i - 0.5*dx_ + step_) + " " + str(2.5e-5) + " " + str(step_*(j + 1)) + ")\n")

cloudFile.write(")")
