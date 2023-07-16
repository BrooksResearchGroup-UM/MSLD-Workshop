#! /usr/bin/env python

import sys

#! ! Set nonbonded options
#! nbonds atom -
#! fswitch -                !! (use "switch" for openMM) Use force switch for elec-nb interactions
#! cdie eps 1 -             !! Set constant dielectric value of 1
#! vdw vfswitch -            !! Set nb VDW using switching function
#! cutnb 15.0 -             !! Set nb list cut off to 15A
#! cutim 15.0 -             !! Set image nb list cut off to 15A
#! ctonnb 10.0 ctofnb 12.0   !! Switching function used from 10A to 12A for VDW interaction

#! Set nonbonded options
#nbonds atom -
# - ! fswitch -                !! Use force switch for elec-nb interactions
#cdie eps 1 -             !! Set constant dielectric value of 1
#vdw vfswitch -            !! Set nb VDW using switching function
#cutnb 12.0 -             !! Was 15.0 Set nb list cut off to 15A
#cutim 12.0 -             !! Was 15.0 Set image nb list cut off to 15A
#ctonnb 9.0 ctofnb 10.0 - !! Were 10.0 and 12.0 Switching function used from 10A to 12A for VDW interaction
#Ewald -                  !! Use Ewald electrostatics
#kappa 0.320 -            !! 0.320 makes direct space decay to 10^-5 at 10 A
#pmEwald -                !! Keyword
#order 6 -                !! Interpolation order
#fftx @pmegrid ffty @pmegrid fftz @pmegrid  !! Choose approximate box length

x=float(sys.argv[1])
allowed=[16,18,20,24,27]
x=x/1.1
factor=1
while x>2*allowed[0]:
  x=x/2
  factor=2*factor
base=2*allowed[0]
for i in range(4,-1,-1):
  if x<allowed[i]:
    base=allowed[i]

print("set pmegrid = %d" % (base*factor,))
print("nbonds atom -")
print("cdie eps 1 -") # Set constant dielectric value of 1
print("vdw vfswitch -") # Set nb VDW using switching function
print("cutnb 12.0 -") # Was 15.0 Set nb list cut off to 15A
print("cutim 12.0 -") # Was 15.0 Set image nb list cut off to 15A
print("ctonnb 9.0 ctofnb 10.0 -") # Were 10.0 and 12.0 Switching function used from 10A to 12A for VDW interaction
print("Ewald -") # Use Ewald electrostatics
print("kappa 0.320 -") # 0.320 makes direct space decay to 10^-5 at 10 A
print("pmEwald -") # Keyword
print("order 6 -") # Interpolation order
print("fftx @pmegrid ffty @pmegrid fftz @pmegrid") # Choose approximate box length
