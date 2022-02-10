import hudson_model
import numpy as np

print("\n\nRun model for multiple porosities and corresponding aspect ratios")
print("#######################################################################")
print("\n\n\n")
#set global variables for anisotropic solid material
Cagg = np.array([[114.966597967223, 29.5984716245580, 25.3827852191623, 0, 0, 0],
                 [29.5984716245580, 114.966597967223, 25.3827852191623, 0, 0, 0],
                 [25.3827852191623, 25.3827852191623, 99.4844041646496, 0, 0, 0],
                 [0, 0, 0, 34.3549850499644, 0, 0],
                 [0, 0, 0, 0, 34.3549850499644, 0],
                 [0, 0, 0, 0, 0, 42.6840631713325]])

#Aggregate bulk density, bulk modulus and shear modulus
Ragg = 2644.22030828962
Kiso, Giso = 54.276991245977946*1e9, 38.642892204430370*1e9

#Mean measured values from triax experiments. Peak vp & vs at 90 MPa
#data = np.array([5608.50,5013.53,3032.33,3066.00,2643.38])
data = np.array([5400,4550,2840,2180])

#Knobs to tweek
fdist = 0.3494
Kfl = 2.2e9                        #source: Cheng 1993 (bulk modulus of water)
Rfl = 1000                         #density of water
fdist = 20.0                     #distribution of cracks in degrees
increment = 0.5                   #degree increment (1.0 would be 1 degree increments)

#Give porosity (as a percentage %) and aspect ratios
phiar = np.array([ 0.5, 1.0, 1.0])
aspar = np.array([ 1.0, 1.0, 3])

#need the dot (.) opperator after hudson to specify the class of functions used
#Run model
model = hudson_model.TIcrackmodel()
model.loadvar(aspar,Kfl,Rfl,Kiso,Giso,Ragg,Cagg,phiar,fdist,increment,data)
model.computeROT()
model.run_model()
model.plotfigure()
