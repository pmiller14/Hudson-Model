import hudson_model
import numpy as np

print "Run Single rock example"
print "#######################################################################"
#set global variables for simple test
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
data = np.array([5400,4550,2840,2180])

Kfl = 2.2e9                        #source: Cheng 1993 (water)
Rfl = 1000                         #density of water
fdist = 20.0
increment = 0.5

#Set porosity and aspect ratios
phiar = np.array([0.007])
aspar = np.array([0.005])

#need the dot (.) opperator after hudson to specify the class of functions used
model = hudson_model.TIcrackmodel()
model.loadvar(aspar,Kfl,Rfl,Kiso,Giso,Ragg,Cagg,phiar,fdist,increment,data)
model.computeROT()
model.run_model()
model.plotfigure()
