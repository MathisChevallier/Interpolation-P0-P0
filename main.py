import medcoupling as mc
import math
#from MEDCouplingRemapper import MEDCouplingRemapper
print("Import are done")

#creation array
d = mc.DataArrayDouble(6,2)
d[:,0] = 3.
d[:,1] = list(range(6))
d[:,1] *= math.pi/3.
d = d.fromPolarToCart() #perd instance initiale
d.setInfoOnComponents(["X [m]","Y [m]"])

print(d.getValues()) #affiche valeur sous forme de liste python
print("Uniform array?", d.magnitude().isUniform(3.,1e-12))

#manipulation
d2 = mc.DataArrayDouble.Aggregate(d)
oldNbOfTuples = d2.getNumberOfTuples()
c,cI = d2.findCommonTuples(1e-12) #retourne tab de tuple et tab d'index des tuples



print("The end")