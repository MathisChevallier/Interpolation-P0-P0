import medcoupling as mc
from MEDCouplingRemapper import MEDCouplingRemapper

##Target mesh
#Array creation (11 nodes)
arr = mc.DataArrayDouble(11) 
#fill from 0.0 to 10.0
arr.iota(0) 
#Mesh creation
trgMesh = mc.MEDCouplingCMesh()
trgMesh.setCoords(arr,arr) #square 10*10 cells
trgMesh = trgMesh.buildUnstructured()

##Source mesh
#Array creation (21 nodes)
arr = mc.DataArrayDouble(21)
#from 0.0 to 10.0
arr.iota(0)
arr *= 0.5
#Mesh creation
srcMesh = mc.MEDCouplingCMesh()
srcMesh.setCoords(arr,arr)
srcMesh = srcMesh.buildUnstructured()
#Triangulation of the first 20 cells
tmp = srcMesh[:20]
tmp.simplexize(0)
srcMesh = mc.MEDCouplingUMesh.MergeUMeshes([tmp, srcMesh[20:]])

##Interpolate
remap = MEDCouplingRemapper()
remap.prepare(srcMesh,trgMesh,"P0P0") #tricky / CPU consuming, match the cells

myMatrix = remap.getCrudeMatrix() #gives the matching srcMesh cell ID for each trgMesh cell
print(myMatrix)

#Check if each matched cells got an area of 1
sumByRow = mc.DataArrayDouble(len(myMatrix))
for i, d in enumerate(myMatrix):
    #i=0
    #d={0: 0.125, 1: 0.125, 2: 0.125, 3: 0.125, 40: 0.25, 41: 0.25}
    sum_values = 0.0
    for key in d:
        sum_values += d[key]
    sumByRow[i] = sum_values
print("Is interpolation well prepared? (Area = 1)", sumByRow.isUniform(1.,1e-12))

#construct a field
srcField = mc.MEDCouplingFieldDouble(mc.ON_CELLS, mc.ONE_TIME)
srcField.setMesh(srcMesh)
srcField.fillFromAnalytic(1,"7-sqrt((x-5.)*(x-5.)+(y-5.)*(y-5.))")
srcField.getArray().setInfoOnComponent(0, "powercell [W]")

srcField.setNature(mc.IntensiveMaximum)

remap.transferField(srcField, 1e300)

trgFieldCV = remap.transferField(srcField,1e300)

integSource = srcField.integral(True)[0]
integTarget =  trgFieldCV.integral(True)[0]
print("IntensiveMaximum -- integrals: %lf == %lf" % (integSource, integTarget))
accSource = srcField.getArray().accumulate()[0]
accTarget = trgFieldCV.getArray().accumulate()[0]
print("IntensiveMaximum -- sums: %lf != %lf" % (accSource, accTarget))

srcField.setNature(mc.ExtensiveConservation)
trgFieldI = remap.transferField(srcField,1e300)

integSource = srcField.integral(True)[0]
integTarget =  trgFieldI.integral(True)[0]
print("ExtensiveConservation -- integrals: %lf != %lf" % (integSource, integTarget))

accSource = srcField.getArray().accumulate()[0]
accTarget = trgFieldI.getArray().accumulate()[0]
print("ExtensiveConservation -- sums: %lf == %lf" % (accSource, accTarget))

