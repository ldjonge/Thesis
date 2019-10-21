import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pypet.utils.explore import cartesian_product
from pypet.trajectory import Trajectory; traj = Trajectory

def sigmoSym(x, Lambda):
    f= 2*np.abs(x)-1
    f=(np.sign(f)*(np.abs(f)**Lambda)+1)/2
    f=(1-np.sin(np.pi*f/2)**2)
    return f

def sigmoCurve(x, stiffness, inflexion):
    aInfl = 1/inflexion
    t=aInfl * x / (1+aInfl*x)
    f=sigmoSym(t, stiffness)
    return f

def singleSim(pars, start):
    for i in range(100):
        new = [x*sigmoCurve(x, pars[0], pars[1]) for x in start]
        newSum = sum(new)
        start = [x/newSum for x in new]
    return start

start = [1/3, 1/3, 1/3]
parSpace = cartesian_product({"stiffness":list(np.arange(0.01,1,0.01)), "inflexion":list(np.arange(0.01,1,0.01))})
parSpace = list(zip(parSpace["stiffness"], parSpace["inflexion"]))
print(parSpace[:10])

output = []
for pars in parSpace:
    out = singleSim(pars, start)
    morphs = 0
    for m in out:
        if m > 0.01:
            morphs += 1
    output.append(morphs)

counts = {1:0, 2:0, 3:0}
for i in output:
    counts[i] += 1

print(counts)
