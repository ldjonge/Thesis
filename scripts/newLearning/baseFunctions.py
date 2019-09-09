import random
import math
import numpy as np
from numpy.random import choice

def randomRound(num):
    if random.random() >= num-math.floor(num):
        return(math.floor(num))
    else:
        return(math.ceil(num))

def wavg(group, avgName, weightName):
    d = group[avgName]
    w = group[weightName]
    try:
        return (d*w).sum() / w.sum()
    except ZeroDivisionError:
        return d.mean()
