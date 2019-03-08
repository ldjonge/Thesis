#!/usr/bin/env python
import math

def alleleFreq(A,I,O):  #Floats only
    try:
        A = float(A)
        I = float(I)
        O = float(O)
    except ValueError:
        print("Please enter valid numbers")
        return
    if A+I+O >= 0.999 and A+I+O <= 1.001:
        pass
    elif A+I+O >= 99.9 and A+I+O <= 100.1:
        A = A/100
        I = I/100
        O = O/100
    else:
        print("Frequencies do not add up to 1")
        return
    #Calculate r from O frequency
    r = math.sqrt(O)
    #Calculate q from r and I frequency
    d = 4*r**2 + 4*I
    q = (-2*r + math.sqrt(d))/2
    #Calculate p from q,r, and A frequency
    d = 4*(q+r)**2+4*A
    p = (-2*(q+r)+math.sqrt(d))/2

    return({"p":p,"q":q,"r":r})

def phenoFreq(p,q,r):
    try:
        p = float(p)
        q = float(q)
        r = float(r)
    except ValueError:
        print("Please enter valid numbers")
        return
    if p+q+r >= 0.999 and p+q+r <= 1.001:
        pass
    elif p+q+r >= 99.9 and p+q+r <= 100.1:
        p = p/100
        q = q/100
        r = r/100
    else:
        print("Frequencies do not add up to 1")
        return
    A = p**2 + 2*p*q + 2*p*r
    I = q**2 + 2*q*r
    O = r**2
    return({"A":A,"I":I,"O":O})
