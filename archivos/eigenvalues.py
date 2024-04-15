
import random
import math
import sys
import shutil

import numpy as np
from numpy  import *
from numpy import linalg

### CALCULATE THE EIGENVALUES ###

def eigenvalues (A):
    
    I=[[1,0,0],[0,1,0],[0,0,1]]
    B=[]
    pi=3.141516
    for i in range(0,3,1):    
        B.append([0,0,0])

    p = A[0][1]*A[0][1] + A[1][2]*A[1][2] + A[0][2]*A[0][2]

    if (p == 0) :
        eig1 = A[0][0]
        eig2 = A[1][1]
        eig3 = A[2][2]
    else :

        q = (A[0][0]+A[1][1]+A[2][2])/3
        p = (A[0][0] - q)*(A[0][0] - q) + (A[1][1] - q)*(A[1][1] - q) + (A[2][2] - q)*(A[2][2] - q) + 2 * p
        
        p = pow(p / 6, 0.5)
        for i in range(0,3,1):
            for j in range(0,3,1) :
                B[i][j]=(1 / p) * (A[i][j] - q * I[i][j])
               
        t1 = B[0][0]*B[1][1]*B[2][2]
        t2=B[0][1]*B[1][2]*B[0][2]
        t3=B[0][1]*B[1][2]*B[0][2]
        t4=B[0][2]*B[1][1]*B[0][2]
        t5=B[1][2]*B[1][2]*B[0][0]
        t6=B[0][1]*B[0][1]*B[2][2]

        r=(t1+t2+t3-t4-t5-t6)/2

        if (r <= -1) :
            phi = pi / 3

        else :
            if (r>=1) :
                phi = 0
            else :
                phi=math.acos(r)/3
   
        eig1 = q + 2 * p * math.cos(phi)
        eig3 = q + 2 * p * math.cos(phi + pi * (2/3))
        eig2 = 3 * q - eig1 - eig3

    return (eig1, eig2, eig3);
