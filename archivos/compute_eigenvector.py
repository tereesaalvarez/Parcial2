
import random
import math
import sys
import shutil


### THIS IS TO CALCULATE THE EIGENVECTORS ###

def compute_eigenvector (eigref, stra) :
    
    eigref+=0.00000000001
    s11=stra[0][0]-eigref
    s12=stra[0][1]+0.00000001
    s13=stra[0][2]+0.00000001
    s22=stra[1][1]-eigref
    s23=stra[1][2]+0.00000001
    s33=stra[2][2]-eigref

    t1=(s13*s12/(-s11)+s23)/(-(s12*s12/(-s11)+s22))
    A=(s12*t1/(-s11))+s13/(-s11)
    B=t1
    C=1
    
    uni=pow(A*A+B*B+C*C,0.5)
    A*=1/uni
    B*=1/uni
    C*=1/uni

    return (A, B, C);
