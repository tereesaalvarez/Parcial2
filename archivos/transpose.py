
import random
import math
import sys
import shutil


### this is to multiply matrices
def transpose (B) :
    
    row_number=len(B)
    colum_number=len(B[0])
    
    Bt= [ [ 0.0 for ig in range(row_number) ] for jg in range(colum_number) ]
    
    for i in range(0,row_number,1):
        for j in range(0,colum_number,1):
            Bt[j][i]=1.0*B[i][j]
    
    return (Bt);
