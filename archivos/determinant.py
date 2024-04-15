

import random
import math
import sys
import shutil

### this is to multiply matrices
def determinant (B) :

    t1=1.0*B[0][0]*B[1][1]*B[2][2]
    t2=1.0*B[0][1]*B[1][2]*B[2][0]
    t3=1.0*B[0][2]*B[1][0]*B[2][1]
    t4=1.0*B[0][2]*B[1][1]*B[2][0]
    t5=1.0*B[0][1]*B[1][0]*B[2][2]
    t6=1.0*B[0][0]*B[1][2]*B[2][1]

    jaco=1.0*(t1+t2+t3-t4-t5-t6)
    
    return (jaco);
