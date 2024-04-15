
import random
import math
import sys
import shutil

from determinant import*

### this is only for tetrahedros
def compute_volume(xev,yev,zev,xn1,yn1,zn1,xn2,yn2,zn2,xn3,yn3,zn3) :

    x1=1.0*(xn1-xev)
    y1=1.0*(yn1-yev)
    z1=1.0*(zn1-zev)

    x2=1.0*(xn2-xev)
    y2=1.0*(yn2-yev)
    z2=1.0*(zn2-zev)

    x3=1.0*(xn3-xev)
    y3=1.0*(yn3-yev)
    z3=1.0*(zn3-zev)
    
    mat=[[x1,x2,x3],[y1,y2,y3],[z1,z2,z3]]

    volume= (1.0/6.0)*abs(determinant (mat))

    return (volume);






















