import random
import math
import sys
import shutil

import numpy as np
from numpy  import *
from numpy import linalg
from matrix_operations import *
from compute_volume import*

def kloc_tet (aux_vec_coord,Eloc):
    #vc es el vector con las coordenadas aux_vec_coord[0]
    voltet=compute_volume(aux_vec_coord[0],aux_vec_coord[1],aux_vec_coord[2],aux_vec_coord[3],aux_vec_coord[4],aux_vec_coord[5],aux_vec_coord[6],aux_vec_coord[7],aux_vec_coord[8],aux_vec_coord[9],aux_vec_coord[10],aux_vec_coord[11])
    vc=[[aux_vec_coord[0],aux_vec_coord[1],aux_vec_coord[2]],[aux_vec_coord[3],aux_vec_coord[4],aux_vec_coord[5]],[aux_vec_coord[6],aux_vec_coord[7],aux_vec_coord[8]],[aux_vec_coord[9],aux_vec_coord[10],aux_vec_coord[11]]]

    vol0= np.array([[vc[0][0],vc[0][1],vc[0][2],1],[vc[1][0],vc[1][1],vc[1][2],1],[vc[2][0],vc[2][1],vc[2][2],1],[vc[3][0],vc[3][1],vc[3][2],1]])
    vol=np.linalg.det(vol0)/6.0
    
    # b1
    b10= np.array([[vc[1][1],vc[1][2],1], [vc[2][1],vc[2][2],1],[vc[3][1],vc[3][2],1]])
    b1=np.linalg.det(b10)*1/(6*vol)
    
    # b2
    b20= np.array([[vc[0][1],vc[0][2],1], [vc[2][1],vc[2][2],1],[vc[3][1],vc[3][2],1]])
    b2=-1.0*np.linalg.det(b20)*1/(6*vol)

    # b3
    b30= np.array([[vc[0][1],vc[0][2],1], [vc[1][1],vc[1][2],1],[vc[3][1],vc[3][2],1]])
    b3=1.0*np.linalg.det(b30)*1/(6*vol)

    # b4
    b40= np.array([[vc[0][1],vc[0][2],1], [vc[1][1],vc[1][2],1],[vc[2][1],vc[2][2],1]])
    b4=-1.0*np.linalg.det(b40)*1/(6*vol)

    # c1
    c10= np.array([[vc[1][0],vc[1][2],1], [vc[2][0],vc[2][2],1],[vc[3][0],vc[3][2],1]])
    c1=-1.0*np.linalg.det(c10)*1/(6*vol)

    # c2
    c20= np.array([[vc[0][0],vc[0][2],1], [vc[2][0],vc[2][2],1],[vc[3][0],vc[3][2],1]])
    c2=1.0*np.linalg.det(c20)*1/(6*vol)

    # c3
    c30= np.array([[vc[0][0],vc[0][2],1], [vc[1][0],vc[1][2],1],[vc[3][0],vc[3][2],1]])
    c3=-1.0*np.linalg.det(c30)*1/(6*vol)

    # c4
    c40= np.array([[vc[0][0],vc[0][2],1], [vc[1][0],vc[1][2],1],[vc[2][0],vc[2][2],1]])
    c4=1.0*np.linalg.det(c40)*1/(6*vol)

    # h1
    h10= np.array([[vc[1][0],vc[1][1],1], [vc[2][0],vc[2][1],1],[vc[3][0],vc[3][1],1]])
    h1=1.0*np.linalg.det(h10)*1/(6*vol)

    # h2
    h20= np.array([[vc[0][0],vc[0][1],1], [vc[2][0],vc[2][1],1],[vc[3][0],vc[3][1],1]])
    h2=-1.0*np.linalg.det(h20)*1/(6*vol)

    # h3
    h30= np.array([[vc[0][0],vc[0][1],1], [vc[1][0],vc[1][1],1],[vc[3][0],vc[3][1],1]])
    h3=1.0*np.linalg.det(h30)*1/(6*vol)

    # h4
    h40= np.array([[vc[0][0],vc[0][1],1], [vc[1][0],vc[1][1],1],[vc[2][0],vc[2][1],1]])
    h4=-1.0*np.linalg.det(h40)*1/(6*vol)
    
    Bmat0=[]
    Bmat0.append([b1, 0, 0, b2, 0, 0, b3, 0, 0, b4, 0, 0 ])
    Bmat0.append([0, c1, 0, 0, c2, 0, 0, c3, 0, 0, c4, 0 ])
    Bmat0.append([0, 0, h1, 0, 0, h2, 0, 0, h3, 0, 0, h4 ])
    Bmat0.append([c1, b1, 0, c2, b2, 0, c3, b3, 0, c4, b4, 0 ])
    Bmat0.append([0, h1, c1, 0, h2, c2, 0, h3, c3, 0, h4, c4 ])
    Bmat0.append([h1, 0, b1, h2, 0, b2, h3, 0, b3, h4, 0, b4 ])

    Bmat=np.array(Bmat0)
    
    Bt=transpose (Bmat)
    Epart=Eloc
    pois=0.2

    lame1=pois*Epart/((1+pois)*(1-2*pois))
    lame2=Epart/(2*(1+pois))

    d1=1.0*lame1+2.0*lame2
    d2=1.0*lame1
    d3=1.0*lame2

    # Compute D
    D=[[d1,d2,d2,0,0,0],[d2,d1,d2,0,0,0],[d2,d2,d1,0,0,0],[0,0,0,d3,0,0],[0,0,0,0,d3,0],[0,0,0,0,0,d3]]
     
    # Compute K
    M1=multiply_mat(Bt,D)
    Kloc0=multiply_mat (M1,Bmat)
    Kloc=np.multiply(Kloc0,abs(vol))
    D1=np.array(D)
    straMat=D1 #np.dot(D,Bmat) #multiply_mat(D,Bmat)

    return (Kloc,straMat,Bmat,voltet);






    
