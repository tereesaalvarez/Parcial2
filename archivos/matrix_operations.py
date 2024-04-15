import random
import math
import sys
import shutil


def multiply_mat (A,B) :
    
    row_number=len(A)
    colum_number=len(B[0])
    
    C= [ [ 0 for ig in range(colum_number) ] for jg in range(row_number) ]
    
    for i in range(0,row_number,1):
        for j in range(0,colum_number,1):
            termij=0
            for k in range(0,len(A[i]),1):
                termij+=A[i][k]*B[k][j]

            C[i][j]=termij
    
    return (C);

### this is to calculate matrices(3x3) determinant
def det3(Matrix_dim3):
    M = []
    M = Matrix_dim3

    d = M[0][0]*(M[1][1]*M[2][2]-M[1][2]*M[2][1]) - M[1][0]*(M[0][1]*M[2][2]-M[0][2]*M[2][1]) + M[2][0]*(M[0][1]*M[1][2]-M[0][2]*M[1][1])

    return(d);

### this is to calculate matrices(4x4) determinant
def det4(M):
    d = 0.0

    M0=[[ M[1][1],M[1][2],M[1][3] ],[ M[2][1],M[2][2],M[2][3] ],[ M[3][1],M[3][2],M[3][3] ]]
    M1=[[ M[0][1],M[0][2],M[0][3] ],[ M[2][1],M[2][2],M[2][3] ],[ M[3][1],M[3][2],M[3][3] ]]
    M2=[[ M[0][1],M[0][2],M[0][3] ],[ M[1][1],M[1][2],M[1][3] ],[ M[3][1],M[3][2],M[3][3] ]]
    M3=[[ M[0][1],M[0][2],M[0][3] ],[ M[1][1],M[1][2],M[1][3] ],[ M[2][1],M[2][2],M[2][3] ]]

    d=det3(M0)+det3(M1)+det3(M2)+det3(M3)

    return(d);


### this is to remove matrices column
def remove_column(Matrix, Column):
        M = Matrix
        l = len(M)

        for i in range (0,l,1):
                M[i].pop(Column)

        return(M);

### this is to copy matrices
def copy (A):
    
    B = []
    for i in range(0,len(A),1):
        B.append([])
        for j in range(0,len(A[0]),1):
            B[i].append(A[i][j])
    
    return (B);

### this is to subtract matrices
def substract (A,B):
    C = [ [ 0 for ig in range(2) ] for jg in range(2) ]
    
    for i in range(0,len(A),1):
        for j in range(0,len(A[0]),1):
            C[i][j] = A[i][j]-B[i][j]
    
    return (C);

    
