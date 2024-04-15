
import random
import math
import sys
import shutil


def print_paraview_micro_stra_stre_disp (v_tetra, stra, stre, disp, lablet, lab):
    
    num_pts=int(len(v_tetra)/3) # 3D
    dim_tet=int(num_pts/4) # Tetrahedro

    b=str(lab)
    c='.vtu'
    nam=lablet+b+c
    
    # v_tetra should be a vector, defined previously
    f=open(nam,'w')
    f.write('<?xml version="1.0"?> \n')
    f.write('<VTKFile type="UnstructuredGrid"> \n')
    f.write('   <UnstructuredGrid> \n')
    f.write('      <Piece NumberOfPoints="')
    f.write(str(int(num_pts)))
    f.write('"  NumberOfCells="')
    f.write(str(int(dim_tet)))
    f.write('">  \n')
    f.write('          <Points>  \n')
    f.write('             <DataArray type="Float64" NumberOfComponents="3" format="ascii">   \n')

# Data array --> Coordinates

    for i in range(0,num_pts,1):
        aux_i=i*3
        f.write( str(v_tetra[int(aux_i+0)]))
        f.write('   ')
        f.write(str( v_tetra[int(aux_i+1)]))
        f.write('   ')
        f.write(str( v_tetra[int(aux_i+2)]))
        f.write('\n')

    f.write('             </DataArray>   \n')
    f.write('          </Points>  \n')
    f.write('          <Cells>  \n')
    f.write('             <DataArray type="Int32" Name="connectivity" format="ascii">   \n')

# Data array --> Connectivity

    for i in range(0,dim_tet,1):
        aux_p=i*4
        f.write( str(aux_p+0))
        f.write('   ')
        f.write(str(aux_p+1))
        f.write('   ')
        f.write(str(aux_p+2))
        f.write('   ')
        f.write(str(aux_p+3))
        f.write('\n')
    
    f.write('              </DataArray>   \n')
    f.write('              <DataArray type="Int32" Name="offsets" format="ascii">   \n')

     # Data array --> Offset

    for i in range(0,dim_tet,1):
        aux_p=(i+1)*4
        f.write( str(aux_p))
        f.write('\n')

    f.write('              </DataArray>   \n')
    f.write('              <DataArray type="UInt8" Name="types" format="ascii">   \n')
    
    for i in range(0,dim_tet,1):
        f.write( str(10))
        f.write('\n')

    f.write('              </DataArray>   \n')
    f.write('          </Cells>  \n')
    f.write('          <PointData>  \n')



    f.write('             <DataArray type="Float64" NumberOfComponents="3" Name="Strains" format="ascii">   \n')
    dim_stra=int(len(stra)/3)
    for i in range(0,dim_stra,1):
        for j in range(0,4,1 ):
            f.write( str(stra[i*3+0]))
            f.write('   ')
            f.write(str(stra[i*3+1]))
            f.write('   ')
            f.write(str(stra[i*3+2]))
            f.write('\n')

    f.write('              </DataArray>   \n')

    f.write('             <DataArray type="Float64" NumberOfComponents="3" Name="Stress" format="ascii">   \n')
    dim_stra=int(len(stra)/3)
    for i in range(0,dim_stra,1):
        for j in range(0,4,1 ):
            f.write( str(stre[i*3+0]))
            f.write('   ')
            f.write(str(stre[i*3+1]))
            f.write('   ')
            f.write(str(stre[i*3+2]))
            f.write('\n')

    f.write('              </DataArray>   \n')

    f.write('             <DataArray type="Float64" NumberOfComponents="3" Name="Displacements" format="ascii">   \n')
    dim_disp=int(len(disp)/3)
    for i in range(0,dim_disp,1):
        
        f.write( str(disp[i*3+0]))
        f.write('   ')
        f.write(str(disp[i*3+1]))
        f.write('   ')
        f.write(str(disp[i*3+2]))
        f.write('\n')

    f.write('              </DataArray>   \n')
    
    f.write('          </PointData>  \n')
    f.write('       </Piece>  \n')
    f.write('   </UnstructuredGrid> \n')
    f.write('</VTKFile> \n')
    f.close()

    return;
