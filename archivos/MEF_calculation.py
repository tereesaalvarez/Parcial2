
from time import time
import random
import math
import sys
import shutil
import numpy as np
from numpy  import *
from numpy import linalg

from compute_eigenvector import*
from determinant import*
from eigenvalues import*
from transpose import*
from print_paraview import*
from transpose import*
from determinant import*
from kloc_tet import*


def calculate_FEM (steps, disp_vert, label_test):


    ##########---------
    # Class
    ##########---------

    class mat_proper :

        Prop_young=[]
        Prop_shear=[]
        Prop_poison=[]
        Prop_dens=[]

        def __init__(self, Prop_young, Prop_shear, Prop_poison, Prop_dens) :
            self.Prop_young=Prop_young
            self.Prop_shear=Prop_shear
            self.Prop_poison=Prop_poison
            self.Prop_dens=Prop_dens

    #####--------------------------------------#####
    #####         READ the new mesh            #####
    #####--------------------------------------#####
    timestart=time()
    # read the nodes and the mesh (connectivity)

    v_nodes=[]
    fcell= open('nodes.txt','r') # x, y, z
    for line in fcell.readlines() :
        vecaux=[float(x) for x in line.split()]
        vecaux2=[vecaux[0],vecaux[1],vecaux[2]]
        v_nodes.append(vecaux2)
    fcell.close()

    v_mesh=[]
    fcell= open('mesh.txt','r') # n1, n2, n3, n4
    for line in fcell.readlines() :
        vecaux=[int(x) for x in line.split()]
        vecaux2=[int(vecaux[0]),int(vecaux[1]),int(vecaux[2]),int(vecaux[3])]
        
        v_mesh.append(vecaux2)

    fcell.close()

    v_mechanical=[] 
    fcell= open('mechanical_mat.txt','r') # n1
    iteraux=0
    for line in fcell.readlines() :
        vecaux=[float(x) for x in line.split()]
        iteraux+=1
        vecaux2=[vecaux[0],vecaux[1],vecaux[2]]
        v_mechanical.append(vecaux2)
    fcell.close()
    
    nodal_connectivity=[]
    for i in range(0,len(v_nodes),1):
        nodal_connectivity.append([])
    for i in range(0,len(v_mesh),1):
        for j in range(0,len(v_mesh[i]),1):
            nodal_connectivity[v_mesh[i][j]].append(i)

    elements=len(v_mesh)
    nodes= len(v_nodes)

    #-------------------------------------------------------------------------#
    ######### END MODULE I ##########

    volume_FEM=[]
    energy_FEM=[]
    stress_FEM=[]
    strains_FEM=[]
    dispnode=[]
    coordsnode=[]
    displacements_FEM=[]
    conect=[]
    for i in range(0,elements,1):
        displacements_FEM.append([])
        conect.append([])
        coordsnode.append([])
        dispnode.append([])
        strains_FEM.append([0,0,0,[]])
        stress_FEM.append([0,0,0])
        energy_FEM.append([0,0])
        volume_FEM.append(0)

    ##### END INITIATE THE VECTORS

    time_loc_ini=time()

    mat_properties=[]

    for i in range(0, elements, 1) :

        E1=v_mechanical[i][0]
        V1=v_mechanical[i][1]
        density=v_mechanical[i][2]
        G12=E1/(2*(1+V1))
        
        ind1=mat_proper([E1, E1, E1],[G12, G12, G12],[V1,V1,V1],[density])
        mat_properties.append(ind1)
   
    fload= open('load_reaction_force_'+label_test+'.txt','w')
    fdisp= open('vertical_disp_'+label_test+'.txt','w')
    ftime= open('timeinfo_'+label_test+'.txt','w')
        
    fload.close()
    fdisp.close()
    ftime.close()

    for iste in range(0,int(steps),1) :

        disp_vert_stepi= -(disp_vert/steps)*(iste+1)
        timeini=time()
        totforce=0

        #-------------------------------------------------------------------------#
                ######### END MODULE II ##########

                #######################################
                #################  Solver #################
                #######################################
        v_vec_sol=[]
        disp_el=[]
        for i in range(0,elements,1):
            disp_el.append([])

        Kglob=[ [ 0 for ig in range(len(v_nodes)*3) ] for jg in range(len(v_nodes)*3) ]
        b_fw= [ 0 for ig in range(len(v_nodes)*3) ] # forces weight

        for l in range(0, elements, 1):
            
            el_loc=l
            Eloc=mat_properties[l].Prop_young[0]
            dens=mat_properties[l].Prop_dens[0]

            # --- Mechanical model --- #

            aux_vec_coord=[]
            for j2 in range(0,4,1):
                naux=v_mesh[el_loc][j2]
                for k2 in range(0,3,1):
                    aux_vec_coord.append(v_nodes[naux][k2])

            (Kloc,M1,B,voltet)=kloc_tet (aux_vec_coord, Eloc)

            vec_loc_ind=[]
            for j2 in range(0,4,1):
                for j3 in range(0,3,1):
                    vec_loc_ind.append(v_mesh[el_loc][j2]*3+j3)

            for iloc in range(0,len(Kloc),1):
                
                indaux=vec_loc_ind[iloc]-3*int(vec_loc_ind[iloc]/3)
                if indaux==1.0:
                    f_peso=(voltet/pow(1000.0,3))*dens*9.8 # en dirección -y
                    b_fw[vec_loc_ind[iloc]]+=-f_peso/4.0
                
                for jloc in range(0,len(Kloc[iloc]),1):
                    Kglob[vec_loc_ind[iloc]][vec_loc_ind[jloc]]+=Kloc[iloc][jloc]

        ####################################################################
        #####            Boundary Conditions   displacement  BC       ######
        ####################################################################

        sup_BC=[]
        inf_BC=[]       
        disp_b=[]
        remove_BC=[]
        dx=[]
        for i in range(0,len(v_nodes),1):
            for j in range(0,3,1):
                disp_b.append(0.0)
                dx.append(0.0)

        # inf BC

        for i in range(0,len(v_nodes),1):
            if ((v_nodes[i][0]>495) and (v_nodes[i][1]<5)) :
                
                nodei=i#int(BC[k][i])
                if inf_BC.count(nodei)==0:
                    inf_BC.append(nodei)
                    for k in range(0,3,1):
                        disp_b[nodei*3+k]=0.0
                        dx[nodei*3+k]=0.0
                        if remove_BC.count(nodei*3+k)==0:
                            remove_BC.append(nodei*3+k)
                    

        for i in range(0,len(v_nodes),1):
            if (v_nodes[i][0]<-495 ) :#and (v_nodes[i][1]<5) :
                
                nodei=i#int(BC[k][i])
                if inf_BC.count(nodei)==0:
                    inf_BC.append(nodei)
                    for k in range(0,3,1):
                        disp_b[nodei*3+k]=0.0
                        dx[nodei*3+k]=0.0
                        if remove_BC.count(nodei*3+k)==0:
                            remove_BC.append(nodei*3+k)

        # sup BC
        
        for i in range(0,len(v_nodes),1):
            if (v_nodes[i][0]<2 and v_nodes[i][0]>-2) and (v_nodes[i][1]>97) :
                
                nodei=i#int(BC[k][i])
                if sup_BC.count(nodei)==0:
                    sup_BC.append(nodei)
                    disp_b[nodei*3+1]=disp_vert_stepi
                    dx[nodei*3+1]=disp_vert_stepi
                if remove_BC.count(nodei*3+1)==0:
                    remove_BC.append(nodei*3+1)

        #####################################################
        #####                Solve                     ######
        #####################################################

        disp_b_mat=np.array(disp_b)
        Kglob_solve0=np.array(Kglob)

        b_bc0=-dot(Kglob_solve0,disp_b_mat) # here imponemos desplazamientos, meter tambien otro vector con fuerzas
        b_fw0=np.array(b_fw)
        b0=b_bc0+b_fw0

        Kglob_solve1=np.delete(Kglob_solve0,remove_BC,1)
        Kglob_solve=np.delete(Kglob_solve1,remove_BC,0)

        b=np.delete(b0,remove_BC,0)

        dx1=linalg.solve(Kglob_solve,b)

        indloc=0
        for iloc in range(0,len(dx),1):
            if remove_BC.count(iloc)==0:
                dx[iloc]=1.0*dx1[indloc]
                indloc+=1

        ###############################################
        #########   Compute Reaction Force   ##########
        ###############################################

        RFx=0
        RFy=0
        RFz=0
        for i in range(0,len(sup_BC),1):    
            nodei=sup_BC[i]
            for j in range(0,len(Kglob[nodei*3+2]),1):
                RFx+=Kglob[nodei*3+0][j]*dx[j]
                RFy+=Kglob[nodei*3+1][j]*dx[j]
                RFz+=Kglob[nodei*3+2][j]*dx[j]

            totforce=abs(RFx*RFx+RFy*RFy+RFz*RFz)
            totforcematrix=totforce
                    
        print('RF(x,y,z)=  ', round(RFx,2), round(RFy,2), round(RFz,2)) # reaction force in x, y and z

        ###############################################
        #########  Compute stress and strain ##########
        ###############################################

        # compute the strain and stress values with the computed E

        strains_el=[]
        stresses_el=[]
        strains_el_eig=[]
        stresses_el_eig=[]
        displacements_el=[]
        coord_nodes_el=[]
        for i in range(0,len(v_mesh),1):
            el_loc=i
            Eloc=mat_properties[l].Prop_young[0]

            # --- Mechanical model --- #
            aux_vec_coord=[]
            disp_loc=[]
            for j2 in range(0,4,1):
                naux=v_mesh[el_loc][j2]
                for k2 in range(0,3,1):
                    aux_vec_coord.append(v_nodes[naux][k2])
                    disp_loc.append(dx[naux*3+k2])
                    
            displacements_el.append(disp_loc)
            coord_nodes_el.append(aux_vec_coord)
            
            (Kloc,straM,B,voltet)=kloc_tet (aux_vec_coord, Eloc)

            stra_mat=B
            stress_mat=straM

            local_disp=np.array(disp_loc)

            matloc=np.array(stra_mat)
            vec_strain_loc=np.dot(matloc,local_disp)
                    
            matloc2=np.array(stress_mat)
            vec_stress_loc=np.dot(matloc2,vec_strain_loc)#local_disp)

            stress=[[vec_stress_loc[0],vec_stress_loc[3],vec_stress_loc[4]],[vec_stress_loc[3],vec_stress_loc[1],vec_stress_loc[5]],[vec_stress_loc[4],vec_stress_loc[5],vec_stress_loc[2]]]
            (eig1, eig2, eig3)=eigenvalues(stress)
            stresses_el_eig.append([eig1,eig2,eig3])
            stresses_el.append(stress)

            strain =[[vec_strain_loc[0],vec_strain_loc[3],vec_strain_loc[4]],[vec_strain_loc[3],vec_strain_loc[1],vec_strain_loc[5]],[vec_strain_loc[4],vec_strain_loc[5],vec_strain_loc[2]]]
            (eig1, eig2, eig3)=eigenvalues(strain)
            strains_el_eig.append([eig1,eig2,eig3])
            strains_el.append(strain)

########### MODULE IV ############
#-------------------------------------------------------------------------#

        for j in range(0, elements, 1) :
            elem_stress=stresses_el_eig[j][0]
            elem_strain=strains_el_eig[j][0]
            strains_FEM[j][0]=strains_el[j][0][0]
            strains_FEM[j][1]=strains_el[j][1][1]
            strains_FEM[j][2]=strains_el[j][2][2]
            
            stress_FEM[j][0]=stresses_el[j][0][0]
            stress_FEM[j][1]=stresses_el[j][1][1]
            stress_FEM[j][2]=stresses_el[j][2][2]

            # displacements of the FE mesh
            displacements_FEM[j]=[]
            for g in range(0,4,1):
                nodeg=v_mesh[j][g]
                for k in range(0,3,1):
                    displacements_FEM[j].append(dx[nodeg*3+k]) 
            vc=[]
            for k in range(0,4,1):
                vc.append(v_nodes[v_mesh[j][k]][0])
                vc.append(v_nodes[v_mesh[j][k]][1])
                vc.append(v_nodes[v_mesh[j][k]][2])
                
            vol_el=compute_volume(vc[0],vc[1],vc[2],vc[3],vc[4],vc[5],vc[6],vc[7],vc[8],vc[9],vc[10],vc[11])  

            s11=strains_el[j][0][0] 
            s22=strains_el[j][1][1] 
            s33=strains_el[j][2][2] 
            s12=strains_el[j][0][1] 
            s13=strains_el[j][0][2] 
            s23=strains_el[j][1][2] 
            stra=[[s11,s12,s13],[s12,s22,s23],[s13,s23,s33]]
            
            strain_calc=[s11,s22,s33]
            (eig1, eig2, eig3)=eigenvalues(stra)

            strains_FEM[j][3]=[s11,s22,s12]
            
            (v1,v2,v3)=compute_eigenvector (strains_FEM[j][0], stra)
           
            eigv=[v1,v2,v3] # extract this: extract the components of the stress

            Eloc=mat_properties[j].Prop_young[0]
            
            volume_FEM[j]=vol_el
            energy_FEM[j][0]=vol_el*(strains_FEM[j][0]*stress_FEM[j][0]+strains_FEM[j][1]*stress_FEM[j][1]+strains_FEM[j][2]*stress_FEM[j][2])*0.5

                
        fload= open('load_reaction_force_'+label_test+'.txt','a')
        fdisp= open('vertical_disp_'+label_test+'.txt','a')
        ftime= open('timeinfo_'+label_test+'.txt','a')
        
        fload.write(str(totforce))
        fload.write('\n')

        fdisp.write(str(disp_vert_stepi))
        fdisp.write('\n')

        ftime.write(str(time()-timeini))
        ftime.write('\n')

        fload.close()
        fdisp.close()
        ftime.close()

        print('Step time  ', time()-timeini)
        
        stepi=int(iste)+1

        v_print_elements_0=[]
        v_print_strains_0=[]
        v_print_stress_0=[]
        v_print_displacement_0=[]

        for j in range(0, elements, 1) :
            for p in range(0, 4, 1) :
                for q in range(0,3,1):
                    coord_loc=v_nodes[v_mesh[j][p]][q]
                    v_print_elements_0.append(coord_loc)
                    v_print_displacement_0.append(displacements_FEM[j][p*3+q]) 
                    
            for p in range(0,3,1) :
                v_print_strains_0.append(strains_FEM[j][3][p])
                v_print_stress_0.append(stress_FEM[j][p])

        
        print_paraview_micro_stra_stre_disp  (v_print_elements_0, v_print_strains_0, v_print_stress_0, v_print_displacement_0, label_test, stepi)

    print('end of step')
    print('Time',time()-timestart)

    load_disp=[]
    fload= open('load_reaction_force_'+label_test+'.txt','r')
    for line in fload.readlines() :
        vecaux=[float(x) for x in line.split()]
        load_disp.append([vecaux[0]])
    fload.close()

    fdisp= open('vertical_disp_'+label_test+'.txt','r')
    ind=0
    for line in fdisp.readlines() :
        vecaux=[float(x) for x in line.split()]
        load_disp[ind].append(vecaux[0])
        ind+=1
    fdisp.close()

    return (disp_vert_stepi, load_disp);





