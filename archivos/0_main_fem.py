
from MEF_calculation import calculate_FEM


############## PARAMETERS ################
disp_vert=0.1
steps=10
##########################################


(disp_vert_stepi, load_disp)=calculate_FEM (steps, disp_vert, 'MEF_beam_bulk_')


print('end of all the calculation')







