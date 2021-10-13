#A function for taking input parameters as a numpy array, then saves the matrix in inputs

import numpy
from scipy import io 

def convert(save_name,input_array):

    #Write the dict from the input array
    ia=input_array
    input_values={    
        'sop_lineouts_data':ia[0],
        'reflectivity_data':ia[1],
        'xw':ia[2],
        'aeta':ia[3],
        'a0':ia[4],
        't0a':ia[5],
        'd':ia[6],
        'c':ia[7],
        't0':ia[8],
        'w':ia[9],
        'km':ia[10],
        'dm':ia[11],
        'cm':ia[12],
        't0m':ia[13],
        'peak_temp':ia[14]
}

    #Now let's save this to a matlab array
    io.savemat('inputs/'+save_name+'.mat',mdict=input_values)