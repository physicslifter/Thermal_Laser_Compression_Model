from scipy import interpolate
from run_model import main as run
import numpy as np
import json

def least_squares_interp(x0,y0,x1,y1):
    #x0: the independent variable data for the model output
    #y0: the dependent variable data for the model output
    #x1: independent variable data for results
    #y1: dependent variable data for results
    f=interpolate.interp1d(x0,y0, bounds_error=False, fill_value='extrapolate')
    interpolated_data=f(x1)

    least_squares=np.sum(np.square(y1-interpolated_data))

    return least_squares

def load_json(fname):
    with open(fname) as json_file:
        data=json.load(json_file)
        
    return data

def get_chi_sq(model_dict, data_dict):
    shots=list(model_dict.keys())
    counter0=0
    val_storage_array=np.zeros((len(shots),3))
    for shot in shots:
        model_shot_data=model_dict[shot]
        faces=list(model_shot_data.keys())
        face_storage_array=np.zeros(len(faces))
        counter1=0
        for face in faces:
            model_face_data=model_shot_data[face]
            experimental_face_data=data_dict[shot][face]
            chi_2=least_squares_interp(model_face_data[0],model_face_data[1], experimental_face_data[0], experimental_face_data[1])
            #face_storage_array[counter1]=chi_2
            val_storage_array[counter0, counter1]=chi_2
            counter1=counter1+1
        #val_storage_array[counter0]=face_storage_array
        counter0=counter0+1
        
    return val_storage_array  
