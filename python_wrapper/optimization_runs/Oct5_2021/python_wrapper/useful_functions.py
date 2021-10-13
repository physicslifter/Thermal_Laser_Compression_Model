#import necessary packages
from scipy import io
from scipy import interpolate
import numpy as np

def chi_sq_interp(x0,y0,x1, y1):
    #x0: the independent variable data for the model output
    #y0: the original dependent variable data for the model output
    #x1: independent variable data for results
    #y1: dependent variable data for results

    #Step 1 - interpolate the model output to match the dataset we have
    f=interpolate.interp1d(x0,y0, fill_value="extrapolate")
    interpolated_data=f(x1)

    #Now for the chi^2 assessment
    chi_sq=np.sum( np.divide(np.square(y0-y1) , y1))

    return chi_sq

def get_chi_sq_from_output_mat(filepath):

    #Start by reading in the data
    d=io.loadmat(filepath)

    #write to variables so chi^sq_interp will accept the arguments
    #variables are named as follows:
    #[input name to chi_sq_interp (x0,y0,x1 or y1)]_[face # (1, 2, or 3)]_[run #: 73, 76 or 80]
    #x0 is only defined once, since it will always be tlist
    #the t_data arrays are only defined once per experiment run
    
    x0=d['tlist'][0] #tlist data for all runs @ all faces

    #run s88773
    x1_73=d['t_data_73'] #time data
    y0_1_73=d['T11'][0] #face 1 model output
    y0_2_73=d['T12'][0] # face 2 model output
    y0_3_73=d['T13'][0] # face 3 model output
    y1_1_73=d['temp1_73'] # face 1 experimental data
    y1_2_73=d['temp2_73'] # face 2 experimental data
    y1_3_73=d['temp3_73'] # face 3 experimental data

    #run s88776
    x1_76=d['t_data_76'] #time data
    y0_1_76=d['T21'][0] #face 1 model output
    y0_2_76=d['T22'][0] # face 2 model output
    y0_3_76=d['T23'][0] # face 3 model output
    y1_1_76=d['temp1_76'] # face 1 experimental data
    y1_2_76=d['temp2_76'] # face 2 experimental data
    y1_3_76=d['temp3_76'] # face 3 experimental data

    #run s88780
    x1_80=d['t_data_80'] #time data
    y0_1_80=d['T31'][0] #face 1 model output
    y0_2_80=d['T32'][0] # face 2 model output
    y0_3_80=d['T33'][0] # face 3 model output
    y1_1_80=d['temp1_80'] # face 1 experimental data
    y1_2_80=d['temp2_80'] # face 2 experimental data
    y1_3_80=d['temp3_80'] # face 3 experimental data

    #Get chi^2 values for each run
    cs731=chi_sq_interp(x0,y0_1_73,x1_73,y1_1_73) # s88773 face 1
    cs732=chi_sq_interp(x0,y0_2_73,x1_73,y1_2_73) # s88773 face 2
    cs733=chi_sq_interp(x0,y0_3_73,x1_73,y1_3_73) # s88773 face 3
    cs761=chi_sq_interp(x0,y0_1_76,x1_76,y1_1_76) # s88773 face 1
    cs762=chi_sq_interp(x0,y0_2_76,x1_76,y1_2_76) # s88773 face 2
    cs763=chi_sq_interp(x0,y0_3_76,x1_76,y1_3_76) # s88773 face 3
    cs801=chi_sq_interp(x0,y0_1_80,x1_80,y1_1_80) # s88773 face 1
    cs802=chi_sq_interp(x0,y0_2_80,x1_80,y1_2_80) # s88773 face 2
    cs803=chi_sq_interp(x0,y0_3_80,x1_80,y1_3_80) # s88773 face 3

    #finally, return the values as an array
    return [cs731,cs732,cs733,cs761,cs762,cs763,cs801,cs802,cs803]