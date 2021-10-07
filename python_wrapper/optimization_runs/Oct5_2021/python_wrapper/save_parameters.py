#A python function for saving any/all parameters required for the transient2DHeatedBlock in a .mat file so that they can be read by the function

#import scipy.io
from scipy import io

#write our function
def write_parameters(peak_temp):
    variables={
        'peak_temp':peak_temp
    }

    #Now let's save the variable
    io.savemat('inputs/function_parameters.mat',mdict=variables)