#Python script for writing the parameters for the BC file
from scipy.io import savemat as s

s('inputs/Default_BC_parameters.mat',mdict={'peak_temp':28000})