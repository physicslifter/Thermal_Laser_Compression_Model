# A simple function for converting the boundary function file back to .m format after editing in python
import os 

#Write command for conversion
filename='BC_external_exp_SampleCopy'
os.system('cmd /c "ren '+filename+'.txt '+filename+'.m"')