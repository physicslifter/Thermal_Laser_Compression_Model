# A simple function for converting the boundary function file to .txt for editing
import os 

#Write command for conversion
filename='BC_external_exp_SampleCopy'
os.system('cmd /c "ren '+filename+'.m '+filename+'.txt"')