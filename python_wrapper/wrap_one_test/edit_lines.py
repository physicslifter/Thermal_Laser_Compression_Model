#Read in lines from the standard file
from scipy import io
import os

#Defining quantities to be changed
#==============================================================
new_input_file='standard_input.mat'
bcfile='BC_standard'
output_file='standard_output.mat'
#==============================================================

#Opening the file contents
#==============================================================
fname='standard_eq' #standard file name
os.system('cmd /c "ren '+fname+'.m '+fname+'.txt"') #change the file from .m to .txt so that it is editable
f=open(fname+'.txt') #Now open the file as a .txt file
s=f.readlines()
f.close() #close the file
#==============================================================

#Rewriting file contents
#==============================================================
#Rewriting the bcfile lines in the matlab file
s[6]='input_parameter_file="'+new_input_file+'";\n' #rewrite line for reading in the input parameter file
s[282]="    thermalBC(thermalmodelT,'Edge',51,'Temperature',@"+bcfile+"); %8, 12\n"
s[283]="    thermalBC(thermalmodelT,'Edge',52,'Temperature',@"+bcfile+"); %9, 13\n"
s[284]="    thermalBC(thermalmodelT,'Edge',53,'Temperature',@"+bcfile+"); %10, 14\n"
s[286]="    thermalBC(thermalmodelT,'Edge',54,'Temperature',@"+bcfile+"); %8, 12\n"
s[287]="    thermalBC(thermalmodelT,'Edge',55,'Temperature',@"+bcfile+"); %9, 13\n"
s[288]="    thermalBC(thermalmodelT,'Edge',56,'Temperature',@"+bcfile+"); %10, 14\n"
s[290]="    thermalBC(thermalmodelT,'Edge',57,'Temperature',@"+bcfile+"); %8, 12\n"
s[291]="    thermalBC(thermalmodelT,'Edge',58,'Temperature',@"+bcfile+"); %9, 13\n"
s[292]="    thermalBC(thermalmodelT,'Edge',59,'Temperature',@"+bcfile+"); %10, 14\n"
s[294]="    thermalBC(thermalmodelT,'Edge',60,'Temperature',@"+bcfile+"); %8, 12\n"
s[295]="    thermalBC(thermalmodelT,'Edge',61,'Temperature',@"+bcfile+"); %9, 13\n"
s[296]="    thermalBC(thermalmodelT,'Edge',62,'Temperature',@"+bcfile+"); %10, 14\n"

s[565]="save('"+output_file+"','tlist','T11','T12','T13','T21','T22','T23','T31','T32','T33','T41','T42','T43','t_data1_73','t_data2_73','t_data3_73','temp1_73','temp2_73','temp3_73','t_data1_76','t_data2_76','t_data3_76','temp1_76','temp2_76','temp3_76','t_data1_80','t_data2_80','t_data3_80','temp1_80','temp2_80','temp3_80','t_data1_83','t_data2_83','t_data3_83','temp1_83','temp2_83','temp3_83')"
#==============================================================

#Rewriting lines & saving...
#==============================================================
with open(fname+'.txt','w') as filehandle:
    filehandle.writelines("%s" % line for line in s)

os.system('cmd /c "ren '+fname+'.txt '+fname+'.m"')
#==============================================================
