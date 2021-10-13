#function(s) for manipulating boundary condition file
import os 

def update(filename,peak_temp,diffusivity,time_shift):
    #Changing the name of the file to .txt so we can edit
    os.system('cmd /c "ren '+filename+'.m '+filename+'.txt"')

    f=open(filename+'.txt') #open file
    s=f.readlines() #read file contents into a list
    f.close() #close the file

    #change lines accordingly
    s[27]='  '+'peak='+str(peak_temp)+'; \n'#peak_temp
    s[28]='  '+'a='+str(diffusivity)+'; \n'#diffusivity
    s[29]='  '+'b='+str(time_shift)+'; \n'#time_shift

    #Now write the new lines to the .txt BC file
    with open(filename+'.txt','w') as filehandle:
        filehandle.writelines("%s" % line for line in s)

    #And finally, change the file back to .m
    os.system('cmd /c "ren '+filename+'.txt '+filename+'.m"')