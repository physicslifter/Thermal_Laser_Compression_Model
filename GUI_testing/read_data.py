# import packages to read from the excel file
import pandas as pd
import time

    #read in file
def write_to_file():
    myfile='../python_wrapper/optimization_runs/Aug17_2021/Run2.xlsx'
    df=pd.read_excel(myfile)

    #write file to optimization_run.txt in a realistic way
    for run_num in range(len(df)):

        #write values to variables
        sls=df['SLS'][run_num]
        peak_temp=df['peak_temp'][run_num]
        a=df['a'][run_num]
        b=df['b'][run_num]
        time_shift=df['time_shift'][run_num]
        diffusivity=df['diffusivity'][run_num]

        #write the values to the txt file
        with open('optimization_data.csv', 'a+') as file_object:
            if run_num>0:
                file_object.write('\n')#newline unless it is the first
            file_object.write(str(run_num)+','+str(sls)+', '+str(peak_temp)+', '+str(a)+', '+str(b)+', '+str(time_shift)+', '+str(diffusivity))

        #wait around for some time before reading in the next line
        time.sleep(1)

if __name__ == '__main__':
    #read_data.py executed as script
    write_to_file()
