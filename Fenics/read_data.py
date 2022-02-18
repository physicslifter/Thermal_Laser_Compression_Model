#%%
#A script for writing our data to a JSON file
import pandas as pd
import json
import numpy as np

#Filenames
fp73='../Data/s88773_sop_lineouts_TP.xlsx'
fp76='../Data/s88776_SOP_TP.xlsx'
fp80='../Data/s88780_SOP_TP.xlsx'
fp83='../Data/s86483_SOP_TP.xlsx'

#Reading in data as pandas dataframes
s88773=pd.read_excel(fp73, sheet_name='June2021')
s88776=pd.read_excel(fp76)
s88780=pd.read_excel(fp80)
s86483=pd.read_excel(fp83)

#====================================================
#Defining start & stop times

#s88773
start1_73=179
start2_73=259 
start3_73=310 
stop1_73=265
stop2_73=415
stop3_73=706

#s88776
start1_76=88
start2_76=93
start3_76=246
stop1_76= 172
stop2_76= 418
stop3_76= 468

#s88780
start1_80=342
start2_80=582
start3_80=377
stop1_80= 412
stop2_80= 732
stop3_80= 800

#s86483
start1_83=342
start2_83=582
start3_83=377
stop1_83= 412
stop2_83= 732
stop3_83= 800

dict_73={
    'face1':[list((s88773['Time'][start1_73:stop1_73].to_numpy()*10**-9).astype(float)), list((s88773['step1_corrected'][start1_73:stop1_73].to_numpy()).astype(float))],
    'face2':[list((s88773['Time'][start2_73:stop2_73].to_numpy()*10**-9).astype(float)), list((s88773['step2_corrected'][start2_73:stop2_73].to_numpy()).astype(float))],
    'face3':[list((s88773['Time'][start3_73:stop3_73].to_numpy()*10**-9).astype(float)), list((s88773['step3_corrected'][start3_73:stop3_73].to_numpy()).astype(float))]
}

dict_76={
    'face1':[list((s88776['time'][start1_76:stop1_76].to_numpy()*10**-9).astype(float)), list((s88776['step1_corrected'][start1_76:stop1_76].to_numpy()).astype(float))],
    'face2':[list((s88776['time'][start2_76:stop2_76].to_numpy()*10**-9).astype(float)), list((s88776['step2_corrected'][start2_76:stop2_76].to_numpy()).astype(float))],
    'face3':[list((s88776['time'][start3_76:stop3_76].to_numpy()*10**-9).astype(float)), list((s88776['step3_corrected'][start3_76:stop3_76].to_numpy()).astype(float))]
}

dict_80={
    'face1':[list((s88780['time'][start1_80:stop1_80].to_numpy()*10**-9).astype(float)), list((s88780['step1_corrected'][start1_80:stop1_80].to_numpy()).astype(float))],
    'face2':[list((s88780['time'][start2_80:stop2_80].to_numpy()*10**-9).astype(float)), list((s88780['step2_corrected'][start2_80:stop2_80].to_numpy()).astype(float))],
    'face3':[list((s88780['time'][start3_80:stop3_80].to_numpy()*10**-9).astype(float)), list((s88780['step3_corrected'][start3_80:stop3_80].to_numpy()).astype(float))]
}

dict_83={
    'face1':[list((s86483['time'][start1_83:stop1_83].to_numpy()*10**-9).astype(float)), list((s86483['step1_corrected'][start1_83:stop1_83].to_numpy()).astype(float))],
    'face2':[list((s86483['time'][start2_83:stop2_83].to_numpy()*10**-9).astype(float)), list((s86483['step2_corrected'][start2_83:stop2_83].to_numpy()).astype(float))],
    'face3':[list((s86483['time'][start3_83:stop3_83].to_numpy()*10**-9).astype(float)), list((s86483['step3_corrected'][start3_83:stop3_83].to_numpy()).astype(float))]
}


data={}
data['s88773']=dict_73
data['s88776']=dict_76
data['s88780']=dict_80
data['s86483']=dict_83

#save our new data dictionary to data_dit.json
data_json=json.dumps(data)
f=open('data_dict.json', 'w')
f.write(data_json)
f.close()