#%%
from fenics import *
import pandas as pd
import os
import numpy as np

class OptimizationGroup:
    pass

class OptimizationData:
    def __init__(self, fp:str):
        self.fp = fp
        self.name = self.fp.split('/')[-1].replace('.xlsx','')
        folder_path = ''
        for k in self.fp.split('/')[:-1]:
            folder_path += k
            folder_path += '/'
        self.folder_path = folder_path
        try:
            self.data = pd.read_excel(fp, 'Data')
        except:
            raise Exception(f"File found, or no Data worksheet found: {fp}")
        try:
            self.metadata = pd.read_excel(fp, 'Metadata')
        except:
            raise Exception("No Metadata worksheet found")
        #=============================================================
        #Checking to ensure that the data file is in an acceptable format
        #Values/keywords that we're looking for are:
        #time
        #lb
        #ub
        #avg_b
        #step1
        #step2
        #step3
        
        #Acceptable keys:
        time_keys = ['time', 'Time']
        lb_keys = ['lb', 'LB']
        ub_keys = ['ub', 'UB']
        avg_b_keys = ['avg_b', 'avgB', 'Average_B', 'AvgB']
        step1_keys = ['step1', 'Step1', 'step 1', 'Step 1','Step_1', 'step_1']
        step2_keys = ['step2', 'Step2', 'step 2', 'Step 2','Step_2', 'step_2']
        step3_keys = ['step3', 'Step3', 'step 3', 'Step 3','Step_3', 'step_3']
        step_keys = [step1_keys, step2_keys, step3_keys]
        #data_labels = ['time', 'lb', 'ub', 'avg_b', 'step data']
        data_labels = ['time', 'lb', 'ub', 'step data']
        #data_labels = ['time', 'step data']
        #data_keys = [time_keys, lb_keys, ub_keys, avg_b_keys, step_keys]
        data_keys = [time_keys, lb_keys, ub_keys, step_keys]
        #data_keys = [time_keys, step_keys]
        for label, key in zip(data_labels, data_keys):
            has_category = False
            for eligible_name in key:
                if type(eligible_name) == str:
                    for data_key in self.data.keys():
                        if data_key == eligible_name:
                            has_category = True
                elif type(eligible_name) == list:
                    for n in eligible_name: #Get each subname
                        for data_key in self.data.keys():
                            #print(data_key)
                            if data_key == n:
                                has_category = True
            #has_category = True
            if not has_category:
                raise Exception(f'Error: no {label} found in the excel file')
        time_key = None
        
        for key in self.data.keys():
            is_time = False
            keep_value = False
            for step in step_keys:
                for k in step:
                    if k == key:
                        keep_value = True
            #for a in avg_b_keys:
            #    if a == key:
            #        keep_value = True
            for u in ub_keys:
                if u == key:
                    keep_value = True
            for l in lb_keys:
                if l == key:
                    keep_value = True
            
            
            for t in time_keys:
                if t == key:
                    keep_value = True
                    is_time = True
                    time_key = t
            
            if keep_value:
                new_key = key.lower().replace(' ','').replace('_','')
                if new_key != key:
                    self.data[new_key] = self.data[key]
                    del self.data[key]
                    
            if not keep_value:
                del self.data[key]
        #print(f'time key: {time_key}')
        self.data['time'] = self.data['time']*10**-9
        self.data['time2'] = self.data.time
        #self.data = self.data.set_index('time')
        #==============================================================
        
        #=============================================================
        #Checking that metadata file is in the right format
        metadata_keys = [ 'xw', 'aeta', 'a0', 't0a', 'reference' ]
        has_xw, has_aeta, has_a0, has_t0a, has_ref = False, False, False, False, False
        for key in metadata_keys:
            has_key = False
            for k in self.metadata.keys():
                if k == key:
                    has_key = True
            if not has_key:
                error_string = f"Data Error: {key} not found in Metadata wksht in {self.fp}"
                raise Exception(error_string)
            
        
        self.xw = self.metadata.xw.item()
        self.aeta = self.metadata.aeta.item()
        self.a0 = self.metadata.a0.item()
        self.t0a = self.metadata.t0a.item()
        self.reference = self.metadata.reference.item()
                    

    def correct_data(self):
        keys = self.data.keys()
        step_keys = []
        corrected_step_keys = []
        face_keys = []
        self.data['avgb'] = (self.data['ub'] + self.data['lb'])/2
        for key in keys:
            if key.lower()[0:4] == 'step':
                step_keys.append(key)
                corrected_step_keys.append(f'{key}_corrected')
        for k, ck in zip(step_keys, corrected_step_keys):
            self.data[ck] = self.data[k] - self.data['avgb']
            #print(ck)
            
        for key in self.data.keys():
            if key[0:4] == 'step':
                if key[-10:] == '_corrected':
                    face_keys.append(key)

        xw = self.metadata.xw.item()
        aeta = self.metadata.aeta.item()
        a0 = self.metadata.a0.item()
        t0a = self.metadata.t0a.item()
        reference = self.metadata.reference.item()
        print('Correcting')
        for key in face_keys:
            print(key)
            SOP_data = self.data[key]
            newSOP=np.real(11605*t0a)/(np.log(1+((1-reference)*a0/(aeta*np.abs(SOP_data)))))
            self.data[key] = newSOP       

try:
    a = OptimizationData('../CleanData/s88773.xlsx')
    print(a.name)
except:
    print('File Not eligible')
    

a = OptimizationData('../CleanData/s88773.xlsx')
# %%
rootdir = '../CleanData'
choices = []
for file in os.listdir(rootdir):
    fp = f'{rootdir}/{file}'
    try:
        a = OptimizationData(fp)
        print(file, 'Format Works')
        choices.append(a.name)
    except:
        print(file, 'Data Format Wrong!')
# %%
