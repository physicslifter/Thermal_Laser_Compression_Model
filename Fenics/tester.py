#%%
import main

run_dict = {'s88773': [1,2,3], #Define which runs & faces to assess
            's88776': [1,2,3], 
            's88780': [1,2,3],
            's86483': [1,2,3]
            }

params=[0.03, # a
        30, # b
        1.8*10**-8, # start_time
        30000 # peak_temp
        ]
#%%
my_group=main.OptimizationGroup(120, run_dict) #define optimization group

#Run 20 times
count=0
while count<20:
    my_group.run(params) # Run Finite Element Model
    count=count+1
print('Done') # Print done when finished
# %%
import opt_func

#Rewrite the global parameters to be the same
vals=main.load_json('default_values_dict.json') # load in default values from the dictionary of default values
tf, num_steps, MgO_length, init_temp, k_1, rho, c, = vals['tf'], vals['num_steps'], vals['MgO_length'], vals['init_temp'], vals['k_1'], vals['rho'], vals['c']

shots=['s88773', 's88776', 's88780', 's86483']
params=[0.03, # a
        30, # b
        1.8*10**-8, # start_time
        30000 # peak_temp
        ]
count=0

#Run 20 times
while count<20:
    bounds=((0,0.05),(0,100),(1.75*10**-8, 2.25*10**-8), (30000,100000)) #define bounds (required for changing optimization params)
    opt_func.rewrite_global_params(num_steps, shots, 'opt_func_time_test', bounds, 10) #rewrite the global parameters
    opt_func.optimizable_function(params)
    count=count+1


# %%
import main
run_dict = {'s88773': [1,2,3], #Define which runs & faces to assess
            's88776': [1,2,3], 
            's88780': [1,2,3],
            's86483': [1,2,3]
            }
       #(min, max)
bounds=((0.01, 0.05), # a
        (0, 100), # b
        (1.75*10**-8, 2.25*10**-8), # start_time
        (30000, 80000)) # peak_temp

opt_data=main.OptimizationData('inter_data_stor.optimization_data.txt')
optimization_group=main.OptimizationGroup(120, run_dict)

main.optimize(optimization_group, opt_data, bounds=bounds, popsize=3, maxiter=3)
# %%
import main
main.optimize(optimization_group, opt_data, bounds=bounds, popsize=10, maxiter=10)


# %%
from pprint import pprint
import datetime
run_dict = {'s88773': [1,2,3], #Define which runs & faces to assess
            's88776': [1,2,3], 
            's88780': [1,2,3],
            's86483': [1,2,3]
            }

       #(min, max)
bounds=((0.01, 0.05), # a
        (0, 100), # b
        (1.75*10**-8, 2.25*10**-8), # start_time
        (30000, 80000)) # peak_temp

start_stop_list=[               #s88773
                 [[179, 240],       # face 1
                  [259, 390],       # face 2
                  [310, 700]],      # face 3 
                                #s88776
                 [[88, 155],        # face 1
                  [93, 400],        # face 2
                  [246, 468]],      # face 3
                                #s88780
                 [[145, 190],       # face 1
                  [240, 350],       # face 2
                  [310, 502]],      # face 3
                                #s86483
                 [[355, 400],       # face 1
                  [582, 732],       # face 2
                  [377, 800]]]      # face 3

name='Test 1'

with open('../pprint_test', 'wt') as out:
        out.write(name+': Metadata')
        out.write('\n')
        out.write('Optimization initialized @: ')
        out.write('\n')
        out.write(str(datetime.datetime.now()))
        out.write('\n')
        out.write(' ')
        out.write('\n')
        out.write('===========================')
        out.write('\n')
        out.write('Runs & Faces: ')
        out.write('\n')
        pprint(run_dict, stream=out)
        out.write(' ')
        out.write('\n')
        out.write('==========================')
        out.write('\n')
        out.write('Bounds: ')
        out.write('\n')
        count=0
        bounds_names=['a', 'b', 'start_time', 'peak_temp']
        for item in bounds:
                out.write(bounds_names[count]+': '+str(item))
                out.write('\n')
                count=count+1
                
        out.write('\n')
        out.write('===========================')
        out.write('\n')
        out.write('Data Chopping: ')
        out.write('\n')
        
# %%
import main
from pprint import pprint
import datetime
run_dict = {'s88773': [1,2,3], #Define which runs & faces to assess
            's88776': [1,2,3], 
            's88780': [1,2,3],
            's86483': [1,2,3]
            }

bounds=((0.01, 0.05), # a
        (0, 100), # b
        (1.75*10**-8, 2.25*10**-8), # start_time
        (30000, 80000)) # peak_temp

opt_data=main.OptimizationData('inter_data_stor.optimization_data.txt')
optimization_group=main.OptimizationGroup(120, run_dict, 'test1')
with open('pprint_test2', 'wt') as out:
        out.write(optimization_group.name+': Metadata')
        out.write('\n')
        out.write('Optimization initialized @: ')
        out.write('\n')
        out.write(str(datetime.datetime.now()))
        out.write('\n')
        out.write(' ')
        out.write('\n')
        out.write('===========================')
        out.write('\n')
        out.write('Runs & Faces: ')
        out.write('\n')
        pprint(optimization_group.run_dictionary, stream=out, width=1)
        out.write(' ')
        out.write('\n')
        out.write('==========================')
        out.write('\n')
        out.write('Bounds: ')
        out.write('\n')
        count=0
        bounds_names=['a', 'b', 'start_time', 'peak_temp']
        for item in bounds:
                out.write(bounds_names[count]+': '+str(item))
                out.write('\n')
                count=count+1
        out.write('\n')
        out.write('===========================')
        out.write('\n')
        out.write('Data Chopping: ')
        out.write('\n')
        out.write('\n')
        for real_data_shot_name in optimization_group.real_data.keys(): # iterate through the shot names
                for run_dict_shot_name in optimization_group.run_dictionary.keys(): # iterate through the shots specified in the run dictionary
                        if real_data_shot_name==run_dict_shot_name: #If the shot from the real data is in the run dictionary
                                out.write(str(run_dict_shot_name)+': ')
                                out.write('\n')
                                data=optimization_group.full_shot_data.start_stop[run_dict_shot_name]
                                for face in optimization_group.run_dictionary[run_dict_shot_name]: # Get the faces for each shot
                                        pretty_name='face'+str(face)
                                        out.write('     '+pretty_name+': '+str(data[face]))
                                        out.write('\n')
                                out.write('\n')
# %%
print(optimization_group.real_data.keys())
# %%
from pprint import pprint
start_stop_list=[               #s88773
                 [[179, 240],       # face 1
                  [259, 390],       # face 2
                  [310, 700]],      # face 3 
                                #s88776
                 [[88, 155],        # face 1
                  [93, 400],        # face 2
                  [246, 468]],      # face 3
                                #s88780
                 [[145, 190],       # face 1
                  [240, 350],       # face 2
                  [310, 502]],      # face 3
                                #s86483
                 [[355, 400],       # face 1
                  [582, 732],       # face 2
                  [377, 800]]]      # face 3

start_stop_dict={}
names=['s88773', 's88776', 's88780', 's86483']
count=0
for starts_and_stops in start_stop_list:
        facecount=0
        #print(starts_and_stops)
        face_dict={}
        for face in starts_and_stops:
                face_dict[facecount]
                facecount=facecount+1
        start_stop_dict[names[count]]=starts_and_stops
        count=count+1
        
pprint(start_stop_dict)
# %%
start_stop_list=[               #s88773
                 [[179, 240],       # face 1
                  [259, 390],       # face 2
                  [310, 700]],      # face 3 
                                #s88776
                 [[88, 155],        # face 1
                  [93, 400],        # face 2
                  [246, 468]],      # face 3
                                #s88780
                 [[145, 190],       # face 1
                  [240, 350],       # face 2
                  [310, 502]],      # face 3
                                #s86483
                 [[355, 400],       # face 1
                  [582, 732],       # face 2
                  [377, 800]]]      # face 3

start_stop_dict={
        's88773':{
                1:[179, 240],
                2:[259, 390],       
                3:[310, 700]
        },
        's88776':{
                1:[88, 155],        
                2:[93, 400],        
                3:[246, 468]
        },
        's88780':{
                1:[145, 190],       
                2:[240, 350],       
                3:[310, 502]
        },
        's86483':{
                1:[355, 400],       
                2:[582, 732],       
                3:[377, 800]
        }
}
        

pprint(start_stop_dict)
# %%
import main
from pprint import pprint
import datetime

run_dict = {'s88773': [1,2], #Define which runs & faces to assess
            }
bounds=((0, 0.06), # a
        (0, 90), # b
        (1.5*10**-8, 2.5*10**-8), # start_time
        (34000, 100000)) # peak_temp
opt_data=main.OptimizationData('inter_data_stor.optimization_data.txt')
group=main.OptimizationGroup(120, run_dict, 's88773_solo_tester_face1&2_V4')
#main.optimize(group)
main.optimize(group, opt_data, bounds=bounds, popsize=3)
# %%
print(group.folder_path)
# %%
run_dict = {'s88773': [1,2], #Define which runs & faces to assess
            's88776': [1], 
            's88780': [1,2],
            's86483': [1,2,3]
            }
# %%
import main
import numpy as np
def get_initial_plot(data_fp:str):
        fp=data_fp+'/optimization_output.csv'
        data=open(fp, "r").read()
        dataArray=data.split('\n')
        item=dataArray[0]
        
        vals=[]
        for k in item[1:-1].split(','):
                vals.append(float(k))
        
        initial_params=vals
                        
        #Now run the FEM model
        #Have to also get the optimization group somehow so we know how to create the model that we want to run
        #Go into metadata file and look for the relevant information
        

a=get_initial_plot('../../winhome/Desktop/optimization_data/20220315/FittingAll_NewFunction_V1')
n=a[0]
# %%
vals=[]
for item in n[1:-1].split(','):
        vals.append(float(item))
iter, sls, a_val, b_val, start_time, peak_temp=vals
type(iter)
# %%
import numpy as np

np.asarray(vals)
# %%
