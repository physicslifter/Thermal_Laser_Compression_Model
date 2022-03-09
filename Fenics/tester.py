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
while count<20:
    bounds=((0,0.05),(0,100),(1.75*10**-8, 2.25*10**-8), (30000,100000)) #define bounds (required for changing optimization params)
    opt_func.rewrite_global_params(num_steps, shots, 'opt_func_time_test', bounds, 10) #rewrite the global parameters
    opt_func.optimizable_function(params)
    count=count+1
# %%
