#%%
import main

#a = 0.025
#b = 100
#pt = 30kK
#st = 2.1e_8

run_dict = {'s88773': [1,2], #Define which runs & faces to assess
            's88776': [1], 
            's88780': [1,2, 3],
            's86483': [1,2,3],
            }

#run_dict = {'s88773': [1,2], #Define which runs & faces to assess
#            }

          #(min, max)
eq1_bounds=((0.01, 0.05), # a
        (0, 100), # b
        (1.75*10**-8, 2.25*10**-8), # start_time
        (30000, 80000)) # peak_temp

my_data = main.Data()
my_data.get_data()
opt_output_data=main.OptimizationData('hello')
my_group=main.OptimizationGroup(60, run_dict, 'Eq1_FittingAll_Popsize10', my_data, equation=1)

params = [
    0.0245, #a
    80, #b
    2.16*10**-8, #start time
    30300 #peak temp
]

my_group.run(params)
my_group.plot_optimization_group("../../winhome/Desktop/optimization_data/20220426/Eq1_FittingAll_Popsize10/params_output2.png")
# %%
