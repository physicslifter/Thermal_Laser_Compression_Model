#%%
import main

run_dict = {'s88773': [1,2], #Define which runs & faces to assess
            's88776': [1], 
            's88780': [1,2,3],
            #'s86483': [1,2,3],
            #'s86480': [1,2,3],
            #'s86484': [1,2,3],
            #'s88774': [1,2,3]
            }

#run_dict = {'s88773': [1,2], #Define which runs & faces to assess
            #}

          #(min, max)
eq1_bounds=((0.01, 0.05), # a
        (0, 100), # b
        (1.75*10**-8, 2.25*10**-8), # start_time
        (30000, 80000)) # peak_temp

            #(min, max)
eq2_bounds=((0.01, 0.05), # a
        (0, 100), # b
        (0,2000), #c
        (1.75*10**-8, 2.25*10**-8), # start_time
        (35000, 80000)) # peak_temp



my_data=main.Data()
my_data.get_data()
opt_output_data=main.OptimizationData('hello')
my_group=main.OptimizationGroup(60, run_dict, 'Eq2_Running73&76&80_Popsize=10', my_data, equation=2)
a=main.optimize(my_group, opt_output_data, eq2_bounds, 10, equation=2)

#%%
import main
#Run end of optimization operations
q=main.PostOptOperations('../../winhome/Desktop/optimization_data/20220415/Eq2_Running73&76&80_Popsize=10')
q.run_end_operations()
# %%

# %%
import main

run_dict = {'s88773': [1,2], #Define which runs & faces to assess
            #'s88776': [1,2,3], 
            #'s88780': [1,2,3],
            #'s86483': [1,2,3],
            #'s86480': [1,2,3],
            #'s86484': [1,2,3],
            #'s88774': [1,2,3]
            }

eq1_default_params=[0.03, # a
        30, # b
        1.8*10**-8, # start_time
        30000 # peak_temp
        ]
eq2_default_params=[0.03, # a
        30, # b
        100, #c
        1.8*10**-8, # start_time
        30000 # peak_temp
        ]

my_data=main.Data()
group1=main.OptimizationGroup(120, run_dict, '', my_data, equation=1)
group2=main.OptimizationGroup(120, run_dict, '', my_data, equation=2)

# %%
count=0
while count<20:
    group1.run(eq1_default_params)
    count+=1
# %%
count=0
while count<20:
    group2.run(eq2_default_params)
    count+=1
# %%
