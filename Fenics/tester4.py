#%%
import main

run_dict = {'s88773': [1,2], #Define which runs & faces to assess
            's88776': [1], 
            's88780': [1,2],
            #'s86483': [1,2,3],
            #'s86480': [1,2,3],
            #'s86484': [1,2,3],
            #'s88774': [1,2,3]
            }


          #(min, max)
eq1_bounds=((0.01, 0.05), # a
        (0, 100), # b
        (1.75*10**-8, 2.25*10**-8), # start_time
        (30000, 80000)) # peak_temp

            #(min, max)
eq2_bounds=((0.01, 0.05), # a
        (0, 100), # b
        (0,10**-5), #c
        (1.75*10**-8, 2.25*10**-8), # start_time
        (35000, 80000)) # peak_temp

my_data=main.Data()
my_data.get_data()
opt_output_data=main.OptimizationData('hello')
my_group=main.OptimizationGroup(120, run_dict, 'Eq3_FittingAll_Popsize10', my_data, equation=3)
main.optimize(my_group, opt_output_data, eq2_bounds, 20, equation=3)

#Run end of optimization operations
q=main.PostOptOperations(my_group.folder_path)
q.run_end_operations()
# %%
