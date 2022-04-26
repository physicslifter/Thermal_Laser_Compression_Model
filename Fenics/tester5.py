
# %%
import main

run_dict = {#'s88773': [1,2], #Define which runs & faces to assess
            #'s88776': [1], 
            #'s88780': [1,2],
            's86483': [1,2,3],
            's86480': [1,2,3],
            's86484': [1,2,3],
            's88774': [1,2,3]
            }

#run_dict = {'s88773': [1,2], #Define which runs & faces to assess
#            }

          #(min, max)
eq1_bounds=((0.01, 0.05), # a
        (0, 100), # b
        (1.75*10**-8, 2.25*10**-8), # start_time
        (30000, 80000)) # peak_temp

            #(min, max)go to figure out what the bounds are
            
eq2_bounds=((0.01, 0.05), # a
        (0, 100), # b
        (0,2000), #c
        (1.75*10**-8, 2.25*10**-8), # start_time
        (35000, 80000)) # peak_temp



my_data=main.Data()
my_data.get_data()
opt_output_data=main.OptimizationData('hello')
my_group=main.OptimizationGroup(60, run_dict, 'Eq1_Run_new', my_data, equation=1)
a=main.optimize(my_group, opt_output_data, eq1_bounds, 10, equation=1)

#Run end of optimization operations
q=main.PostOptOperations(my_group.folder_path)
q.run_end_operations()
# %%
