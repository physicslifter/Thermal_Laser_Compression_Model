import json
from tokenize import Double
import numpy as np

def eval_to_pop_array(input_file, num_params, popsize):
    pullData=open(input_file, 'r').read()
    dataArray = pullData.split('\n')
    data_length=len(dataArray)
    generation_size=num_params*popsize
    num_generations=int((data_length-1)/generation_size)+1
    pop_array=np.zeros((num_generations, generation_size, num_params+2))
    
    count=1
    for line in dataArray:
        if line !='':
            
            it, sls, a, b, pt, ts = line.split(',')
            it=float(it)
            sls=float(sls)
            a=float(a)
            b=float(b)
            pt=float(pt)
            ts=float(ts)
            if it==-1:
                it=0
            else:
                it=it
            eval_no_in_pop=int(it%generation_size)
            generation_num=int(it/generation_size)
            pop_array[generation_num][eval_no_in_pop][0]=it
            pop_array[generation_num][eval_no_in_pop][1]=sls
            pop_array[generation_num][eval_no_in_pop][2]=a
            pop_array[generation_num][eval_no_in_pop][3]=b
            pop_array[generation_num][eval_no_in_pop][4]=pt
            pop_array[generation_num][eval_no_in_pop][5]=ts
        
    return pop_array

def pop_to_gen(pop_array):
    generation_array=np.zeros((np.shape(pop_array)[0],np.shape(pop_array)[2])) #create generation array as an array of zeros with length equal to the number of generations
    count1=0
    for generation in pop_array:
        count2=0
        for x in generation:
            it, sls, a, b, peak_temp, start_time = x
            if count2==0:
                    index_lowest=count2
            else:
                    if generation[index_lowest][0]>sls:
                            index_lowest=count2
            count2+=1
        generation_array[count1]=generation[index_lowest]
        generation_array[count1][0]=count1
        
        count1+=1
        
    return generation_array[:-1]

def eval_to_gen(input_file, num_params, popsize):
    pop=eval_to_pop_array(input_file, num_params, popsize)
    gen=pop_to_gen(pop)
    np.savetxt('gen.csv', gen, delimiter= ",")
    