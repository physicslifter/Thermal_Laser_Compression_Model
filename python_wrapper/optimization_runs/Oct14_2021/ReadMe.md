# Optimization Runs on October 14, 2021

## Run 1

### Notes:  
    - Fitting face 1 for s88776, using the 2-input function (all 4 shots, only the k-parameters as input)
    - failed because I started with an a-parameter that was too high by about an order of magnitude
    - fixed Matlab problem. I tried optimizing this same function yesterday (Run 3 from October 13). Changing a & b parameters had no result on the function, and so I temporarily abandoned the technique. It was later revealed that 

### Initial Params:
    - params0=[0.5,10]


## Run 2

### NotesL
    - same as the previous run, but with the parameters fixed to new values

### Initial Params:
    - params0=[0.05,10]

## Run 3

### Notes 
    - now fitting all runs holding peak_temp and start_time constant

## Run 4 

### Notes
    - Fitting all runs, but with higher peak_temp and slightly lower start_time. 
    - run aborted because it tried sending a value of a that was too low and negative, and the run quit

### Starting parameters

### Parameters:
  ![Parameters, run 5](run5_params.PNG)


## Run 5

### Notes
    - back to the simple_four_run technique
    - I think there is a reasonable likelihood of this approach working, because I've now fixed the problem that was causing peak_temp and start_time to be the only parameters that change (a and b held constant in the .m file)
        - this is releant, because

# Starting parameters:
    initial params: params0=[20000.0,0.05,10,2.5*10**-8]

### Parameters:
  ![Parameters, run 5](run5_params.PNG)

