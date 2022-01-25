# Optimization runs on October 12, 2021

## Run 1  

### Notes
Since I got the chopped data to converge for the new run (83), I am now trying fitting all data runs at once to see if this converges to a decent value. 

### Initial parameters:  
     params0=[36000.0,0.017,30,2.0*10**-8]  
  
### Initial fit:  
  ![Initial fit, run 1](run1_initial.png)  

### Final Fit:  
  ![Final fit, run 1](run1_final.png)  
  
### Parameters:
  ![Parameters, run 1](run1_parameters.png)


## Run 2

### Notes
Fitting all runs starting from totally random parameters (as opposed to the "good fit" parameters that we have been using)

- it turns out that I messed up the a parameter by making it too large by an order of magnitude. The fit did not converge towards a good/reasonable 

### Initial parameters:  
      params1=[12000.0,0.38,48,1.2*10**-8] 
  
### Initial fit:  
  ![Initial fit, run 2](run2_initial.png)  

### Final Fit:  
  ![Final fit, run 2](run2_final.png)  
  
### Parameters:
  ![Parameters, run 2](run2_parameters.png)

## Run 3

### Notes
Fitting all runs starting from totally random parameters (as opposed to the "good fit" parameters that we have been using)

- same as run 2, except the starting a parameter is 10x the value
- still does not converge to a good value..
    - start_time is way too low, peak_temp is too low
    - need to do next run, with higher start_time as an initial parameter

### Initial parameters:  
      params1=[12000.0,0.038,48,1.2*10**-8] 
  
### Initial fit:  
  ![Initial fit, run 3](run3_initial.png)  

### Final Fit:  
  ![Final fit, run 3](run3_final.png)  
  
### Parameters:
  ![Parameters, run 3](run3_parameters.png)

## Run 4