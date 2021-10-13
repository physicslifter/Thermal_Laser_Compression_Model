# Optimization runs on October 13, 2021

## Run 1
### Parameters  
    - initial parameters:  params0=[55000.0, 0.03,50.0, 3.0*10**-8]

### Initial fit:  
  ![Initial fit, run 1](run1_initial.PNG)  

### Final Fit:  
  ![Final fit, run 1](run1_final.PNG)  
  
### Parameters:
  ![Parameters, run 1](run1_params.PNG)

## Run 2
    - similar to run 1, but starting from different initial parameters

### Starting parameters:
    - params1=[45000, 0.03, 50.0, 2e-08]

### Initial fit:  
  ![Initial fit, run 2](run2_initial.PNG)  

### Final Fit:  
  ![Final fit, run 2](run2_final.PNG)  
  
### Parameters:
  ![Parameters, run 2](run2_params.png)

## Run 3

### Notes:
    - now trying to optimize for the first two faces on the first shot, shot s88773
    - using new function that holds peak_temp and start_time constant while optimizing for the two k parameters

### Parameters  
    - initial parameters:  params0=[0.03,30]

### Initial fit:  
  ![Initial fit, run 3](run3_initial.png)  

### Final Fit:  
  ![Final fit, run 3](run3_final.png)  
  
### Parameters:
  ![Parameters, run 3](run3_params.png)


## Run 4

### Notes:
  
### Parameters  
    - initial parameters:  params0=[30000,0.025,30.0,2.0*10**-8]
    - fitting with somewhat random initial parameters, and trying to get the optimization to converge

### Initial fit:  
  ![Initial fit, run 4](run4_initial.png)  

### Final Fit:  
  ![Final fit, run 4](run4_final.png)  
  
### Parameters:
  ![Parameters, run 4](run4_params.png)


# Run 5
### Notes:
  
### Parameters  
    - initial parameters:  params0=[36000.0,0.017,30.0,2.0*10**-8]
    - fitting with optimal params from earlier (specified above)

### Initial fit:  
  ![Initial fit, run 5](run5_initial.png)  

### Final Fit:  
  ![Final fit, run 5](run5_final.png)  
  
### Parameters:
  ![Parameters, run 5](run5_params.png)


# Run 6
### Notes:
  
### Parameters  
    - initial parameters:  params0=[36000.0,0.017,30.0,2.0*10**-8]
    - same as the fit above, except we are only trying to fit the first face of the second run only
    - fit works well, started with a pretty good fit, so I really only marginally improved from there

### Initial fit:  
  ![Initial fit, run 6](run6_initial.png)  

### Final Fit:  
  ![Final fit, run 6](run6_final.png)  
  
### Parameters:
  ![Parameters, run 6](run6_params.png)


# Run 7
### Notes:
Fitting the s88780 data with the same initial parameters as runs 5 & 6

### Parameters  
    - initial parameters:  params0=[36000.0,0.017,30.0,2.0*10**-8]


### Initial fit:  
  ![Initial fit, run 7](run7_initial.png)  

### Final Fit:  
  ![Final fit, run 7](run7_final.png)  
  
### Parameters:
  ![Parameters, run 7](run7_params.png)


# Run 8
### Notes:
Fitting the s86483 data with the same initial parameters as runs 5, 6 & 7

### Parameters  
    - initial parameters:  params0=[36000.0,0.017,30.0,2.0*10**-8]


### Initial fit:  
  ![Initial fit, run 8](run8_initial.png)  

### Final Fit:  
  ![Final fit, run 8](run8_final.png)  
  
### Parameters:
  ![Parameters, run 8](run8_params.png)