# Thermal_Laser_Compression_Model

Finite Element Model (FEM) high temperature thermal laser compression for deriving conductivity and diffusivity.


## Dependencies
Python 3.8 is required. Currently, this model uses the python engine for matlab to run the MATLAB FEM in python. The engine is only compatible with these versions of Python.

## Tutorial
once you have installed dependencies and have a current version of the repository:
1. Open command line
2. navigate to python_wrapper subfolder
3. to run the model:
    from scipy import io
    from matlab import engine
    import runFEM2 as f
    from matplotlib import pyplot as plt
    
    #set default parameters
    #currently, the only parameter is the peak temperature, and so we will set it accordingly
    my_peak_temp=28000 #default from .m script
    
    #run the FEM, saving the output as mydata
    mydata=f.run_model(my_peak_temp)
    
    #plot the results
    f.plot_results(mydata)
    plt.show()
    
