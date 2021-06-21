def say_hello(name):
    return 'hello, '+str(name)+'!'
    print('hello, '+str(name)+'!')

#import matlab engine package
def fem_run(input_matrix,run_name):
    eng=matlab.engine.start_matlab()

    #Run the matlab function with the associated inputs
    eng.RunFEM(input_matrix,run_name,nargout=0)

    #exit/quit the engine
    eng.quit()

    #show user the run is complete
    print('run complete')

